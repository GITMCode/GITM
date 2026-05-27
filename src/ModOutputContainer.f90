! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE
!
! ModOutputContainer.f90
!
! Containerized outputs for GITM.  Each output type defines its
! variables once (define_var, schema phase), then fills them per timestep
! (put, data phase).  Backends consume the container.
!
! Two-phase design:
!   Phase 1 — SCHEMA: call define_var() once per file open, all ranks.
!   Phase 2 — DATA:   call prepare(iBlock), then put() on participating ranks.
!
! Container precision = kind(1.0) — tracks compile-time `real`.  A double
! build (real = real*8) produces 8-byte output; a single build produces 4.
! This matches mpiio's MPI_REAL convention and avoids silent precision loss.

module ModOutputContainer

  implicit none

  ! Container storage precision tracks the model's compile-time precision
  integer, parameter :: output_kind = kind(1.0)

  ! Grid kind constants — controls how backends interpret dims and participation.
  integer, parameter :: GRID_GEO_3D = 1  ! geographic 3D block-partitioned
  integer, parameter :: GRID_GEO_2D = 2  ! geographic 2D (integrated, e.g. TEC)
  integer, parameter :: GRID_MAG_2D = 3  ! magnetic 2D (e.g. 2DMEL, only iBlock==1 writes)
  integer, parameter :: GRID_HIME   = 4  ! regional (only DoSaveHIMEPlot ranks write)

  integer, parameter :: MaxContainerVars = 300

  type :: Variable
    character(len=64)  :: name     = ''
    character(len=128) :: longName = ''
    character(len=32)  :: units    = ''
    logical            :: is_axis  = .false.  ! netcdf: write as coordinate variable
    integer            :: rank     = 0        ! 0..3
    integer            :: shape3(3) = [1, 1, 1]  ! padded with 1s; drives allocation
    logical            :: defined  = .false.
    real(output_kind), allocatable :: data(:, :, :)
  end type Variable

  type :: OutputContainer
    character(len=5)            :: cType     = '     '
    integer                     :: gridKind  = GRID_GEO_3D
    type(Variable), allocatable :: vars(:)
    integer                     :: nVars     = 0
    logical                     :: this_rank_writes = .false.
    logical                     :: any_rank_skips   = .false.  ! set at file open via allreduce
    logical                     :: schema_locked    = .false.  ! true after open_file
  contains
    procedure :: define_var
    procedure :: prepare
    procedure :: reset
    procedure :: find_var
    procedure :: put_0d
    procedure :: put_1d
    procedure :: put_2d
    procedure :: put_3d
    generic   :: put => put_0d, put_1d, put_2d, put_3d
  end type OutputContainer

contains

  ! ---------------------------------------------------------------------------
  ! define_var — register a variable in the schema.
  !
  ! Must be called before open_file (schema phase).  All ranks call identically.
  ! Allocates data storage (reused across timesteps).
  ! ---------------------------------------------------------------------------
  subroutine define_var(this, name, units, shape3, longName, is_axis)
    class(OutputContainer), intent(inout) :: this
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: units
    integer,          intent(in) :: shape3(3)
    character(len=*), intent(in), optional :: longName
    logical,          intent(in), optional :: is_axis

    type(Variable), allocatable :: tmp(:)
    integer :: n

    if (this%schema_locked) then
      write(*, '(a,a,a)') 'ModOutputContainer: WARNING define_var called after schema locked for ', &
        trim(this%cType), ' — ignored'
      return
    end if

    ! Grow vars array by one
    n = this%nVars + 1
    if (.not. allocated(this%vars)) then
      allocate(this%vars(MaxContainerVars))
    end if

    this%vars(n)%name   = trim(name)
    this%vars(n)%units  = trim(units)
    this%vars(n)%shape3 = shape3
    this%vars(n)%rank   = count(shape3 > 1)  ! non-degenerate dims
    this%vars(n)%defined = .true.

    if (present(longName)) then
      this%vars(n)%longName = trim(longName)
    else
      this%vars(n)%longName = trim(name)
    end if

    if (present(is_axis)) then
      this%vars(n)%is_axis = is_axis
    else
      this%vars(n)%is_axis = .false.
    end if

    ! Allocate storage once; reused across timesteps.
    if (.not. allocated(this%vars(n)%data)) then
      allocate(this%vars(n)%data(shape3(1), shape3(2), shape3(3)))
    end if
    this%vars(n)%data = 0.0_output_kind

    this%nVars = n
  end subroutine define_var

  ! ---------------------------------------------------------------------------
  ! prepare — set this_rank_writes for this block, per container's gridKind.
  !
  ! Called by the producer at the start of each fill_<type> call.
  ! ---------------------------------------------------------------------------
  subroutine prepare(this, iBlock)
    use ModGITM, only: iProc
    class(OutputContainer), intent(inout) :: this
    integer, intent(in) :: iBlock

    select case (this%gridKind)
    case (GRID_GEO_3D, GRID_GEO_2D)
      this%this_rank_writes = .true.
    case (GRID_MAG_2D)
      ! Magnetic grid data is gathered to rank 0; iBlock is always 1 per rank.
      this%this_rank_writes = (iProc == 0)
    case (GRID_HIME)
      ! HIME participation is a per-rank flag set externally; do not override here.
      ! Caller must set this_rank_writes before calling put().
    case default
      this%this_rank_writes = .true.
    end select
  end subroutine prepare

  ! ---------------------------------------------------------------------------
  ! reset — clear participation flag between timesteps.  Keeps allocations.
  ! ---------------------------------------------------------------------------
  subroutine reset(this)
    class(OutputContainer), intent(inout) :: this
    this%this_rank_writes = .false.
    this%any_rank_skips   = .false.
  end subroutine reset

  ! ---------------------------------------------------------------------------
  ! find_var — return index of named variable, or 0 if not found.
  ! ---------------------------------------------------------------------------
  integer function find_var(this, name)
    class(OutputContainer), intent(in) :: this
    character(len=*), intent(in) :: name
    integer :: i

    find_var = 0
    do i = 1, this%nVars
      if (trim(this%vars(i)%name) == trim(name)) then
        find_var = i
        return
      end if
    end do
  end function find_var

  ! ---------------------------------------------------------------------------
  ! put_0d — store a scalar value.
  ! ---------------------------------------------------------------------------
  subroutine put_0d(this, name, scalar)
    class(OutputContainer), intent(inout) :: this
    character(len=*), intent(in) :: name
    real(output_kind), intent(in) :: scalar
    integer :: idx

    if (.not. this%this_rank_writes) return
    idx = this%find_var(name)
    if (idx == 0) then
      write(*, '(a,a,a,a)') 'ModOutputContainer: put_0d unknown variable "', &
        trim(name), '" in type ', trim(this%cType)
      return
    end if
    this%vars(idx)%data(1, 1, 1) = scalar
  end subroutine put_0d

  ! ---------------------------------------------------------------------------
  ! put_1d — store a 1D array.
  ! ---------------------------------------------------------------------------
  subroutine put_1d(this, name, arr)
    class(OutputContainer), intent(inout) :: this
    character(len=*), intent(in) :: name
    real(output_kind), intent(in) :: arr(:)
    integer :: idx, n

    if (.not. this%this_rank_writes) return
    idx = this%find_var(name)
    if (idx == 0) then
      write(*, '(a,a,a,a)') 'ModOutputContainer: put_1d unknown variable "', &
        trim(name), '" in type ', trim(this%cType)
      return
    end if
    n = size(arr)
    if (n /= this%vars(idx)%shape3(1)) then
      write(*, '(a,a,a,i0,a,i0)') 'ModOutputContainer: put_1d shape mismatch for "', &
        trim(name), '": expected ', this%vars(idx)%shape3(1), ', got ', n
      return
    end if
    this%vars(idx)%data(:, 1, 1) = arr
  end subroutine put_1d

  ! ---------------------------------------------------------------------------
  ! put_2d — store a 2D array.
  ! ---------------------------------------------------------------------------
  subroutine put_2d(this, name, arr)
    class(OutputContainer), intent(inout) :: this
    character(len=*), intent(in) :: name
    real(output_kind), intent(in) :: arr(:, :)
    integer :: idx, s1, s2

    if (.not. this%this_rank_writes) return
    idx = this%find_var(name)
    if (idx == 0) then
      write(*, '(a,a,a,a)') 'ModOutputContainer: put_2d unknown variable "', &
        trim(name), '" in type ', trim(this%cType)
      return
    end if
    s1 = size(arr, 1); s2 = size(arr, 2)
    if (s1 /= this%vars(idx)%shape3(1) .or. s2 /= this%vars(idx)%shape3(2)) then
      write(*, '(a,a,a,2(i0,a),a,2(i0,a))') &
        'ModOutputContainer: put_2d shape mismatch for "', trim(name), '": expected ', &
        this%vars(idx)%shape3(1), 'x', this%vars(idx)%shape3(2), ', got ', &
        s1, 'x', s2, ''
      return
    end if
    this%vars(idx)%data(:, :, 1) = arr
  end subroutine put_2d

  ! ---------------------------------------------------------------------------
  ! put_3d — store a 3D array.
  ! ---------------------------------------------------------------------------
  subroutine put_3d(this, name, arr)
    class(OutputContainer), intent(inout) :: this
    character(len=*), intent(in) :: name
    real(output_kind), intent(in) :: arr(:, :, :)
    integer :: idx, s1, s2, s3

    if (.not. this%this_rank_writes) return
    idx = this%find_var(name)
    if (idx == 0) then
      write(*, '(a,a,a,a)') 'ModOutputContainer: put_3d unknown variable "', &
        trim(name), '" in type ', trim(this%cType)
      return
    end if
    s1 = size(arr, 1); s2 = size(arr, 2); s3 = size(arr, 3)
    if (s1 /= this%vars(idx)%shape3(1) .or. &
        s2 /= this%vars(idx)%shape3(2) .or. &
        s3 /= this%vars(idx)%shape3(3)) then
      write(*, '(a,a,a,3(i0,a),a,3(i0,a))') &
        'ModOutputContainer: put_3d shape mismatch for "', trim(name), '": expected ', &
        this%vars(idx)%shape3(1), 'x', this%vars(idx)%shape3(2), 'x', &
        this%vars(idx)%shape3(3), ', got ', s1, 'x', s2, 'x', s3, ''
      return
    end if
    this%vars(idx)%data = arr
  end subroutine put_3d

end module ModOutputContainer
