! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE
!
! ModOutputContainer.f90
!
! Containerized outputs for GITM.  Backends consume the container; both APIs
! below leave it in the same state, so no backend can tell which built it.
!
! ONE-SHOT (put_field/put_axis) is preferred: a single fill_<type> both defines
!   and fills.  The first put of a name defines the variable (shape = size(arr)
!   VERBATIM, metadata from ModOutputVarRegistry); later puts must match.  A
!   priming pass (init_output_containers runs each fill once, then locks the
!   schema) makes the schema complete before open_file, where netcdf issues
!   collective nfmpi_def_var for the whole file.
!   The CALLER owns the slice: put stores exactly what it is handed, with no
!   ghost-cell or grid-kind inference (3DALL puts 2-ghost full blocks; 2DANC
!   puts 0-ghost interior slices plus one altitude LEVEL of a 3-D array).
!   nGhostCells is metadata (what the netcdf backend strips), decoupled from
!   the stored slice.
!
! EXPLICIT TWO-PHASE (define_var once per open on all ranks, then prepare +
!   put on participating ranks): for magnetic grids (2DMEL), HME/HIME regionals,
!   the rank-0 probes 1DGLO/1DTHM/1DCHM/1DEMI, and USR.  USR is permanent by design:
!   users put from anywhere, so a put-derived schema would be empty at the
!   first output; registration IS its schema.
!
! Container precision = kind(1.0), tracking compile-time `real`, matching
! mpiio's MPI_REAL convention.

module ModOutputContainer

  use ModOutputVarRegistry, only: lookup_var_meta

  implicit none

  ! Container storage precision tracks the model's compile-time precision
  integer, parameter :: output_kind = kind(1.0)

  ! Grid kind constants: controls how backends interpret dims and participation.
  ! Participation rules are set in prepare(): GEO_3D/GEO_2D = all ranks write,
  ! MAG_2D = only iProc==0, HIME = caller sets per-rank flag externally.
  ! GRID_GEO_2D doubles as the "all geographic ranks participate" default for
  ! 1D and 0D outputs (no dedicated GRID_GEO_1D/_0D; participation rule is
  ! identical, and the dim count is carried by each variable's shape3).
  integer, parameter :: GRID_GEO_3D = 1  ! geographic 3D block-partitioned
  integer, parameter :: GRID_GEO_2D = 2  ! geographic 2D / 1D / 0D (all ranks participate)
  integer, parameter :: GRID_MAG_2D = 3  ! magnetic 2D (e.g. 2DMEL); iProc==0 only
  integer, parameter :: GRID_HIME = 4  ! regional (only DoSaveHIMEPlot ranks write)

  integer, parameter :: MaxContainerVars = 300

  type :: Variable
    character(len=64)  :: name = ''
    character(len=128) :: longName = ''
    character(len=32)  :: units = ''
    logical            :: is_axis = .false.  ! netcdf: write as coordinate variable
    integer            :: rank = 0        ! 0..3
    integer            :: shape3(3) = [1, 1, 1]  ! padded with 1s; drives allocation
    logical            :: defined = .false.
    real(output_kind), allocatable :: data(:, :, :)
  end type Variable

  type :: OutputContainer
    character(len=5)            :: cType = '     '
    integer                     :: gridKind = GRID_GEO_3D
    type(Variable), allocatable :: vars(:)
    integer                     :: nVars = 0
    logical                     :: this_rank_writes = .false.
    logical                     :: any_rank_skips = .false.  ! set at file open via allreduce
    ! Once primed, a put of an unknown name is a hard error: a per-rank or
    ! per-timestep schema change would desync the collective def_var.
    logical                     :: schema_locked = .false.
    ! Drives the per-fill wipe and put-count guards. Off for unconverted types;
    ! USR in particular needs its values to persist between puts.
    logical                     :: is_one_shot = .false.
    ! Puts since the last prepare(); -1 when no one-shot fill is active. The
    ! sentinel lets finish_fill run from one place after the dispatch chain.
    integer                     :: nPuts = -1
    ! Per-type ghost cell width (per side); 0 = interior-only. Set to 2 to keep
    ! GITM's ghosts in legacy/mpiio; netcdf strips them for GEO_* gridKinds.
    integer                     :: nGhostCells = 0
    ! Per-type fill routine, registered in init_output_containers. fill_iface is
    ! declared after this type: an interface importing OutputContainer cannot
    ! precede it, and gfortran resolves the forward reference within one spec part.
    procedure(fill_iface), pointer, nopass :: fill => null()
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
    ! One-shot API.  A generic works here (unlike a shape selector) because the
    ! array's RANK is part of the signature, so the compiler dispatches on it.
    procedure :: put_field_1d
    procedure :: put_field_2d
    procedure :: put_field_3d
    generic   :: put_field => put_field_1d, put_field_2d, put_field_3d
    procedure :: put_axis => put_axis_1d
    procedure :: finish_fill
    procedure, private :: ensure_var
  end type OutputContainer

  abstract interface
    subroutine fill_iface(c, iBlock)
      import :: OutputContainer
      type(OutputContainer), intent(inout) :: c
      integer, intent(in) :: iBlock
    end subroutine fill_iface
  end interface

contains

  ! ---------------------------------------------------------------------------
  ! define_var: register a variable in the schema.
  !
  ! Must be called before open_file (schema phase).  All ranks call identically.
  ! Allocates data storage (reused across timesteps).
  ! ---------------------------------------------------------------------------
  subroutine define_var(this, name, units, shape3, longName, is_axis)
    class(OutputContainer), intent(inout) :: this
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: units
    integer, intent(in) :: shape3(3)
    character(len=*), intent(in), optional :: longName
    logical, intent(in), optional :: is_axis

    integer :: n
    character(len=32)  :: regUnits
    character(len=128) :: regLongName
    logical :: inRegistry

    if (this%schema_locked) then
      write(*, '(a,a,a)') 'ModOutputContainer: WARNING define_var called after schema locked for ', &
        trim(this%cType), ': ignored'
      return
    endif

    if (this%nVars >= MaxContainerVars) &
      call stop_gitm('ModOutputContainer: too many variables in type '//trim(this%cType)// &
                     '; raise MaxContainerVars')

    n = this%nVars + 1
    if (.not. allocated(this%vars)) then
      allocate(this%vars(MaxContainerVars))
    endif

    this%vars(n)%name = trim(name)
    this%vars(n)%units = trim(units)
    this%vars(n)%shape3 = shape3
    this%vars(n)%rank = count(shape3 > 1)  ! non-degenerate dims
    this%vars(n)%defined = .true.

    if (present(longName)) then
      this%vars(n)%longName = trim(longName)
    else
      ! No explicit longName: adopt the registry's canonical text when the name
      ! is registered, matching ensure_var, so define-path types cannot drift
      ! from the central registry.  Registry units are discarded here; units is a
      ! required argument and is never overridden.
      regUnits = ''
      regLongName = ''
      inRegistry = lookup_var_meta(name, regUnits, regLongName)
      if (inRegistry) then
        this%vars(n)%longName = trim(regLongName)
      else
        this%vars(n)%longName = trim(name)
      endif
    endif

    if (present(is_axis)) then
      this%vars(n)%is_axis = is_axis
    else
      this%vars(n)%is_axis = .false.
    endif

    ! Allocate storage once; reused across timesteps.
    if (.not. allocated(this%vars(n)%data)) then
      allocate(this%vars(n)%data(shape3(1), shape3(2), shape3(3)))
    endif
    this%vars(n)%data = 0.0_output_kind

    this%nVars = n
  end subroutine define_var

  ! ---------------------------------------------------------------------------
  ! ensure_var: resolve name to a variable index, defining it on first touch.
  !
  ! First touch: append the variable, taking its shape verbatim from the caller
  ! and its metadata from the registry (units=/longName= override it, which is
  ! how the computed species/ion names carry metadata a static table cannot).
  ! Later touches: return the existing index after checking the shape still
  ! matches; a fill whose slice changes size mid-run would silently corrupt
  ! every subsequent file.
  !
  ! Leaves exactly the state define_var would, so backends cannot distinguish
  ! the two paths.
  ! ---------------------------------------------------------------------------
  subroutine ensure_var(this, name, shape3, is_axis, idx, units, longName)
    class(OutputContainer), intent(inout) :: this
    character(len=*), intent(in) :: name
    integer, intent(in) :: shape3(3)
    logical, intent(in) :: is_axis
    integer, intent(out) :: idx
    character(len=*), intent(in), optional :: units
    character(len=*), intent(in), optional :: longName

    character(len=32)  :: varUnits
    character(len=128) :: varLongName
    logical :: inRegistry
    character(len=160) :: errmsg

    if (this%nPuts >= 0) this%nPuts = this%nPuts + 1  ! -1 = not counting (see prepare)
    idx = this%find_var(name)

    if (idx > 0) then
      ! Overwrite path: shape is fixed at first touch.
      if (any(this%vars(idx)%shape3 /= shape3)) then
        write(errmsg, '(a,a,a,a,a,3(i0,a),3(i0,a))') &
          'ModOutputContainer: put shape mismatch for "', trim(name), '" in type ', trim(this%cType), &
          ': defined as ', this%vars(idx)%shape3(1), 'x', this%vars(idx)%shape3(2), 'x', &
          this%vars(idx)%shape3(3), ', got ', shape3(1), 'x', shape3(2), 'x', shape3(3), ''
        call stop_gitm(trim(errmsg))
      endif
      return
    endif

    ! ---- first touch: define ----
    if (this%schema_locked) &
      call stop_gitm('ModOutputContainer: put of unknown variable "'//trim(name)//'" in type '// &
                     trim(this%cType)//' after the schema was primed.  A converted fill_<type> must '// &
                     'put the same variables on every rank and every timestep; no conditional puts.')

    if (this%nVars >= MaxContainerVars) &
      call stop_gitm('ModOutputContainer: too many variables in type '//trim(this%cType)// &
                     '; raise MaxContainerVars')

    if (.not. allocated(this%vars)) allocate(this%vars(MaxContainerVars))

    varUnits = ''
    varLongName = ''
    inRegistry = lookup_var_meta(name, varUnits, varLongName)
    if (.not. inRegistry) then
      ! Not a registered name (a computed species/ion name): the call site owns
      ! the metadata.  longName falls back to the name, matching define_var, but
      ! units cannot be guessed; refuse rather than emit empty units.
      varLongName = trim(name)
      if (.not. present(units)) &
        call stop_gitm('ModOutputContainer: variable "'//trim(name)//'" in type '//trim(this%cType)// &
                       ' has no ModOutputVarRegistry entry; pass units= at the put site or register it')
    endif
    if (present(units)) varUnits = trim(units)
    if (present(longName)) varLongName = trim(longName)

    idx = this%nVars + 1
    this%vars(idx)%name = trim(name)
    this%vars(idx)%units = varUnits
    this%vars(idx)%longName = varLongName
    this%vars(idx)%shape3 = shape3
    this%vars(idx)%rank = count(shape3 > 1)  ! non-degenerate dims, as define_var
    this%vars(idx)%is_axis = is_axis
    this%vars(idx)%defined = .true.

    if (.not. allocated(this%vars(idx)%data)) then
      allocate(this%vars(idx)%data(shape3(1), shape3(2), shape3(3)))
    endif
    this%vars(idx)%data = 0.0_output_kind

    this%nVars = idx
  end subroutine ensure_var

  ! ---------------------------------------------------------------------------
  ! finish_fill: assert a one-shot fill put EVERY variable exactly once.
  !
  ! This is the primary guard against a stale variable: without it, a fill that
  ! quietly stops putting some name keeps shipping whatever was in %data (with
  ! prepare's wipe, zeros), and nothing in the file says so.  It equally catches
  ! a double put, and a fill that puts nothing at all.
  !
  ! Called from ONE place (output_common, after the fill/write dispatch chain)
  ! so that adding an output type cannot forget it.  It therefore runs after
  ! reset() has cleared this_rank_writes, which is why participation is read
  ! from the nPuts sentinel (armed in prepare) rather than re-tested here.
  ! Unconverted types and non-writing ranks sit at -1 and are skipped.
  ! ---------------------------------------------------------------------------
  subroutine finish_fill(this)
    class(OutputContainer), intent(inout) :: this
    character(len=160) :: errmsg

    if (this%nPuts < 0) return  ! no one-shot fill ran on this rank

    if (this%nPuts /= this%nVars) then
      write(errmsg, '(a,a,a,i0,a,i0,a)') &
        'ModOutputContainer: fill_', trim(this%cType), ' made ', this%nPuts, &
        ' put(s) but the type has ', this%nVars, &
        ' variable(s); every variable must be put exactly once, unconditionally.'
      call stop_gitm(trim(errmsg))
    endif

    ! Disarm so a finish_fill reached without an intervening prepare() cannot
    ! pass on this fill's stale count.
    this%nPuts = -1
  end subroutine finish_fill

  ! ---------------------------------------------------------------------------
  ! put_axis_1d: store a 1D coordinate variable (netcdf writes it as a
  ! coordinate).  Otherwise identical to put_field_1d.
  ! ---------------------------------------------------------------------------
  subroutine put_axis_1d(this, name, arr, units, longName)
    class(OutputContainer), intent(inout) :: this
    character(len=*), intent(in) :: name
    real, intent(in) :: arr(:)
    character(len=*), intent(in), optional :: units
    character(len=*), intent(in), optional :: longName
    integer :: idx

    if (.not. this%this_rank_writes) return
    call this%ensure_var(name, [size(arr), 1, 1], .true., idx, units, longName)
    this%vars(idx)%data(:, 1, 1) = real(arr, output_kind)
  end subroutine put_axis_1d

  ! ---------------------------------------------------------------------------
  ! put_field_1d/2d/3d: define-on-first-touch + store, verbatim.
  !
  ! The dummy is default real (the model's kind), NOT output_kind: taking model
  ! reals lets the r8 -> output_kind conversion happen in the same assignment
  ! that stores, so there is exactly ONE copy and call sites need no
  ! real(..., output_kind) wrapper (which would itself materialize a temp).
  ! ---------------------------------------------------------------------------
  subroutine put_field_1d(this, name, arr, units, longName)
    class(OutputContainer), intent(inout) :: this
    character(len=*), intent(in) :: name
    real, intent(in) :: arr(:)
    character(len=*), intent(in), optional :: units
    character(len=*), intent(in), optional :: longName
    integer :: idx

    if (.not. this%this_rank_writes) return
    call this%ensure_var(name, [size(arr), 1, 1], .false., idx, units, longName)
    this%vars(idx)%data(:, 1, 1) = real(arr, output_kind)
  end subroutine put_field_1d

  subroutine put_field_2d(this, name, arr, units, longName)
    class(OutputContainer), intent(inout) :: this
    character(len=*), intent(in) :: name
    real, intent(in) :: arr(:, :)
    character(len=*), intent(in), optional :: units
    character(len=*), intent(in), optional :: longName
    integer :: idx

    if (.not. this%this_rank_writes) return
    call this%ensure_var(name, [size(arr, 1), size(arr, 2), 1], .false., idx, units, longName)
    this%vars(idx)%data(:, :, 1) = real(arr, output_kind)
  end subroutine put_field_2d

  subroutine put_field_3d(this, name, arr, units, longName)
    class(OutputContainer), intent(inout) :: this
    character(len=*), intent(in) :: name
    real, intent(in) :: arr(:, :, :)
    character(len=*), intent(in), optional :: units
    character(len=*), intent(in), optional :: longName
    integer :: idx

    if (.not. this%this_rank_writes) return
    call this%ensure_var(name, [size(arr, 1), size(arr, 2), size(arr, 3)], .false., idx, units, longName)
    this%vars(idx)%data = real(arr, output_kind)
  end subroutine put_field_3d

  ! ---------------------------------------------------------------------------
  ! prepare: set this_rank_writes for this block, per container's gridKind.
  !
  ! Called by the producer at the start of each fill_<type> call.
  ! ---------------------------------------------------------------------------
  subroutine prepare(this, iBlock)
    use ModGITM, only: iProc
    class(OutputContainer), intent(inout) :: this
    integer, intent(in) :: iBlock
    integer :: i

    select case (this%gridKind)
    case (GRID_GEO_3D, GRID_GEO_2D)
      this%this_rank_writes = .true.
    case (GRID_MAG_2D)
      this%this_rank_writes = (iProc == 0)
    case (GRID_HIME)
      ! HIME participation is a per-rank flag set externally; do not override here.
      ! Caller must set this_rank_writes before calling put().
    case default
      this%this_rank_writes = .true.
    end select

    ! One-shot fill about to run on this rank: arm the put counter and wipe the
    ! data, so what lands in the file depends only on THIS fill's puts and never
    ! on a previous block or timestep. Unconverted types and non-writing ranks
    ! keep nPuts = -1 and are left alone: USR in particular relies on
    ! last-known values persisting between puts.
    if (this%is_one_shot .and. this%this_rank_writes) then
      this%nPuts = 0
      do i = 1, this%nVars
        if (allocated(this%vars(i)%data)) this%vars(i)%data = 0.0_output_kind
      enddo
    else
      this%nPuts = -1
    endif
  end subroutine prepare

  ! ---------------------------------------------------------------------------
  ! reset: clear participation flag between timesteps.  Keeps allocations.
  ! ---------------------------------------------------------------------------
  subroutine reset(this)
    class(OutputContainer), intent(inout) :: this
    this%this_rank_writes = .false.
    this%any_rank_skips = .false.
  end subroutine reset

  ! ---------------------------------------------------------------------------
  ! find_var: return index of named variable, or 0 if not found.
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
      endif
    enddo
  end function find_var

  ! ---------------------------------------------------------------------------
  ! put_0d: store a scalar value.
  ! ---------------------------------------------------------------------------
  subroutine put_0d(this, name, scalar)
    class(OutputContainer), intent(inout) :: this
    character(len=*), intent(in) :: name
    real(output_kind), intent(in) :: scalar
    integer :: idx

    if (.not. this%this_rank_writes) return
    idx = this%find_var(name)
    if (idx == 0) &
      call stop_gitm('ModOutputContainer: put_0d unknown variable "'// &
                     trim(name)//'" in type '//trim(this%cType))
    this%vars(idx)%data(1, 1, 1) = scalar
  end subroutine put_0d

  ! ---------------------------------------------------------------------------
  ! put_1d: store a 1D array.
  ! ---------------------------------------------------------------------------
  subroutine put_1d(this, name, arr)
    class(OutputContainer), intent(inout) :: this
    character(len=*), intent(in) :: name
    real(output_kind), intent(in) :: arr(:)
    integer :: idx, n
    character(len=128) :: errmsg

    if (.not. this%this_rank_writes) return
    idx = this%find_var(name)
    if (idx == 0) &
      call stop_gitm('ModOutputContainer: put_1d unknown variable "'// &
                     trim(name)//'" in type '//trim(this%cType))
    n = size(arr)
    if (n /= this%vars(idx)%shape3(1)) then
      write(errmsg, '(a,a,a,i0,a,i0)') 'ModOutputContainer: put_1d shape mismatch for "', &
        trim(name), '": expected ', this%vars(idx)%shape3(1), ', got ', n
      call stop_gitm(trim(errmsg))
    endif
    this%vars(idx)%data(:, 1, 1) = arr
  end subroutine put_1d

  ! ---------------------------------------------------------------------------
  ! put_2d: store a 2D array.
  ! ---------------------------------------------------------------------------
  subroutine put_2d(this, name, arr)
    class(OutputContainer), intent(inout) :: this
    character(len=*), intent(in) :: name
    real(output_kind), intent(in) :: arr(:, :)
    integer :: idx, s1, s2
    character(len=128) :: errmsg

    if (.not. this%this_rank_writes) return
    idx = this%find_var(name)
    if (idx == 0) &
      call stop_gitm('ModOutputContainer: put_2d unknown variable "'// &
                     trim(name)//'" in type '//trim(this%cType))
    s1 = size(arr, 1); s2 = size(arr, 2)
    if (s1 /= this%vars(idx)%shape3(1) .or. s2 /= this%vars(idx)%shape3(2)) then
      write(errmsg, '(a,a,a,2(i0,a),a,2(i0,a))') &
        'ModOutputContainer: put_2d shape mismatch for "', trim(name), '": expected ', &
        this%vars(idx)%shape3(1), 'x', this%vars(idx)%shape3(2), ', got ', &
        s1, 'x', s2, ''
      call stop_gitm(trim(errmsg))
    endif
    this%vars(idx)%data(:, :, 1) = arr
  end subroutine put_2d

  ! ---------------------------------------------------------------------------
  ! put_3d: store a 3D array.
  ! ---------------------------------------------------------------------------
  subroutine put_3d(this, name, arr)
    class(OutputContainer), intent(inout) :: this
    character(len=*), intent(in) :: name
    real(output_kind), intent(in) :: arr(:, :, :)
    integer :: idx, s1, s2, s3
    character(len=128) :: errmsg

    if (.not. this%this_rank_writes) return
    idx = this%find_var(name)
    if (idx == 0) &
      call stop_gitm('ModOutputContainer: put_3d unknown variable "'// &
                     trim(name)//'" in type '//trim(this%cType))
    s1 = size(arr, 1); s2 = size(arr, 2); s3 = size(arr, 3)
    if (s1 /= this%vars(idx)%shape3(1) .or. &
        s2 /= this%vars(idx)%shape3(2) .or. &
        s3 /= this%vars(idx)%shape3(3)) then
      write(errmsg, '(a,a,a,3(i0,a),3(i0,a))') &
        'ModOutputContainer: put_3d shape mismatch for "', trim(name), '": expected ', &
        this%vars(idx)%shape3(1), 'x', this%vars(idx)%shape3(2), 'x', &
        this%vars(idx)%shape3(3), ', got ', s1, 'x', s2, 'x', s3, ''
      call stop_gitm(trim(errmsg))
    endif
    this%vars(idx)%data = arr
  end subroutine put_3d

end module ModOutputContainer
