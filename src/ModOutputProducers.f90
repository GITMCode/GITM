! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE
!
! ModOutputProducers.f90
!
! Container-based output producers.  Each output type has two procedures:
!   define_schema_<type>(c)   — register variable names/units/shapes once
!   fill_<type>(c, iBlock)    — put data into the container each timestep
!
! Module-level containers(:) array lives for the whole run; allocated once in
! init_output_containers, reused across timesteps.
!
! Ghost cells (per-type, opt-in):
!   Containers default to interior-only outputs (`c%nGhostCells = 0`).  Any
!   define_schema_<type> may opt in by setting `c%nGhostCells = 2` and padding
!   every variable's spatial dims by 2*nGhostCells; the matching fill_<type>
!   then `put`s ghost-inclusive slices (e.g. -1-nGC:nLons+nGC+2).  Legacy and mpiio
!   backends write the padded container as-is; the netcdf backend strips ghost
!   cells automatically for geographic gridKinds (mag-grid types are never
!   stripped).  Magnetic, 1D, and 0D types have no ghost concept — leave the
!   default 0.

module ModOutputProducers

  use ModOutputContainer

  implicit none

  type(OutputContainer), allocatable, target :: containers(:)

  ! Index constants — one per output type handled here.
  integer, parameter :: iCont_2DMEL = 1
  integer, parameter :: iCont_2DTEC = 2
  integer, parameter :: iCont_2DGEL = 3
  integer, parameter :: iCont_3DNEU = 4
  integer, parameter :: iCont_3DALL = 5
  integer, parameter :: iCont_3DION = 6
  integer, parameter :: iCont_3DTHM = 7
  integer, parameter :: iCont_3DCHM = 8
  integer, parameter :: iCont_3DLST = 9
  integer, parameter :: iCont_3DGLO = 10
  integer, parameter :: iCont_3DMAG = 11
  integer, parameter :: iCont_3DHME = 12
  integer, parameter :: iCont_3DMOH = 13
  integer, parameter :: iCont_3DMOV = 14
  integer, parameter :: iCont_2DANC = 15
  integer, parameter :: iCont_2DHME = 16
  integer, parameter :: iCont_1DALL = 17
  integer, parameter :: iCont_0DALL = 18
  integer, parameter :: iCont_1DGLO = 19
  integer, parameter :: iCont_1DTHM = 20
  integer, parameter :: iCont_1DCHM = 21
  integer, parameter :: iCont_1DNEW = 22
  integer, parameter :: iCont_3DUSR = 23
  integer, parameter :: iCont_2DUSR = 24
  integer, parameter :: iCont_1DUSR = 25
  integer, parameter :: iCont_0DUSR = 26
  integer, parameter :: nContainers = 26

  logical :: containers_initialized = .false.

  ! ---------------------------------------------------------------------------
  ! USR variable registry — populated by register_usr_var() from user.f90 before
  ! the first output timestep.  Coord vars (Longitude/Latitude/Altitude) are
  ! always defined unconditionally in define_schema_*dusr and are NOT listed here.
  ! ---------------------------------------------------------------------------
  integer, parameter :: MaxUsrVars = 100

  type :: UsrVarEntry
    character(len=60) :: name     = ''
    character(len=40) :: units    = ''
    character(len=60) :: longName = ''
  end type UsrVarEntry

  type(UsrVarEntry) :: UsrVars3D(MaxUsrVars)
  integer           :: nUsrVars3D = 0
  type(UsrVarEntry) :: UsrVars2D(MaxUsrVars)
  integer           :: nUsrVars2D = 0
  type(UsrVarEntry) :: UsrVars1D(MaxUsrVars)
  integer           :: nUsrVars1D = 0
  type(UsrVarEntry) :: UsrVars0D(MaxUsrVars)
  integer           :: nUsrVars0D = 0

  ! Generic interface for storing user-defined output data.
  ! Dispatches on array rank: 3D (lon×lat×alt), 2D (lon×lat), 1D (alt profile).
  ! See add_usr_output_3d for full documentation.
  interface add_usr_output
    module procedure add_usr_output_3d
    module procedure add_usr_output_2d
    module procedure add_usr_output_1d
  end interface add_usr_output

contains

  ! ---------------------------------------------------------------------------
  ! is_known_output_type — return .true. if cType is a registered output code.
  ! Pure string check, no container state required (safe to call before init).
  ! ---------------------------------------------------------------------------
  logical function is_known_output_type(cType)
    character(len=5), intent(in) :: cType
    select case (trim(cType))
    case ('2DMEL', '2DTEC', '2DGEL', '3DNEU', '3DALL', '3DION', '3DTHM', '3DCHM', &
          '3DLST', '3DGLO', '3DMAG', '3DHME', '3DMOH', '3DMOV', '2DANC', '2DHME', &
          '1DALL', '0DALL', '1DGLO', '1DTHM', '1DCHM', '1DNEW', &
          '3DUSR', '2DUSR', '1DUSR', '0DUSR')
      is_known_output_type = .true.
    case default
      is_known_output_type = .false.
    end select
  end function is_known_output_type

  ! ---------------------------------------------------------------------------
  ! get_container_idx — return the containers(:) index for cType, or -1.
  ! Requires init_output_containers() to have been called.
  ! ---------------------------------------------------------------------------
  integer function get_container_idx(cType)
    character(len=5), intent(in) :: cType
    select case (trim(cType))
    case ('2DMEL'); get_container_idx = iCont_2DMEL
    case ('2DTEC'); get_container_idx = iCont_2DTEC
    case ('2DGEL'); get_container_idx = iCont_2DGEL
    case ('3DNEU'); get_container_idx = iCont_3DNEU
    case ('3DALL'); get_container_idx = iCont_3DALL
    case ('3DION'); get_container_idx = iCont_3DION
    case ('3DTHM'); get_container_idx = iCont_3DTHM
    case ('3DCHM'); get_container_idx = iCont_3DCHM
    case ('3DLST'); get_container_idx = iCont_3DLST
    case ('3DGLO'); get_container_idx = iCont_3DGLO
    case ('3DMAG'); get_container_idx = iCont_3DMAG
    case ('3DHME'); get_container_idx = iCont_3DHME
    case ('3DMOH'); get_container_idx = iCont_3DMOH
    case ('3DMOV'); get_container_idx = iCont_3DMOV
    case ('2DANC'); get_container_idx = iCont_2DANC
    case ('2DHME'); get_container_idx = iCont_2DHME
    case ('1DALL'); get_container_idx = iCont_1DALL
    case ('0DALL'); get_container_idx = iCont_0DALL
    case ('1DGLO'); get_container_idx = iCont_1DGLO
    case ('1DTHM'); get_container_idx = iCont_1DTHM
    case ('1DCHM'); get_container_idx = iCont_1DCHM
    case ('1DNEW'); get_container_idx = iCont_1DNEW
    case ('3DUSR'); get_container_idx = iCont_3DUSR
    case ('2DUSR'); get_container_idx = iCont_2DUSR
    case ('1DUSR'); get_container_idx = iCont_1DUSR
    case ('0DUSR'); get_container_idx = iCont_0DUSR
    case default;   get_container_idx = -1
    end select
  end function get_container_idx

  ! ---------------------------------------------------------------------------
  ! register_usr_var — add one user-defined variable to the appropriate USR list.
  ! Called from init_usr_output_registry in user.f90 before the first output.
  ! dims: 3, 2, 1, or 0 (spatial dimensionality of the output type).
  ! Coord vars (Longitude, Latitude, Altitude) must NOT be registered here;
  ! they are defined unconditionally in each define_schema_*dusr.
  ! ---------------------------------------------------------------------------
  subroutine register_usr_var(dims, name, units, longName)
    use ModUserGITM, only: nUserOutputs
    integer, intent(in)          :: dims
    character(len=*), intent(in) :: name, units
    character(len=*), intent(in), optional :: longName
    character(len=60) :: lName

    lName = name
    if (present(longName)) lName = longName

    select case (dims)
    case (3)
      nUsrVars3D = nUsrVars3D + 1
      if (nUsrVars3D > nUserOutputs) call stop_gitm("register_usr_var: too many 3DUSR vars — increase nUserOutputs")
      if (nUsrVars3D > MaxUsrVars)   call stop_gitm("register_usr_var: nUsrVars3D exceeds MaxUsrVars")
      UsrVars3D(nUsrVars3D) = UsrVarEntry(name, units, lName)
    case (2)
      nUsrVars2D = nUsrVars2D + 1
      if (nUsrVars2D > nUserOutputs) call stop_gitm("register_usr_var: too many 2DUSR vars — increase nUserOutputs")
      if (nUsrVars2D > MaxUsrVars)   call stop_gitm("register_usr_var: nUsrVars2D exceeds MaxUsrVars")
      UsrVars2D(nUsrVars2D) = UsrVarEntry(name, units, lName)
    case (1)
      nUsrVars1D = nUsrVars1D + 1
      if (nUsrVars1D > nUserOutputs) call stop_gitm("register_usr_var: too many 1DUSR vars — increase nUserOutputs")
      if (nUsrVars1D > MaxUsrVars)   call stop_gitm("register_usr_var: nUsrVars1D exceeds MaxUsrVars")
      UsrVars1D(nUsrVars1D) = UsrVarEntry(name, units, lName)
    case (0)
      nUsrVars0D = nUsrVars0D + 1
      if (nUsrVars0D > nUserOutputs) call stop_gitm("register_usr_var: too many 0DUSR vars — increase nUserOutputs")
      if (nUsrVars0D > MaxUsrVars)   call stop_gitm("register_usr_var: nUsrVars0D exceeds MaxUsrVars")
      UsrVars0D(nUsrVars0D) = UsrVarEntry(name, units, lName)
    end select
  end subroutine register_usr_var

  ! ---------------------------------------------------------------------------
  ! find_usr_var_idx — return 1-based index of name in UsrVars*D, or -1.
  ! ---------------------------------------------------------------------------
  integer function find_usr_var_idx(dim, name) result(idx)
    integer, intent(in)          :: dim
    character(len=*), intent(in) :: name
    integer :: i
    idx = -1
    select case (dim)
    case (3)
      do i = 1, nUsrVars3D
        if (trim(UsrVars3D(i)%name) == trim(name)) then; idx = i; return; end if
      end do
    case (2)
      do i = 1, nUsrVars2D
        if (trim(UsrVars2D(i)%name) == trim(name)) then; idx = i; return; end if
      end do
    case (1)
      do i = 1, nUsrVars1D
        if (trim(UsrVars1D(i)%name) == trim(name)) then; idx = i; return; end if
      end do
    case (0)
      do i = 1, nUsrVars0D
        if (trim(UsrVars0D(i)%name) == trim(name)) then; idx = i; return; end if
      end do
    end select
  end function find_usr_var_idx

  ! ---------------------------------------------------------------------------
  ! add_usr_output_3d — store a 3D (lon×lat×alt) user-defined field.
  !
  ! Call from physics code each timestep to provide data for 3DUSR output.
  ! First call with a new name registers the variable; subsequent calls update
  ! the stored data.  All new-name calls must occur before the first output
  ! write — see init_usr_output_registry in user.f90.
  !
  ! Ghost cell detection (size(arr,1)):
  !   nLons     — interior array, stored as-is
  !   nLons+4   — ghost-inclusive array (2 ghosts each side), stripped to interior
  !   anything else (incl. magnetic grid sizes) — stop_gitm
  !
  ! iBlock defaults to 1 when not supplied.
  ! ---------------------------------------------------------------------------
  subroutine add_usr_output_3d(arr, name, units, iBlock, longName)
    use ModSizeGitm, only: nLons, nLats, nAlts
    use ModUserGITM, only: UserData3D
    use ModElectrodynamics, only: nMagLons, nMagLats
    real, intent(in)                       :: arr(:,:,:)
    character(len=*), intent(in)           :: name, units
    integer, intent(in), optional          :: iBlock
    character(len=*), intent(in), optional :: longName
    integer :: iBlk, idx, n1

    iBlk = 1; if (present(iBlock)) iBlk = iBlock
    n1 = size(arr, 1)

    if (n1 == nMagLons .or. n1 == nMagLons + 1) &
      call stop_gitm("add_usr_output: magnetic grid arrays not supported ('"//trim(name)//"')")
    if (n1 /= nLons .and. n1 /= nLons + 4) &
      call stop_gitm("add_usr_output: unrecognized shape for '"//trim(name)//"' — pass (nLons,nLats,nAlts) or ghost-inclusive")

    idx = find_usr_var_idx(3, name)
    if (idx < 0) then
      call register_usr_var(3, name, units, longName)
      idx = find_usr_var_idx(3, name)
    end if

    if (n1 == nLons + 4) then
      UserData3D(1:nLons, 1:nLats, 1:nAlts, idx, iBlk) = arr(3:nLons+2, 3:nLats+2, 3:nAlts+2)
    else
      UserData3D(1:nLons, 1:nLats, 1:nAlts, idx, iBlk) = arr
    end if
  end subroutine add_usr_output_3d

  ! ---------------------------------------------------------------------------
  ! add_usr_output_2d — store a 2D (lon×lat) user-defined field.
  ! Ghost detection: nLons → interior; nLons+4 → ghost-inclusive, stripped.
  ! ---------------------------------------------------------------------------
  subroutine add_usr_output_2d(arr, name, units, iBlock, longName)
    use ModSizeGitm, only: nLons, nLats
    use ModUserGITM, only: UserData2D
    use ModElectrodynamics, only: nMagLons, nMagLats
    real, intent(in)                       :: arr(:,:)
    character(len=*), intent(in)           :: name, units
    integer, intent(in), optional          :: iBlock
    character(len=*), intent(in), optional :: longName
    integer :: iBlk, idx, n1

    iBlk = 1; if (present(iBlock)) iBlk = iBlock
    n1 = size(arr, 1)

    if (n1 == nMagLons .or. n1 == nMagLons + 1) &
      call stop_gitm("add_usr_output: magnetic grid arrays not supported ('"//trim(name)//"')")
    if (n1 /= nLons .and. n1 /= nLons + 4) &
      call stop_gitm("add_usr_output: unrecognized shape for '"//trim(name)//"' — pass (nLons,nLats) or ghost-inclusive")

    idx = find_usr_var_idx(2, name)
    if (idx < 0) then
      call register_usr_var(2, name, units, longName)
      idx = find_usr_var_idx(2, name)
    end if

    if (n1 == nLons + 4) then
      UserData2D(1:nLons, 1:nLats, 1, idx, iBlk) = arr(3:nLons+2, 3:nLats+2)
    else
      UserData2D(1:nLons, 1:nLats, 1, idx, iBlk) = arr
    end if
  end subroutine add_usr_output_2d

  ! ---------------------------------------------------------------------------
  ! add_usr_output_1d — store a 1D (altitude profile) user-defined field.
  ! Data is stored globally (not per-block); iBlock is accepted but ignored.
  ! arr must have size nAlts.
  ! ---------------------------------------------------------------------------
  subroutine add_usr_output_1d(arr, name, units, iBlock, longName)
    use ModSizeGitm, only: nAlts
    use ModUserGITM, only: UserData1D
    real, intent(in)                       :: arr(:)
    character(len=*), intent(in)           :: name, units
    integer, intent(in), optional          :: iBlock
    character(len=*), intent(in), optional :: longName
    integer :: idx

    if (size(arr) /= nAlts) &
      call stop_gitm("add_usr_output: 1D array must have size nAlts ('"//trim(name)//"')")

    idx = find_usr_var_idx(1, name)
    if (idx < 0) then
      call register_usr_var(1, name, units, longName)
      idx = find_usr_var_idx(1, name)
    end if

    UserData1D(1, 1, 1:nAlts, idx) = arr
  end subroutine add_usr_output_1d

  ! ---------------------------------------------------------------------------
  ! init_output_containers — allocate and define schemas.  Called once.
  ! ---------------------------------------------------------------------------
  subroutine init_output_containers()
    if (containers_initialized) return
    allocate(containers(nContainers))
    call define_schema_2dmel(containers(iCont_2DMEL))
    call define_schema_2dtec(containers(iCont_2DTEC))
    call define_schema_2dgel(containers(iCont_2DGEL))
    call define_schema_3dneu(containers(iCont_3DNEU))
    call define_schema_3dall(containers(iCont_3DALL))
    call define_schema_3dion(containers(iCont_3DION))
    call define_schema_3dthm(containers(iCont_3DTHM))
    call define_schema_3dchm(containers(iCont_3DCHM))
    call define_schema_3dlst(containers(iCont_3DLST))
    call define_schema_3dglo(containers(iCont_3DGLO))
    call define_schema_3dmag(containers(iCont_3DMAG))
    call define_schema_3dhme(containers(iCont_3DHME))
    call define_schema_3dmoh(containers(iCont_3DMOH))
    call define_schema_3dmov(containers(iCont_3DMOV))
    call define_schema_2danc(containers(iCont_2DANC))
    call define_schema_2dhme(containers(iCont_2DHME))
    call define_schema_1dall(containers(iCont_1DALL))
    call define_schema_0dall(containers(iCont_0DALL))
    call define_schema_1dglo(containers(iCont_1DGLO))
    call define_schema_1dthm(containers(iCont_1DTHM))
    call define_schema_1dchm(containers(iCont_1DCHM))
    call define_schema_1dnew(containers(iCont_1DNEW))
    call define_schema_3dusr(containers(iCont_3DUSR))
    call define_schema_2dusr(containers(iCont_2DUSR))
    call define_schema_1dusr(containers(iCont_1DUSR))
    call define_schema_0dusr(containers(iCont_0DUSR))
    containers_initialized = .true.
  end subroutine init_output_containers

  ! ---------------------------------------------------------------------------
  ! define_schema_2dmel — register all 2DMEL variables.
  !
  ! Grid: magnetic (nMagLons+1) x nMagLats, no altitude.
  ! Axis variables are 1D; all others are 2D (degenerate 3rd dim = 1).
  ! Units: degrees for lon/lat (deg, not rad), hours for MLT.
  ! GeoLat/GeoLon are stored in radians in the model; conversion at put time.
  ! ---------------------------------------------------------------------------
  subroutine define_schema_2dmel(c)
    use ModElectrodynamics, only: nMagLats, nMagLons
    type(OutputContainer), intent(inout) :: c

    c%cType    = '2DMEL'
    ! gridKind drives BOTH participation (prepare: only iProc==0 writes) AND
    ! file layout (mpiio_write_container: MAG_2D writes at offset 0, not
    ! striped by global block index). Changing the participation rule for
    ! MAG_2D requires updating the offset logic in lockstep.
    c%gridKind = GRID_MAG_2D

    ! Axis variables (1D coordinate arrays).
    call c%define_var('mlon', units='degrees_east', &
                      shape3=[nMagLons + 1, 1, 1], is_axis=.true., &
                      longName='Magnetic longitude')
    call c%define_var('mlat', units='degrees_north', &
                      shape3=[nMagLats, 1, 1], is_axis=.true., &
                      longName='Magnetic latitude')

    ! 2D data variables (nMagLons+1 x nMagLats x 1).
    call c%define_var('MLT', units='hr', &
                      shape3=[nMagLons + 1, nMagLats, 1], &
                      longName='Magnetic local time')
    call c%define_var('GeoLat', units='degrees_north', &
                      shape3=[nMagLons + 1, nMagLats, 1], &
                      longName='Geographic latitude')
    call c%define_var('GeoLon', units='degrees_east', &
                      shape3=[nMagLons + 1, nMagLats, 1], &
                      longName='Geographic longitude')
    call c%define_var('PedCond', units='S', &
                      shape3=[nMagLons + 1, nMagLats, 1], &
                      longName='Height-integrated Pedersen conductance')
    call c%define_var('HalCond', units='S', &
                      shape3=[nMagLons + 1, nMagLats, 1], &
                      longName='Height-integrated Hall conductance')
    call c%define_var('DivJuAlt', units='A/m2', &
                      shape3=[nMagLons + 1, nMagLats, 1], &
                      longName='Height-integrated divergence of horizontal current')
    call c%define_var('FL_length', units='m', &
                      shape3=[nMagLons + 1, nMagLats, 1], &
                      longName='Field line length')
    call c%define_var('Sigma_PP', units='S', &
                      shape3=[nMagLons + 1, nMagLats, 1], &
                      longName='Sigma PP conductance component')
    call c%define_var('Sigma_LL', units='S', &
                      shape3=[nMagLons + 1, nMagLats, 1], &
                      longName='Sigma LL conductance component')
    call c%define_var('Sigma_H', units='S', &
                      shape3=[nMagLons + 1, nMagLats, 1], &
                      longName='Hall conductance component')
    call c%define_var('Sigma_C', units='S', &
                      shape3=[nMagLons + 1, nMagLats, 1], &
                      longName='Cowling conductance component')
    call c%define_var('Sigma_PL', units='S', &
                      shape3=[nMagLons + 1, nMagLats, 1], &
                      longName='Sigma PL conductance component')
    call c%define_var('Sigma_LP', units='S', &
                      shape3=[nMagLons + 1, nMagLats, 1], &
                      longName='Sigma LP conductance component')
    call c%define_var('K_mphi', units='', &
                      shape3=[nMagLons + 1, nMagLats, 1], &
                      longName='Dynamo K_mphi coefficient')
    call c%define_var('K_mlambda', units='', &
                      shape3=[nMagLons + 1, nMagLats, 1], &
                      longName='Dynamo K_mlambda coefficient')
    call c%define_var('Solver_A', units='', &
                      shape3=[nMagLons + 1, nMagLats, 1], &
                      longName='Dynamo solver coefficient A')
    call c%define_var('Solver_B', units='', &
                      shape3=[nMagLons + 1, nMagLats, 1], &
                      longName='Dynamo solver coefficient B')
    call c%define_var('Solver_C', units='', &
                      shape3=[nMagLons + 1, nMagLats, 1], &
                      longName='Dynamo solver coefficient C')
    call c%define_var('Solver_D', units='', &
                      shape3=[nMagLons + 1, nMagLats, 1], &
                      longName='Dynamo solver coefficient D')
    call c%define_var('Solver_E', units='', &
                      shape3=[nMagLons + 1, nMagLats, 1], &
                      longName='Dynamo solver coefficient E')
    call c%define_var('Solver_S', units='', &
                      shape3=[nMagLons + 1, nMagLats, 1], &
                      longName='Dynamo solver coefficient S')
    call c%define_var('Phi_dynamo', units='V', &
                      shape3=[nMagLons + 1, nMagLats, 1], &
                      longName='Dynamo electric potential')
    call c%define_var('Ed1new', units='V/m', &
                      shape3=[nMagLons + 1, nMagLats, 1], &
                      longName='Dynamo electric field component 1')
    call c%define_var('Ed2new', units='V/m', &
                      shape3=[nMagLons + 1, nMagLats, 1], &
                      longName='Dynamo electric field component 2')
    call c%define_var('Kphi', units='', &
                      shape3=[nMagLons + 1, nMagLats, 1], &
                      longName='Dynamo Kphi coefficient')
    call c%define_var('Klamda', units='', &
                      shape3=[nMagLons + 1, nMagLats, 1], &
                      longName='Dynamo Klambda coefficient')
  end subroutine define_schema_2dmel

  ! ---------------------------------------------------------------------------
  ! fill_2dmel — put 2DMEL data into the container for iBlock.
  !
  ! Only iBlock == 1 writes (GRID_MAG_2D); all others return after prepare().
  ! GeoLat/GeoLon are in radians in the model; multiplied by cRadToDeg here.
  ! MagLonMC/MagLatMC are already in degrees.
  ! ---------------------------------------------------------------------------
  subroutine fill_2dmel(c, iBlock)
    use ModElectrodynamics
    use ModConst, only: cRadToDeg
    type(OutputContainer), intent(inout) :: c
    integer, intent(in) :: iBlock

    call c%prepare(iBlock)
    if (.not. c%this_rank_writes) return

    ! Guard: arrays may not be allocated when dynamo is disabled.
    if (.not. (allocated(MagLonMC) .and. allocated(MagLatMC))) return

    ! Axis arrays: 1D slices from the 2D grid (all rows are identical for axes).
    call c%put('mlon', real(MagLonMC(1:nMagLons + 1, 1), output_kind))
    call c%put('mlat', real(MagLatMC(1, 1:nMagLats), output_kind))

    ! 2D data variables.
    call c%put('MLT',      real(MagLocTimeMC, output_kind))
    call c%put('GeoLat',   real(GeoLatMC * cRadToDeg, output_kind))
    call c%put('GeoLon',   real(GeoLonMC * cRadToDeg, output_kind))
    call c%put('PedCond',  real(SigmaPedersenMC, output_kind))
    call c%put('HalCond',  real(SigmaHallMC, output_kind))
    call c%put('DivJuAlt', real(DivJuAltMC, output_kind))
    call c%put('FL_length',real(LengthMC, output_kind))
    call c%put('Sigma_PP', real(SigmaPPMC, output_kind))
    call c%put('Sigma_LL', real(SigmaLLMC, output_kind))
    call c%put('Sigma_H',  real(SigmaHHMC, output_kind))
    call c%put('Sigma_C',  real(SigmaCCMC, output_kind))
    call c%put('Sigma_PL', real(SigmaPLMC, output_kind))
    call c%put('Sigma_LP', real(SigmaLPMC, output_kind))
    call c%put('K_mphi',    real(KDpmMC, output_kind))
    call c%put('K_mlambda', real(KDlmMC, output_kind))
    call c%put('Solver_A',  real(solver_a_mc, output_kind))
    call c%put('Solver_B',  real(solver_b_mc, output_kind))
    call c%put('Solver_C',  real(solver_c_mc, output_kind))
    call c%put('Solver_D',  real(solver_d_mc, output_kind))
    call c%put('Solver_E',  real(solver_e_mc, output_kind))
    call c%put('Solver_S',  real(solver_s_mc, output_kind))
    call c%put('Phi_dynamo',real(DynamoPotentialMC, output_kind))
    call c%put('Ed1new',    real(Ed1new, output_kind))
    call c%put('Ed2new',    real(Ed2new, output_kind))
    call c%put('Kphi',      real(kpmMC, output_kind))
    call c%put('Klamda',    real(klmMC, output_kind))
  end subroutine fill_2dmel

  ! ---------------------------------------------------------------------------
  ! define_schema_2dtec — register all 2DTEC variables.
  !
  ! Grid: geographic 2D, per-block, nLons x nLats.
  ! Longitude/Latitude stored in radians in model; converted to degrees at put time.
  ! ---------------------------------------------------------------------------
  subroutine define_schema_2dtec(c)
    use ModGITM, only: nLons, nLats
    type(OutputContainer), intent(inout) :: c
    integer :: nGC

    c%cType    = '2DTEC'
    c%gridKind = GRID_GEO_2D
    c%nGhostCells = 0
    nGC = c%nGhostCells

    ! Axis variables (1D coordinate arrays).
    call c%define_var('Longitude', units='degrees_east', &
                      shape3=[nLons + 2*nGC, 1, 1], is_axis=.true., &
                      longName='Geographic longitude')
    call c%define_var('Latitude', units='degrees_north', &
                      shape3=[nLats + 2*nGC, 1, 1], is_axis=.true., &
                      longName='Geographic latitude')

    ! 2D data variables (nLons x nLats x 1).
    call c%define_var('Altitude', units='m', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, 1], &
                      longName='Altitude above surface')
    call c%define_var('SZA', units='rad', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, 1], &
                      longName='Solar zenith angle')
    call c%define_var('TEC', units='TECU', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, 1], &
                      longName='Vertical total electron content')
  end subroutine define_schema_2dtec

  ! ---------------------------------------------------------------------------
  ! fill_2dtec — put 2DTEC data into the container for iBlock.
  !
  ! All blocks participate (GRID_GEO_2D).
  ! Longitude/Latitude are in radians in the model; converted to degrees here.
  ! ---------------------------------------------------------------------------
  subroutine fill_2dtec(c, iBlock)
    use ModGITM
    use ModEUV, only: Sza
    use ModConst, only: cRadToDeg
    type(OutputContainer), intent(inout) :: c
    integer :: nGC
    integer, intent(in) :: iBlock

    call c%prepare(iBlock)
    if (.not. c%this_rank_writes) return
    nGC = c%nGhostCells

    call calc_vtec(iBlock)

    call c%put('Longitude', real(Longitude(1-nGC:nLons+nGC, iBlock) * cRadToDeg, output_kind))
    call c%put('Latitude',  real(Latitude(1-nGC:nLats+nGC, iBlock) * cRadToDeg, output_kind))
    call c%put('Altitude',  real(Altitude_GB(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1, iBlock), output_kind))
    call c%put('SZA',       real(Sza(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, iBlock), output_kind))
    call c%put('TEC',       real(VTEC(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, iBlock), output_kind))
  end subroutine fill_2dtec

  ! ---------------------------------------------------------------------------
  ! define_schema_2dgel — register all 2DGEL variables.
  !
  ! Grid: geographic 2D, per-block, nLons x nLats.
  ! Longitude/Latitude stored in radians in model; converted to degrees at put time.
  ! ---------------------------------------------------------------------------
  subroutine define_schema_2dgel(c)
    use ModGITM, only: nLons, nLats
    type(OutputContainer), intent(inout) :: c
    integer :: nGC

    c%cType    = '2DGEL'
    c%gridKind = GRID_GEO_2D
    c%nGhostCells = 0
    nGC = c%nGhostCells

    ! Axis variables (1D coordinate arrays).
    call c%define_var('Longitude', units='degrees_east', &
                      shape3=[nLons + 2*nGC, 1, 1], is_axis=.true., &
                      longName='Geographic longitude')
    call c%define_var('Latitude', units='degrees_north', &
                      shape3=[nLats + 2*nGC, 1, 1], is_axis=.true., &
                      longName='Geographic latitude')

    ! 2D data variables (nLons x nLats x 1).
    call c%define_var('Altitude', units='m', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, 1], &
                      longName='Altitude above surface')
    call c%define_var('pot', units='V', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, 1], &
                      longName='Electric potential')
    call c%define_var('PedCond', units='S', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, 1], &
                      longName='Height-integrated Pedersen conductance')
    call c%define_var('HalCond', units='S', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, 1], &
                      longName='Height-integrated Hall conductance')
    call c%define_var('AveE', units='eV', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, 1], &
                      longName='Diffuse electron precipitation average energy')
    call c%define_var('eFlux', units='ergs/cm2/s', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, 1], &
                      longName='Diffuse electron precipitation energy flux')
    call c%define_var('AveE_W', units='eV', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, 1], &
                      longName='Wave-driven electron precipitation average energy')
    call c%define_var('eFlux_W', units='ergs/cm2/s', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, 1], &
                      longName='Wave-driven electron precipitation energy flux')
    call c%define_var('AveE_M', units='eV', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, 1], &
                      longName='Monoenergetic electron precipitation average energy')
    call c%define_var('eFlux_M', units='ergs/cm2/s', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, 1], &
                      longName='Monoenergetic electron precipitation energy flux')
    call c%define_var('AveE_I', units='eV', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, 1], &
                      longName='Ion precipitation average energy')
    call c%define_var('eFlux_I', units='ergs/cm2/s', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, 1], &
                      longName='Ion precipitation energy flux')
    call c%define_var('DivJuAlt', units='A/m2', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, 1], &
                      longName='Height-integrated divergence of horizontal current')
    call c%define_var('PedFLCond', units='S', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, 1], &
                      longName='Field-line integrated Pedersen conductance')
    call c%define_var('HalFLCond', units='S', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, 1], &
                      longName='Field-line integrated Hall conductance')
    call c%define_var('DivJuFL', units='A/m2', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, 1], &
                      longName='Field-line divergence of current')
    call c%define_var('FL_length', units='m', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, 1], &
                      longName='Field line length')
  end subroutine define_schema_2dgel

  ! ---------------------------------------------------------------------------
  ! fill_2dgel — put 2DGEL data into the container for iBlock.
  !
  ! All blocks participate (GRID_GEO_2D).
  ! Longitude/Latitude are in radians in the model; converted to degrees here.
  ! ModElectrodynamics 2D arrays (ElectronAverageEnergyDiffuse etc.) are
  ! block-local: no iBlock index, hold the current block's data at call time.
  ! ---------------------------------------------------------------------------
  subroutine fill_2dgel(c, iBlock)
    use ModGITM
    use ModElectrodynamics
    use ModConst, only: cRadToDeg
    type(OutputContainer), intent(inout) :: c
    integer :: nGC
    integer, intent(in) :: iBlock

    call c%prepare(iBlock)
    if (.not. c%this_rank_writes) return
    nGC = c%nGhostCells

    call c%put('Longitude', real(Longitude(1-nGC:nLons+nGC, iBlock) * cRadToDeg, output_kind))
    call c%put('Latitude',  real(Latitude(1-nGC:nLats+nGC, iBlock) * cRadToDeg, output_kind))
    call c%put('Altitude',  real(Altitude_GB(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1, iBlock), output_kind))
    call c%put('pot',       real(Potential(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1, iBlock), output_kind))
    call c%put('PedCond',   real(PedersenConductance(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, iBlock), output_kind))
    call c%put('HalCond',   real(HallConductance(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, iBlock), output_kind))
    call c%put('AveE',      real(ElectronAverageEnergyDiffuse(1-nGC:nLons+nGC, 1-nGC:nLats+nGC), output_kind))
    call c%put('eFlux',     real(ElectronEnergyFluxDiffuse(1-nGC:nLons+nGC, 1-nGC:nLats+nGC), output_kind))
    call c%put('AveE_W',    real(ElectronAverageEnergyWave(1-nGC:nLons+nGC, 1-nGC:nLats+nGC), output_kind))
    call c%put('eFlux_W',   real(ElectronEnergyFluxWave(1-nGC:nLons+nGC, 1-nGC:nLats+nGC), output_kind))
    call c%put('AveE_M',    real(ElectronAverageEnergyMono(1-nGC:nLons+nGC, 1-nGC:nLats+nGC), output_kind))
    call c%put('eFlux_M',   real(ElectronEnergyFluxMono(1-nGC:nLons+nGC, 1-nGC:nLats+nGC), output_kind))
    call c%put('AveE_I',    real(IonAverageEnergy(1-nGC:nLons+nGC, 1-nGC:nLats+nGC), output_kind))
    call c%put('eFlux_I',   real(IonEnergyFlux(1-nGC:nLons+nGC, 1-nGC:nLats+nGC), output_kind))
    call c%put('DivJuAlt',  real(DivJuAlt(1-nGC:nLons+nGC, 1-nGC:nLats+nGC), output_kind))
    call c%put('PedFLCond', real(PedersenFieldLine(1-nGC:nLons+nGC, 1-nGC:nLats+nGC), output_kind))
    call c%put('HalFLCond', real(HallFieldLine(1-nGC:nLons+nGC, 1-nGC:nLats+nGC), output_kind))
    call c%put('DivJuFL',   real(DivJuFieldLine(1-nGC:nLons+nGC, 1-nGC:nLats+nGC), output_kind))
    call c%put('FL_length', real(LengthFieldLine(1-nGC:nLons+nGC, 1-nGC:nLats+nGC), output_kind))
  end subroutine fill_2dgel

  ! ---------------------------------------------------------------------------
  ! define_schema_3dneu — register all 3DNEU variables.
  ! ---------------------------------------------------------------------------
  subroutine define_schema_3dneu(c)
    use ModGITM, only: nLons, nLats, nAlts
    use ModPlanet, only: nSpeciesTotal, nSpecies, cSpecies
    type(OutputContainer), intent(inout) :: c
    integer :: nGC
    integer :: i

    c%cType    = '3DNEU'
    c%gridKind = GRID_GEO_3D
    c%nGhostCells = 2
    nGC = c%nGhostCells

    call c%define_var('Longitude', units='degrees_east', &
                      shape3=[nLons + 2*nGC, 1, 1], is_axis=.true., &
                      longName='Geographic longitude')
    call c%define_var('Latitude', units='degrees_north', &
                      shape3=[nLats + 2*nGC, 1, 1], is_axis=.true., &
                      longName='Geographic latitude')
    call c%define_var('Altitude', units='m', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], &
                      longName='Altitude above surface')
    call c%define_var('Rho', units='kg/m3', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], &
                      longName='Total neutral mass density')
    do i = 1, nSpeciesTotal
      call c%define_var('['//trim(cSpecies(i))//']', units='/m3', &
                        shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], &
                        longName='Number density of '//trim(cSpecies(i)))
    end do
    call c%define_var('Temperature', units='K', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], &
                      longName='Neutral temperature')
    call c%define_var('Vn (east)', units='m/s', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], &
                      longName='Eastward neutral wind')
    call c%define_var('Vn (north)', units='m/s', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], &
                      longName='Northward neutral wind')
    call c%define_var('Vn (up)', units='m/s', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], &
                      longName='Vertical neutral wind')
    do i = 1, nSpecies
      call c%define_var('Vn (up,'//trim(cSpecies(i))//')', units='m/s', &
                        shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], &
                        longName='Vertical neutral wind for '//trim(cSpecies(i)))
    end do
  end subroutine define_schema_3dneu

  ! ---------------------------------------------------------------------------
  ! fill_3dneu — put 3DNEU data into the container for iBlock.
  ! ---------------------------------------------------------------------------
  subroutine fill_3dneu(c, iBlock)
    use ModGITM
    use ModConst, only: cRadToDeg
    use ModPlanet, only: nSpeciesTotal, nSpecies, cSpecies
    type(OutputContainer), intent(inout) :: c
    integer :: nGC
    integer, intent(in) :: iBlock
    integer :: i

    call c%prepare(iBlock)
    if (.not. c%this_rank_writes) return
    nGC = c%nGhostCells

    call c%put('Longitude', real(Longitude(1-nGC:nLons+nGC, iBlock) * cRadToDeg, output_kind))
    call c%put('Latitude',  real(Latitude(1-nGC:nLats+nGC, iBlock) * cRadToDeg, output_kind))
    call c%put('Altitude',  real(Altitude_GB(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, iBlock), output_kind))
    call c%put('Rho',       real(Rho(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, iBlock), output_kind))
    do i = 1, nSpeciesTotal
      call c%put('['//trim(cSpecies(i))//']', &
                 real(NDensityS(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, i, iBlock), output_kind))
    end do
    call c%put('Temperature', &
               real(Temperature(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, iBlock) * &
                    TempUnit(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC), output_kind))
    call c%put('Vn (east)',  real(Velocity(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, 1, iBlock), output_kind))
    call c%put('Vn (north)', real(Velocity(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, 2, iBlock), output_kind))
    call c%put('Vn (up)',    real(Velocity(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, 3, iBlock), output_kind))
    do i = 1, nSpecies
      call c%put('Vn (up,'//trim(cSpecies(i))//')', &
                 real(VerticalVelocity(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, i, iBlock), output_kind))
    end do
  end subroutine fill_3dneu

  ! ---------------------------------------------------------------------------
  ! define_schema_3dall — register all 3DALL variables.
  ! ---------------------------------------------------------------------------
  subroutine define_schema_3dall(c)
    use ModGITM, only: nLons, nLats, nAlts
    use ModPlanet, only: nSpeciesTotal, nSpecies, nIons, cSpecies, cIons
    type(OutputContainer), intent(inout) :: c
    integer :: nGC
    integer :: i

    c%cType    = '3DALL'
    c%gridKind = GRID_GEO_3D
    c%nGhostCells = 2
    nGC = c%nGhostCells

    call c%define_var('Longitude', units='degrees_east', &
                      shape3=[nLons + 2*nGC, 1, 1], is_axis=.true., &
                      longName='Geographic longitude')
    call c%define_var('Latitude', units='degrees_north', &
                      shape3=[nLats + 2*nGC, 1, 1], is_axis=.true., &
                      longName='Geographic latitude')
    call c%define_var('Altitude', units='m', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], &
                      longName='Altitude above surface')
    call c%define_var('Rho', units='kg/m3', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], &
                      longName='Total neutral mass density')
    do i = 1, nSpeciesTotal
      call c%define_var('['//trim(cSpecies(i))//']', units='/m3', &
                        shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], &
                        longName='Number density of '//trim(cSpecies(i)))
    end do
    call c%define_var('Temperature', units='K', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], &
                      longName='Neutral temperature')
    call c%define_var('Vn (east)', units='m/s', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], &
                      longName='Eastward neutral wind')
    call c%define_var('Vn (north)', units='m/s', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], &
                      longName='Northward neutral wind')
    call c%define_var('Vn (up)', units='m/s', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], &
                      longName='Vertical neutral wind')
    do i = 1, nSpecies
      call c%define_var('Vn (up,'//trim(cSpecies(i))//')', units='m/s', &
                        shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], &
                        longName='Vertical neutral wind for '//trim(cSpecies(i)))
    end do
    do i = 1, nIons
      call c%define_var('['//trim(cIons(i))//']', units='/m3', &
                        shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], &
                        longName='Number density of '//trim(cIons(i)))
    end do
    call c%define_var('eTemperature', units='K', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], &
                      longName='Electron temperature')
    call c%define_var('iTemperature', units='K', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], &
                      longName='Ion temperature')
    call c%define_var('Vi (east)', units='m/s', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], &
                      longName='Eastward ion drift')
    call c%define_var('Vi (north)', units='m/s', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], &
                      longName='Northward ion drift')
    call c%define_var('Vi (up)', units='m/s', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], &
                      longName='Vertical ion drift')
  end subroutine define_schema_3dall

  ! ---------------------------------------------------------------------------
  ! fill_3dall — put 3DALL data into the container for iBlock.
  ! ---------------------------------------------------------------------------
  subroutine fill_3dall(c, iBlock)
    use ModGITM
    use ModConst, only: cRadToDeg
    use ModPlanet, only: nSpeciesTotal, nSpecies, nIons, cSpecies, cIons
    type(OutputContainer), intent(inout) :: c
    integer :: nGC
    integer, intent(in) :: iBlock
    integer :: i

    call c%prepare(iBlock)
    if (.not. c%this_rank_writes) return
    nGC = c%nGhostCells

    call c%put('Longitude', real(Longitude(1-nGC:nLons+nGC, iBlock) * cRadToDeg, output_kind))
    call c%put('Latitude',  real(Latitude(1-nGC:nLats+nGC, iBlock) * cRadToDeg, output_kind))
    call c%put('Altitude',  real(Altitude_GB(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, iBlock), output_kind))
    call c%put('Rho',       real(Rho(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, iBlock), output_kind))
    do i = 1, nSpeciesTotal
      call c%put('['//trim(cSpecies(i))//']', &
                 real(NDensityS(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, i, iBlock), output_kind))
    end do
    call c%put('Temperature', &
               real(Temperature(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, iBlock) * &
                    TempUnit(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC), output_kind))
    call c%put('Vn (east)',  real(Velocity(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, 1, iBlock), output_kind))
    call c%put('Vn (north)', real(Velocity(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, 2, iBlock), output_kind))
    call c%put('Vn (up)',    real(Velocity(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, 3, iBlock), output_kind))
    do i = 1, nSpecies
      call c%put('Vn (up,'//trim(cSpecies(i))//')', &
                 real(VerticalVelocity(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, i, iBlock), output_kind))
    end do
    do i = 1, nIons
      call c%put('['//trim(cIons(i))//']', &
                 real(IDensityS(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, i, iBlock), output_kind))
    end do
    call c%put('eTemperature', real(eTemperature(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, iBlock), output_kind))
    call c%put('iTemperature', real(ITemperature(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, iBlock), output_kind))
    call c%put('Vi (east)',    real(IVelocity(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, 1, iBlock), output_kind))
    call c%put('Vi (north)',   real(IVelocity(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, 2, iBlock), output_kind))
    call c%put('Vi (up)',      real(IVelocity(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, 3, iBlock), output_kind))
  end subroutine fill_3dall

  ! ---------------------------------------------------------------------------
  ! define_schema_3dion — register all 3DION variables.
  ! ---------------------------------------------------------------------------
  subroutine define_schema_3dion(c)
    use ModGITM, only: nLons, nLats, nAlts
    use ModPlanet, only: nIons, cIons
    type(OutputContainer), intent(inout) :: c
    integer :: nGC
    integer :: i

    c%cType    = '3DION'
    c%gridKind = GRID_GEO_3D
    c%nGhostCells = 2
    nGC = c%nGhostCells

    call c%define_var('Longitude', units='degrees_east', &
                      shape3=[nLons + 2*nGC, 1, 1], is_axis=.true., &
                      longName='Geographic longitude')
    call c%define_var('Latitude', units='degrees_north', &
                      shape3=[nLats + 2*nGC, 1, 1], is_axis=.true., &
                      longName='Geographic latitude')
    call c%define_var('Altitude', units='m', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], &
                      longName='Altitude above surface')
    do i = 1, nIons
      call c%define_var('['//trim(cIons(i))//']', units='/m3', &
                        shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], &
                        longName='Number density of '//trim(cIons(i)))
    end do
    call c%define_var('eTemperature', units='K', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], &
                      longName='Electron temperature')
    call c%define_var('iTemperature', units='K', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], &
                      longName='Ion temperature')
    call c%define_var('Vi (east)',  units='m/s', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], &
                      longName='Eastward ion drift')
    call c%define_var('Vi (north)', units='m/s', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], &
                      longName='Northward ion drift')
    call c%define_var('Vi (up)',    units='m/s', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], &
                      longName='Vertical ion drift')
    call c%define_var('Ed1', units='V/m', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], &
                      longName='Electric field component 1 (magnetic east)')
    call c%define_var('Ed2', units='V/m', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], &
                      longName='Electric field component 2 (magnetic north)')
    call c%define_var('Je1', units='A/m2', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], &
                      longName='Current density component 1 (magnetic east)')
    call c%define_var('Je2', units='A/m2', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], &
                      longName='Current density component 2 (magnetic north)')
    call c%define_var('Magnetic Latitude',  units='deg', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], &
                      longName='Magnetic latitude')
    call c%define_var('Magnetic Longitude', units='deg', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], &
                      longName='Magnetic longitude')
    call c%define_var('B.F. East',     units='T', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], &
                      longName='Eastward magnetic field component')
    call c%define_var('B.F. North',    units='T', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], &
                      longName='Northward magnetic field component')
    call c%define_var('B.F. Vertical', units='T', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], &
                      longName='Vertical magnetic field component')
    call c%define_var('B.F. Magnitude', units='T', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], &
                      longName='Magnetic field magnitude')
    call c%define_var('Potential', units='V', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], &
                      longName='Electric potential')
    call c%define_var('E.F. East',      units='V/m', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], &
                      longName='Eastward electric field component')
    call c%define_var('E.F. North',     units='V/m', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], &
                      longName='Northward electric field component')
    call c%define_var('E.F. Vertical',  units='V/m', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], &
                      longName='Vertical electric field component')
    call c%define_var('E.F. Magnitude', units='V/m', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], &
                      longName='Electric field magnitude')
    call c%define_var('IN Collision Freq', units='Hz', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], &
                      longName='Ion-neutral collision frequency')
    call c%define_var('PressGrad (east)',  units='m/s2', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], &
                      longName='Pressure gradient acceleration eastward')
    call c%define_var('PressGrad (north)', units='m/s2', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], &
                      longName='Pressure gradient acceleration northward')
    call c%define_var('PressGrad (up)',   units='m/s2', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], &
                      longName='Pressure gradient acceleration upward')
  end subroutine define_schema_3dion

  ! ---------------------------------------------------------------------------
  ! fill_3dion — put 3DION data into the container for iBlock.
  ! ---------------------------------------------------------------------------
  subroutine fill_3dion(c, iBlock)
    use ModGITM
    use ModConst, only: cRadToDeg
    use ModPlanet, only: nIons, cIons
    use ModElectrodynamics
    type(OutputContainer), intent(inout) :: c
    integer :: nGC
    integer, intent(in) :: iBlock
    integer :: i

    call c%prepare(iBlock)
    if (.not. c%this_rank_writes) return
    nGC = c%nGhostCells

    call c%put('Longitude', real(Longitude(1-nGC:nLons+nGC, iBlock) * cRadToDeg, output_kind))
    call c%put('Latitude',  real(Latitude(1-nGC:nLats+nGC, iBlock) * cRadToDeg, output_kind))
    call c%put('Altitude',  real(Altitude_GB(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, iBlock), output_kind))
    do i = 1, nIons
      call c%put('['//trim(cIons(i))//']', &
                 real(IDensityS(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, i, iBlock), output_kind))
    end do
    call c%put('eTemperature', real(eTemperature(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, iBlock), output_kind))
    call c%put('iTemperature', real(ITemperature(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, iBlock), output_kind))
    call c%put('Vi (east)',    real(IVelocity(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, 1, iBlock), output_kind))
    call c%put('Vi (north)',   real(IVelocity(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, 2, iBlock), output_kind))
    call c%put('Vi (up)',      real(IVelocity(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, 3, iBlock), output_kind))
    call c%put('Ed1', real(ed1(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC), output_kind))
    call c%put('Ed2', real(ed2(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC), output_kind))
    call c%put('Je1', real(je1(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC), output_kind))
    call c%put('Je2', real(je2(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC), output_kind))
    call c%put('Magnetic Latitude',  real(MLatitude(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, iBlock), output_kind))
    call c%put('Magnetic Longitude', real(MLongitude(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, iBlock), output_kind))
    call c%put('B.F. East',     real(B0(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, 1, iBlock), output_kind))
    call c%put('B.F. North',    real(B0(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, 2, iBlock), output_kind))
    call c%put('B.F. Vertical', real(B0(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, 3, iBlock), output_kind))
    call c%put('B.F. Magnitude',real(B0(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, 4, iBlock), output_kind))
    call c%put('Potential', real(Potential(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, iBlock), output_kind))
    call c%put('E.F. East',      real(EField(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, 1), output_kind))
    call c%put('E.F. North',     real(EField(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, 2), output_kind))
    call c%put('E.F. Vertical',  real(EField(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, 3), output_kind))
    call c%put('E.F. Magnitude', &
               real(sqrt(sum(EField(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, :)**2, dim=4)), output_kind))
    call c%put('IN Collision Freq', real(Collisions(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, iVIN_), output_kind))
    call c%put('PressGrad (east)',  real(IonPressureGradient(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, 1, iBlock), output_kind))
    call c%put('PressGrad (north)', real(IonPressureGradient(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, 2, iBlock), output_kind))
    call c%put('PressGrad (up)',    real(IonPressureGradient(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, 3, iBlock), output_kind))
  end subroutine fill_3dion

  ! ---------------------------------------------------------------------------
  ! define_schema_3dthm — register all 3DTHM variables.
  ! ---------------------------------------------------------------------------
  subroutine define_schema_3dthm(c)
    use ModGITM, only: nLons, nLats, nAlts
    type(OutputContainer), intent(inout) :: c
    integer :: nGC

    c%cType    = '3DTHM'
    c%gridKind = GRID_GEO_3D
    c%nGhostCells = 0
    nGC = c%nGhostCells

    call c%define_var('Longitude', units='degrees_east', &
                      shape3=[nLons + 2*nGC, 1, 1], is_axis=.true., &
                      longName='Geographic longitude')
    call c%define_var('Latitude', units='degrees_north', &
                      shape3=[nLats + 2*nGC, 1, 1], is_axis=.true., &
                      longName='Geographic latitude')
    call c%define_var('Altitude', units='m', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], &
                      longName='Altitude above surface')
    call c%define_var('EUV Heating (K/s)',              units='K/s', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='EUV heating rate')
    call c%define_var('Conduction (K/s)',               units='K/s', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Thermal conduction heating rate')
    call c%define_var('Molecular Conduction (K/s)',     units='K/s', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Molecular conduction heating rate')
    call c%define_var('Eddy Conduction (K/s)',          units='K/s', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Eddy conduction heating rate')
    call c%define_var('Eddy Adiabatic Conduction (K/s)', units='K/s', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Eddy adiabatic conduction heating rate')
    call c%define_var('Chemical Heating (K/s)',         units='K/s', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Chemical heating rate')
    call c%define_var('Joule Heating (K/s)',            units='K/s', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Joule heating rate')
    call c%define_var('NO Cooling (K/s)',               units='K/s', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='NO infrared cooling rate')
    call c%define_var('O Cooling (K/s)',                units='K/s', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='O infrared cooling rate')
    call c%define_var('Total Abs EUV',                 units='W/m2', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Total absorbed EUV flux')
    call c%define_var('Cp',                            units='J/kg/K', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Specific heat at constant pressure')
    call c%define_var('Rho',                           units='kg/m3', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Total neutral mass density')
    call c%define_var('E-Field Mag',                   units='V/m', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Electric field magnitude')
    call c%define_var('Sigma Ped',                     units='S/m', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Pedersen conductivity')
    call c%define_var('Ionization Rate O_3P',          units='/s',  shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Ionization rate for O(3P)')
    call c%define_var('Ionization Rate O2',            units='/s',  shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Ionization rate for O2')
    call c%define_var('Ionization Rate N2',            units='/s',  shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Ionization rate for N2')
  end subroutine define_schema_3dthm

  ! ---------------------------------------------------------------------------
  ! fill_3dthm — put 3DTHM data into the container for iBlock.
  ! ---------------------------------------------------------------------------
  subroutine fill_3dthm(c, iBlock)
    use ModGITM
    use ModConst, only: cRadToDeg
    use ModSources
    use ModElectrodynamics, only: Sigma_Pedersen
    use ModEUV, only: EuvTotal
    type(OutputContainer), intent(inout) :: c
    integer :: nGC
    integer, intent(in) :: iBlock

    call c%prepare(iBlock)
    if (.not. c%this_rank_writes) return
    nGC = c%nGhostCells

    call c%put('Longitude', real(Longitude(1-nGC:nLons+nGC, iBlock) * cRadToDeg, output_kind))
    call c%put('Latitude',  real(Latitude(1-nGC:nLats+nGC, iBlock) * cRadToDeg, output_kind))
    call c%put('Altitude',  real(Altitude_GB(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, iBlock), output_kind))
    call c%put('EUV Heating (K/s)', &
               real(EuvHeating(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, iBlock) * &
                    TempUnit(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC), output_kind))
    call c%put('Conduction (K/s)', &
               real(Conduction(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC) * &
                    TempUnit(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC) / dt, output_kind))
    call c%put('Molecular Conduction (K/s)',      real(MoleConduction(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC), output_kind))
    call c%put('Eddy Conduction (K/s)',           real(EddyCond(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC), output_kind))
    call c%put('Eddy Adiabatic Conduction (K/s)', real(EddyCondAdia(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC), output_kind))
    call c%put('Chemical Heating (K/s)', &
               real(ChemicalHeatingRate(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC) * &
                    TempUnit(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC) / dt, output_kind))
    call c%put('Joule Heating (K/s)', &
               real(JouleHeating(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC) * &
                    TempUnit(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC), output_kind))
    call c%put('NO Cooling (K/s)', &
               real(-NOCooling(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC) * &
                    TempUnit(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC), output_kind))
    call c%put('O Cooling (K/s)', &
               real(-OCooling(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC) * &
                    TempUnit(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC), output_kind))
    call c%put('Total Abs EUV', real(EuvTotal(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, iBlock), output_kind))
    call c%put('Cp',            real(cp(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, iBlock), output_kind))
    call c%put('Rho',           real(Rho(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, iBlock), output_kind))
    call c%put('E-Field Mag',   &
               real(sqrt(sum(EField(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, :)**2, dim=4)), output_kind))
    call c%put('Sigma Ped',     real(Sigma_Pedersen(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC), output_kind))
    call c%put('Ionization Rate O_3P', real(AuroralIonRateS(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, iO_3P_, iBlock), output_kind))
    call c%put('Ionization Rate O2',   real(AuroralIonRateS(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, iO2_,   iBlock), output_kind))
    call c%put('Ionization Rate N2',   real(AuroralIonRateS(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, iN2_,   iBlock), output_kind))
  end subroutine fill_3dthm

  ! ---------------------------------------------------------------------------
  ! define_schema_3dchm — register all 3DCHM variables.
  ! ---------------------------------------------------------------------------
  subroutine define_schema_3dchm(c)
    use ModGITM, only: nLons, nLats, nAlts
    type(OutputContainer), intent(inout) :: c
    integer :: nGC

    c%cType    = '3DCHM'
    c%gridKind = GRID_GEO_3D
    c%nGhostCells = 0
    nGC = c%nGhostCells

    call c%define_var('Longitude', units='degrees_east', &
                      shape3=[nLons + 2*nGC, 1, 1], is_axis=.true., &
                      longName='Geographic longitude')
    call c%define_var('Latitude', units='degrees_north', &
                      shape3=[nLats + 2*nGC, 1, 1], is_axis=.true., &
                      longName='Geographic latitude')
    call c%define_var('Altitude', units='m', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], &
                      longName='Altitude above surface')
    call c%define_var('N2+ + e',    units='/m3/s', shape3=[nLons,nLats,nAlts], longName='N2+ + e recombination rate')
    call c%define_var('O2+ + e',    units='/m3/s', shape3=[nLons,nLats,nAlts], longName='O2+ + e recombination rate')
    call c%define_var('N2+ + O',    units='/m3/s', shape3=[nLons,nLats,nAlts], longName='N2+ + O reaction rate')
    call c%define_var('NO+ + e',    units='/m3/s', shape3=[nLons,nLats,nAlts], longName='NO+ + e recombination rate')
    call c%define_var('N+ + O2',    units='/m3/s', shape3=[nLons,nLats,nAlts], longName='N+ + O2 reaction rate')
    call c%define_var('NO + N',     units='/m3/s', shape3=[nLons,nLats,nAlts], longName='NO + N reaction rate')
    call c%define_var('O+ + O2',    units='/m3/s', shape3=[nLons,nLats,nAlts], longName='O+ + O2 reaction rate')
    call c%define_var('N + O2',     units='/m3/s', shape3=[nLons,nLats,nAlts], longName='N + O2 reaction rate')
    call c%define_var('O2+ + N',    units='/m3/s', shape3=[nLons,nLats,nAlts], longName='O2+ + N reaction rate')
    call c%define_var('O2+ + NO',   units='/m3/s', shape3=[nLons,nLats,nAlts], longName='O2+ + NO reaction rate')
    call c%define_var('O2+ + N2',   units='/m3/s', shape3=[nLons,nLats,nAlts], longName='O2+ + N2 reaction rate')
    call c%define_var('N2+ + O2',   units='/m3/s', shape3=[nLons,nLats,nAlts], longName='N2+ + O2 reaction rate')
    call c%define_var('N+ + O',     units='/m3/s', shape3=[nLons,nLats,nAlts], longName='N+ + O reaction rate')
    call c%define_var('O+ + N2',    units='/m3/s', shape3=[nLons,nLats,nAlts], longName='O+ + N2 reaction rate')
    call c%define_var('O(1D) + N2', units='/m3/s', shape3=[nLons,nLats,nAlts], longName='O(1D) + N2 quench rate')
    call c%define_var('O(1D) + O2', units='/m3/s', shape3=[nLons,nLats,nAlts], longName='O(1D) + O2 quench rate')
    call c%define_var('O(1D) + O',  units='/m3/s', shape3=[nLons,nLats,nAlts], longName='O(1D) + O quench rate')
    call c%define_var('O(1D) + e',  units='/m3/s', shape3=[nLons,nLats,nAlts], longName='O(1D) + e quench rate')
    call c%define_var('N(2D) + O2', units='/m3/s', shape3=[nLons,nLats,nAlts], longName='N(2D) + O2 reaction rate')
    call c%define_var('O+(2D)+e',   units='/m3/s', shape3=[nLons,nLats,nAlts], longName='O+(2D) + e recombination rate')
    call c%define_var('N(2D) + O',  units='/m3/s', shape3=[nLons,nLats,nAlts], longName='N(2D) + O quench rate')
    call c%define_var('N(2D) + e',  units='/m3/s', shape3=[nLons,nLats,nAlts], longName='N(2D) + e quench rate')
    call c%define_var('O+(2D + N2', units='/m3/s', shape3=[nLons,nLats,nAlts], longName='O+(2D) + N2 reaction rate')
    call c%define_var('O+(2P) + e', units='/m3/s', shape3=[nLons,nLats,nAlts], longName='O+(2P) + e recombination rate')
    call c%define_var('O+(2P) + O', units='/m3/s', shape3=[nLons,nLats,nAlts], longName='O+(2P) + O reaction rate')
    call c%define_var('O+(2P) + N2', units='/m3/s', shape3=[nLons,nLats,nAlts], longName='O+(2P) + N2 reaction rate')
    call c%define_var('Chemical Heating Rate', units='K/s', shape3=[nLons,nLats,nAlts], longName='Chemical heating rate')
  end subroutine define_schema_3dchm

  ! ---------------------------------------------------------------------------
  ! fill_3dchm — put 3DCHM data into the container for iBlock.
  ! ---------------------------------------------------------------------------
  subroutine fill_3dchm(c, iBlock)
    use ModGITM
    use ModConst, only: cRadToDeg
    use ModSources
    use ModConstants, only: Element_Charge
    type(OutputContainer), intent(inout) :: c
    integer :: nGC
    integer, intent(in) :: iBlock
    integer :: i
    character(len=20), parameter :: chem_names(nReactions) = [ &
      'N2+ + e             ', 'O2+ + e             ', 'N2+ + O             ', &
      'NO+ + e             ', 'N+ + O2             ', 'NO + N              ', &
      'O+ + O2             ', 'N + O2              ', 'O2+ + N             ', &
      'O2+ + NO            ', 'O2+ + N2            ', 'N2+ + O2            ', &
      'N+ + O              ', 'O+ + N2             ', 'O(1D) + N2          ', &
      'O(1D) + O2          ', 'O(1D) + O           ', 'O(1D) + e           ', &
      'N(2D) + O2          ', 'O+(2D)+e            ', 'N(2D) + O           ', &
      'N(2D) + e           ', 'O+(2D + N2          ', 'O+(2P) + e          ', &
      'O+(2P) + O          ', 'O+(2P) + N2         ']

    call c%prepare(iBlock)
    if (.not. c%this_rank_writes) return
    nGC = c%nGhostCells

    call c%put('Longitude', real(Longitude(1-nGC:nLons+nGC, iBlock) * cRadToDeg, output_kind))
    call c%put('Latitude',  real(Latitude(1-nGC:nLats+nGC, iBlock) * cRadToDeg, output_kind))
    call c%put('Altitude',  real(Altitude_GB(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, iBlock), output_kind))
    do i = 1, nReactions
      call c%put(trim(chem_names(i)), &
                 real(ChemicalHeatingSpecies(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, i) / Element_Charge, output_kind))
    end do
    call c%put('Chemical Heating Rate', &
               real(ChemicalHeatingRate(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC) * &
                    cp(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, iBlock) * &
                    Rho(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, iBlock) * &
                    TempUnit(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC) / Element_Charge, output_kind))
  end subroutine fill_3dchm

  ! ---------------------------------------------------------------------------
  ! define_schema_3dlst — register variables for 3DLST dynamically.
  ! ---------------------------------------------------------------------------
  subroutine define_schema_3dlst(c)
    use ModGITM, only: nLons, nLats, nAlts
    use ModPlanet, only: nSpeciesTotal, nSpecies, nIons, cSpecies, cIons
    use ModInputs, only: iRhoOutputList, iNeutralDensityOutputList, &
                         iNeutralWindOutputList, iIonDensityOutputList, &
                         iIonWindOutputList, iTemperatureOutputList
    type(OutputContainer), intent(inout) :: c
    integer :: nGC
    integer :: i

    c%cType    = '3DLST'
    c%gridKind = GRID_GEO_3D
    c%nGhostCells = 2
    nGC = c%nGhostCells

    call c%define_var('Longitude', units='degrees_east', &
                      shape3=[nLons + 2*nGC, 1, 1], is_axis=.true., &
                      longName='Geographic longitude')
    call c%define_var('Latitude', units='degrees_north', &
                      shape3=[nLats + 2*nGC, 1, 1], is_axis=.true., &
                      longName='Geographic latitude')
    call c%define_var('Altitude', units='m', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], &
                      longName='Altitude')

    if (iRhoOutputList) then
      call c%define_var('Rho', units='kg/m3', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Total neutral mass density')
    endif
    do i = 1, nSpeciesTotal
      if (iNeutralDensityOutputList(i)) then
        call c%define_var('['//trim(cSpecies(i))//']', units='/m3', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], &
                          longName='Number density of '//trim(cSpecies(i)))
      endif
    enddo
    if (iNeutralWindOutputList(1)) then
      call c%define_var('Vn (east)', units='m/s', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Eastward neutral wind')
    endif
    if (iNeutralWindOutputList(2)) then
      call c%define_var('Vn (north)', units='m/s', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Northward neutral wind')
    endif
    if (iNeutralWindOutputList(3)) then
      call c%define_var('Vn (up)', units='m/s', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Vertical neutral wind')
    endif
    do i = 1, nIons
      if (iIonDensityOutputList(i)) then
        call c%define_var('['//trim(cIons(i))//']', units='/m3', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], &
                          longName='Number density of '//trim(cIons(i)))
      endif
    enddo
    if (iIonWindOutputList(1)) then
      call c%define_var('Vi (east)', units='m/s', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Eastward ion drift')
    endif
    if (iIonWindOutputList(2)) then
      call c%define_var('Vi (north)', units='m/s', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Northward ion drift')
    endif
    if (iIonWindOutputList(3)) then
      call c%define_var('Vi (up)', units='m/s', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Vertical ion drift')
    endif
    if (iTemperatureOutputList(1)) then
      call c%define_var('Neutral Temperature', units='K', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Neutral temperature')
    endif
    if (iTemperatureOutputList(2)) then
      call c%define_var('Ion Temperature', units='K', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Ion temperature')
    endif
    if (iTemperatureOutputList(3)) then
      call c%define_var('Electron Temperature', units='K', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Electron temperature')
    endif
  end subroutine define_schema_3dlst

  ! ---------------------------------------------------------------------------
  ! fill_3dlst — put 3DLST data into the container for iBlock.
  ! ---------------------------------------------------------------------------
  subroutine fill_3dlst(c, iBlock)
    use ModGITM
    use ModConst, only: cRadToDeg
    use ModPlanet, only: nSpeciesTotal, nSpecies, nIons, cSpecies, cIons
    use ModInputs, only: iRhoOutputList, iNeutralDensityOutputList, &
                         iNeutralWindOutputList, iIonDensityOutputList, &
                         iIonWindOutputList, iTemperatureOutputList
    type(OutputContainer), intent(inout) :: c
    integer :: nGC
    integer, intent(in) :: iBlock
    integer :: i

    call c%prepare(iBlock)
    if (.not. c%this_rank_writes) return
    nGC = c%nGhostCells

    call c%put('Longitude', real(Longitude(1-nGC:nLons+nGC, iBlock) * cRadToDeg, output_kind))
    call c%put('Latitude',  real(Latitude(1-nGC:nLats+nGC, iBlock) * cRadToDeg, output_kind))
    call c%put('Altitude',  real(Altitude_GB(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, iBlock), output_kind))

    if (iRhoOutputList) then
      call c%put('Rho', real(Rho(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, iBlock), output_kind))
    endif
    do i = 1, nSpeciesTotal
      if (iNeutralDensityOutputList(i)) then
        call c%put('['//trim(cSpecies(i))//']', &
                   real(NDensityS(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, i, iBlock), output_kind))
      endif
    enddo
    if (iNeutralWindOutputList(1)) then
      call c%put('Vn (east)', real(Velocity(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, 1, iBlock), output_kind))
    endif
    if (iNeutralWindOutputList(2)) then
      call c%put('Vn (north)', real(Velocity(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, 2, iBlock), output_kind))
    endif
    if (iNeutralWindOutputList(3)) then
      call c%put('Vn (up)', real(Velocity(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, 3, iBlock), output_kind))
    endif
    do i = 1, nIons
      if (iIonDensityOutputList(i)) then
        call c%put('['//trim(cIons(i))//']', &
                   real(IDensityS(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, i, iBlock), output_kind))
      endif
    enddo
    if (iIonWindOutputList(1)) then
      call c%put('Vi (east)', real(IVelocity(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, 1, iBlock), output_kind))
    endif
    if (iIonWindOutputList(2)) then
      call c%put('Vi (north)', real(IVelocity(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, 2, iBlock), output_kind))
    endif
    if (iIonWindOutputList(3)) then
      call c%put('Vi (up)', real(IVelocity(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, 3, iBlock), output_kind))
    endif
    if (iTemperatureOutputList(1)) then
      call c%put('Neutral Temperature', &
                 real(Temperature(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, iBlock) * &
                      TempUnit(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC), output_kind))
    endif
    if (iTemperatureOutputList(2)) then
      call c%put('Ion Temperature', real(ITemperature(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, iBlock), output_kind))
    endif
    if (iTemperatureOutputList(3)) then
      call c%put('Electron Temperature', real(eTemperature(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, iBlock), output_kind))
    endif
  end subroutine fill_3dlst

  ! ---------------------------------------------------------------------------
  ! define_schema_3dglo — register variables for 3DGLO (glow emission stub).
  ! ---------------------------------------------------------------------------
  subroutine define_schema_3dglo(c)
    use ModGITM, only: nLons, nLats, nAlts
    type(OutputContainer), intent(inout) :: c
    integer :: nGC

    c%cType    = '3DGLO'
    c%gridKind = GRID_GEO_3D
    c%nGhostCells = 2
    nGC = c%nGhostCells

    call c%define_var('Longitude', units='degrees_east', shape3=[nLons + 2*nGC, 1, 1], is_axis=.true., &
                      longName='Geographic longitude')
    call c%define_var('Latitude', units='degrees_north', shape3=[nLats + 2*nGC, 1, 1], is_axis=.true., &
                      longName='Geographic latitude')
    call c%define_var('Altitude', units='m', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Altitude')
    call c%define_var('6300 A Emission', units='', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='6300 Angstrom emission')
    call c%define_var('PhotoElectronUp', units='', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Photoelectron flux up')
    call c%define_var('PhotoElectronDown', units='', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Photoelectron flux down')
  end subroutine define_schema_3dglo

  ! ---------------------------------------------------------------------------
  ! fill_3dglo — put 3DGLO stub data into the container.
  ! ---------------------------------------------------------------------------
  subroutine fill_3dglo(c, iBlock)
    use ModGITM
    use ModConst, only: cRadToDeg
    type(OutputContainer), intent(inout) :: c
    integer :: nGC
    integer, intent(in) :: iBlock
    real(output_kind), allocatable :: zero_buf(:,:,:)

    call c%prepare(iBlock)
    if (.not. c%this_rank_writes) return
    nGC = c%nGhostCells

    call c%put('Longitude', real(Longitude(1-nGC:nLons+nGC, iBlock) * cRadToDeg, output_kind))
    call c%put('Latitude',  real(Latitude(1-nGC:nLats+nGC, iBlock) * cRadToDeg, output_kind))
    call c%put('Altitude',  real(Altitude_GB(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, iBlock), output_kind))

    allocate(zero_buf(nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC))
    zero_buf = 0.0_output_kind
    call c%put('6300 A Emission', zero_buf)
    call c%put('PhotoElectronUp', zero_buf)
    call c%put('PhotoElectronDown', zero_buf)
    deallocate(zero_buf)
  end subroutine fill_3dglo

  ! ---------------------------------------------------------------------------
  ! define_schema_3dmag — register variables for 3DMAG.
  ! ---------------------------------------------------------------------------
  subroutine define_schema_3dmag(c)
    use ModGITM, only: nLons, nLats, nAlts
    type(OutputContainer), intent(inout) :: c
    integer :: nGC

    c%cType    = '3DMAG'
    c%gridKind = GRID_GEO_3D
    c%nGhostCells = 2
    nGC = c%nGhostCells

    call c%define_var('Longitude', units='degrees_east', shape3=[nLons + 2*nGC, 1, 1], is_axis=.true., &
                      longName='Geographic longitude')
    call c%define_var('Latitude', units='degrees_north', shape3=[nLats + 2*nGC, 1, 1], is_axis=.true., &
                      longName='Geographic latitude')
    call c%define_var('Altitude', units='m', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Altitude')
    call c%define_var('Magnetic Latitude', units='deg', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Magnetic latitude')
    call c%define_var('Magnetic Longitude', units='deg', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Magnetic longitude')
    call c%define_var('B.F. East', units='T', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Magnetic field eastward')
    call c%define_var('B.F. North', units='T', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Magnetic field northward')
    call c%define_var('B.F. Vertical', units='T', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Magnetic field vertical')
    call c%define_var('B.F. Magnitude', units='T', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Magnetic field magnitude')
  end subroutine define_schema_3dmag

  ! ---------------------------------------------------------------------------
  ! fill_3dmag — put 3DMAG variables into the container.
  ! ---------------------------------------------------------------------------
  subroutine fill_3dmag(c, iBlock)
    use ModGITM
    use ModConst, only: cRadToDeg
    type(OutputContainer), intent(inout) :: c
    integer :: nGC
    integer, intent(in) :: iBlock

    call c%prepare(iBlock)
    if (.not. c%this_rank_writes) return
    nGC = c%nGhostCells

    call c%put('Longitude', real(Longitude(1-nGC:nLons+nGC, iBlock) * cRadToDeg, output_kind))
    call c%put('Latitude',  real(Latitude(1-nGC:nLats+nGC, iBlock) * cRadToDeg, output_kind))
    call c%put('Altitude',  real(Altitude_GB(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, iBlock), output_kind))
    call c%put('Magnetic Latitude', real(mLatitude(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, iBlock), output_kind))
    call c%put('Magnetic Longitude', real(mLongitude(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, iBlock), output_kind))
    call c%put('B.F. East', real(B0(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, 1, iBlock), output_kind))
    call c%put('B.F. North', real(B0(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, 2, iBlock), output_kind))
    call c%put('B.F. Vertical', real(B0(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, 3, iBlock), output_kind))
    call c%put('B.F. Magnitude', real(B0(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, 4, iBlock), output_kind))
  end subroutine fill_3dmag

  ! ---------------------------------------------------------------------------
  ! define_schema_3dhme — register variables for 3DHME.
  ! ---------------------------------------------------------------------------
  subroutine define_schema_3dhme(c)
    use ModGITM, only: nLons, nLats, nAlts
    use ModPlanet, only: nSpeciesTotal, nSpecies, nIons, cSpecies, cIons
    type(OutputContainer), intent(inout) :: c
    integer :: nGC
    integer :: i

    c%cType    = '3DHME'
    c%gridKind = GRID_HIME
    c%nGhostCells = 0
    nGC = c%nGhostCells

    call c%define_var('Longitude', units='degrees_east', shape3=[nLons + 2*nGC, 1, 1], is_axis=.true., &
                      longName='Geographic longitude')
    call c%define_var('Latitude', units='degrees_north', shape3=[nLats + 2*nGC, 1, 1], is_axis=.true., &
                      longName='Geographic latitude')
    call c%define_var('Altitude', units='m', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Altitude')
    call c%define_var('Rho', units='kg/m3', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Total neutral mass density')
    do i = 1, nSpeciesTotal
      call c%define_var('['//trim(cSpecies(i))//']', units='/m3', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], &
                        longName='Number density of '//trim(cSpecies(i)))
    enddo
    call c%define_var('Temperature', units='K', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Neutral temperature')
    call c%define_var('Vn (east)', units='m/s', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Eastward neutral wind')
    call c%define_var('Vn (north)', units='m/s', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Northward neutral wind')
    call c%define_var('Vn (up)', units='m/s', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Vertical neutral wind')
    do i = 1, nSpecies
      call c%define_var('Vn (up,'//trim(cSpecies(i))//')', units='m/s', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], &
                        longName='Vertical neutral wind for '//trim(cSpecies(i)))
    enddo
    do i = 1, nIons
      call c%define_var('['//trim(cIons(i))//']', units='/m3', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], &
                        longName='Number density of '//trim(cIons(i)))
    enddo
    call c%define_var('eTemperature', units='K', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Electron temperature')
    call c%define_var('iTemperature', units='K', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Ion temperature')
    call c%define_var('Vi (east)', units='m/s', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Eastward ion drift')
    call c%define_var('Vi (north)', units='m/s', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Northward ion drift')
    call c%define_var('Vi (up)', units='m/s', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Vertical ion drift')
    call c%define_var('PhotoElectron Heating', units='K', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Photoelectron heating')
    call c%define_var('Joule Heating', units='K', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Joule heating')
    call c%define_var('Specific Heat', units='J/kg/K', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Specific heat Cp')
    call c%define_var('Magnetic Latitude', units='deg', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Magnetic latitude')
    call c%define_var('Magnetic Longitude', units='deg', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Magnetic longitude')
    call c%define_var('B.F. East', units='T', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Magnetic field eastward')
    call c%define_var('B.F. North', units='T', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Magnetic field northward')
    call c%define_var('B.F. Vertical', units='T', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Magnetic field vertical')
    call c%define_var('B.F. Magnitude', units='T', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Magnetic field magnitude')
    call c%define_var('Potential', units='V', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Electric potential')
    call c%define_var('PotentialY', units='V', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Potential Y component')
    call c%define_var('E.F. East', units='V/m', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Electric field eastward')
    call c%define_var('E.F. North', units='V/m', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Electric field northward')
    call c%define_var('E.F. Vertical', units='V/m', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Electric field vertical')
  end subroutine define_schema_3dhme

  ! ---------------------------------------------------------------------------
  ! fill_3dhme — put 3DHME regional variables into the container.
  ! ---------------------------------------------------------------------------
  subroutine fill_3dhme(c, iBlock)
    use ModGITM
    use ModConst, only: cRadToDeg
    use ModConstants, only: pi
    use ModInputs, only: HIMEPlotLonStart, HIMEPlotLonEnd, HIMEPlotLatStart, HIMEPlotLatEnd
    use ModPlanet, only: nSpeciesTotal, nSpecies, nIons, cSpecies, cIons
    use ModSources, only: PhotoElectronHeating, JouleHeating
    type(OutputContainer), intent(inout) :: c
    integer :: nGC
    integer, intent(in) :: iBlock
    integer :: i, iLat, iLon
    logical :: DoSaveHIMEPlot

    DoSaveHIMEPlot = .false.
    do iLat = -1, nLats + 2
      do iLon = -1, nLons + 2
        if (Longitude(iLon, iBlock) >= HIMEPlotLonStart*pi/180.0_output_kind &
            .and. Longitude(iLon, iBlock) <= HIMEPlotLonEnd*pi/180.0_output_kind &
            .and. Latitude(iLat, iBlock) >= HIMEPlotLatStart*pi/180.0_output_kind &
            .and. Latitude(iLat, iBlock) <= HIMEPlotLatEnd*pi/180.0_output_kind) then
          DoSaveHIMEPlot = .true.
          exit
        endif
      enddo
    enddo

    c%this_rank_writes = DoSaveHIMEPlot
    if (.not. c%this_rank_writes) return
    nGC = c%nGhostCells

    call c%put('Longitude', real(Longitude(1-nGC:nLons+nGC, iBlock) * cRadToDeg, output_kind))
    call c%put('Latitude',  real(Latitude(1-nGC:nLats+nGC, iBlock) * cRadToDeg, output_kind))
    call c%put('Altitude',  real(Altitude_GB(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, iBlock), output_kind))
    call c%put('Rho',       real(Rho(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, iBlock), output_kind))
    do i = 1, nSpeciesTotal
      call c%put('['//trim(cSpecies(i))//']', &
                 real(NDensityS(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, i, iBlock), output_kind))
    end do
    call c%put('Temperature', &
               real(Temperature(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, iBlock) * &
                    TempUnit(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC), output_kind))
    call c%put('Vn (east)',  real(Velocity(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, 1, iBlock), output_kind))
    call c%put('Vn (north)', real(Velocity(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, 2, iBlock), output_kind))
    call c%put('Vn (up)',    real(Velocity(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, 3, iBlock), output_kind))
    do i = 1, nSpecies
      call c%put('Vn (up,'//trim(cSpecies(i))//')', &
                 real(VerticalVelocity(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, i, iBlock), output_kind))
    end do
    do i = 1, nIons
      call c%put('['//trim(cIons(i))//']', &
                 real(IDensityS(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, i, iBlock), output_kind))
    end do
    call c%put('eTemperature', real(eTemperature(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, iBlock), output_kind))
    call c%put('iTemperature', real(ITemperature(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, iBlock), output_kind))
    call c%put('Vi (east)',    real(IVelocity(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, 1, iBlock), output_kind))
    call c%put('Vi (north)',   real(IVelocity(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, 2, iBlock), output_kind))
    call c%put('Vi (up)',      real(IVelocity(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, 3, iBlock), output_kind))
    call c%put('PhotoElectron Heating', &
               real(PhotoElectronHeating(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, iBlock) * &
                    dt * TempUnit(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC), output_kind))
    call c%put('Joule Heating', &
               real(JouleHeating(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC) * &
                    dt * TempUnit(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC), output_kind))
    call c%put('Specific Heat', real(cp(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, iBlock), output_kind))
    call c%put('Magnetic Latitude', real(mLatitude(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, iBlock), output_kind))
    call c%put('Magnetic Longitude', real(mLongitude(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, iBlock), output_kind))
    call c%put('B.F. East', real(B0(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, 1, iBlock), output_kind))
    call c%put('B.F. North', real(B0(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, 2, iBlock), output_kind))
    call c%put('B.F. Vertical', real(B0(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, 3, iBlock), output_kind))
    call c%put('B.F. Magnitude', real(B0(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, 4, iBlock), output_kind))
    call c%put('Potential', real(Potential(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, iBlock), output_kind))
    call c%put('PotentialY', real(PotentialY(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, iBlock), output_kind))
    call c%put('E.F. East', real(EField(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, 1), output_kind))
    call c%put('E.F. North', real(EField(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, 2), output_kind))
    call c%put('E.F. Vertical', real(EField(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, 3), output_kind))
  end subroutine fill_3dhme

  ! ---------------------------------------------------------------------------
  ! define_schema_3dmoh — register variables for 3DMOH.
  ! ---------------------------------------------------------------------------
  subroutine define_schema_3dmoh(c)
    use ModGITM, only: nLons, nLats, nAlts
    type(OutputContainer), intent(inout) :: c
    integer :: nGC

    c%cType    = '3DMOH'
    c%gridKind = GRID_GEO_3D
    c%nGhostCells = 0
    nGC = c%nGhostCells

    call c%define_var('Longitude', units='degrees_east', shape3=[nLons + 2*nGC, 1, 1], is_axis=.true., &
                      longName='Geographic longitude')
    call c%define_var('Latitude', units='degrees_north', shape3=[nLats + 2*nGC, 1, 1], is_axis=.true., &
                      longName='Geographic latitude')
    call c%define_var('Altitude', units='m', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Altitude')
    call c%define_var('Rho', units='kg/m3', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Mass density')
    call c%define_var('Vn (east)', units='m/s', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Eastward wind')
    call c%define_var('Vn (north)', units='m/s', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Northward wind')
    call c%define_var('Visc_Ve_rshear', units='m/s2', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Viscous east shear accel')
    call c%define_var('Visc_Vn_rshear', units='m/s2', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Viscous north shear accel')
    call c%define_var('IonDrag_east', units='m/s2', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Ion drag east accel')
    call c%define_var('IonDrag_north', units='m/s2', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Ion drag north accel')
    call c%define_var('Horiz_Adv_Ve', units='m/s2', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Horiz advection east accel')
    call c%define_var('Horiz_Adv_Vn', units='m/s2', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Horiz advection north accel')
    call c%define_var('PressGrad (east)', units='m/s2', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Press gradient east accel')
    call c%define_var('PressGrad (north)', units='m/s2', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Press gradient north accel')
    call c%define_var('Coriolis_east', units='m/s2', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Coriolis east accel')
    call c%define_var('Coriolis_north', units='m/s2', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Coriolis north accel')
    call c%define_var('Centrifugal_north', units='m/s2', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Centrifugal north accel')
  end subroutine define_schema_3dmoh

  ! ---------------------------------------------------------------------------
  ! fill_3dmoh — put 3DMOH variables into the container.
  ! ---------------------------------------------------------------------------
  subroutine fill_3dmoh(c, iBlock)
    use ModGITM
    use ModConst, only: cRadToDeg
    use ModSources
    type(OutputContainer), intent(inout) :: c
    integer :: nGC
    integer, intent(in) :: iBlock

    call c%prepare(iBlock)
    if (.not. c%this_rank_writes) return
    nGC = c%nGhostCells

    call c%put('Longitude', real(Longitude(1-nGC:nLons+nGC, iBlock) * cRadToDeg, output_kind))
    call c%put('Latitude',  real(Latitude(1-nGC:nLats+nGC, iBlock) * cRadToDeg, output_kind))
    call c%put('Altitude',  real(Altitude_GB(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, iBlock), output_kind))
    call c%put('Rho',       real(Rho(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, iBlock), output_kind))
    call c%put('Vn (east)',  real(Velocity(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, iEast_, iBlock), output_kind))
    call c%put('Vn (north)', real(Velocity(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, iNorth_, iBlock), output_kind))
    call c%put('Visc_Ve_rshear', real(Viscosity(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, iEast_)/dt, output_kind))
    call c%put('Visc_Vn_rshear', real(Viscosity(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, iNorth_)/dt, output_kind))
    call c%put('IonDrag_east', real(IonDrag(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, iEast_), output_kind))
    call c%put('IonDrag_north', real(IonDrag(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, iNorth_), output_kind))
    call c%put('Horiz_Adv_Ve', real(HorizAdvection(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, 1), output_kind))
    call c%put('Horiz_Adv_Vn', real(HorizAdvection(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, 2), output_kind))
    call c%put('PressGrad (east)', real(HorizPressureGrad(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, 1), output_kind))
    call c%put('PressGrad (north)', real(HorizPressureGrad(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, 2), output_kind))
    call c%put('Coriolis_east', real(HorizCoriolis(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, 1), output_kind))
    call c%put('Coriolis_north', real(HorizCoriolis(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, 2), output_kind))
    call c%put('Centrifugal_north', real(Centrifugal(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, 2), output_kind))
  end subroutine fill_3dmoh

  ! ---------------------------------------------------------------------------
  ! define_schema_3dmov — register variables for 3DMOV.
  ! ---------------------------------------------------------------------------
  subroutine define_schema_3dmov(c)
    use ModGITM, only: nLons, nLats, nAlts
    use ModPlanet, only: nSpecies, cSpecies
    type(OutputContainer), intent(inout) :: c
    integer :: nGC
    integer :: i

    c%cType    = '3DMOV'
    c%gridKind = GRID_GEO_3D
    c%nGhostCells = 0
    nGC = c%nGhostCells

    call c%define_var('Longitude', units='degrees_east', shape3=[nLons + 2*nGC, 1, 1], is_axis=.true., &
                      longName='Geographic longitude')
    call c%define_var('Latitude', units='degrees_north', shape3=[nLats + 2*nGC, 1, 1], is_axis=.true., &
                      longName='Geographic latitude')
    call c%define_var('Altitude', units='m', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Altitude')
    do i = 1, nSpecies
      call c%define_var('Vn (up,'//trim(cSpecies(i))//')', units='m/s', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], &
                        longName='Vertical wind for '//trim(cSpecies(i)))
      call c%define_var('IonDrag_up_'//trim(cSpecies(i)), units='m/s2', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], &
                        longName='Ion drag upward for '//trim(cSpecies(i)))
      call c%define_var('Horiz_Adv_Vn (up,'//trim(cSpecies(i))//')', units='m/s2', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], &
                        longName='Horiz advection upward for '//trim(cSpecies(i)))
    enddo
    call c%define_var('Coriolis_up', units='m/s2', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Coriolis upward accel')
    call c%define_var('Centrif_up', units='m/s2', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Centrifugal upward accel')
    call c%define_var('EffectiveGravity', units='m/s2', shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC], longName='Effective gravity accel')
  end subroutine define_schema_3dmov

  ! ---------------------------------------------------------------------------
  ! fill_3dmov — put 3DMOV variables into the container.
  ! ---------------------------------------------------------------------------
  subroutine fill_3dmov(c, iBlock)
    use ModGITM
    use ModConst, only: cRadToDeg
    use ModSources
    type(OutputContainer), intent(inout) :: c
    integer :: nGC
    integer, intent(in) :: iBlock
    integer :: iSpecies

    call c%prepare(iBlock)
    if (.not. c%this_rank_writes) return
    nGC = c%nGhostCells

    call c%put('Longitude', real(Longitude(1-nGC:nLons+nGC, iBlock) * cRadToDeg, output_kind))
    call c%put('Latitude',  real(Latitude(1-nGC:nLats+nGC, iBlock) * cRadToDeg, output_kind))
    call c%put('Altitude',  real(Altitude_GB(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, iBlock), output_kind))
    do iSpecies = 1, nSpecies
      call c%put('Vn (up,'//trim(cSpecies(iSpecies))//')', &
                 real(VerticalVelocity(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, iSpecies, iBlock), output_kind))
      call c%put('IonDrag_up_'//trim(cSpecies(iSpecies)), &
                 real(VerticalIonDrag(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, iSpecies), output_kind))
      call c%put('Horiz_Adv_Vn (up,'//trim(cSpecies(iSpecies))//')', &
                 real(VertAdvection(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC, iSpecies), output_kind))
    enddo
    call c%put('Coriolis_up', real(VertCoriolis(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC), output_kind))
    call c%put('Centrif_up', real(VertCentrifugal(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC), output_kind))
    call c%put('EffectiveGravity', real(EffectiveGravity(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1-nGC:nAlts+nGC), output_kind))
  end subroutine fill_3dmov

  ! ---------------------------------------------------------------------------
  ! define_schema_2danc — register variables for 2DANC.
  ! ---------------------------------------------------------------------------
  subroutine define_schema_2danc(c)
    use ModGITM, only: nLons, nLats
    type(OutputContainer), intent(inout) :: c
    integer :: nGC

    c%cType    = '2DANC'
    c%gridKind = GRID_GEO_2D
    c%nGhostCells = 0
    nGC = c%nGhostCells

    call c%define_var('Longitude', units='degrees_east', shape3=[nLons + 2*nGC, 1, 1], is_axis=.true., &
                      longName='Geographic longitude')
    call c%define_var('Latitude', units='degrees_north', shape3=[nLats + 2*nGC, 1, 1], is_axis=.true., &
                      longName='Geographic latitude')
    call c%define_var('Altitude', units='m', shape3=[nLons + 2*nGC, nLats + 2*nGC, 1], longName='Altitude')
    call c%define_var('LT', units='hr', shape3=[nLons + 2*nGC, nLats + 2*nGC, 1], longName='Local Time')
    call c%define_var('SZA', units='rad', shape3=[nLons + 2*nGC, nLats + 2*nGC, 1], longName='Solar zenith angle')
    call c%define_var('TEC', units='TECU', shape3=[nLons + 2*nGC, nLats + 2*nGC, 1], longName='Total electron content')
    call c%define_var('JH_int', units='W/m2', shape3=[nLons + 2*nGC, nLats + 2*nGC, 1], longName='Alt-integrated Joule heating')
    call c%define_var('HT_int', units='W/m2', shape3=[nLons + 2*nGC, nLats + 2*nGC, 1], longName='Alt-integrated Heat transfer')
    call c%define_var('EUV_int', units='W/m2', shape3=[nLons + 2*nGC, nLats + 2*nGC, 1], longName='Alt-integrated EUV heating')
    call c%define_var('PE_int', units='W/m2', shape3=[nLons + 2*nGC, nLats + 2*nGC, 1], longName='Alt-integrated PE heating')
    call c%define_var('Chem_int', units='W/m2', shape3=[nLons + 2*nGC, nLats + 2*nGC, 1], longName='Alt-integrated Chem heating')
    call c%define_var('RadCool_int', units='W/m2', shape3=[nLons + 2*nGC, nLats + 2*nGC, 1], longName='Alt-integrated Rad cooling')
    call c%define_var('CO2Cool_int', units='W/m2', shape3=[nLons + 2*nGC, nLats + 2*nGC, 1], longName='Alt-integrated CO2 cooling')
    call c%define_var('NOCool_int', units='W/m2', shape3=[nLons + 2*nGC, nLats + 2*nGC, 1], longName='Alt-integrated NO cooling')
    call c%define_var('OCool_int', units='W/m2', shape3=[nLons + 2*nGC, nLats + 2*nGC, 1], longName='Alt-integrated O cooling')
  end subroutine define_schema_2danc

  ! ---------------------------------------------------------------------------
  ! fill_2danc — put 2DANC integrated variables into the container.
  ! ---------------------------------------------------------------------------
  subroutine fill_2danc(c, iBlock)
    use ModGITM
    use ModConst, only: cRadToDeg
    use ModElectrodynamics
    use ModEUV, only: Sza
    use ModSources
    type(OutputContainer), intent(inout) :: c
    integer :: nGC
    integer, intent(in) :: iBlock
    real(output_kind) :: lt_buf(nLons, nLats)
    integer :: iLat

    call c%prepare(iBlock)
    if (.not. c%this_rank_writes) return
    nGC = c%nGhostCells

    call calc_vtec(iBlock)

    call c%put('Longitude', real(Longitude(1-nGC:nLons+nGC, iBlock) * cRadToDeg, output_kind))
    call c%put('Latitude',  real(Latitude(1-nGC:nLats+nGC, iBlock) * cRadToDeg, output_kind))
    call c%put('Altitude',  real(Altitude_GB(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, nAlts, iBlock), output_kind))

    ! Replicate LocalTime(1-nGC:nLons+nGC) along Latitude dimension
    do iLat = 1, nLats
      lt_buf(:, iLat) = real(LocalTime(1-nGC:nLons+nGC), output_kind)
    enddo
    call c%put('LT', lt_buf)

    call c%put('SZA', real(Sza(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, iBlock), output_kind))
    call c%put('TEC', real(VTEC(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, iBlock), output_kind))
    call c%put('JH_int', real(JouleHeating2d(1-nGC:nLons+nGC, 1-nGC:nLats+nGC), output_kind))
    call c%put('HT_int', real(HeatTransfer2d(1-nGC:nLons+nGC, 1-nGC:nLats+nGC), output_kind))
    call c%put('EUV_int', real(EuvHeating2d(1-nGC:nLons+nGC, 1-nGC:nLats+nGC), output_kind))
    call c%put('PE_int', real(PhotoElectronHeating2d(1-nGC:nLons+nGC, 1-nGC:nLats+nGC), output_kind))
    call c%put('Chem_int', real(ChemicalHeating2d(1-nGC:nLons+nGC, 1-nGC:nLats+nGC), output_kind))
    call c%put('RadCool_int', real(RadiativeCooling2d(1-nGC:nLons+nGC, 1-nGC:nLats+nGC), output_kind))
    call c%put('CO2Cool_int', real(CO2Cooling2d(1-nGC:nLons+nGC, 1-nGC:nLats+nGC), output_kind))
    call c%put('NOCool_int', real(NOCooling2d(1-nGC:nLons+nGC, 1-nGC:nLats+nGC), output_kind))
    call c%put('OCool_int', real(OCooling2d(1-nGC:nLons+nGC, 1-nGC:nLats+nGC), output_kind))
  end subroutine fill_2danc

  ! ---------------------------------------------------------------------------
  ! define_schema_2dhme — register variables for 2DHME.
  ! ---------------------------------------------------------------------------
  subroutine define_schema_2dhme(c)
    use ModGITM, only: nLons, nLats
    type(OutputContainer), intent(inout) :: c
    integer :: nGC

    c%cType    = '2DHME'
    c%gridKind = GRID_HIME
    c%nGhostCells = 0
    nGC = c%nGhostCells

    call c%define_var('Longitude', units='degrees_east', shape3=[nLons + 2*nGC, 1, 1], is_axis=.true., &
                      longName='Geographic longitude')
    call c%define_var('Latitude', units='degrees_north', shape3=[nLats + 2*nGC, 1, 1], is_axis=.true., &
                      longName='Geographic latitude')
    call c%define_var('Altitude', units='m', shape3=[nLons + 2*nGC, nLats + 2*nGC, 1], longName='Altitude')
    call c%define_var('LT', units='hr', shape3=[nLons + 2*nGC, nLats + 2*nGC, 1], longName='Local Time')
    call c%define_var('TEC', units='TECU', shape3=[nLons + 2*nGC, nLats + 2*nGC, 1], longName='Total electron content')
  end subroutine define_schema_2dhme

  ! ---------------------------------------------------------------------------
  ! fill_2dhme — put 2DHME regional variables into the container.
  ! ---------------------------------------------------------------------------
  subroutine fill_2dhme(c, iBlock)
    use ModGITM
    use ModConst, only: cRadToDeg
    use ModConstants, only: pi
    use ModInputs, only: HIMEPlotLonStart, HIMEPlotLonEnd, HIMEPlotLatStart, HIMEPlotLatEnd
    use ModElectrodynamics
    type(OutputContainer), intent(inout) :: c
    integer :: nGC
    integer, intent(in) :: iBlock
    real(output_kind) :: lt_buf(nLons, nLats)
    integer :: iLat, iLon
    logical :: DoSaveHIMEPlot

    DoSaveHIMEPlot = .false.
    do iLat = -1, nLats + 2
      do iLon = -1, nLons + 2
        if (Longitude(iLon, iBlock) >= HIMEPlotLonStart*pi/180.0_output_kind &
            .and. Longitude(iLon, iBlock) <= HIMEPlotLonEnd*pi/180.0_output_kind &
            .and. Latitude(iLat, iBlock) >= HIMEPlotLatStart*pi/180.0_output_kind &
            .and. Latitude(iLat, iBlock) <= HIMEPlotLatEnd*pi/180.0_output_kind) then
          DoSaveHIMEPlot = .true.
          exit
        endif
      enddo
    enddo

    c%this_rank_writes = DoSaveHIMEPlot
    if (.not. c%this_rank_writes) return
    nGC = c%nGhostCells

    call calc_vtec(iBlock)

    call c%put('Longitude', real(Longitude(1-nGC:nLons+nGC, iBlock) * cRadToDeg, output_kind))
    call c%put('Latitude',  real(Latitude(1-nGC:nLats+nGC, iBlock) * cRadToDeg, output_kind))
    call c%put('Altitude',  real(Altitude_GB(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, 1, iBlock), output_kind))

    do iLat = 1, nLats
      lt_buf(:, iLat) = real(LocalTime(1-nGC:nLons+nGC), output_kind)
    enddo
    call c%put('LT', lt_buf)
    call c%put('TEC', real(VTEC(1-nGC:nLons+nGC, 1-nGC:nLats+nGC, iBlock), output_kind))
  end subroutine fill_2dhme

  ! ---------------------------------------------------------------------------
  ! define_schema_1dall — register variables for 1DALL.
  ! ---------------------------------------------------------------------------
  subroutine define_schema_1dall(c)
    use ModGITM, only: nAlts
    use ModPlanet, only: nSpeciesTotal, nSpecies, nIons, cSpecies, cIons
    type(OutputContainer), intent(inout) :: c
    integer :: i

    c%cType    = '1DALL'
    c%gridKind = GRID_GEO_2D

    call c%define_var('Longitude', units='degrees_east', shape3=[nAlts, 1, 1], is_axis=.true.)
    call c%define_var('Latitude', units='degrees_north', shape3=[nAlts, 1, 1], is_axis=.true.)
    call c%define_var('Altitude', units='m', shape3=[nAlts, 1, 1])
    call c%define_var('Rho', units='kg/m3', shape3=[nAlts, 1, 1])
    do i = 1, nSpeciesTotal
      call c%define_var('['//trim(cSpecies(i))//']', units='/m3', shape3=[nAlts, 1, 1])
    enddo
    call c%define_var('Temperature', units='K', shape3=[nAlts, 1, 1])
    call c%define_var('Vn (east)', units='m/s', shape3=[nAlts, 1, 1])
    call c%define_var('Vn (north)', units='m/s', shape3=[nAlts, 1, 1])
    call c%define_var('Vn (up)', units='m/s', shape3=[nAlts, 1, 1])
    do i = 1, nSpecies
      call c%define_var('Vn (up,'//trim(cSpecies(i))//')', units='m/s', shape3=[nAlts, 1, 1])
    enddo
    do i = 1, nIons
      call c%define_var('['//trim(cIons(i))//']', units='/m3', shape3=[nAlts, 1, 1])
    enddo
    call c%define_var('eTemperature', units='K', shape3=[nAlts, 1, 1])
    call c%define_var('iTemperature', units='K', shape3=[nAlts, 1, 1])
    call c%define_var('Vi (east)', units='m/s', shape3=[nAlts, 1, 1])
    call c%define_var('Vi (north)', units='m/s', shape3=[nAlts, 1, 1])
    call c%define_var('Vi (up)', units='m/s', shape3=[nAlts, 1, 1])
  end subroutine define_schema_1dall

  ! ---------------------------------------------------------------------------
  ! fill_1dall — put 1DALL interpolated vertical column data.
  ! ---------------------------------------------------------------------------
  subroutine fill_1dall(c, iBlock, iiLon, iiLat, rLon, rLat)
    use ModGITM
    use ModConst, only: cRadToDeg
    use ModPlanet, only: nSpeciesTotal, nSpecies, nIons, cSpecies, cIons
    type(OutputContainer), intent(inout) :: c
    integer, intent(in) :: iBlock, iiLon, iiLat
    real, intent(in) :: rLon, rLat
    real(output_kind) :: tmp_lon(nAlts), tmp_lat(nAlts), tmp_data(nAlts)
    real :: slice2d(0:nLons + 1, 0:nLats + 1)
    integer :: iAlt, i

    call c%prepare(iBlock)
    if (.not. c%this_rank_writes) return

    tmp_lon = real((rLon * Longitude(iiLon, iBlock) + (1.0 - rLon) * Longitude(iiLon + 1, iBlock)) * cRadToDeg, output_kind)
    tmp_lat = real((rLat * Latitude(iiLat, iBlock) + (1.0 - rLat) * Latitude(iiLat + 1, iBlock)) * cRadToDeg, output_kind)
    call c%put('Longitude', tmp_lon)
    call c%put('Latitude',  tmp_lat)

    do iAlt = 1, nAlts
      tmp_data(iAlt) = real(Altitude_GB(iiLon, iiLat, iAlt, iBlock), output_kind)
    enddo
    call c%put('Altitude', tmp_data)

    do iAlt = 1, nAlts
      slice2d = Rho(0:nLons + 1, 0:nLats + 1, iAlt, iBlock)
      tmp_data(iAlt) = real(inter2d(slice2d, iiLon, iiLat, rLon, rLat), output_kind)
    enddo
    call c%put('Rho', tmp_data)

    do i = 1, nSpeciesTotal
      do iAlt = 1, nAlts
        slice2d = NDensityS(0:nLons + 1, 0:nLats + 1, iAlt, i, iBlock)
        tmp_data(iAlt) = real(inter2d(slice2d, iiLon, iiLat, rLon, rLat), output_kind)
      enddo
      call c%put('['//trim(cSpecies(i))//']', tmp_data)
    enddo

    do iAlt = 1, nAlts
      slice2d = Temperature(0:nLons + 1, 0:nLats + 1, iAlt, iBlock) * &
                TempUnit(0:nLons + 1, 0:nLats + 1, iAlt)
      tmp_data(iAlt) = real(inter2d(slice2d, iiLon, iiLat, rLon, rLat), output_kind)
    enddo
    call c%put('Temperature', tmp_data)

    do i = 1, 3
      do iAlt = 1, nAlts
        slice2d = Velocity(0:nLons + 1, 0:nLats + 1, iAlt, i, iBlock)
        tmp_data(iAlt) = real(inter2d(slice2d, iiLon, iiLat, rLon, rLat), output_kind)
      enddo
      if (i == 1) call c%put('Vn (east)', tmp_data)
      if (i == 2) call c%put('Vn (north)', tmp_data)
      if (i == 3) call c%put('Vn (up)', tmp_data)
    enddo

    do i = 1, nSpecies
      do iAlt = 1, nAlts
        slice2d = VerticalVelocity(0:nLons + 1, 0:nLats + 1, iAlt, i, iBlock)
        tmp_data(iAlt) = real(inter2d(slice2d, iiLon, iiLat, rLon, rLat), output_kind)
      enddo
      call c%put('Vn (up,'//trim(cSpecies(i))//')', tmp_data)
    enddo

    do i = 1, nIons
      do iAlt = 1, nAlts
        slice2d = IDensityS(0:nLons + 1, 0:nLats + 1, iAlt, i, iBlock)
        tmp_data(iAlt) = real(inter2d(slice2d, iiLon, iiLat, rLon, rLat), output_kind)
      enddo
      call c%put('['//trim(cIons(i))//']', tmp_data)
    enddo

    do iAlt = 1, nAlts
      slice2d = eTemperature(0:nLons + 1, 0:nLats + 1, iAlt, iBlock)
      tmp_data(iAlt) = real(inter2d(slice2d, iiLon, iiLat, rLon, rLat), output_kind)
    enddo
    call c%put('eTemperature', tmp_data)

    do iAlt = 1, nAlts
      slice2d = ITemperature(0:nLons + 1, 0:nLats + 1, iAlt, iBlock)
      tmp_data(iAlt) = real(inter2d(slice2d, iiLon, iiLat, rLon, rLat), output_kind)
    enddo
    call c%put('iTemperature', tmp_data)

    do i = 1, 3
      do iAlt = 1, nAlts
        slice2d = IVelocity(0:nLons + 1, 0:nLats + 1, iAlt, i, iBlock)
        tmp_data(iAlt) = real(inter2d(slice2d, iiLon, iiLat, rLon, rLat), output_kind)
      enddo
      if (i == 1) call c%put('Vi (east)', tmp_data)
      if (i == 2) call c%put('Vi (north)', tmp_data)
      if (i == 3) call c%put('Vi (up)', tmp_data)
    enddo
  end subroutine fill_1dall

  ! ---------------------------------------------------------------------------
  ! define_schema_0dall — register variables for 0DALL.
  ! ---------------------------------------------------------------------------
  subroutine define_schema_0dall(c)
    use ModPlanet, only: nSpeciesTotal, nSpecies, nIons, cSpecies, cIons
    type(OutputContainer), intent(inout) :: c
    integer :: i

    c%cType    = '0DALL'
    c%gridKind = GRID_GEO_2D

    call c%define_var('Longitude', units='degrees_east', shape3=[1, 1, 1], is_axis=.true.)
    call c%define_var('Latitude', units='degrees_north', shape3=[1, 1, 1], is_axis=.true.)
    call c%define_var('Altitude', units='m', shape3=[1, 1, 1])
    call c%define_var('Rho', units='kg/m3', shape3=[1, 1, 1])
    do i = 1, nSpeciesTotal
      call c%define_var('['//trim(cSpecies(i))//']', units='/m3', shape3=[1, 1, 1])
    enddo
    call c%define_var('Temperature', units='K', shape3=[1, 1, 1])
    call c%define_var('Vn (east)', units='m/s', shape3=[1, 1, 1])
    call c%define_var('Vn (north)', units='m/s', shape3=[1, 1, 1])
    call c%define_var('Vn (up)', units='m/s', shape3=[1, 1, 1])
    do i = 1, nSpecies
      call c%define_var('Vn (up,'//trim(cSpecies(i))//')', units='m/s', shape3=[1, 1, 1])
    enddo
    do i = 1, nIons
      call c%define_var('['//trim(cIons(i))//']', units='/m3', shape3=[1, 1, 1])
    enddo
    call c%define_var('eTemperature', units='K', shape3=[1, 1, 1])
    call c%define_var('iTemperature', units='K', shape3=[1, 1, 1])
    call c%define_var('Vi (east)', units='m/s', shape3=[1, 1, 1])
    call c%define_var('Vi (north)', units='m/s', shape3=[1, 1, 1])
    call c%define_var('Vi (up)', units='m/s', shape3=[1, 1, 1])
    do i = 1, nSpecies
      call c%define_var(trim(cSpecies(i))//' Mixing Ratio', units='', shape3=[1, 1, 1])
    enddo
    call c%define_var('RadCooling', units='K/s', shape3=[1, 1, 1])
    call c%define_var('EuvHeating', units='K/s', shape3=[1, 1, 1])
    call c%define_var('Conduction', units='K/s', shape3=[1, 1, 1])
    call c%define_var('Heat Balance Total', units='K/s', shape3=[1, 1, 1])
    call c%define_var('Heating Efficiency', units='', shape3=[1, 1, 1])
  end subroutine define_schema_0dall

  ! ---------------------------------------------------------------------------
  ! fill_0dall — put 0DALL 3D-interpolated point data.
  ! ---------------------------------------------------------------------------
  subroutine fill_0dall(c, iBlock, iiLon, iiLat, iiAlt, rLon, rLat, rAlt)
    use ModGITM
    use ModConst, only: cRadToDeg
    use ModEUV, only: HeatingEfficiency_CB
    use ModSources, only: JouleHeating, RadCooling, EuvHeating, Conduction
    use ModPlanet, only: nSpeciesTotal, nSpecies, nIons, cSpecies, cIons
    type(OutputContainer), intent(inout) :: c
    integer, intent(in) :: iBlock, iiLon, iiLat, iiAlt
    real, intent(in) :: rLon, rLat, rAlt
    real :: Tmp3d(0:nLons + 1, 0:nLats + 1, 0:nAlts + 1)
    integer :: jAlt, i

    call c%prepare(iBlock)
    if (.not. c%this_rank_writes) return

    jAlt = max(min(iiAlt, nAlts), 1)

    call c%put('Longitude', real((rLon * Longitude(iiLon, iBlock) + (1.0 - rLon) * Longitude(iiLon + 1, iBlock)) * cRadToDeg, output_kind))
    call c%put('Latitude',  real((rLat * Latitude(iiLat, iBlock) + (1.0 - rLat) * Latitude(iiLat + 1, iBlock)) * cRadToDeg, output_kind))
    call c%put('Altitude',  real(rAlt * Altitude_GB(iiLon, iiLat, iiAlt, iBlock) + &
                                 (1.0 - rAlt) * Altitude_GB(iiLon + 1, iiLat + 1, iiAlt + 1, iBlock), output_kind))

    Tmp3d = Rho(0:nLons + 1, 0:nLats + 1, 0:nAlts + 1, iBlock)
    call c%put('Rho', real(inter3d(Tmp3d, iiLon, iiLat, iiAlt, rLon, rLat, rAlt), output_kind))

    do i = 1, nSpeciesTotal
      Tmp3d = NDensityS(0:nLons + 1, 0:nLats + 1, 0:nAlts + 1, i, iBlock)
      call c%put('['//trim(cSpecies(i))//']', real(inter3d(Tmp3d, iiLon, iiLat, iiAlt, rLon, rLat, rAlt), output_kind))
    enddo

    Tmp3d = Temperature(0:nLons + 1, 0:nLats + 1, 0:nAlts + 1, iBlock) * &
            TempUnit(0:nLons + 1, 0:nLats + 1, 0:nAlts + 1)
    call c%put('Temperature', real(inter3d(Tmp3d, iiLon, iiLat, iiAlt, rLon, rLat, rAlt), output_kind))

    do i = 1, 3
      Tmp3d = Velocity(0:nLons + 1, 0:nLats + 1, 0:nAlts + 1, i, iBlock)
      if (i == 1) call c%put('Vn (east)', real(inter3d(Tmp3d, iiLon, iiLat, iiAlt, rLon, rLat, rAlt), output_kind))
      if (i == 2) call c%put('Vn (north)', real(inter3d(Tmp3d, iiLon, iiLat, iiAlt, rLon, rLat, rAlt), output_kind))
      if (i == 3) call c%put('Vn (up)', real(inter3d(Tmp3d, iiLon, iiLat, iiAlt, rLon, rLat, rAlt), output_kind))
    enddo

    do i = 1, nSpecies
      Tmp3d = VerticalVelocity(0:nLons + 1, 0:nLats + 1, 0:nAlts + 1, i, iBlock)
      call c%put('Vn (up,'//trim(cSpecies(i))//')', real(inter3d(Tmp3d, iiLon, iiLat, iiAlt, rLon, rLat, rAlt), output_kind))
    enddo

    do i = 1, nIons
      Tmp3d = IDensityS(0:nLons + 1, 0:nLats + 1, 0:nAlts + 1, i, iBlock)
      call c%put('['//trim(cIons(i))//']', real(inter3d(Tmp3d, iiLon, iiLat, iiAlt, rLon, rLat, rAlt), output_kind))
    enddo

    Tmp3d = eTemperature(0:nLons + 1, 0:nLats + 1, 0:nAlts + 1, iBlock)
    call c%put('eTemperature', real(inter3d(Tmp3d, iiLon, iiLat, iiAlt, rLon, rLat, rAlt), output_kind))
    Tmp3d = iTemperature(0:nLons + 1, 0:nLats + 1, 0:nAlts + 1, iBlock)
    call c%put('iTemperature', real(inter3d(Tmp3d, iiLon, iiLat, iiAlt, rLon, rLat, rAlt), output_kind))

    Tmp3d = IVelocity(0:nLons + 1, 0:nLats + 1, 0:nAlts + 1, iEast_, iBlock)
    call c%put('Vi (east)', real(inter3d(Tmp3d, iiLon, iiLat, iiAlt, rLon, rLat, rAlt), output_kind))
    Tmp3d = IVelocity(0:nLons + 1, 0:nLats + 1, 0:nAlts + 1, iNorth_, iBlock)
    call c%put('Vi (north)', real(inter3d(Tmp3d, iiLon, iiLat, iiAlt, rLon, rLat, rAlt), output_kind))
    Tmp3d = IVelocity(0:nLons + 1, 0:nLats + 1, 0:nAlts + 1, iUp_, iBlock)
    call c%put('Vi (up)', real(inter3d(Tmp3d, iiLon, iiLat, iiAlt, rLon, rLat, rAlt), output_kind))

    do i = 1, nSpecies
      Tmp3d = NDensityS(0:nLons + 1, 0:nLats + 1, 0:nAlts + 1, i, iBlock) / &
              NDensity(0:nLons + 1, 0:nLats + 1, 0:nAlts + 1, iBlock)
      call c%put(trim(cSpecies(i))//' Mixing Ratio', real(inter3d(Tmp3d, iiLon, iiLat, iiAlt, rLon, rLat, rAlt), output_kind))
    enddo

    call c%put('RadCooling', real(dt*RadCooling(1, 1, jAlt, iBlock)*TempUnit(1, 1, jAlt), output_kind))
    call c%put('EuvHeating', real(dt*EuvHeating(1, 1, jAlt, iBlock)*TempUnit(1, 1, jAlt), output_kind))
    call c%put('Conduction', real(Conduction(1, 1, jAlt)*TempUnit(1, 1, jAlt), output_kind))
    call c%put('Heat Balance Total', real((dt*EuvHeating(1, 1, jAlt, iBlock)*TempUnit(1, 1, jAlt) - &
                                          dt*RadCooling(1, 1, jAlt, iBlock)*TempUnit(1, 1, jAlt) + &
                                          Conduction(1, 1, jAlt)*TempUnit(1, 1, jAlt)), output_kind))
    call c%put('Heating Efficiency', real(HeatingEfficiency_CB(1, 1, jAlt, iBlock), output_kind))
  end subroutine fill_0dall

  ! ---------------------------------------------------------------------------
  ! define_schema_1dglo — register variables for 1DGLO.
  ! ---------------------------------------------------------------------------
  subroutine define_schema_1dglo(c)
    use ModGITM, only: nAlts
    type(OutputContainer), intent(inout) :: c

    c%cType    = '1DGLO'
    c%gridKind = GRID_MAG_2D

    call c%define_var('Longitude', units='degrees_east', shape3=[nAlts, 1, 1], is_axis=.true.)
    call c%define_var('Latitude', units='degrees_north', shape3=[nAlts, 1, 1], is_axis=.true.)
    call c%define_var('Altitude', units='m', shape3=[nAlts, 1, 1])
    call c%define_var('6300 A Emission', units='', shape3=[nAlts, 1, 1])
    call c%define_var('PhotoElectronUp', units='', shape3=[nAlts, 1, 1])
    call c%define_var('PhotoElectronDown', units='', shape3=[nAlts, 1, 1])
  end subroutine define_schema_1dglo

  ! ---------------------------------------------------------------------------
  ! fill_1dglo — put 1DGLO stub variables.
  ! ---------------------------------------------------------------------------
  subroutine fill_1dglo(c, iBlock)
    use ModGITM
    use ModConst, only: cRadToDeg
    type(OutputContainer), intent(inout) :: c
    integer, intent(in) :: iBlock
    real(output_kind) :: tmp_lon(nAlts), tmp_lat(nAlts), tmp_alt(nAlts), zero_arr(nAlts)
    integer :: iAlt

    call c%prepare(iBlock)
    if (.not. c%this_rank_writes) return

    tmp_lon = real(Longitude(1, iBlock) * cRadToDeg, output_kind)
    tmp_lat = real(Latitude(1, iBlock) * cRadToDeg, output_kind)
    call c%put('Longitude', tmp_lon)
    call c%put('Latitude',  tmp_lat)

    do iAlt = 1, nAlts
      tmp_alt(iAlt) = real(Altitude_GB(1, 1, iAlt, iBlock), output_kind)
    enddo
    call c%put('Altitude', tmp_alt)

    zero_arr = 0.0_output_kind
    call c%put('6300 A Emission', zero_arr)
    call c%put('PhotoElectronUp', zero_arr)
    call c%put('PhotoElectronDown', zero_arr)
  end subroutine fill_1dglo

  ! ---------------------------------------------------------------------------
  ! define_schema_1dthm — register variables for 1DTHM.
  ! ---------------------------------------------------------------------------
  subroutine define_schema_1dthm(c)
    use ModGITM, only: nAlts
    use ModPlanet, only: nSpeciesTotal, cSpecies
    type(OutputContainer), intent(inout) :: c
    integer :: i

    c%cType    = '1DTHM'
    c%gridKind = GRID_MAG_2D

    call c%define_var('Longitude', units='degrees_east', shape3=[nAlts, 1, 1], is_axis=.true.)
    call c%define_var('Latitude', units='degrees_north', shape3=[nAlts, 1, 1], is_axis=.true.)
    call c%define_var('Altitude', units='m', shape3=[nAlts, 1, 1])
    call c%define_var('EUV Heating (K/s)', units='K/s', shape3=[nAlts, 1, 1])
    call c%define_var('Conduction (K/s)', units='K/s', shape3=[nAlts, 1, 1])
    call c%define_var('Molecular Conduction (K/s)', units='K/s', shape3=[nAlts, 1, 1])
    call c%define_var('Eddy Conduction (K/s)', units='K/s', shape3=[nAlts, 1, 1])
    call c%define_var('Eddy Adiabatic Conduction (K/s)', units='K/s', shape3=[nAlts, 1, 1])
    call c%define_var('Chemical Heating (K/s)', units='K/s', shape3=[nAlts, 1, 1])
    call c%define_var('Joule Heating (K/s)', units='K/s', shape3=[nAlts, 1, 1])
    call c%define_var('NO Cooling (K/s)', units='K/s', shape3=[nAlts, 1, 1])
    call c%define_var('O Cooling (K/s)', units='K/s', shape3=[nAlts, 1, 1])
    call c%define_var('Total Abs EUV', units='W/m2', shape3=[nAlts, 1, 1])
    do i = 1, nSpeciesTotal
      call c%define_var('Production Rate '//trim(cSpecies(i)), units='/m3/s', shape3=[nAlts, 1, 1])
    enddo
    do i = 1, nSpeciesTotal
      call c%define_var('Loss Rate '//trim(cSpecies(i)), units='/m3/s', shape3=[nAlts, 1, 1])
    enddo
  end subroutine define_schema_1dthm

  ! ---------------------------------------------------------------------------
  ! fill_1dthm — put 1DTHM variables at point (1,1).
  ! ---------------------------------------------------------------------------
  subroutine fill_1dthm(c, iBlock)
    use ModGITM
    use ModConst, only: cRadToDeg
    use ModSources
    use ModEUV, only: EuvTotal
    use ModPlanet, only: nSpeciesTotal, cSpecies
    type(OutputContainer), intent(inout) :: c
    integer, intent(in) :: iBlock
    real(output_kind) :: tmp_lon(nAlts), tmp_lat(nAlts), tmp_data(nAlts)
    integer :: iAlt, iiAlt, i

    call c%prepare(iBlock)
    if (.not. c%this_rank_writes) return

    tmp_lon = real(Longitude(1, 1) * cRadToDeg, output_kind)
    tmp_lat = real(Latitude(1, 1) * cRadToDeg, output_kind)
    call c%put('Longitude', tmp_lon)
    call c%put('Latitude',  tmp_lat)

    do iAlt = 1, nAlts
      tmp_data(iAlt) = real(Altitude_GB(1, 1, iAlt, 1), output_kind)
    enddo
    call c%put('Altitude', tmp_data)

    do iAlt = 1, nAlts
      iiAlt = max(min(iAlt, nAlts), 1)
      tmp_data(iAlt) = real(EuvHeating(1, 1, iiAlt, 1)*dt*TempUnit(1, 1, iiAlt), output_kind)
    enddo
    call c%put('EUV Heating (K/s)', tmp_data)

    do iAlt = 1, nAlts
      iiAlt = max(min(iAlt, nAlts), 1)
      tmp_data(iAlt) = real(Conduction(1, 1, iiAlt)*TempUnit(1, 1, iiAlt), output_kind)
    enddo
    call c%put('Conduction (K/s)', tmp_data)

    do iAlt = 1, nAlts
      iiAlt = max(min(iAlt, nAlts), 1)
      tmp_data(iAlt) = real(MoleConduction(1, 1, iiAlt), output_kind)
    enddo
    call c%put('Molecular Conduction (K/s)', tmp_data)

    do iAlt = 1, nAlts
      iiAlt = max(min(iAlt, nAlts), 1)
      tmp_data(iAlt) = real(EddyCond(1, 1, iiAlt), output_kind)
    enddo
    call c%put('Eddy Conduction (K/s)', tmp_data)

    do iAlt = 1, nAlts
      iiAlt = max(min(iAlt, nAlts), 1)
      tmp_data(iAlt) = real(EddyCondAdia(1, 1, iiAlt), output_kind)
    enddo
    call c%put('Eddy Adiabatic Conduction (K/s)', tmp_data)

    do iAlt = 1, nAlts
      iiAlt = max(min(iAlt, nAlts), 1)
      tmp_data(iAlt) = real(ChemicalHeatingRate(1, 1, iiAlt)*TempUnit(1, 1, iiAlt), output_kind)
    enddo
    call c%put('Chemical Heating (K/s)', tmp_data)

    do iAlt = 1, nAlts
      iiAlt = max(min(iAlt, nAlts), 1)
      tmp_data(iAlt) = real(JouleHeating(1, 1, iiAlt)*dt*TempUnit(1, 1, iiAlt), output_kind)
    enddo
    call c%put('Joule Heating (K/s)', tmp_data)

    do iAlt = 1, nAlts
      iiAlt = max(min(iAlt, nAlts), 1)
      tmp_data(iAlt) = real(-RadCooling(1, 1, iiAlt, 1)*dt*TempUnit(1, 1, iiAlt), output_kind)
    enddo
    call c%put('NO Cooling (K/s)', tmp_data)

    do iAlt = 1, nAlts
      iiAlt = max(min(iAlt, nAlts), 1)
      tmp_data(iAlt) = real(-OCooling(1, 1, iiAlt)*dt*TempUnit(1, 1, iiAlt), output_kind)
    enddo
    call c%put('O Cooling (K/s)', tmp_data)

    do iAlt = 1, nAlts
      iiAlt = max(min(iAlt, nAlts), 1)
      tmp_data(iAlt) = real(EuvTotal(1, 1, iiAlt, 1)*dt, output_kind)
    enddo
    call c%put('Total Abs EUV', tmp_data)

    do i = 1, nSpeciesTotal
      do iAlt = 1, nAlts
        iiAlt = max(min(iAlt, nAlts), 1)
        tmp_data(iAlt) = real(NeutralSourcesTotal(iiAlt, i), output_kind)
      enddo
      call c%put('Production Rate '//trim(cSpecies(i)), tmp_data)
    enddo

    do i = 1, nSpeciesTotal
      do iAlt = 1, nAlts
        iiAlt = max(min(iAlt, nAlts), 1)
        tmp_data(iAlt) = real(NeutralLossesTotal(iiAlt, i), output_kind)
      enddo
      call c%put('Loss Rate '//trim(cSpecies(i)), tmp_data)
    enddo
  end subroutine fill_1dthm

  ! ---------------------------------------------------------------------------
  ! define_schema_1dchm — register variables for 1DCHM.
  ! ---------------------------------------------------------------------------
  subroutine define_schema_1dchm(c)
    use ModGITM, only: nAlts
    type(OutputContainer), intent(inout) :: c

    c%cType    = '1DCHM'
    c%gridKind = GRID_MAG_2D

    call c%define_var('Longitude', units='degrees_east', shape3=[nAlts, 1, 1], is_axis=.true.)
    call c%define_var('Latitude', units='degrees_north', shape3=[nAlts, 1, 1], is_axis=.true.)
    call c%define_var('Altitude', units='m', shape3=[nAlts, 1, 1])
    call c%define_var('N2+ + e', units='/m3/s', shape3=[nAlts, 1, 1])
    call c%define_var('O2+ + e', units='/m3/s', shape3=[nAlts, 1, 1])
    call c%define_var('N2+ + O', units='/m3/s', shape3=[nAlts, 1, 1])
    call c%define_var('NO+ + e', units='/m3/s', shape3=[nAlts, 1, 1])
    call c%define_var('N+ + O2', units='/m3/s', shape3=[nAlts, 1, 1])
    call c%define_var('NO + N', units='/m3/s', shape3=[nAlts, 1, 1])
    call c%define_var('O+ + O2', units='/m3/s', shape3=[nAlts, 1, 1])
    call c%define_var('N + O2', units='/m3/s', shape3=[nAlts, 1, 1])
    call c%define_var('O2+ + N', units='/m3/s', shape3=[nAlts, 1, 1])
    call c%define_var('O2+ + NO', units='/m3/s', shape3=[nAlts, 1, 1])
    call c%define_var('O2+ + N2', units='/m3/s', shape3=[nAlts, 1, 1])
    call c%define_var('N2+ + O2', units='/m3/s', shape3=[nAlts, 1, 1])
    call c%define_var('N+ + O', units='/m3/s', shape3=[nAlts, 1, 1])
    call c%define_var('O!+ + N2', units='/m3/s', shape3=[nAlts, 1, 1])
    call c%define_var('O(1D) + N2', units='/m3/s', shape3=[nAlts, 1, 1])
    call c%define_var('O(1D) + O2', units='/m3/s', shape3=[nAlts, 1, 1])
    call c%define_var('O(1D) + O', units='/m3/s', shape3=[nAlts, 1, 1])
    call c%define_var('O(1D) + e', units='/m3/s', shape3=[nAlts, 1, 1])
    call c%define_var('N(2D) + O2', units='/m3/s', shape3=[nAlts, 1, 1])
    call c%define_var('O+(2D)+e', units='/m3/s', shape3=[nAlts, 1, 1])
    call c%define_var('N(2D) + O', units='/m3/s', shape3=[nAlts, 1, 1])
    call c%define_var('N(2D) + e', units='/m3/s', shape3=[nAlts, 1, 1])
    call c%define_var('O+(2D + N2', units='/m3/s', shape3=[nAlts, 1, 1])
    call c%define_var('O+(2P) + e', units='/m3/s', shape3=[nAlts, 1, 1])
    call c%define_var('O+(2P) + O', units='/m3/s', shape3=[nAlts, 1, 1])
    call c%define_var('O+(2P) + N2', units='/m3/s', shape3=[nAlts, 1, 1])
    call c%define_var('Chemical Heating Rate', units='K/s', shape3=[nAlts, 1, 1])
  end subroutine define_schema_1dchm

  ! ---------------------------------------------------------------------------
  ! fill_1dchm — put 1DCHM variables at point (1,1).
  ! ---------------------------------------------------------------------------
  subroutine fill_1dchm(c, iBlock)
    use ModGITM
    use ModConst, only: cRadToDeg
    use ModSources
    use ModConstants, only: Element_Charge
    type(OutputContainer), intent(inout) :: c
    integer, intent(in) :: iBlock
    real(output_kind) :: tmp_lon(nAlts), tmp_lat(nAlts), tmp_data(nAlts)
    integer :: iAlt, iiAlt, i

    call c%prepare(iBlock)
    if (.not. c%this_rank_writes) return

    tmp_lon = real(Longitude(1, 1) * cRadToDeg, output_kind)
    tmp_lat = real(Latitude(1, 1) * cRadToDeg, output_kind)
    call c%put('Longitude', tmp_lon)
    call c%put('Latitude',  tmp_lat)

    do iAlt = 1, nAlts
      tmp_data(iAlt) = real(Altitude_GB(1, 1, iAlt, 1), output_kind)
    enddo
    call c%put('Altitude', tmp_data)

    do i = 1, nReactions
      do iAlt = 1, nAlts
        iiAlt = max(min(iAlt, nAlts), 1)
        tmp_data(iAlt) = real(ChemicalHeatingSpecies(1, 1, iiAlt, i) / Element_Charge, output_kind)
      enddo
      if (i == 1) call c%put('N2+ + e', tmp_data)
      if (i == 2) call c%put('O2+ + e', tmp_data)
      if (i == 3) call c%put('N2+ + O', tmp_data)
      if (i == 4) call c%put('NO+ + e', tmp_data)
      if (i == 5) call c%put('N+ + O2', tmp_data)
      if (i == 6) call c%put('NO + N', tmp_data)
      if (i == 7) call c%put('O+ + O2', tmp_data)
      if (i == 8) call c%put('N + O2', tmp_data)
      if (i == 9) call c%put('O2+ + N', tmp_data)
      if (i == 10) call c%put('O2+ + NO', tmp_data)
      if (i == 11) call c%put('O2+ + N2', tmp_data)
      if (i == 12) call c%put('N2+ + O2', tmp_data)
      if (i == 13) call c%put('N+ + O', tmp_data)
      if (i == 14) call c%put('O!+ + N2', tmp_data)
      if (i == 15) call c%put('O(1D) + N2', tmp_data)
      if (i == 16) call c%put('O(1D) + O2', tmp_data)
      if (i == 17) call c%put('O(1D) + O', tmp_data)
      if (i == 18) call c%put('O(1D) + e', tmp_data)
      if (i == 19) call c%put('N(2D) + O2', tmp_data)
      if (i == 20) call c%put('O+(2D)+e', tmp_data)
      if (i == 21) call c%put('N(2D) + O', tmp_data)
      if (i == 22) call c%put('N(2D) + e', tmp_data)
      if (i == 23) call c%put('O+(2D + N2', tmp_data)
      if (i == 24) call c%put('O+(2P) + e', tmp_data)
      if (i == 25) call c%put('O+(2P) + O', tmp_data)
      if (i == 26) call c%put('O+(2P) + N2', tmp_data)
    enddo

    do iAlt = 1, nAlts
      iiAlt = max(min(iAlt, nAlts), 1)
      tmp_data(iAlt) = real(ChemicalHeatingRate(1, 1, iiAlt) * &
                            cp(1, 1, iiAlt, 1) * &
                            Rho(1, 1, iiAlt, 1) * &
                            TempUnit(1, 1, iiAlt) / Element_Charge, output_kind)
    enddo
    call c%put('Chemical Heating Rate', tmp_data)
  end subroutine fill_1dchm

  ! ---------------------------------------------------------------------------
  ! define_schema_1dnew — register variables for 1DNEW.
  ! ---------------------------------------------------------------------------
  subroutine define_schema_1dnew(c)
    use ModGITM, only: nAlts
    use ModPlanet, only: nSpeciesTotal, nSpecies, nIons, cSpecies, cIons
    type(OutputContainer), intent(inout) :: c
    integer :: i

    c%cType    = '1DNEW'
    c%gridKind = GRID_GEO_2D

    call c%define_var('Longitude', units='degrees_east', shape3=[nAlts, 1, 1], is_axis=.true.)
    call c%define_var('Local Time', units='hr', shape3=[nAlts, 1, 1])
    call c%define_var('Latitude', units='degrees_north', shape3=[nAlts, 1, 1])
    call c%define_var('Solar Zenith Angle', units='rad', shape3=[nAlts, 1, 1])
    call c%define_var('Altitude', units='m', shape3=[nAlts, 1, 1])
    call c%define_var('Rho', units='kg/m3', shape3=[nAlts, 1, 1])
    do i = 1, nSpeciesTotal
      call c%define_var('['//trim(cSpecies(i))//']', units='/m3', shape3=[nAlts, 1, 1])
    enddo
    call c%define_var('Temperature', units='K', shape3=[nAlts, 1, 1])
    call c%define_var('Vn (east)', units='m/s', shape3=[nAlts, 1, 1])
    call c%define_var('Vn (north)', units='m/s', shape3=[nAlts, 1, 1])
    call c%define_var('Vn (up)', units='m/s', shape3=[nAlts, 1, 1])
    do i = 1, nSpecies
      call c%define_var('Vn (up,'//trim(cSpecies(i))//')', units='m/s', shape3=[nAlts, 1, 1])
    enddo
    do i = 1, nIons
      call c%define_var('['//trim(cIons(i))//']', units='/m3', shape3=[nAlts, 1, 1])
    enddo
    call c%define_var('eTemperature', units='K', shape3=[nAlts, 1, 1])
    call c%define_var('iTemperature', units='K', shape3=[nAlts, 1, 1])
    call c%define_var('Vi (east)', units='m/s', shape3=[nAlts, 1, 1])
    call c%define_var('Vi (north)', units='m/s', shape3=[nAlts, 1, 1])
    call c%define_var('Vi (up)', units='m/s', shape3=[nAlts, 1, 1])
    do i = 1, nSpecies
      call c%define_var(trim(cSpecies(i))//' Mixing Ratio', units='', shape3=[nAlts, 1, 1])
    enddo
  end subroutine define_schema_1dnew

  ! ---------------------------------------------------------------------------
  ! fill_1dnew — put 1DNEW interpolated vertical column variables.
  ! ---------------------------------------------------------------------------
  subroutine fill_1dnew(c, iBlock, iiLon, iiLat, rLon, rLat)
    use ModGITM
    use ModConst, only: cRadToDeg
    use ModEUV, only: Sza
    use ModPlanet, only: nSpeciesTotal, nSpecies, nIons, cSpecies, cIons
    type(OutputContainer), intent(inout) :: c
    integer, intent(in) :: iBlock, iiLon, iiLat
    real, intent(in) :: rLon, rLat
    real(output_kind) :: tmp_lon(nAlts), tmp_lat(nAlts), tmp_lt(nAlts), tmp_sza(nAlts), tmp_data(nAlts)
    real :: slice2d(0:nLons + 1, 0:nLats + 1)
    integer :: iAlt, i

    call c%prepare(iBlock)
    if (.not. c%this_rank_writes) return

    tmp_lon = real((rLon * Longitude(iiLon, iBlock) + (1.0 - rLon) * Longitude(iiLon + 1, iBlock)) * cRadToDeg, output_kind)
    tmp_lat = real((rLat * Latitude(iiLat, iBlock) + (1.0 - rLat) * Latitude(iiLat + 1, iBlock)) * cRadToDeg, output_kind)
    tmp_lt  = real(LocalTime(iiLon), output_kind)
    tmp_sza = real(rLon*rLat*Sza(iiLon, iiLat, iBlock) + &
                   (1.0 - rLon)*rLat*Sza(iiLon + 1, iiLat, iBlock) + &
                   rLon*(1.0 - rLat)*Sza(iiLon, iiLat + 1, iBlock) + &
                   (1.0 - rLon)*(1.0 - rLat)*Sza(iiLon + 1, iiLat + 1, iBlock), output_kind)

    call c%put('Longitude', tmp_lon)
    call c%put('Local Time', tmp_lt)
    call c%put('Latitude',  tmp_lat)
    call c%put('Solar Zenith Angle', tmp_sza)

    do iAlt = 1, nAlts
      tmp_data(iAlt) = real(Altitude_GB(iiLon, iiLat, iAlt, iBlock), output_kind)
    enddo
    call c%put('Altitude', tmp_data)

    do iAlt = 1, nAlts
      slice2d = Rho(0:nLons + 1, 0:nLats + 1, iAlt, iBlock)
      tmp_data(iAlt) = real(inter2d(slice2d, iiLon, iiLat, rLon, rLat), output_kind)
    enddo
    call c%put('Rho', tmp_data)

    do i = 1, nSpeciesTotal
      do iAlt = 1, nAlts
        slice2d = NDensityS(0:nLons + 1, 0:nLats + 1, iAlt, i, iBlock)
        tmp_data(iAlt) = real(inter2d(slice2d, iiLon, iiLat, rLon, rLat), output_kind)
      enddo
      call c%put('['//trim(cSpecies(i))//']', tmp_data)
    enddo

    do iAlt = 1, nAlts
      slice2d = Temperature(0:nLons + 1, 0:nLats + 1, iAlt, iBlock) * &
                TempUnit(0:nLons + 1, 0:nLats + 1, iAlt)
      tmp_data(iAlt) = real(inter2d(slice2d, iiLon, iiLat, rLon, rLat), output_kind)
    enddo
    call c%put('Temperature', tmp_data)

    do i = 1, 3
      do iAlt = 1, nAlts
        slice2d = Velocity(0:nLons + 1, 0:nLats + 1, iAlt, i, iBlock)
        tmp_data(iAlt) = real(inter2d(slice2d, iiLon, iiLat, rLon, rLat), output_kind)
      enddo
      if (i == 1) call c%put('Vn (east)', tmp_data)
      if (i == 2) call c%put('Vn (north)', tmp_data)
      if (i == 3) call c%put('Vn (up)', tmp_data)
    enddo

    do i = 1, nSpecies
      do iAlt = 1, nAlts
        slice2d = VerticalVelocity(0:nLons + 1, 0:nLats + 1, iAlt, i, iBlock)
        tmp_data(iAlt) = real(inter2d(slice2d, iiLon, iiLat, rLon, rLat), output_kind)
      enddo
      call c%put('Vn (up,'//trim(cSpecies(i))//')', tmp_data)
    enddo

    do i = 1, nIons
      do iAlt = 1, nAlts
        slice2d = IDensityS(0:nLons + 1, 0:nLats + 1, iAlt, i, iBlock)
        tmp_data(iAlt) = real(inter2d(slice2d, iiLon, iiLat, rLon, rLat), output_kind)
      enddo
      call c%put('['//trim(cIons(i))//']', tmp_data)
    enddo

    do iAlt = 1, nAlts
      slice2d = eTemperature(0:nLons + 1, 0:nLats + 1, iAlt, iBlock)
      tmp_data(iAlt) = real(inter2d(slice2d, iiLon, iiLat, rLon, rLat), output_kind)
    enddo
    call c%put('eTemperature', tmp_data)

    do iAlt = 1, nAlts
      slice2d = ITemperature(0:nLons + 1, 0:nLats + 1, iAlt, iBlock)
      tmp_data(iAlt) = real(inter2d(slice2d, iiLon, iiLat, rLon, rLat), output_kind)
    enddo
    call c%put('iTemperature', tmp_data)

    do i = 1, 3
      do iAlt = 1, nAlts
        slice2d = IVelocity(0:nLons + 1, 0:nLats + 1, iAlt, i, iBlock)
        tmp_data(iAlt) = real(inter2d(slice2d, iiLon, iiLat, rLon, rLat), output_kind)
      enddo
      if (i == 1) call c%put('Vi (east)', tmp_data)
      if (i == 2) call c%put('Vi (north)', tmp_data)
      if (i == 3) call c%put('Vi (up)', tmp_data)
    enddo

    do i = 1, nSpecies
      do iAlt = 1, nAlts
        slice2d = NDensityS(0:nLons + 1, 0:nLats + 1, iAlt, i, iBlock)/ &
                  NDensity(0:nLons + 1, 0:nLats + 1, iAlt, iBlock)
        tmp_data(iAlt) = real(inter2d(slice2d, iiLon, iiLat, rLon, rLat), output_kind)
      enddo
      call c%put(trim(cSpecies(i))//' Mixing Ratio', tmp_data)
    enddo
  end subroutine fill_1dnew

  ! ---------------------------------------------------------------------------
  ! define_schema_3dusr — define variables for the 3DUSR container.
  ! ---------------------------------------------------------------------------
  subroutine define_schema_3dusr(c)
    use ModSizeGitm, only: nLons, nLats, nAlts
    type(OutputContainer), intent(inout) :: c
    integer :: i

    c%cType       = '3DUSR'
    c%gridKind    = GRID_GEO_3D
    c%nGhostCells = 0

    call c%define_var('Longitude', units='degrees_east',  shape3=[nLons, 1, 1],           is_axis=.true., longName='Geographic longitude')
    call c%define_var('Latitude',  units='degrees_north', shape3=[nLats, 1, 1],           is_axis=.true., longName='Geographic latitude')
    call c%define_var('Altitude',  units='m',             shape3=[nLons, nLats, nAlts],   longName='Altitude')
    do i = 1, nUsrVars3D
      call c%define_var(trim(UsrVars3D(i)%name), units=trim(UsrVars3D(i)%units), &
                        shape3=[nLons, nLats, nAlts], longName=trim(UsrVars3D(i)%longName))
    end do
  end subroutine define_schema_3dusr

  ! ---------------------------------------------------------------------------
  ! fill_3dusr — put 3DUSR variables into the container.
  ! ---------------------------------------------------------------------------
  subroutine fill_3dusr(c, iBlock)
    use ModGITM
    use ModConst, only: cRadToDeg
    use ModUserGITM, only: UserData3D
    type(OutputContainer), intent(inout) :: c
    integer, intent(in) :: iBlock
    integer :: i

    call c%prepare(iBlock)
    if (.not. c%this_rank_writes) return

    call c%put('Longitude', real(Longitude(1:nLons, iBlock) * cRadToDeg, output_kind))
    call c%put('Latitude',  real(Latitude(1:nLats, iBlock) * cRadToDeg, output_kind))
    call c%put('Altitude',  real(Altitude_GB(1:nLons, 1:nLats, 1:nAlts, iBlock), output_kind))
    do i = 1, nUsrVars3D
      call c%put(trim(UsrVars3D(i)%name), real(UserData3D(1:nLons, 1:nLats, 1:nAlts, i, iBlock), output_kind))
    end do
  end subroutine fill_3dusr

  ! ---------------------------------------------------------------------------
  ! define_schema_2dusr — define variables for the 2DUSR container.
  ! ---------------------------------------------------------------------------
  subroutine define_schema_2dusr(c)
    use ModGITM, only: nLons, nLats
    type(OutputContainer), intent(inout) :: c
    integer :: i

    c%cType       = '2DUSR'
    c%gridKind    = GRID_GEO_2D
    c%nGhostCells = 0

    call c%define_var('Longitude', units='degrees_east',  shape3=[nLons, 1, 1],    is_axis=.true., longName='Geographic longitude')
    call c%define_var('Latitude',  units='degrees_north', shape3=[nLats, 1, 1],    is_axis=.true., longName='Geographic latitude')
    call c%define_var('Altitude',  units='m',             shape3=[nLons, nLats, 1], longName='Altitude')
    do i = 1, nUsrVars2D
      call c%define_var(trim(UsrVars2D(i)%name), units=trim(UsrVars2D(i)%units), &
                        shape3=[nLons, nLats, 1], longName=trim(UsrVars2D(i)%longName))
    end do
  end subroutine define_schema_2dusr

  ! ---------------------------------------------------------------------------
  ! fill_2dusr — put 2DUSR variables into the container.
  ! ---------------------------------------------------------------------------
  subroutine fill_2dusr(c, iBlock)
    use ModGITM
    use ModConst, only: cRadToDeg
    use ModUserGITM, only: UserData2D
    type(OutputContainer), intent(inout) :: c
    integer, intent(in) :: iBlock
    integer :: i

    call c%prepare(iBlock)
    if (.not. c%this_rank_writes) return

    call c%put('Longitude', real(Longitude(1:nLons, iBlock) * cRadToDeg, output_kind))
    call c%put('Latitude',  real(Latitude(1:nLats, iBlock) * cRadToDeg, output_kind))
    call c%put('Altitude',  real(Altitude_GB(1:nLons, 1:nLats, 1, iBlock), output_kind))
    do i = 1, nUsrVars2D
      call c%put(trim(UsrVars2D(i)%name), real(UserData2D(1:nLons, 1:nLats, 1, i, iBlock), output_kind))
    end do
  end subroutine fill_2dusr

  ! ---------------------------------------------------------------------------
  ! define_schema_1dusr — define variables for the 1DUSR container.
  ! ---------------------------------------------------------------------------
  subroutine define_schema_1dusr(c)
    use ModGITM, only: nAlts
    type(OutputContainer), intent(inout) :: c
    integer :: i

    c%cType       = '1DUSR'
    c%gridKind    = GRID_GEO_2D
    c%nGhostCells = 0

    call c%define_var('Longitude', units='degrees_east',  shape3=[nAlts, 1, 1], is_axis=.true.)
    call c%define_var('Latitude',  units='degrees_north', shape3=[nAlts, 1, 1], is_axis=.true.)
    call c%define_var('Altitude',  units='m',             shape3=[nAlts, 1, 1])
    do i = 1, nUsrVars1D
      call c%define_var(trim(UsrVars1D(i)%name), units=trim(UsrVars1D(i)%units), &
                        shape3=[nAlts, 1, 1], longName=trim(UsrVars1D(i)%longName))
    end do
  end subroutine define_schema_1dusr

  ! ---------------------------------------------------------------------------
  ! fill_1dusr — put 1DUSR variables into the container.
  ! ---------------------------------------------------------------------------
  subroutine fill_1dusr(c, iBlock, iiLon, iiLat, rLon, rLat)
    use ModGITM
    use ModConst, only: cRadToDeg
    use ModUserGITM, only: UserData1D
    type(OutputContainer), intent(inout) :: c
    integer, intent(in) :: iBlock, iiLon, iiLat
    real, intent(in) :: rLon, rLat
    real(output_kind) :: tmp(nAlts)
    integer :: i, iAlt

    call c%prepare(iBlock)
    if (.not. c%this_rank_writes) return

    tmp = real((rLon * Longitude(iiLon, iBlock) + (1.0 - rLon) * Longitude(iiLon+1, iBlock)) * cRadToDeg, output_kind)
    call c%put('Longitude', tmp)
    tmp = real((rLat * Latitude(iiLat, iBlock) + (1.0 - rLat) * Latitude(iiLat+1, iBlock)) * cRadToDeg, output_kind)
    call c%put('Latitude', tmp)
    do iAlt = 1, nAlts
      tmp(iAlt) = real(Altitude_GB(iiLon, iiLat, iAlt, iBlock), output_kind)
    end do
    call c%put('Altitude', tmp)
    do i = 1, nUsrVars1D
      do iAlt = 1, nAlts
        tmp(iAlt) = real(UserData1D(iiLon, iiLat, iAlt, i), output_kind)
      end do
      call c%put(trim(UsrVars1D(i)%name), tmp)
    end do
  end subroutine fill_1dusr

  ! ---------------------------------------------------------------------------
  ! define_schema_0dusr — define variables for the 0DUSR container.
  ! ---------------------------------------------------------------------------
  subroutine define_schema_0dusr(c)
    type(OutputContainer), intent(inout) :: c
    integer :: i

    c%cType       = '0DUSR'
    c%gridKind    = GRID_GEO_2D
    c%nGhostCells = 0

    call c%define_var('Longitude', units='degrees_east',  shape3=[1, 1, 1], is_axis=.true.)
    call c%define_var('Latitude',  units='degrees_north', shape3=[1, 1, 1], is_axis=.true.)
    call c%define_var('Altitude',  units='m',             shape3=[1, 1, 1])
    do i = 1, nUsrVars0D
      call c%define_var(trim(UsrVars0D(i)%name), units=trim(UsrVars0D(i)%units), &
                        shape3=[1, 1, 1], longName=trim(UsrVars0D(i)%longName))
    end do
  end subroutine define_schema_0dusr

  ! ---------------------------------------------------------------------------
  ! fill_0dusr — put 0DUSR variables into the container.
  ! ---------------------------------------------------------------------------
  subroutine fill_0dusr(c, iBlock, iiLon, iiLat, iiAlt, rLon, rLat, rAlt)
    use ModGITM
    use ModConst, only: cRadToDeg
    use ModUserGITM, only: UserData1D
    type(OutputContainer), intent(inout) :: c
    integer, intent(in) :: iBlock, iiLon, iiLat, iiAlt
    real, intent(in) :: rLon, rLat, rAlt
    integer :: i

    call c%prepare(iBlock)
    if (.not. c%this_rank_writes) return

    call c%put('Longitude', real((rLon * Longitude(iiLon, iBlock) + (1.0 - rLon) * Longitude(iiLon+1, iBlock)) * cRadToDeg, output_kind))
    call c%put('Latitude',  real((rLat * Latitude(iiLat, iBlock)  + (1.0 - rLat) * Latitude(iiLat+1, iBlock))  * cRadToDeg, output_kind))
    call c%put('Altitude',  real(rAlt * Altitude_GB(iiLon, iiLat, iiAlt, iBlock) + &
                                 (1.0 - rAlt) * Altitude_GB(iiLon+1, iiLat+1, iiAlt+1, iBlock), output_kind))
    do i = 1, nUsrVars0D
      call c%put(trim(UsrVars0D(i)%name), real(UserData1D(iiLon, iiLat, iiAlt, i), output_kind))
    end do
  end subroutine fill_0dusr

  ! ---------------------------------------------------------------------------
  ! inter2d — Bilinear interpolation on a 2D grid slice.
  ! ---------------------------------------------------------------------------
  real function inter2d(variable, iiLon, iiLat, rLon, rLat)
    real, intent(in) :: variable(:, :), rLon, rLat
    integer, intent(in) :: iiLon, iiLat

    inter2d = &
      rLon*rLat*variable(iiLon, iiLat) + &
      (1.0 - rLon)*rLat*variable(iiLon + 1, iiLat) + &
      rLon*(1.0 - rLat)*variable(iiLon, iiLat + 1) + &
      (1.0 - rLon)*(1.0 - rLat)*variable(iiLon + 1, iiLat + 1)
  end function inter2d

  ! ---------------------------------------------------------------------------
  ! inter3d — Trilinear interpolation on a 3D grid.
  ! ---------------------------------------------------------------------------
  real function inter3d(variable, iiLon, iiLat, iiAlt, rLon, rLat, rAlt)
    real, intent(in) :: variable(:, :, :)
    real, intent(in) :: rLon, rLat, rAlt
    integer, intent(in) :: iiLon, iiLat, iiAlt

    inter3d = &
      rLon*rLat*rAlt*variable(iiLon, iiLat, iiAlt) + &
      (1.0 - rLon)*rLat*rAlt*variable(iiLon + 1, iiLat, iiAlt) + &
      rLon*(1.0 - rLat)*rAlt*variable(iiLon, iiLat + 1, iiAlt) + &
      (1.0 - rLon)*(1.0 - rLat)*rAlt*variable(iiLon + 1, iiLat + 1, iiAlt) + &
      rLon*rLat*(1.0 - rAlt)*variable(iiLon, iiLat, iiAlt + 1) + &
      (1.0 - rLon)*rLat*(1.0 - rAlt)*variable(iiLon + 1, iiLat, iiAlt + 1) + &
      rLon*(1.0 - rLat)*(1.0 - rAlt)*variable(iiLon, iiLat + 1, iiAlt + 1) + &
      (1.0 - rLon)*(1.0 - rLat)*(1.0 - rAlt)*variable(iiLon + 1, iiLat + 1, iiAlt + 1)
  end function inter3d

end module ModOutputProducers
