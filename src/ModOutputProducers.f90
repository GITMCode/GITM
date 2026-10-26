! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE
!
! ModOutputProducers.f90
!
! Container-based output producers.  Each output type has two procedures:
!   define_schema_<type>(c)   register variable names/units/shapes once
!   fill_<type>(c, iBlock)    put data into the container each timestep
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
!   stripped).  Magnetic, 1D, and 0D types have no ghost concept, so leave the
!   default 0.

module ModOutputProducers

  use ModOutputContainer
  use ModUserGITM, only: nUserOutputs

  implicit none

  type(OutputContainer), allocatable, target :: containers(:)

  ! Canonical type list: the single source for container identity and order.
  ! Container i holds type cContainerTypes(i); get_container_idx and
  ! is_known_output_type derive from this array (no container state involved,
  ! so both are safe before init, as bad_outputtype calls them at input parse time),
  ! and assert_types_fully_wired verifies every container landed in its slot.
  character(len=5), parameter :: cContainerTypes(26) = [ &
    '2DMEL', '2DTEC', '2DGEL', '3DNEU', '3DALL', '3DION', '3DTHM', '3DCHM', &
    '3DLST', '3DGLO', '3DMAG', '3DHME', '3DMOH', '3DMOV', '2DANC', '2DHME', &
    '1DALL', '0DALL', '1DGLO', '1DTHM', '1DCHM', '1DNEW', &
    '3DUSR', '2DUSR', '1DUSR', '0DUSR']
  integer, parameter :: nContainers = size(cContainerTypes)

  logical :: containers_initialized = .false.

  ! Probe location for the 1D/0D point-output fills (Route A1). Set by the
  ! output driver (output_common) before dispatch; read by the probe fills.
  ! Module state instead of fill arguments so every fill shares the uniform
  ! (c, iBlock) signature. Defaults are the dummy-index priming values: safe
  ! on every rank, and what init-time priming of the probe types uses.
  integer :: iProbeLon = 1, iProbeLat = 1, iProbeAlt = 1
  real :: rProbeLon = 0.5, rProbeLat = 0.5, rProbeAlt = 0.5

  ! ---------------------------------------------------------------------------
  ! USR variable registry: populated by register_usr_var() from user.f90 before
  ! the first output timestep.  Coord vars (Longitude/Latitude/Altitude) are
  ! always defined unconditionally in define_schema_*dusr and are NOT listed here.
  !
  ! Storage strategy: allocatable, grows on demand (doubling) starting from
  ! UsrInitCap.  The registry itself is uncapped; the real ceiling is
  ! nUserOutputs (from ModUserGITM), which sizes the actual UserData*D arrays.
  ! Decoupling them means raising nUserOutputs is the only knob users ever touch
  ! to support more USR variables.
  ! ---------------------------------------------------------------------------
  integer, parameter :: UsrInitCap = 16

  type :: UsrVarEntry
    character(len=60) :: name = ''
    character(len=40) :: units = ''
    character(len=60) :: longName = ''
  end type UsrVarEntry

  type(UsrVarEntry), allocatable :: UsrVars3D(:)
  integer                        :: nUsrVars3D = 0
  type(UsrVarEntry), allocatable :: UsrVars2D(:)
  integer                        :: nUsrVars2D = 0
  type(UsrVarEntry), allocatable :: UsrVars1D(:)
  integer                        :: nUsrVars1D = 0
  type(UsrVarEntry), allocatable :: UsrVars0D(:)
  integer                        :: nUsrVars0D = 0

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
  ! is_known_output_type: return .true. if cType is a registered output code.
  ! Pure string check, no container state required (safe to call before init).
  ! ---------------------------------------------------------------------------
  logical function is_known_output_type(cType)
    character(len=5), intent(in) :: cType
    is_known_output_type = get_container_idx(cType) > 0
  end function is_known_output_type

  ! ---------------------------------------------------------------------------
  ! get_container_idx: return the containers(:) index for cType, or -1.
  ! Position in cContainerTypes; no container state required (safe before init).
  ! ---------------------------------------------------------------------------
  integer function get_container_idx(cType)
    character(len=5), intent(in) :: cType
    integer :: i
    do i = 1, nContainers
      if (cType == cContainerTypes(i)) then
        get_container_idx = i
        return
      endif
    enddo
    get_container_idx = -1
  end function get_container_idx

  ! ---------------------------------------------------------------------------
  ! register_usr_var: add one user-defined variable to the appropriate USR list.
  ! Called from init_usr_output_registry in user.f90 before the first output.
  ! dims: 3, 2, 1, or 0 (spatial dimensionality of the output type).
  ! Coord vars (Longitude, Latitude, Altitude) must NOT be registered here;
  ! they are defined unconditionally in each define_schema_*dusr.
  ! ---------------------------------------------------------------------------
  subroutine register_usr_var(dims, name, units, longName)
    integer, intent(in)          :: dims
    character(len=*), intent(in) :: name, units
    character(len=*), intent(in), optional :: longName
    character(len=60) :: lName

    lName = name
    if (present(longName)) lName = longName

    ! The cap check guards the UserData*D data arrays (sized nUserOutputs).
    ! The registry itself is uncapped; grow_usr_registry handles allocation.
    select case (dims)
    case (3)
      nUsrVars3D = nUsrVars3D + 1
      if (nUsrVars3D > nUserOutputs) &
        call stop_gitm("register_usr_var: too many 3DUSR vars; increase nUserOutputs in ModUser.f90")
      call grow_usr_registry(UsrVars3D, nUsrVars3D)
      UsrVars3D(nUsrVars3D) = UsrVarEntry(name, units, lName)
    case (2)
      nUsrVars2D = nUsrVars2D + 1
      if (nUsrVars2D > nUserOutputs) call &
        stop_gitm("register_usr_var: too many 2DUSR vars; increase nUserOutputs in ModUser.f90")
      call grow_usr_registry(UsrVars2D, nUsrVars2D)
      UsrVars2D(nUsrVars2D) = UsrVarEntry(name, units, lName)
    case (1)
      nUsrVars1D = nUsrVars1D + 1
      if (nUsrVars1D > nUserOutputs) &
        call stop_gitm("register_usr_var: too many 1DUSR vars; increase nUserOutputs in ModUser.f90")
      call grow_usr_registry(UsrVars1D, nUsrVars1D)
      UsrVars1D(nUsrVars1D) = UsrVarEntry(name, units, lName)
    case (0)
      nUsrVars0D = nUsrVars0D + 1
      if (nUsrVars0D > nUserOutputs) &
        call stop_gitm("register_usr_var: too many 0DUSR vars; increase nUserOutputs in ModUser.f90")
      call grow_usr_registry(UsrVars0D, nUsrVars0D)
      UsrVars0D(nUsrVars0D) = UsrVarEntry(name, units, lName)
    end select
  end subroutine register_usr_var

  ! ---------------------------------------------------------------------------
  ! grow_usr_registry: ensure arr is allocated and has capacity >= n.
  ! Doubling growth from UsrInitCap.  Registration is one-time at startup so
  ! the copy cost is negligible (tiny string records, few growth events).
  ! ---------------------------------------------------------------------------
  subroutine grow_usr_registry(arr, n)
    type(UsrVarEntry), allocatable, intent(inout) :: arr(:)
    integer, intent(in) :: n
    type(UsrVarEntry), allocatable :: tmp(:)
    integer :: cap

    if (.not. allocated(arr)) then
      allocate(arr(max(UsrInitCap, n)))
      return
    endif
    if (size(arr) >= n) return
    cap = size(arr)
    do while (cap < n)
      cap = cap*2
    enddo
    allocate(tmp(cap))
    tmp(1:size(arr)) = arr
    call move_alloc(tmp, arr)
  end subroutine grow_usr_registry

  ! ---------------------------------------------------------------------------
  ! find_usr_var_idx: return 1-based index of name in UsrVars*D, or -1.
  ! ---------------------------------------------------------------------------
  integer function find_usr_var_idx(dim, name) result(idx)
    integer, intent(in)          :: dim
    character(len=*), intent(in) :: name
    integer :: i
    idx = -1
    select case (dim)
    case (3)
      do i = 1, nUsrVars3D
        if (trim(UsrVars3D(i)%name) == trim(name)) then; idx = i; return; endif
      enddo
    case (2)
      do i = 1, nUsrVars2D
        if (trim(UsrVars2D(i)%name) == trim(name)) then; idx = i; return; endif
      enddo
    case (1)
      do i = 1, nUsrVars1D
        if (trim(UsrVars1D(i)%name) == trim(name)) then; idx = i; return; endif
      enddo
    case (0)
      do i = 1, nUsrVars0D
        if (trim(UsrVars0D(i)%name) == trim(name)) then; idx = i; return; endif
      enddo
    end select
  end function find_usr_var_idx

  ! ---------------------------------------------------------------------------
  ! add_usr_output_3d: store a 3D (lon×lat×alt) user-defined field.
  !
  ! Call from physics code each timestep to provide data for 3DUSR output.
  ! First call with a new name registers the variable; subsequent calls update
  ! the stored data.  All new-name calls must occur before the first output
  ! write; see init_usr_output_registry in user.f90.
  !
  ! Ghost cell detection (size(arr,1)):
  !   nLons     : interior array, stored as-is
  !   nLons+4   : ghost-inclusive array (2 ghosts each side), stripped to interior
  !   anything else (incl. magnetic grid sizes): stop_gitm
  !
  ! iBlock defaults to 1 when not supplied.
  ! ---------------------------------------------------------------------------
  subroutine add_usr_output_3d(arr, name, units, iBlock, longName)
    use ModSizeGitm, only: nLons, nLats, nAlts
    use ModUserGITM, only: UserData3D
    use ModElectrodynamics, only: nMagLons, nMagLats
    real, intent(in)                       :: arr(:, :, :)
    character(len=*), intent(in)           :: name, units
    integer, intent(in), optional          :: iBlock
    character(len=*), intent(in), optional :: longName
    integer :: iBlk, idx, n1

    iBlk = 1; if (present(iBlock)) iBlk = iBlock
    n1 = size(arr, 1)

    if (n1 == nMagLons .or. n1 == nMagLons + 1) &
      call stop_gitm("add_usr_output: magnetic grid arrays not supported ('"//trim(name)//"')")
    if (n1 /= nLons .and. n1 /= nLons + 4) &
      call stop_gitm("add_usr_output: unrecognized shape for '" &
                     //trim(name)//"'; pass (nLons,nLats,nAlts) or ghost-inclusive")

    idx = find_usr_var_idx(3, name)
    if (idx < 0) then
      call register_usr_var(3, name, units, longName)
      idx = find_usr_var_idx(3, name)
    endif

    if (n1 == nLons + 4) then
      UserData3D(1:nLons, 1:nLats, 1:nAlts, idx, iBlk) = arr(3:nLons + 2, 3:nLats + 2, 3:nAlts + 2)
    else
      UserData3D(1:nLons, 1:nLats, 1:nAlts, idx, iBlk) = arr
    endif
  end subroutine add_usr_output_3d

  ! ---------------------------------------------------------------------------
  ! add_usr_output_2d: store a 2D (lon×lat) user-defined field.
  ! Ghost detection: nLons -> interior; nLons+4 -> ghost-inclusive, stripped.
  ! ---------------------------------------------------------------------------
  subroutine add_usr_output_2d(arr, name, units, iBlock, longName)
    use ModSizeGitm, only: nLons, nLats
    use ModUserGITM, only: UserData2D
    use ModElectrodynamics, only: nMagLons, nMagLats
    real, intent(in)                       :: arr(:, :)
    character(len=*), intent(in)           :: name, units
    integer, intent(in), optional          :: iBlock
    character(len=*), intent(in), optional :: longName
    integer :: iBlk, idx, n1

    iBlk = 1; if (present(iBlock)) iBlk = iBlock
    n1 = size(arr, 1)

    if (n1 == nMagLons .or. n1 == nMagLons + 1) &
      call stop_gitm("add_usr_output: magnetic grid arrays not supported ('"//trim(name)//"')")
    if (n1 /= nLons .and. n1 /= nLons + 4) &
      call stop_gitm("add_usr_output: unrecognized shape for '" &
                     //trim(name)//"'; pass (nLons,nLats) or ghost-inclusive")

    idx = find_usr_var_idx(2, name)
    if (idx < 0) then
      call register_usr_var(2, name, units, longName)
      idx = find_usr_var_idx(2, name)
    endif

    if (n1 == nLons + 4) then
      UserData2D(1:nLons, 1:nLats, 1, idx, iBlk) = arr(3:nLons + 2, 3:nLats + 2)
    else
      UserData2D(1:nLons, 1:nLats, 1, idx, iBlk) = arr
    endif
  end subroutine add_usr_output_2d

  ! ---------------------------------------------------------------------------
  ! add_usr_output_1d: store a 1D (altitude profile) user-defined field.
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
    endif

    UserData1D(1, 1, 1:nAlts, idx) = arr
  end subroutine add_usr_output_1d

  ! ---------------------------------------------------------------------------
  ! init_output_containers: allocate and populate every container's schema.
  ! Called once (idempotent); both output_common and netcdf_open_file call it.
  !
  ! Two kinds of type coexist here and the backends cannot tell them apart:
  !
  !   * Unconverted types (mag grids, HME/HIME regionals, USR) declare their
  !     schema explicitly via define_schema_<type>.
  !
  !   * Converted one-shot types have no define side at all; their schema is
  !     materialized by PRIMING: call fill_<type>(c, 1) once so the first put of
  !     each name defines it, then lock the schema.  This is what lets the
  !     schema be DERIVED from the fill code while still being complete before
  !     open_file, where netcdf issues collective nfmpi_def_var for the whole
  !     file.  (Deferring definition to the block loop instead would make
  !     def_var rank-dependent, the same bug class as the 2DMEL write deadlock.)
  !
  ! Priming is safe because it runs at the first output, after the model is
  ! initialized: every array a converted fill touches has been allocated by
  ! then.  That allocation-timing precondition must hold for EVERY converted
  ! fill in EVERY run (priming runs them all, whether or not the type was
  ! requested) and is what keeps the lazily-allocated mag-grid types out.
  ! Priming's VALUES are thrown away (every real fill overwrites them), so only
  ! the names/shapes/metadata it establishes matter.
  ! ---------------------------------------------------------------------------
  subroutine init_output_containers()
    integer :: i
    if (containers_initialized) return
    allocate(containers(nContainers))

    ! ---- Unconverted types: schemas declared explicitly -------------------
    call define_schema_2dmel(containers(get_container_idx('2DMEL')))
    call define_schema_3dhme(containers(get_container_idx('3DHME')))
    call define_schema_2dhme(containers(get_container_idx('2DHME')))
    call define_schema_1dglo(containers(get_container_idx('1DGLO')))
    call define_schema_1dthm(containers(get_container_idx('1DTHM')))
    call define_schema_1dchm(containers(get_container_idx('1DCHM')))
    call define_schema_3dusr(containers(get_container_idx('3DUSR')))
    call define_schema_2dusr(containers(get_container_idx('2DUSR')))
    call define_schema_1dusr(containers(get_container_idx('1DUSR')))
    call define_schema_0dusr(containers(get_container_idx('0DUSR')))

    ! ---- Fill registration: one line per type -----------------------------
    ! Adding a type = write fill_<type> + one line here (+ registry rows for
    ! fixed names). The dispatch in output_common runs whatever is registered.
    call set_fill('3DALL', fill_3dall)
    call set_fill('2DANC', fill_2danc)
    call set_fill('3DNEU', fill_3dneu)
    call set_fill('3DION', fill_3dion)
    call set_fill('3DMAG', fill_3dmag)
    call set_fill('2DTEC', fill_2dtec)
    call set_fill('2DGEL', fill_2dgel)
    call set_fill('3DTHM', fill_3dthm)
    call set_fill('3DCHM', fill_3dchm)
    call set_fill('3DMOH', fill_3dmoh)
    call set_fill('3DMOV', fill_3dmov)
    call set_fill('3DLST', fill_3dlst)
    call set_fill('3DGLO', fill_3dglo)
    call set_fill('2DMEL', fill_2dmel)
    call set_fill('3DHME', fill_3dhme)
    call set_fill('2DHME', fill_2dhme)
    call set_fill('1DGLO', fill_1dglo)
    call set_fill('1DTHM', fill_1dthm)
    call set_fill('1DCHM', fill_1dchm)
    call set_fill('3DUSR', fill_3dusr)
    call set_fill('2DUSR', fill_2dusr)
    ! Probe fills read the iProbe*/rProbe* module state (set by the driver
    ! before dispatch) instead of taking interpolation arguments.
    call set_fill('1DALL', fill_1dall)
    call set_fill('0DALL', fill_0dall)
    call set_fill('1DNEW', fill_1dnew)
    call set_fill('1DUSR', fill_1dusr)
    call set_fill('0DUSR', fill_0dusr)

    ! Set before priming: the flag's job is idempotency, and priming runs real
    ! fill code, so it must not be able to re-enter this routine.
    containers_initialized = .true.

    ! ---- Converted one-shot types: the schema IS the fill routine ----------
    ! Any container still without a schema is a converted type: run its fill
    ! once against block 1 purely for its define-on-first-touch side effect
    ! (the data it stores is garbage that the first real fill overwrites),
    ! then finish_priming locks the schema. Unconverted types were populated
    ! above (nVars > 0), so the loop never runs their fills at init, which
    ! is what keeps the lazily-allocated mag-grid types out.
    do i = 1, nContainers
      if (containers(i)%nVars == 0 .and. associated(containers(i)%fill)) then
        call containers(i)%fill(containers(i), 1)
        call finish_priming(containers(i))
      endif
    enddo

    call assert_types_fully_wired()
  end subroutine init_output_containers

  ! ---------------------------------------------------------------------------
  ! set_fill: register a type's fill routine on its container.
  ! ---------------------------------------------------------------------------
  subroutine set_fill(cType, f)
    character(len=*), intent(in) :: cType
    procedure(fill_iface) :: f
    containers(get_container_idx(cType))%fill => f
  end subroutine set_fill

  ! ---------------------------------------------------------------------------
  ! assert_types_fully_wired: init-time self-check over containers(:).
  ! Catches the silent add-a-type mistakes at startup instead of at analysis
  ! time: a slot whose define/prime call was forgotten (no cType or empty
  ! schema), a cType missing from the lookup lists, or an index/cType mismatch.
  ! ---------------------------------------------------------------------------
  subroutine assert_types_fully_wired()
    integer :: i
    do i = 1, nContainers
      if (len_trim(containers(i)%cType) == 0) &
        call stop_gitm('output container slot has no cType'// &
                       ' (forgot its define/prime call in init_output_containers?)')
      if (get_container_idx(containers(i)%cType) /= i) &
        call stop_gitm('output container index/cType mismatch: '// &
                       trim(containers(i)%cType))
      if (.not. is_known_output_type(containers(i)%cType)) &
        call stop_gitm('output type missing from is_known_output_type: '// &
                       trim(containers(i)%cType))
      if (containers(i)%nVars == 0) &
        call stop_gitm('output container has an empty schema'// &
                       ' (forgot define/prime?): '//trim(containers(i)%cType))
      if (.not. associated(containers(i)%fill)) &
        call stop_gitm('output container has no fill registered'// &
                       ' (forgot its set_fill line?): '//trim(containers(i)%cType))
    enddo
  end subroutine assert_types_fully_wired

  ! ---------------------------------------------------------------------------
  ! finish_priming: close out a converted type's priming pass.
  !
  ! Clears the participation flag that the priming fill's prepare() set (the
  ! real per-block fills set it again), then locks the schema so any later put
  ! of a name the priming pass did not establish aborts instead of quietly
  ! adding a variable on some ranks or some timesteps.
  ! ---------------------------------------------------------------------------
  subroutine finish_priming(c)
    type(OutputContainer), intent(inout) :: c
    call c%reset()
    c%schema_locked = .true.
    c%is_one_shot = .true.
  end subroutine finish_priming

  ! ---------------------------------------------------------------------------
  ! define_schema_2dmel: register all 2DMEL variables.
  !
  ! Grid: magnetic (nMagLons+1) x nMagLats, no altitude.
  ! Axis variables are 1D; all others are 2D (degenerate 3rd dim = 1).
  ! Units: MagLon/MagLat are already degrees in the model; MLT is hours.
  ! GeoLat/GeoLon are PUT RAW, in radians -- the declared 'degrees_east' /
  ! 'degrees_north' below describe what the NETCDF backend publishes, since it
  ! is the only backend that converts (see nc_is_geo_lonlat in
  ! ModOutputNetCDF.F90).  legacy/mpiio .bin carry radians, matching develop.
  ! ---------------------------------------------------------------------------
  subroutine define_schema_2dmel(c)
    use ModElectrodynamics, only: nMagLats, nMagLons
    type(OutputContainer), intent(inout) :: c

    c%cType = '2DMEL'
    ! gridKind drives BOTH participation (prepare: only iProc==0 writes) AND
    ! file layout (mpiio_write_container: MAG_2D writes at offset 0, not
    ! striped by global block index). Changing the participation rule for
    ! MAG_2D requires updating the offset logic in lockstep.
    c%gridKind = GRID_MAG_2D

    ! Axis variables (1D coordinate arrays).
    call c%define_var('mlon', units='degrees', &
                      shape3=[nMagLons + 1, 1, 1], is_axis=.true., &
                      longName='Magnetic longitude')
    call c%define_var('mlat', units='degrees', &
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
                      shape3=[nMagLons + 1, nMagLats, 1])
    call c%define_var('HalCond', units='S', &
                      shape3=[nMagLons + 1, nMagLats, 1])
    call c%define_var('DivJuAlt', units='A/m2', &
                      shape3=[nMagLons + 1, nMagLats, 1])
    call c%define_var('FL_length', units='m', &
                      shape3=[nMagLons + 1, nMagLats, 1])
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
    call c%define_var('Klambda', units='', &
                      shape3=[nMagLons + 1, nMagLats, 1], &
                      longName='Dynamo Klambda coefficient')
  end subroutine define_schema_2dmel

  ! ---------------------------------------------------------------------------
  ! fill_2dmel: put 2DMEL data into the container for iBlock.
  !
  ! Only iBlock == 1 writes (GRID_MAG_2D); all others return after prepare().
  ! GeoLat/GeoLon are written RAW, in radians, matching develop's .bin contract
  ! (develop's output_common.f90 has no cRadToDeg anywhere).  The netcdf backend
  ! converts on write; legacy/mpiio ship radians.  MagLonMC/MagLatMC are already
  ! in degrees, so this type mixes the two on purpose -- as develop does.
  ! ---------------------------------------------------------------------------
  subroutine fill_2dmel(c, iBlock)
    use ModElectrodynamics
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
    call c%put('MLT', real(MagLocTimeMC, output_kind))
    call c%put('GeoLat', real(GeoLatMC, output_kind))
    call c%put('GeoLon', real(GeoLonMC, output_kind))
    call c%put('PedCond', real(SigmaPedersenMC, output_kind))
    call c%put('HalCond', real(SigmaHallMC, output_kind))
    call c%put('DivJuAlt', real(DivJuFieldLineMC, output_kind))
    call c%put('FL_length', real(LengthMC, output_kind))
    call c%put('Sigma_PP', real(SigmaPPMC, output_kind))
    call c%put('Sigma_LL', real(SigmaLLMC, output_kind))
    call c%put('Sigma_H', real(SigmaHHMC, output_kind))
    call c%put('Sigma_C', real(SigmaCCMC, output_kind))
    call c%put('Sigma_PL', real(SigmaPLMC, output_kind))
    call c%put('Sigma_LP', real(SigmaLPMC, output_kind))
    call c%put('K_mphi', real(KDpmMC, output_kind))
    call c%put('K_mlambda', real(KDlmMC, output_kind))
    call c%put('Solver_A', real(solver_a_mc, output_kind))
    call c%put('Solver_B', real(solver_b_mc, output_kind))
    call c%put('Solver_C', real(solver_c_mc, output_kind))
    call c%put('Solver_D', real(solver_d_mc, output_kind))
    call c%put('Solver_E', real(solver_e_mc, output_kind))
    call c%put('Solver_S', real(solver_s_mc, output_kind))
    call c%put('Phi_dynamo', real(DynamoPotentialMC, output_kind))
    call c%put('Ed1new', real(Ed1new, output_kind))
    call c%put('Ed2new', real(Ed2new, output_kind))
    call c%put('Kphi', real(kpmMC, output_kind))
    call c%put('Klambda', real(klmMC, output_kind))
  end subroutine fill_2dmel

  ! ---------------------------------------------------------------------------
  ! fill_2dtec: 2DTEC schema AND data, one-shot (no define_schema_2dtec).
  !
  ! Interior-writing (nGhostCells = 0), like 2DANC. Longitude/Latitude/
  ! Altitude_GB/VTEC are ghosted (-1:nLons+2, ...) and sliced to the interior;
  ! Sza is (nLons, nLats, nBlocksMax) with NO ghosts, so the block slice already
  ! IS the interior.
  ! ---------------------------------------------------------------------------
  subroutine fill_2dtec(c, iBlock)
    use ModGITM
    use ModEUV, only: Sza
    type(OutputContainer), intent(inout) :: c
    integer, intent(in) :: iBlock

    c%cType = '2DTEC'
    c%gridKind = GRID_GEO_2D
    c%nGhostCells = 0

    call c%prepare(iBlock)
    if (.not. c%this_rank_writes) return

    call calc_vtec(iBlock)

    ! NOTE: put order is the file's variable order for legacy/mpiio, so keep it.
    call c%put_axis('Longitude', Longitude(1:nLons, iBlock))
    call c%put_axis('Latitude', Latitude(1:nLats, iBlock))
    call c%put_field('Altitude', Altitude_GB(1:nLons, 1:nLats, 1, iBlock))  ! ghosted, one level
    call c%put_field('SZA', Sza(:, :, iBlock))
    call c%put_field('TEC', VTEC(1:nLons, 1:nLats, iBlock))
  end subroutine fill_2dtec

  ! ---------------------------------------------------------------------------
  ! fill_2dgel: 2DGEL schema AND data, one-shot (no define_schema_2dgel).
  !
  ! Interior-writing (nGhostCells = 0). Unlike 2DTEC, every array here is
  ! ghosted (-1:nLons+2, -1:nLats+2, ...): Potential/PedersenConductance/
  ! HallConductance carry a block dim, the precipitation and field-line arrays
  ! do not; every put still slices to the interior explicitly.
  ! ModGITM/ModElectrodynamics 2D arrays (ElectronAverageEnergyDiffuse etc.)
  ! are block-local: no iBlock index, hold the current block's data at call
  ! time.
  ! ---------------------------------------------------------------------------
  subroutine fill_2dgel(c, iBlock)
    use ModGITM
    use ModElectrodynamics
    type(OutputContainer), intent(inout) :: c
    integer, intent(in) :: iBlock

    c%cType = '2DGEL'
    c%gridKind = GRID_GEO_2D
    c%nGhostCells = 0

    call c%prepare(iBlock)
    if (.not. c%this_rank_writes) return

    call c%put_axis('Longitude', Longitude(1:nLons, iBlock))
    call c%put_axis('Latitude', Latitude(1:nLats, iBlock))
    call c%put_field('Altitude', Altitude_GB(1:nLons, 1:nLats, 1, iBlock))  ! ghosted, one level
    call c%put_field('Potential', Potential(1:nLons, 1:nLats, 1, iBlock))       ! ghosted, one level
    call c%put_field('PedCond', PedersenConductance(1:nLons, 1:nLats, iBlock))
    call c%put_field('HalCond', HallConductance(1:nLons, 1:nLats, iBlock))
    call c%put_field('AveE', ElectronAverageEnergyDiffuse(1:nLons, 1:nLats))
    call c%put_field('eFlux', ElectronEnergyFluxDiffuse(1:nLons, 1:nLats))
    call c%put_field('AveE_W', ElectronAverageEnergyWave(1:nLons, 1:nLats))
    call c%put_field('eFlux_W', ElectronEnergyFluxWave(1:nLons, 1:nLats))
    call c%put_field('AveE_M', ElectronAverageEnergyMono(1:nLons, 1:nLats))
    call c%put_field('eFlux_M', ElectronEnergyFluxMono(1:nLons, 1:nLats))
    call c%put_field('AveE_I', IonAverageEnergy(1:nLons, 1:nLats))
    call c%put_field('eFlux_I', IonEnergyFlux(1:nLons, 1:nLats))
    call c%put_field('DivJuAlt', DivJuAlt(1:nLons, 1:nLats))
    call c%put_field('PedFLCond', PedersenFieldLine(1:nLons, 1:nLats))
    call c%put_field('HalFLCond', HallFieldLine(1:nLons, 1:nLats))
    call c%put_field('DivJuFL', DivJuFieldLine(1:nLons, 1:nLats))
    call c%put_field('FL_length', LengthFieldLine(1:nLons, 1:nLats))
  end subroutine fill_2dgel

  ! ---------------------------------------------------------------------------
  ! fill_3dneu: 3DNEU schema AND data, one-shot (no define_schema_3dneu).
  ! 2-ghosted like 3DALL: the whole array is the slice.
  ! ---------------------------------------------------------------------------
  subroutine fill_3dneu(c, iBlock)
    use ModGITM
    use ModPlanet, only: nSpeciesTotal, nSpecies, cSpecies
    type(OutputContainer), intent(inout) :: c
    integer, intent(in) :: iBlock
    integer :: i

    c%cType = '3DNEU'
    c%gridKind = GRID_GEO_3D
    c%nGhostCells = 2

    call c%prepare(iBlock)
    if (.not. c%this_rank_writes) return

    call c%put_axis('Longitude', Longitude(:, iBlock))
    call c%put_axis('Latitude', Latitude(:, iBlock))
    call c%put_field('Altitude', Altitude_GB(:, :, :, iBlock))
    call c%put_field('Rho', Rho(:, :, :, iBlock))
    do i = 1, nSpeciesTotal
      call c%put_field('['//trim(cSpecies(i))//']', NDensityS(:, :, :, i, iBlock), &
                       units='/m3', longName='Number density of '//trim(cSpecies(i)))
    enddo
    call c%put_field('Temperature', Temperature(:, :, :, iBlock)*TempUnit)
    call c%put_field('Vn (east)', Velocity(:, :, :, 1, iBlock))
    call c%put_field('Vn (north)', Velocity(:, :, :, 2, iBlock))
    call c%put_field('Vn (up)', Velocity(:, :, :, 3, iBlock))
    do i = 1, nSpecies
      call c%put_field('Vn (up,'//trim(cSpecies(i))//')', VerticalVelocity(:, :, :, i, iBlock), &
                       units='m/s', longName='Vertical neutral wind for '//trim(cSpecies(i)))
    enddo
  end subroutine fill_3dneu

  ! ---------------------------------------------------------------------------
  ! fill_3dall: 3DALL schema AND data, one-shot (no define_schema_3dall).
  !
  ! Priming (init_output_containers) defines every variable as a side effect of
  ! the first put; later calls just refill, so schema and data can't disagree.
  ! This is the canonical 2-ghost 3D type: every array is declared with exactly
  ! the ghosts written (-1:nLons+2, -1:nLats+2, -1:nAlts+2), so the whole array
  ! is the slice and is passed plainly. nGhostCells records that width;
  ! legacy/mpiio write it as-is, netcdf strips that many cells per side. Types
  ! that write an interior subset (e.g. 2DANC) slice explicitly instead.
  ! ---------------------------------------------------------------------------
  subroutine fill_3dall(c, iBlock)
    use ModGITM
    use ModPlanet, only: nSpeciesTotal, nSpecies, nIons, cSpecies, cIons
    type(OutputContainer), intent(inout) :: c
    integer, intent(in) :: iBlock
    integer :: i

    ! Set before prepare(): it reads gridKind, and output_common reads gridKind
    ! and nGhostCells before it ever calls this routine (priming establishes
    ! them for that first read).
    c%cType = '3DALL'
    c%gridKind = GRID_GEO_3D
    c%nGhostCells = 2

    call c%prepare(iBlock)
    if (.not. c%this_rank_writes) return

    call c%put_axis('Longitude', Longitude(:, iBlock))
    call c%put_axis('Latitude', Latitude(:, iBlock))
    call c%put_field('Altitude', Altitude_GB(:, :, :, iBlock))
    call c%put_field('Rho', Rho(:, :, :, iBlock))
    do i = 1, nSpeciesTotal
      call c%put_field('['//trim(cSpecies(i))//']', NDensityS(:, :, :, i, iBlock), &
                       units='/m3', longName='Number density of '//trim(cSpecies(i)))
    enddo
    call c%put_field('Temperature', Temperature(:, :, :, iBlock)*TempUnit)
    call c%put_field('Vn (east)', Velocity(:, :, :, 1, iBlock))
    call c%put_field('Vn (north)', Velocity(:, :, :, 2, iBlock))
    call c%put_field('Vn (up)', Velocity(:, :, :, 3, iBlock))
    do i = 1, nSpecies
      call c%put_field('Vn (up,'//trim(cSpecies(i))//')', VerticalVelocity(:, :, :, i, iBlock), &
                       units='m/s', longName='Vertical neutral wind for '//trim(cSpecies(i)))
    enddo
    do i = 1, nIons
      call c%put_field('['//trim(cIons(i))//']', IDensityS(:, :, :, i, iBlock), &
                       units='/m3', longName='Number density of '//trim(cIons(i)))
    enddo
    call c%put_field('eTemperature', eTemperature(:, :, :, iBlock))
    call c%put_field('iTemperature', ITemperature(:, :, :, iBlock))
    call c%put_field('Vi (east)', IVelocity(:, :, :, 1, iBlock))
    call c%put_field('Vi (north)', IVelocity(:, :, :, 2, iBlock))
    call c%put_field('Vi (up)', IVelocity(:, :, :, 3, iBlock))
  end subroutine fill_3dall

  ! ---------------------------------------------------------------------------
  ! fill_3dion: 3DION schema AND data, one-shot (no define_schema_3dion).
  ! 2-ghosted like 3DALL, including Collisions, MLatitude, MLongitude, B0,
  ! EField, Potential, IonPressureGradient: the whole array is the slice.
  ! ---------------------------------------------------------------------------
  subroutine fill_3dion(c, iBlock)
    use ModGITM
    use ModPlanet, only: nIons, cIons
    use ModElectrodynamics
    type(OutputContainer), intent(inout) :: c
    integer, intent(in) :: iBlock
    integer :: i

    c%cType = '3DION'
    c%gridKind = GRID_GEO_3D
    c%nGhostCells = 2

    call c%prepare(iBlock)
    if (.not. c%this_rank_writes) return

    call c%put_axis('Longitude', Longitude(:, iBlock))
    call c%put_axis('Latitude', Latitude(:, iBlock))
    call c%put_field('Altitude', Altitude_GB(:, :, :, iBlock))
    do i = 1, nIons
      call c%put_field('['//trim(cIons(i))//']', IDensityS(:, :, :, i, iBlock), &
                       units='/m3', longName='Number density of '//trim(cIons(i)))
    enddo
    call c%put_field('eTemperature', eTemperature(:, :, :, iBlock))
    call c%put_field('iTemperature', ITemperature(:, :, :, iBlock))
    call c%put_field('Vi (east)', IVelocity(:, :, :, 1, iBlock))
    call c%put_field('Vi (north)', IVelocity(:, :, :, 2, iBlock))
    call c%put_field('Vi (up)', IVelocity(:, :, :, 3, iBlock))
    call c%put_field('Ed1', ed1)
    call c%put_field('Ed2', ed2)
    call c%put_field('Je1', je1)
    call c%put_field('Je2', je2)
    call c%put_field('Magnetic Latitude', MLatitude(:, :, :, iBlock))
    call c%put_field('Magnetic Longitude', MLongitude(:, :, :, iBlock))
    call c%put_field('B.F. East', B0(:, :, :, 1, iBlock))
    call c%put_field('B.F. North', B0(:, :, :, 2, iBlock))
    call c%put_field('B.F. Vertical', B0(:, :, :, 3, iBlock))
    call c%put_field('B.F. Magnitude', B0(:, :, :, 4, iBlock))
    call c%put_field('Potential', Potential(:, :, :, iBlock))
    call c%put_field('E.F. East', EField(:, :, :, 1))
    call c%put_field('E.F. North', EField(:, :, :, 2))
    call c%put_field('E.F. Vertical', EField(:, :, :, 3))
    call c%put_field('E.F. Magnitude', sqrt(sum(EField**2, dim=4)))
    call c%put_field('IN Collision Freq', Collisions(:, :, :, iVIN_))
    call c%put_field('IonPressGrad (east)', IonPressureGradient(:, :, :, 1, iBlock))
    call c%put_field('IonPressGrad (north)', IonPressureGradient(:, :, :, 2, iBlock))
    call c%put_field('IonPressGrad (up)', IonPressureGradient(:, :, :, 3, iBlock))
  end subroutine fill_3dion

  ! ---------------------------------------------------------------------------
  ! fill_3dthm: 3DTHM schema AND data, one-shot (no define_schema_3dthm).
  !
  ! Interior-writing (nGhostCells = 0), MIXED ghosting: Longitude/Latitude/
  ! Altitude_GB/Rho/TempUnit/EField/Sigma_Pedersen are 2-ghost and sliced to
  ! the interior; Conduction/MoleConduction carry ghosts only in altitude
  ! (sliced to 1:nAlts); cp is per-dimension (whole lon/lat, sliced altitude);
  ! EuvHeating/EuvTotal/AuroralIonRateS/JouleHeating/ChemicalHeatingRate/
  ! NOCooling/OCooling are ghost-less and pass whole.
  ! ---------------------------------------------------------------------------
  subroutine fill_3dthm(c, iBlock)
    use ModGITM
    use ModSources
    use ModElectrodynamics, only: Sigma_Pedersen
    use ModEUV, only: EuvTotal
    type(OutputContainer), intent(inout) :: c
    integer, intent(in) :: iBlock

    c%cType = '3DTHM'
    c%gridKind = GRID_GEO_3D
    c%nGhostCells = 0

    call c%prepare(iBlock)
    if (.not. c%this_rank_writes) return

    call c%put_axis('Longitude', Longitude(1:nLons, iBlock))
    call c%put_axis('Latitude', Latitude(1:nLats, iBlock))
    call c%put_field('Altitude', Altitude_GB(1:nLons, 1:nLats, 1:nAlts, iBlock))
    call c%put_field('EUV Heating Rate', &                                  ! no ghosts * ghosted
                     EuvHeating(:, :, :, iBlock)*TempUnit(1:nLons, 1:nLats, 1:nAlts))
    call c%put_field('Conduction', &                                  ! alt ghosts only * ghosted
                     per_dt(Conduction(:, :, 1:nAlts)*TempUnit(1:nLons, 1:nLats, 1:nAlts)))
    ! MoleConduction from the last (full-Dt) calc_conduction_1d call is a
    ! Kelvin temperature change over the step, so /dt alone gives K/s.
    call c%put_field('Molecular Conduction', per_dt(MoleConduction(:, :, 1:nAlts)))  ! alt ghosts only
    call c%put_field('Chemical Heating Rate', &                            ! no ghosts * ghosted
                     per_dt(ChemicalHeatingRate*TempUnit(1:nLons, 1:nLats, 1:nAlts)))
    call c%put_field('Joule Heating Rate', &                               ! no ghosts * ghosted
                     JouleHeating*TempUnit(1:nLons, 1:nLats, 1:nAlts))
    call c%put_field('NO Cooling Rate', &                                  ! no ghosts * ghosted
                     -NOCooling*TempUnit(1:nLons, 1:nLats, 1:nAlts))
    call c%put_field('O Cooling Rate', &                                   ! no ghosts * ghosted
                     -OCooling*TempUnit(1:nLons, 1:nLats, 1:nAlts))
    call c%put_field('Total Abs EUV', EuvTotal(:, :, :, iBlock))
    call c%put_field('Cp', cp(:, :, 1:nAlts, iBlock))                       ! per-dim: lon/lat whole, alt sliced
    call c%put_field('Rho', Rho(1:nLons, 1:nLats, 1:nAlts, iBlock))
    call c%put_field('E.F. Magnitude', &
                     sqrt(sum(EField(1:nLons, 1:nLats, 1:nAlts, :)**2, dim=4)))
    call c%put_field('Sigma Ped', Sigma_Pedersen(1:nLons, 1:nLats, 1:nAlts))
    call c%put_field('Ionization Rate O_3P', AuroralIonRateS(:, :, :, iO_3P_, iBlock))
    call c%put_field('Ionization Rate O2', AuroralIonRateS(:, :, :, iO2_, iBlock))
    call c%put_field('Ionization Rate N2', AuroralIonRateS(:, :, :, iN2_, iBlock))
  end subroutine fill_3dthm

  ! ---------------------------------------------------------------------------
  ! fill_3dchm: 3DCHM schema AND data, one-shot (no define_schema_3dchm).
  !
  ! Interior-writing (nGhostCells = 0), MIXED ghosting: Longitude/Latitude/
  ! Altitude_GB/Rho/TempUnit are 2-ghost and sliced to the interior; cp is
  ! per-dimension (whole lon/lat, sliced altitude); ChemicalHeatingSpecies and
  ! ChemicalHeatingRate are ghost-less and pass whole.
  ! ---------------------------------------------------------------------------
  subroutine fill_3dchm(c, iBlock)
    use ModGITM
    use ModSources
    use ModConstants, only: Element_Charge
    type(OutputContainer), intent(inout) :: c
    integer, intent(in) :: iBlock
    integer :: i
    character(len=20), parameter :: &
      chem_names(nReactions) = [ &
      'N2+ + e             ', 'O2+ + e             ', 'N2+ + O             ', &
      'NO+ + e             ', 'N+ + O2             ', 'NO + N              ', &
      'O+ + O2             ', 'N + O2              ', 'O2+ + N             ', &
      'O2+ + NO            ', 'O2+ + N2            ', 'N2+ + O2            ', &
      'N+ + O              ', 'O+ + N2             ', 'O(1D) + N2          ', &
      'O(1D) + O2          ', 'O(1D) + O           ', 'O(1D) + e           ', &
      'N(2D) + O2          ', 'O+(2D) + e          ', 'N(2D) + O           ', &
      'N(2D) + e           ', 'O+(2D) + N2         ', 'O+(2P) + e          ', &
      'O+(2P) + O          ', 'O+(2P) + N2         ']
    ! Parallels chem_names index-for-index.
    character(len=40), parameter :: &
      chem_longnames(nReactions) = &
      [character(len=40) :: &
       'N2+ + e recombination heating rate', 'O2+ + e recombination heating rate', &
       'N2+ + O reaction heating rate', 'NO+ + e recombination heating rate', &
       'N+ + O2 reaction heating rate', 'NO + N reaction heating rate', &
       'O+ + O2 reaction heating rate', 'N + O2 reaction heating rate', &
       'O2+ + N reaction heating rate', 'O2+ + NO reaction heating rate', &
       'O2+ + N2 reaction heating rate', 'N2+ + O2 reaction heating rate', &
       'N+ + O reaction heating rate', 'O+ + N2 reaction heating rate', &
       'O(1D) + N2 quench heating rate', 'O(1D) + O2 quench heating rate', &
       'O(1D) + O quench heating rate', 'O(1D) + e quench heating rate', &
       'N(2D) + O2 reaction heating rate', 'O+(2D) + e recombination heating rate', &
       'N(2D) + O quench heating rate', 'N(2D) + e quench heating rate', &
       'O+(2D) + N2 reaction heating rate', 'O+(2P) + e recombination heating rate', &
       'O+(2P) + O reaction heating rate', 'O+(2P) + N2 reaction heating rate']

    c%cType = '3DCHM'
    c%gridKind = GRID_GEO_3D
    c%nGhostCells = 0

    call c%prepare(iBlock)
    if (.not. c%this_rank_writes) return

    call c%put_axis('Longitude', Longitude(1:nLons, iBlock))
    call c%put_axis('Latitude', Latitude(1:nLats, iBlock))
    call c%put_field('Altitude', Altitude_GB(1:nLons, 1:nLats, 1:nAlts, iBlock))
    do i = 1, nReactions
      ! calc_chemistry accumulates ChemicalHeatingSpecies as <term>*DtSub over
      ! one physics step, so /Dt recovers the mean heating rate over the step
      ! (Dt at write time is still the step the accumulation integrated over).
      call c%put_field(trim(chem_names(i)), per_dt(ChemicalHeatingSpecies(:, :, :, i)/Element_Charge), &
                       units='eV/m3/s', longName=trim(chem_longnames(i)))
    enddo
    call c%put_field('Chemical Heating', &                                  ! no ghosts * per-dim * ghosted * ghosted
                     per_dt(ChemicalHeatingRate*cp(:, :, 1:nAlts, iBlock)* &
                     Rho(1:nLons, 1:nLats, 1:nAlts, iBlock)* &
                     TempUnit(1:nLons, 1:nLats, 1:nAlts)/Element_Charge))
  end subroutine fill_3dchm

  ! ---------------------------------------------------------------------------
  ! fill_3dlst: 3DLST schema AND data, one-shot (no define_schema_3dlst).
  ! 2-ghosted like 3DALL: the whole array is the slice.
  !
  ! Puts are gated by the UAM.in i*OutputList flags, which are uniform across
  ! ranks and time-invariant, so priming sees the same conditions as every fill.
  ! That is what makes the conditional puts safe; a rank- or time-dependent
  ! condition would desync the schema. All list flags default .true.
  ! HARD INVARIANT: the flags are set once (input parse / init_planet) before
  ! the first output and never written again. Anything that mutates them
  ! mid-run, or moves priming earlier, breaks the locked schema and aborts.
  ! ---------------------------------------------------------------------------
  subroutine fill_3dlst(c, iBlock)
    use ModGITM
    use ModPlanet, only: nSpeciesTotal, nIons, cSpecies, cIons
    use ModInputs, only: iRhoOutputList, iNeutralDensityOutputList, &
                         iNeutralWindOutputList, iIonDensityOutputList, &
                         iIonWindOutputList, iTemperatureOutputList
    type(OutputContainer), intent(inout) :: c
    integer, intent(in) :: iBlock
    integer :: i

    c%cType = '3DLST'
    c%gridKind = GRID_GEO_3D
    c%nGhostCells = 2

    call c%prepare(iBlock)
    if (.not. c%this_rank_writes) return

    call c%put_axis('Longitude', Longitude(:, iBlock))
    call c%put_axis('Latitude', Latitude(:, iBlock))
    call c%put_field('Altitude', Altitude_GB(:, :, :, iBlock))

    if (iRhoOutputList) call c%put_field('Rho', Rho(:, :, :, iBlock))
    do i = 1, nSpeciesTotal
      if (iNeutralDensityOutputList(i)) then
        call c%put_field('['//trim(cSpecies(i))//']', NDensityS(:, :, :, i, iBlock), &
                         units='/m3', longName='Number density of '//trim(cSpecies(i)))
      endif
    enddo
    if (iNeutralWindOutputList(1)) call c%put_field('Vn (east)', Velocity(:, :, :, 1, iBlock))
    if (iNeutralWindOutputList(2)) call c%put_field('Vn (north)', Velocity(:, :, :, 2, iBlock))
    if (iNeutralWindOutputList(3)) call c%put_field('Vn (up)', Velocity(:, :, :, 3, iBlock))
    do i = 1, nIons
      if (iIonDensityOutputList(i)) then
        call c%put_field('['//trim(cIons(i))//']', IDensityS(:, :, :, i, iBlock), &
                         units='/m3', longName='Number density of '//trim(cIons(i)))
      endif
    enddo
    if (iIonWindOutputList(1)) call c%put_field('Vi (east)', IVelocity(:, :, :, 1, iBlock))
    if (iIonWindOutputList(2)) call c%put_field('Vi (north)', IVelocity(:, :, :, 2, iBlock))
    if (iIonWindOutputList(3)) call c%put_field('Vi (up)', IVelocity(:, :, :, 3, iBlock))
    if (iTemperatureOutputList(1)) &
      call c%put_field('Temperature', Temperature(:, :, :, iBlock)*TempUnit)
    if (iTemperatureOutputList(2)) call c%put_field('iTemperature', ITemperature(:, :, :, iBlock))
    if (iTemperatureOutputList(3)) call c%put_field('eTemperature', eTemperature(:, :, :, iBlock))
  end subroutine fill_3dlst

  ! ---------------------------------------------------------------------------
  ! fill_3dglo: 3DGLO schema AND data, one-shot (no define_schema_3dglo).
  ! 2-ghosted like 3DALL: the whole array is the slice.
  !
  ! The three emission/photoelectron fields are placeholders (not yet computed);
  ! zero_buf is a default-real scratch buffer at the full 2-ghost extent. Their
  ! units are empty strings carried through as registry rows, so today's output
  ! is unchanged.
  ! ---------------------------------------------------------------------------
  subroutine fill_3dglo(c, iBlock)
    use ModGITM
    type(OutputContainer), intent(inout) :: c
    integer, intent(in) :: iBlock
    real, allocatable :: zero_buf(:, :, :)  ! model kind, to match the put_field dummy

    ! 3DGLO has never been wired to a source of emission data on any GITM
    ! lineage: develop's output_3dglo returns immediately with its entire write
    ! commented out, and the puts below are placeholder zeros.  A column of
    ! zeros labelled '6300 A Emission' cannot be told apart from real data by
    ! anything downstream, so a run that actually asks for this type stops here
    ! rather than being handed a file that looks like science.
    !
    ! Guarded on schema_locked because init_output_containers PRIMES this fill
    ! once at startup -- that priming pass is where the schema below comes from,
    ! and it runs whether or not the type was requested.  finish_priming sets
    ! schema_locked afterwards, so only a real fill can reach the halt.
    if (c%schema_locked) &
      call stop_gitm('3DGLO writes placeholder zeros, not emission data. '// &
                     'Remove 3DGLO from #SAVEPLOT.')

    c%cType = '3DGLO'
    c%gridKind = GRID_GEO_3D
    c%nGhostCells = 2

    call c%prepare(iBlock)
    if (.not. c%this_rank_writes) return

    call c%put_axis('Longitude', Longitude(:, iBlock))
    call c%put_axis('Latitude', Latitude(:, iBlock))
    call c%put_field('Altitude', Altitude_GB(:, :, :, iBlock))

    allocate(zero_buf(nLons + 2*c%nGhostCells, nLats + 2*c%nGhostCells, nAlts + 2*c%nGhostCells))
    zero_buf = 0.0
    call c%put_field('6300 A Emission', zero_buf)
    call c%put_field('PhotoElectronUp', zero_buf)
    call c%put_field('PhotoElectronDown', zero_buf)
    deallocate(zero_buf)
  end subroutine fill_3dglo

  ! ---------------------------------------------------------------------------
  ! fill_3dmag: 3DMAG schema AND data, one-shot (no define_schema_3dmag).
  ! 2-ghosted like 3DALL: the whole array is the slice.
  ! ---------------------------------------------------------------------------
  subroutine fill_3dmag(c, iBlock)
    use ModGITM
    type(OutputContainer), intent(inout) :: c
    integer, intent(in) :: iBlock

    c%cType = '3DMAG'
    c%gridKind = GRID_GEO_3D
    c%nGhostCells = 2

    call c%prepare(iBlock)
    if (.not. c%this_rank_writes) return

    call c%put_axis('Longitude', Longitude(:, iBlock))
    call c%put_axis('Latitude', Latitude(:, iBlock))
    call c%put_field('Altitude', Altitude_GB(:, :, :, iBlock))
    call c%put_field('Magnetic Latitude', MLatitude(:, :, :, iBlock))
    call c%put_field('Magnetic Longitude', MLongitude(:, :, :, iBlock))
    call c%put_field('B.F. East', B0(:, :, :, 1, iBlock))
    call c%put_field('B.F. North', B0(:, :, :, 2, iBlock))
    call c%put_field('B.F. Vertical', B0(:, :, :, 3, iBlock))
    call c%put_field('B.F. Magnitude', B0(:, :, :, 4, iBlock))
  end subroutine fill_3dmag

  ! ---------------------------------------------------------------------------
  ! define_schema_3dhme: register variables for 3DHME.
  ! ---------------------------------------------------------------------------
  subroutine define_schema_3dhme(c)
    use ModGITM, only: nLons, nLats, nAlts
    use ModPlanet, only: nSpeciesTotal, nSpecies, nIons, cSpecies, cIons
    type(OutputContainer), intent(inout) :: c
    integer :: nGC
    integer :: i

    c%cType = '3DHME'
    c%gridKind = GRID_HIME
    c%nGhostCells = 0
    nGC = c%nGhostCells

    call c%define_var('Longitude', units='degrees_east', shape3=[nLons + 2*nGC, 1, 1], is_axis=.true.)
    call c%define_var('Latitude', units='degrees_north', shape3=[nLats + 2*nGC, 1, 1], is_axis=.true.)
    call c%define_var('Altitude', units='m', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC])
    call c%define_var('Rho', units='kg/m3', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC])
    do i = 1, nSpeciesTotal
      call c%define_var('['//trim(cSpecies(i))//']', units='/m3', &
                        longName='Number density of '//trim(cSpecies(i)), &
                        shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC])
    enddo
    call c%define_var('Temperature', units='K', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC])
    call c%define_var('Vn (east)', units='m/s', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC])
    call c%define_var('Vn (north)', units='m/s', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC])
    call c%define_var('Vn (up)', units='m/s', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC])
    do i = 1, nSpecies
      call c%define_var('Vn (up,'//trim(cSpecies(i))//')', units='m/s', &
                        longName='Vertical neutral wind for '//trim(cSpecies(i)), &
                        shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC])
    enddo
    do i = 1, nIons
      call c%define_var('['//trim(cIons(i))//']', units='/m3', &
                        longName='Number density of '//trim(cIons(i)), &
                        shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC])
    enddo
    call c%define_var('eTemperature', units='K', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC])
    call c%define_var('iTemperature', units='K', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC])
    call c%define_var('Vi (east)', units='m/s', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC])
    call c%define_var('Vi (north)', units='m/s', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC])
    call c%define_var('Vi (up)', units='m/s', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC])
    call c%define_var('PhotoElectron Heating', units='K', longName='Photoelectron heating', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC])
    call c%define_var('Joule Heating', units='K', longName='Joule heating', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC])
    call c%define_var('Specific Heat', units='J/kg/K', longName='Specific heat Cp', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC])
    call c%define_var('Magnetic Latitude', units='deg', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC])
    call c%define_var('Magnetic Longitude', units='deg', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC])
    call c%define_var('B.F. East', units='T', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC])
    call c%define_var('B.F. North', units='T', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC])
    call c%define_var('B.F. Vertical', units='T', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC])
    call c%define_var('B.F. Magnitude', units='T', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC])
    call c%define_var('Potential', units='V', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC])
    call c%define_var('PotentialY', units='V', longName='Potential Y component', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC])
    call c%define_var('E.F. East', units='V/m', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC])
    call c%define_var('E.F. North', units='V/m', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC])
    call c%define_var('E.F. Vertical', units='V/m', &
                      shape3=[nLons + 2*nGC, nLats + 2*nGC, nAlts + 2*nGC])
  end subroutine define_schema_3dhme

  ! ---------------------------------------------------------------------------
  ! fill_3dhme: put 3DHME regional variables into the container.
  ! ---------------------------------------------------------------------------
  subroutine fill_3dhme(c, iBlock)
    use ModGITM
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

    call c%put('Longitude', real(Longitude(1 - nGC:nLons + nGC, iBlock), output_kind))
    call c%put('Latitude', real(Latitude(1 - nGC:nLats + nGC, iBlock), output_kind))
    call c%put('Altitude', &
               real(Altitude_GB(1 - nGC:nLons + nGC, 1 - nGC:nLats + nGC, 1 - nGC:nAlts + nGC, iBlock), &
                    output_kind))
    call c%put('Rho', real(Rho(1 - nGC:nLons + nGC, 1 - nGC:nLats + nGC, 1 - nGC:nAlts + nGC, iBlock), &
                           output_kind))
    do i = 1, nSpeciesTotal
      call c%put('['//trim(cSpecies(i))//']', &
                 real(NDensityS(1 - nGC:nLons + nGC, 1 - nGC:nLats + nGC, 1 - nGC:nAlts + nGC, i, iBlock), &
                      output_kind))
    enddo
    call c%put('Temperature', &
               real(Temperature(1 - nGC:nLons + nGC, 1 - nGC:nLats + nGC, 1 - nGC:nAlts + nGC, iBlock) &
                    *TempUnit(1 - nGC:nLons + nGC, 1 - nGC:nLats + nGC, 1 - nGC:nAlts + nGC), &
                    output_kind))
    call c%put('Vn (east)', &
               real(Velocity(1 - nGC:nLons + nGC, 1 - nGC:nLats + nGC, 1 - nGC:nAlts + nGC, 1, iBlock), &
                    output_kind))
    call c%put('Vn (north)', &
               real(Velocity(1 - nGC:nLons + nGC, 1 - nGC:nLats + nGC, 1 - nGC:nAlts + nGC, 2, iBlock), &
                    output_kind))
    call c%put('Vn (up)', &
               real(Velocity(1 - nGC:nLons + nGC, 1 - nGC:nLats + nGC, 1 - nGC:nAlts + nGC, 3, iBlock), &
                    output_kind))
    do i = 1, nSpecies
      call c%put('Vn (up,'//trim(cSpecies(i))//')', &
                 real(VerticalVelocity(1 - nGC:nLons + nGC, 1 - nGC:nLats + nGC, 1 - nGC:nAlts + nGC, i, iBlock), &
                      output_kind))
    enddo
    do i = 1, nIons
      call c%put('['//trim(cIons(i))//']', &
                 real(IDensityS(1 - nGC:nLons + nGC, 1 - nGC:nLats + nGC, 1 - nGC:nAlts + nGC, i, iBlock), &
                      output_kind))
    enddo
    call c%put('eTemperature', &
               real(eTemperature(1 - nGC:nLons + nGC, 1 - nGC:nLats + nGC, 1 - nGC:nAlts + nGC, iBlock), &
                    output_kind))
    call c%put('iTemperature', &
               real(ITemperature(1 - nGC:nLons + nGC, 1 - nGC:nLats + nGC, 1 - nGC:nAlts + nGC, iBlock), &
                    output_kind))
    call c%put('Vi (east)', &
               real(IVelocity(1 - nGC:nLons + nGC, 1 - nGC:nLats + nGC, 1 - nGC:nAlts + nGC, 1, iBlock), &
                    output_kind))
    call c%put('Vi (north)', &
               real(IVelocity(1 - nGC:nLons + nGC, 1 - nGC:nLats + nGC, 1 - nGC:nAlts + nGC, 2, iBlock), &
                    output_kind))
    call c%put('Vi (up)', &
               real(IVelocity(1 - nGC:nLons + nGC, 1 - nGC:nLats + nGC, 1 - nGC:nAlts + nGC, 3, iBlock), &
                    output_kind))
    call c%put('PhotoElectron Heating', &
               real(PhotoElectronHeating(1 - nGC:nLons + nGC, 1 - nGC:nLats + nGC, 1 - nGC:nAlts + nGC, iBlock) &
                    *dt*TempUnit(1 - nGC:nLons + nGC, 1 - nGC:nLats + nGC, 1 - nGC:nAlts + nGC), &
                    output_kind))
    call c%put('Joule Heating', &
               real(JouleHeating(1 - nGC:nLons + nGC, 1 - nGC:nLats + nGC, 1 - nGC:nAlts + nGC) &
                    *dt*TempUnit(1 - nGC:nLons + nGC, 1 - nGC:nLats + nGC, 1 - nGC:nAlts + nGC), &
                    output_kind))
    call c%put('Specific Heat', &
               real(cp(1 - nGC:nLons + nGC, 1 - nGC:nLats + nGC, 1 - nGC:nAlts + nGC, iBlock), &
                    output_kind))
    call c%put('Magnetic Latitude', &
               real(mLatitude(1 - nGC:nLons + nGC, 1 - nGC:nLats + nGC, 1 - nGC:nAlts + nGC, iBlock), &
                    output_kind))
    call c%put('Magnetic Longitude', &
               real(mLongitude(1 - nGC:nLons + nGC, 1 - nGC:nLats + nGC, 1 - nGC:nAlts + nGC, iBlock), &
                    output_kind))
    call c%put('B.F. East', &
               real(B0(1 - nGC:nLons + nGC, 1 - nGC:nLats + nGC, 1 - nGC:nAlts + nGC, 1, iBlock), &
                    output_kind))
    call c%put('B.F. North', &
               real(B0(1 - nGC:nLons + nGC, 1 - nGC:nLats + nGC, 1 - nGC:nAlts + nGC, 2, iBlock), &
                    output_kind))
    call c%put('B.F. Vertical', &
               real(B0(1 - nGC:nLons + nGC, 1 - nGC:nLats + nGC, 1 - nGC:nAlts + nGC, 3, iBlock), &
                    output_kind))
    call c%put('B.F. Magnitude', &
               real(B0(1 - nGC:nLons + nGC, 1 - nGC:nLats + nGC, 1 - nGC:nAlts + nGC, 4, iBlock), &
                    output_kind))
    call c%put('Potential', &
               real(Potential(1 - nGC:nLons + nGC, 1 - nGC:nLats + nGC, 1 - nGC:nAlts + nGC, iBlock), &
                    output_kind))
    call c%put('PotentialY', &
               real(PotentialY(1 - nGC:nLons + nGC, 1 - nGC:nLats + nGC, 1 - nGC:nAlts + nGC, iBlock), &
                    output_kind))
    call c%put('E.F. East', &
               real(EField(1 - nGC:nLons + nGC, 1 - nGC:nLats + nGC, 1 - nGC:nAlts + nGC, 1), &
                    output_kind))
    call c%put('E.F. North', &
               real(EField(1 - nGC:nLons + nGC, 1 - nGC:nLats + nGC, 1 - nGC:nAlts + nGC, 2), &
                    output_kind))
    call c%put('E.F. Vertical', &
               real(EField(1 - nGC:nLons + nGC, 1 - nGC:nLats + nGC, 1 - nGC:nAlts + nGC, 3), &
                    output_kind))
  end subroutine fill_3dhme

  ! ---------------------------------------------------------------------------
  ! fill_3dmoh: 3DMOH schema AND data, one-shot (no define_schema_3dmoh).
  !
  ! Interior-writing (nGhostCells = 0), MIXED ghosting: Longitude/Latitude/
  ! Altitude_GB/Rho/Velocity are 2-ghost and sliced to the interior; Viscosity
  ! carries ghosts only in altitude (sliced to 1:nAlts); IonDrag/
  ! HorizAdvection/HorizPressureGrad/HorizCoriolis/Centrifugal are ghost-less
  ! and pass whole.
  ! ---------------------------------------------------------------------------
  subroutine fill_3dmoh(c, iBlock)
    use ModGITM
    use ModSources
    type(OutputContainer), intent(inout) :: c
    integer, intent(in) :: iBlock

    c%cType = '3DMOH'
    c%gridKind = GRID_GEO_3D
    c%nGhostCells = 0

    call c%prepare(iBlock)
    if (.not. c%this_rank_writes) return

    call c%put_axis('Longitude', Longitude(1:nLons, iBlock))
    call c%put_axis('Latitude', Latitude(1:nLats, iBlock))
    call c%put_field('Altitude', Altitude_GB(1:nLons, 1:nLats, 1:nAlts, iBlock))
    call c%put_field('Rho', Rho(1:nLons, 1:nLats, 1:nAlts, iBlock))
    call c%put_field('Vn (east)', Velocity(1:nLons, 1:nLats, 1:nAlts, iEast_, iBlock))
    call c%put_field('Vn (north)', Velocity(1:nLons, 1:nLats, 1:nAlts, iNorth_, iBlock))
    call c%put_field('Visc_Ve_rshear', per_dt(Viscosity(:, :, 1:nAlts, iEast_)))     ! alt ghosts only
    call c%put_field('Visc_Vn_rshear', per_dt(Viscosity(:, :, 1:nAlts, iNorth_)))    ! alt ghosts only
    call c%put_field('IonDrag_east', IonDrag(:, :, :, iEast_))
    call c%put_field('IonDrag_north', IonDrag(:, :, :, iNorth_))
    call c%put_field('Horiz_Adv_Ve', HorizAdvection(:, :, :, 1))
    call c%put_field('Horiz_Adv_Vn', HorizAdvection(:, :, :, 2))
    call c%put_field('PressGrad (east)', HorizPressureGrad(:, :, :, 1))
    call c%put_field('PressGrad (north)', HorizPressureGrad(:, :, :, 2))
    call c%put_field('Coriolis_east', HorizCoriolis(:, :, :, 1))
    call c%put_field('Coriolis_north', HorizCoriolis(:, :, :, 2))
    call c%put_field('Centrifugal_north', Centrifugal(:, :, :, 2))
  end subroutine fill_3dmoh

  ! ---------------------------------------------------------------------------
  ! fill_3dmov: 3DMOV schema AND data, one-shot (no define_schema_3dmov).
  !
  ! Interior-writing (nGhostCells = 0), MIXED ghosting: Longitude/Latitude/
  ! Altitude_GB/VerticalVelocity are 2-ghost and sliced to the interior;
  ! VerticalIonDrag/VertAdvection/VertCoriolis/VertCentrifugal/
  ! EffectiveGravity are ghost-less and pass whole.
  ! ---------------------------------------------------------------------------
  subroutine fill_3dmov(c, iBlock)
    use ModGITM
    use ModSources
    type(OutputContainer), intent(inout) :: c
    integer, intent(in) :: iBlock
    integer :: iSpecies

    c%cType = '3DMOV'
    c%gridKind = GRID_GEO_3D
    c%nGhostCells = 0

    call c%prepare(iBlock)
    if (.not. c%this_rank_writes) return

    call c%put_axis('Longitude', Longitude(1:nLons, iBlock))
    call c%put_axis('Latitude', Latitude(1:nLats, iBlock))
    call c%put_field('Altitude', Altitude_GB(1:nLons, 1:nLats, 1:nAlts, iBlock))
    do iSpecies = 1, nSpecies
      call c%put_field('Vn (up,'//trim(cSpecies(iSpecies))//')', &
                       VerticalVelocity(1:nLons, 1:nLats, 1:nAlts, iSpecies, iBlock), &
                       units='m/s', longName='Vertical wind for '//trim(cSpecies(iSpecies)))
      call c%put_field('IonDrag_up_'//trim(cSpecies(iSpecies)), &
                       VerticalIonDrag(:, :, :, iSpecies), &
                       units='m/s2', longName='Ion drag upward for '//trim(cSpecies(iSpecies)))
      call c%put_field('Horiz_Adv_Vn (up,'//trim(cSpecies(iSpecies))//')', &
                       VertAdvection(:, :, :, iSpecies), &
                       units='m/s2', longName='Horiz advection upward for '//trim(cSpecies(iSpecies)))
    enddo
    call c%put_field('Coriolis_up', VertCoriolis)
    call c%put_field('Centrif_up', VertCentrifugal)
    call c%put_field('EffectiveGravity', EffectiveGravity)
  end subroutine fill_3dmov

  ! ---------------------------------------------------------------------------
  ! fill_2danc: 2DANC schema AND data, one-shot (no define_schema_2danc).
  !
  ! The arrays feeding this ONE container disagree about ghost cells:
  !   * VTEC, Altitude_GB, Longitude, Latitude, LocalTime are ghosted
  !     (-1:nLons+2, ...) and must be sliced down to the interior here;
  !   * Sza and the ModSources *2d arrays are declared (nLons, nLats): they
  !     have NO ghosts, so the whole array already IS the interior;
  !   * lt_buf is synthesized at interior size;
  !   * Altitude is a single altitude LEVEL of a 3-D array.
  ! No container-level "strip N ghosts" rule can describe that set, so the
  ! caller owns the slice.
  !
  ! Do NOT reintroduce a 1-nGC:nLons+nGC idiom here: it is only ever correct
  ! at nGC=0, and raising it would send Sza and JouleHeating2d out of bounds
  ! while VTEC kept working.  Explicit slices keep the difference visible at
  ! each call site.
  ! ---------------------------------------------------------------------------
  subroutine fill_2danc(c, iBlock)
    use ModGITM
    use ModElectrodynamics
    use ModEUV, only: Sza
    use ModSources
    type(OutputContainer), intent(inout) :: c
    integer, intent(in) :: iBlock
    real :: lt_buf(nLons, nLats)  ! model kind, to match the put_field dummy
    integer :: iLat

    c%cType = '2DANC'
    c%gridKind = GRID_GEO_2D
    c%nGhostCells = 0

    call c%prepare(iBlock)
    if (.not. c%this_rank_writes) return

    call calc_vtec(iBlock)

    ! NOTE: put order is the file's variable order for legacy/mpiio, so keep it.
    ! Ghosted arrays are sliced 1:nLons/1:nLats; ghost-less ones pass whole.
    call c%put_axis('Longitude', Longitude(1:nLons, iBlock))
    call c%put_axis('Latitude', Latitude(1:nLats, iBlock))
    ! Candidate for removal (along with dropping Altitude from other 2D types) once
    ! PyITM's altitude handling no longer depends on it.
    call c%put_field('Altitude', Altitude_GB(1:nLons, 1:nLats, nAlts, iBlock))  ! ghosted, one level

    ! Replicate the interior of LocalTime (ghosted) along the Latitude dimension
    do iLat = 1, nLats
      lt_buf(:, iLat) = LocalTime(1:nLons)
    enddo
    call c%put_field('LT', lt_buf)                                       ! interior-sized already

    call c%put_field('SZA', Sza(:, :, iBlock))
    call c%put_field('TEC', VTEC(1:nLons, 1:nLats, iBlock))
    call c%put_field('JH_int', JouleHeating2d)                           ! no ghosts, no block dim
    call c%put_field('HT_int', HeatTransfer2d)
    call c%put_field('EUV_int', EuvHeating2d)
    call c%put_field('PE_int', PhotoElectronHeating2d)
    ! ChemicalHeating2d is J/m2 deposited over the step; /Dt gives the W/m2
    ! the registry promises.
    call c%put_field('Chem_int', per_dt(ChemicalHeating2d))
    call c%put_field('RadCool_int', RadiativeCooling2d)
    call c%put_field('CO2Cool_int', CO2Cooling2d)
    call c%put_field('NOCool_int', NOCooling2d)
    call c%put_field('OCool_int', OCooling2d)
  end subroutine fill_2danc

  ! ---------------------------------------------------------------------------
  ! define_schema_2dhme: register variables for 2DHME.
  ! ---------------------------------------------------------------------------
  subroutine define_schema_2dhme(c)
    use ModGITM, only: nLons, nLats
    type(OutputContainer), intent(inout) :: c
    integer :: nGC

    c%cType = '2DHME'
    c%gridKind = GRID_HIME
    c%nGhostCells = 0
    nGC = c%nGhostCells

    call c%define_var('Longitude', units='degrees_east', shape3=[nLons + 2*nGC, 1, 1], is_axis=.true.)
    call c%define_var('Latitude', units='degrees_north', shape3=[nLats + 2*nGC, 1, 1], is_axis=.true.)
    call c%define_var('Altitude', units='m', shape3=[nLons + 2*nGC, nLats + 2*nGC, 1])
    call c%define_var('LT', units='hr', shape3=[nLons + 2*nGC, nLats + 2*nGC, 1])
    call c%define_var('TEC', units='TECU', shape3=[nLons + 2*nGC, nLats + 2*nGC, 1])
  end subroutine define_schema_2dhme

  ! ---------------------------------------------------------------------------
  ! fill_2dhme: put 2DHME regional variables into the container.
  ! ---------------------------------------------------------------------------
  subroutine fill_2dhme(c, iBlock)
    use ModGITM
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

    call c%put('Longitude', real(Longitude(1 - nGC:nLons + nGC, iBlock), output_kind))
    call c%put('Latitude', real(Latitude(1 - nGC:nLats + nGC, iBlock), output_kind))
    call c%put('Altitude', real(Altitude_GB(1 - nGC:nLons + nGC, 1 - nGC:nLats + nGC, 1, iBlock), &
                                output_kind))

    do iLat = 1, nLats
      lt_buf(:, iLat) = real(LocalTime(1 - nGC:nLons + nGC), output_kind)
    enddo
    call c%put('LT', lt_buf)
    call c%put('TEC', real(VTEC(1 - nGC:nLons + nGC, 1 - nGC:nLats + nGC, iBlock), output_kind))
  end subroutine fill_2dhme

  ! ---------------------------------------------------------------------------
  ! fill_1dall: put 1DALL interpolated vertical column data.
  ! ---------------------------------------------------------------------------
  subroutine fill_1dall(c, iBlock)
    use ModGITM
    use ModPlanet, only: nSpeciesTotal, nSpecies, nIons, cSpecies, cIons
    type(OutputContainer), intent(inout) :: c
    integer, intent(in) :: iBlock
    ! Route A1: probe location comes from module state (set by the driver).
    integer :: iiLon, iiLat
    real :: rLon, rLat
    real :: tmp_lon(nAlts), tmp_lat(nAlts), tmp_data(nAlts)
    real :: slice2d(0:nLons + 1, 0:nLats + 1)
    integer :: iAlt, i

    iiLon = iProbeLon; iiLat = iProbeLat
    rLon = rProbeLon; rLat = rProbeLat

    c%cType = '1DALL'
    c%gridKind = GRID_GEO_2D

    call c%prepare(iBlock)
    if (.not. c%this_rank_writes) return

    tmp_lon = (rLon*Longitude(iiLon, iBlock) + (1.0 - rLon)*Longitude(iiLon + 1, iBlock))
    tmp_lat = (rLat*Latitude(iiLat, iBlock) + (1.0 - rLat)*Latitude(iiLat + 1, iBlock))
    call c%put_axis('Longitude', tmp_lon)
    call c%put_axis('Latitude', tmp_lat)

    do iAlt = 1, nAlts
      tmp_data(iAlt) = Altitude_GB(iiLon, iiLat, iAlt, iBlock)
    enddo
    call c%put_field('Altitude', tmp_data)

    do iAlt = 1, nAlts
      slice2d = Rho(0:nLons + 1, 0:nLats + 1, iAlt, iBlock)
      tmp_data(iAlt) = inter2d(slice2d, iiLon, iiLat, rLon, rLat)
    enddo
    call c%put_field('Rho', tmp_data)

    do i = 1, nSpeciesTotal
      do iAlt = 1, nAlts
        slice2d = NDensityS(0:nLons + 1, 0:nLats + 1, iAlt, i, iBlock)
        tmp_data(iAlt) = inter2d(slice2d, iiLon, iiLat, rLon, rLat)
      enddo
      call c%put_field('['//trim(cSpecies(i))//']', tmp_data, units='/m3')
    enddo

    do iAlt = 1, nAlts
      slice2d = Temperature(0:nLons + 1, 0:nLats + 1, iAlt, iBlock)* &
                TempUnit(0:nLons + 1, 0:nLats + 1, iAlt)
      tmp_data(iAlt) = inter2d(slice2d, iiLon, iiLat, rLon, rLat)
    enddo
    call c%put_field('Temperature', tmp_data)

    do i = 1, 3
      do iAlt = 1, nAlts
        slice2d = Velocity(0:nLons + 1, 0:nLats + 1, iAlt, i, iBlock)
        tmp_data(iAlt) = inter2d(slice2d, iiLon, iiLat, rLon, rLat)
      enddo
      if (i == 1) call c%put_field('Vn (east)', tmp_data)
      if (i == 2) call c%put_field('Vn (north)', tmp_data)
      if (i == 3) call c%put_field('Vn (up)', tmp_data)
    enddo

    do i = 1, nSpecies
      do iAlt = 1, nAlts
        slice2d = VerticalVelocity(0:nLons + 1, 0:nLats + 1, iAlt, i, iBlock)
        tmp_data(iAlt) = inter2d(slice2d, iiLon, iiLat, rLon, rLat)
      enddo
      call c%put_field('Vn (up,'//trim(cSpecies(i))//')', tmp_data, units='m/s')
    enddo

    do i = 1, nIons
      do iAlt = 1, nAlts
        slice2d = IDensityS(0:nLons + 1, 0:nLats + 1, iAlt, i, iBlock)
        tmp_data(iAlt) = inter2d(slice2d, iiLon, iiLat, rLon, rLat)
      enddo
      call c%put_field('['//trim(cIons(i))//']', tmp_data, units='/m3')
    enddo

    do iAlt = 1, nAlts
      slice2d = eTemperature(0:nLons + 1, 0:nLats + 1, iAlt, iBlock)
      tmp_data(iAlt) = inter2d(slice2d, iiLon, iiLat, rLon, rLat)
    enddo
    call c%put_field('eTemperature', tmp_data)

    do iAlt = 1, nAlts
      slice2d = ITemperature(0:nLons + 1, 0:nLats + 1, iAlt, iBlock)
      tmp_data(iAlt) = inter2d(slice2d, iiLon, iiLat, rLon, rLat)
    enddo
    call c%put_field('iTemperature', tmp_data)

    do i = 1, 3
      do iAlt = 1, nAlts
        slice2d = IVelocity(0:nLons + 1, 0:nLats + 1, iAlt, i, iBlock)
        tmp_data(iAlt) = inter2d(slice2d, iiLon, iiLat, rLon, rLat)
      enddo
      if (i == 1) call c%put_field('Vi (east)', tmp_data)
      if (i == 2) call c%put_field('Vi (north)', tmp_data)
      if (i == 3) call c%put_field('Vi (up)', tmp_data)
    enddo
  end subroutine fill_1dall

  ! ---------------------------------------------------------------------------
  ! fill_0dall: put 0DALL 3D-interpolated point data.
  ! ---------------------------------------------------------------------------
  subroutine fill_0dall(c, iBlock)
    use ModGITM
    use ModEUV, only: HeatingEfficiency_CB
    use ModSources, only: JouleHeating, RadCooling, EuvHeating, Conduction
    use ModPlanet, only: nSpeciesTotal, nSpecies, nIons, cSpecies, cIons
    type(OutputContainer), intent(inout) :: c
    integer, intent(in) :: iBlock
    ! Route A1: probe location comes from module state (set by the driver).
    integer :: iiLon, iiLat, iiAlt
    real :: rLon, rLat, rAlt
    real :: Tmp3d(0:nLons + 1, 0:nLats + 1, 0:nAlts + 1)
    integer :: jAlt, i

    iiLon = iProbeLon; iiLat = iProbeLat; iiAlt = iProbeAlt
    rLon = rProbeLon; rLat = rProbeLat; rAlt = rProbeAlt

    c%cType = '0DALL'
    c%gridKind = GRID_GEO_2D

    call c%prepare(iBlock)
    if (.not. c%this_rank_writes) return

    jAlt = max(min(iiAlt, nAlts), 1)

    call c%put_axis('Longitude', &
               [(rLon*Longitude(iiLon, iBlock) + (1.0 - rLon)*Longitude(iiLon + 1, iBlock))])
    call c%put_axis('Latitude', &
               [(rLat*Latitude(iiLat, iBlock) + (1.0 - rLat)*Latitude(iiLat + 1, iBlock))])
    call c%put_field('Altitude', [rAlt*Altitude_GB(iiLon, iiLat, iiAlt, iBlock) + &
                                (1.0 - rAlt)*Altitude_GB(iiLon + 1, iiLat + 1, iiAlt + 1, iBlock)])

    Tmp3d = Rho(0:nLons + 1, 0:nLats + 1, 0:nAlts + 1, iBlock)
    call c%put_field('Rho', [inter3d(Tmp3d, iiLon, iiLat, iiAlt, rLon, rLat, rAlt)])

    do i = 1, nSpeciesTotal
      Tmp3d = NDensityS(0:nLons + 1, 0:nLats + 1, 0:nAlts + 1, i, iBlock)
      call c%put_field('['//trim(cSpecies(i))//']', &
                       [inter3d(Tmp3d, iiLon, iiLat, iiAlt, rLon, rLat, rAlt)], units='/m3')
    enddo

    Tmp3d = Temperature(0:nLons + 1, 0:nLats + 1, 0:nAlts + 1, iBlock)* &
            TempUnit(0:nLons + 1, 0:nLats + 1, 0:nAlts + 1)
    call c%put_field('Temperature', [inter3d(Tmp3d, iiLon, iiLat, iiAlt, rLon, rLat, rAlt)])

    do i = 1, 3
      Tmp3d = Velocity(0:nLons + 1, 0:nLats + 1, 0:nAlts + 1, i, iBlock)
      if (i == 1) call c%put_field('Vn (east)', [inter3d(Tmp3d, iiLon, iiLat, iiAlt, rLon, rLat, rAlt)])
      if (i == 2) call c%put_field('Vn (north)', [inter3d(Tmp3d, iiLon, iiLat, iiAlt, rLon, rLat, rAlt)])
      if (i == 3) call c%put_field('Vn (up)', [inter3d(Tmp3d, iiLon, iiLat, iiAlt, rLon, rLat, rAlt)])
    enddo

    do i = 1, nSpecies
      Tmp3d = VerticalVelocity(0:nLons + 1, 0:nLats + 1, 0:nAlts + 1, i, iBlock)
      call c%put_field('Vn (up,'//trim(cSpecies(i))//')', &
                 [inter3d(Tmp3d, iiLon, iiLat, iiAlt, rLon, rLat, rAlt)], units='m/s')
    enddo

    do i = 1, nIons
      Tmp3d = IDensityS(0:nLons + 1, 0:nLats + 1, 0:nAlts + 1, i, iBlock)
      call c%put_field('['//trim(cIons(i))//']', [inter3d(Tmp3d, iiLon, iiLat, iiAlt, rLon, rLat, rAlt)], units='/m3')
    enddo

    Tmp3d = eTemperature(0:nLons + 1, 0:nLats + 1, 0:nAlts + 1, iBlock)
    call c%put_field('eTemperature', [inter3d(Tmp3d, iiLon, iiLat, iiAlt, rLon, rLat, rAlt)])
    Tmp3d = iTemperature(0:nLons + 1, 0:nLats + 1, 0:nAlts + 1, iBlock)
    call c%put_field('iTemperature', [inter3d(Tmp3d, iiLon, iiLat, iiAlt, rLon, rLat, rAlt)])

    Tmp3d = IVelocity(0:nLons + 1, 0:nLats + 1, 0:nAlts + 1, iEast_, iBlock)
    call c%put_field('Vi (east)', [inter3d(Tmp3d, iiLon, iiLat, iiAlt, rLon, rLat, rAlt)])
    Tmp3d = IVelocity(0:nLons + 1, 0:nLats + 1, 0:nAlts + 1, iNorth_, iBlock)
    call c%put_field('Vi (north)', [inter3d(Tmp3d, iiLon, iiLat, iiAlt, rLon, rLat, rAlt)])
    Tmp3d = IVelocity(0:nLons + 1, 0:nLats + 1, 0:nAlts + 1, iUp_, iBlock)
    call c%put_field('Vi (up)', [inter3d(Tmp3d, iiLon, iiLat, iiAlt, rLon, rLat, rAlt)])

    do i = 1, nSpecies
      Tmp3d = NDensityS(0:nLons + 1, 0:nLats + 1, 0:nAlts + 1, i, iBlock)/ &
              NDensity(0:nLons + 1, 0:nLats + 1, 0:nAlts + 1, iBlock)
      call c%put_field(trim(cSpecies(i))//' Mixing Ratio', &
                 [inter3d(Tmp3d, iiLon, iiLat, iiAlt, rLon, rLat, rAlt)], units='')
    enddo

    call c%put_field('RadCooling', [RadCooling(1, 1, jAlt, iBlock)*TempUnit(1, 1, jAlt)], units='K/s')
    call c%put_field('EuvHeating', [EuvHeating(1, 1, jAlt, iBlock)*TempUnit(1, 1, jAlt)], units='K/s')
    call c%put_field('Conduction', [per_dt(Conduction(1, 1, jAlt)*TempUnit(1, 1, jAlt))])
    call c%put_field('Heat Balance Total', [(EuvHeating(1, 1, jAlt, iBlock)*TempUnit(1, 1, jAlt) - &
                                           RadCooling(1, 1, jAlt, iBlock)*TempUnit(1, 1, jAlt) + &
                                           per_dt(Conduction(1, 1, jAlt)*TempUnit(1, 1, jAlt)))], units='K/s')
    call c%put_field('Heating Efficiency', [HeatingEfficiency_CB(1, 1, jAlt, iBlock)], units='')
  end subroutine fill_0dall

  ! ---------------------------------------------------------------------------
  ! define_schema_1dglo: register variables for 1DGLO.
  ! ---------------------------------------------------------------------------
  subroutine define_schema_1dglo(c)
    use ModGITM, only: nAlts
    type(OutputContainer), intent(inout) :: c

    c%cType = '1DGLO'
    c%gridKind = GRID_MAG_2D

    call c%define_var('Longitude', units='degrees_east', shape3=[nAlts, 1, 1], is_axis=.true.)
    call c%define_var('Latitude', units='degrees_north', shape3=[nAlts, 1, 1], is_axis=.true.)
    call c%define_var('Altitude', units='m', shape3=[nAlts, 1, 1])
    call c%define_var('6300 A Emission', units='', shape3=[nAlts, 1, 1])
    call c%define_var('PhotoElectronUp', units='', shape3=[nAlts, 1, 1])
    call c%define_var('PhotoElectronDown', units='', shape3=[nAlts, 1, 1])
  end subroutine define_schema_1dglo

  ! ---------------------------------------------------------------------------
  ! fill_1dglo: put 1DGLO stub variables.
  ! ---------------------------------------------------------------------------
  subroutine fill_1dglo(c, iBlock)
    use ModGITM
    type(OutputContainer), intent(inout) :: c
    integer, intent(in) :: iBlock
    real(output_kind) :: tmp_lon(nAlts), tmp_lat(nAlts), tmp_alt(nAlts), zero_arr(nAlts)
    integer :: iAlt

    ! Same story as 3DGLO -- never wired, the puts below are placeholder zeros.
    ! No schema_locked guard is needed here: 1DGLO declares its schema through
    ! define_schema_1dglo, so its nVars is already non-zero at init and
    ! init_output_containers never primes this fill.
    call stop_gitm('1DGLO writes placeholder zeros, not emission data. '// &
                   'Remove 1DGLO from #SAVEPLOT.')

    call c%prepare(iBlock)
    if (.not. c%this_rank_writes) return

    tmp_lon = real(Longitude(1, iBlock), output_kind)
    tmp_lat = real(Latitude(1, iBlock), output_kind)
    call c%put('Longitude', tmp_lon)
    call c%put('Latitude', tmp_lat)

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
  ! define_schema_1dthm: register variables for 1DTHM.
  ! ---------------------------------------------------------------------------
  subroutine define_schema_1dthm(c)
    use ModGITM, only: nAlts
    use ModPlanet, only: nSpeciesTotal, cSpecies
    type(OutputContainer), intent(inout) :: c
    integer :: i

    c%cType = '1DTHM'
    c%gridKind = GRID_MAG_2D

    call c%define_var('Longitude', units='degrees_east', shape3=[nAlts, 1, 1], is_axis=.true.)
    call c%define_var('Latitude', units='degrees_north', shape3=[nAlts, 1, 1], is_axis=.true.)
    call c%define_var('Altitude', units='m', shape3=[nAlts, 1, 1])
    call c%define_var('EUV Heating (K/s)', units='K/s', shape3=[nAlts, 1, 1])
    call c%define_var('Conduction (K/s)', units='K/s', shape3=[nAlts, 1, 1])
    call c%define_var('Molecular Conduction (K/s)', units='K/s', shape3=[nAlts, 1, 1])
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
  ! fill_1dthm: put 1DTHM variables at point (1,1).
  ! ---------------------------------------------------------------------------
  subroutine fill_1dthm(c, iBlock)
    use ModGITM
    use ModSources
    use ModEUV, only: EuvTotal
    use ModPlanet, only: nSpeciesTotal, cSpecies
    type(OutputContainer), intent(inout) :: c
    integer, intent(in) :: iBlock
    real(output_kind) :: tmp_lon(nAlts), tmp_lat(nAlts), tmp_data(nAlts)
    integer :: iAlt, iiAlt, i

    call c%prepare(iBlock)
    if (.not. c%this_rank_writes) return

    tmp_lon = real(Longitude(1, 1), output_kind)
    tmp_lat = real(Latitude(1, 1), output_kind)
    call c%put('Longitude', tmp_lon)
    call c%put('Latitude', tmp_lat)

    do iAlt = 1, nAlts
      tmp_data(iAlt) = real(Altitude_GB(1, 1, iAlt, 1), output_kind)
    enddo
    call c%put('Altitude', tmp_data)

    do iAlt = 1, nAlts
      iiAlt = max(min(iAlt, nAlts), 1)
      tmp_data(iAlt) = real(EuvHeating(1, 1, iiAlt, 1)*TempUnit(1, 1, iiAlt), output_kind)
    enddo
    call c%put('EUV Heating (K/s)', tmp_data)

    do iAlt = 1, nAlts
      iiAlt = max(min(iAlt, nAlts), 1)
      tmp_data(iAlt) = real(per_dt(Conduction(1, 1, iiAlt)*TempUnit(1, 1, iiAlt)), output_kind)
    enddo
    call c%put('Conduction (K/s)', tmp_data)

    do iAlt = 1, nAlts
      iiAlt = max(min(iAlt, nAlts), 1)
      ! MoleConduction from the last (full-Dt) calc_conduction_1d call is a
      ! Kelvin temperature change over the step, so /dt alone gives K/s.
      tmp_data(iAlt) = real(per_dt(MoleConduction(1, 1, iiAlt)), output_kind)
    enddo
    call c%put('Molecular Conduction (K/s)', tmp_data)

    do iAlt = 1, nAlts
      iiAlt = max(min(iAlt, nAlts), 1)
      tmp_data(iAlt) = real(per_dt(ChemicalHeatingRate(1, 1, iiAlt)*TempUnit(1, 1, iiAlt)), output_kind)
    enddo
    call c%put('Chemical Heating (K/s)', tmp_data)

    do iAlt = 1, nAlts
      iiAlt = max(min(iAlt, nAlts), 1)
      tmp_data(iAlt) = real(JouleHeating(1, 1, iiAlt)*TempUnit(1, 1, iiAlt), output_kind)
    enddo
    call c%put('Joule Heating (K/s)', tmp_data)

    do iAlt = 1, nAlts
      iiAlt = max(min(iAlt, nAlts), 1)
      ! NOCooling, not RadCooling: RadCooling is the combined radiative term,
      ! so writing it under this label made the column total radiative cooling
      ! wearing NO's name.  3DTHM already used NOCooling here.
      tmp_data(iAlt) = real(-NOCooling(1, 1, iiAlt)*TempUnit(1, 1, iiAlt), output_kind)
    enddo
    call c%put('NO Cooling (K/s)', tmp_data)

    do iAlt = 1, nAlts
      iiAlt = max(min(iAlt, nAlts), 1)
      tmp_data(iAlt) = real(-OCooling(1, 1, iiAlt)*TempUnit(1, 1, iiAlt), output_kind)
    enddo
    call c%put('O Cooling (K/s)', tmp_data)

    do iAlt = 1, nAlts
      iiAlt = max(min(iAlt, nAlts), 1)
      tmp_data(iAlt) = real(EuvTotal(1, 1, iiAlt, 1), output_kind)
    enddo
    call c%put('Total Abs EUV', tmp_data)

    ! NeutralSourcesTotal/NeutralLossesTotal hold NeutralSources*DtSub from
    ! calc_chemistry.Earth.f90, i.e. a per-step INTEGRAL -- a density, not the
    ! /m3/s these columns are named for -- so /dt is what makes the name true.
    !
    ! The sub-step loop there (`do while (DtTotal < Dt)`) runs exactly ONCE:
    ! DtSub is only ever assigned `Dt - DtTotal`, and DtTotal is only ever
    ! advanced by DtSub, so the first pass sets DtSub = Dt and ends the loop.
    ! There is no chemistry sub-cycling, so this is not a time-weighted mean
    ! over sub-steps -- /dt recovers NeutralSources EXACTLY, the instantaneous
    ! rate. (If sub-cycling is ever restored, sum(DtSub) is still Dt and /dt
    ! becomes the correct mean, so this does not need revisiting.)
    do i = 1, nSpeciesTotal
      do iAlt = 1, nAlts
        iiAlt = max(min(iAlt, nAlts), 1)
        tmp_data(iAlt) = real(per_dt(NeutralSourcesTotal(iiAlt, i)), output_kind)
      enddo
      call c%put('Production Rate '//trim(cSpecies(i)), tmp_data)
    enddo

    do i = 1, nSpeciesTotal
      do iAlt = 1, nAlts
        iiAlt = max(min(iAlt, nAlts), 1)
        tmp_data(iAlt) = real(per_dt(NeutralLossesTotal(iiAlt, i)), output_kind)
      enddo
      call c%put('Loss Rate '//trim(cSpecies(i)), tmp_data)
    enddo
  end subroutine fill_1dthm

  ! ---------------------------------------------------------------------------
  ! define_schema_1dchm: register variables for 1DCHM.
  ! ---------------------------------------------------------------------------
  subroutine define_schema_1dchm(c)
    use ModGITM, only: nAlts
    type(OutputContainer), intent(inout) :: c

    c%cType = '1DCHM'
    c%gridKind = GRID_MAG_2D

    call c%define_var('Longitude', units='degrees_east', shape3=[nAlts, 1, 1], is_axis=.true.)
    call c%define_var('Latitude', units='degrees_north', shape3=[nAlts, 1, 1], is_axis=.true.)
    call c%define_var('Altitude', units='m', shape3=[nAlts, 1, 1])
    ! Same ChemicalHeatingSpecies/Element_Charge/Dt heating rates that 3DCHM
    ! writes; they sum to the 'Chemical Heating' total below. See fill_1dchm.
    call c%define_var('N2+ + e', units='eV/m3/s', shape3=[nAlts, 1, 1])
    call c%define_var('O2+ + e', units='eV/m3/s', shape3=[nAlts, 1, 1])
    call c%define_var('N2+ + O', units='eV/m3/s', shape3=[nAlts, 1, 1])
    call c%define_var('NO+ + e', units='eV/m3/s', shape3=[nAlts, 1, 1])
    call c%define_var('N+ + O2', units='eV/m3/s', shape3=[nAlts, 1, 1])
    call c%define_var('NO + N', units='eV/m3/s', shape3=[nAlts, 1, 1])
    call c%define_var('O+ + O2', units='eV/m3/s', shape3=[nAlts, 1, 1])
    call c%define_var('N + O2', units='eV/m3/s', shape3=[nAlts, 1, 1])
    call c%define_var('O2+ + N', units='eV/m3/s', shape3=[nAlts, 1, 1])
    call c%define_var('O2+ + NO', units='eV/m3/s', shape3=[nAlts, 1, 1])
    call c%define_var('O2+ + N2', units='eV/m3/s', shape3=[nAlts, 1, 1])
    call c%define_var('N2+ + O2', units='eV/m3/s', shape3=[nAlts, 1, 1])
    call c%define_var('N+ + O', units='eV/m3/s', shape3=[nAlts, 1, 1])
    call c%define_var('O+ + N2', units='eV/m3/s', shape3=[nAlts, 1, 1])
    call c%define_var('O(1D) + N2', units='eV/m3/s', shape3=[nAlts, 1, 1])
    call c%define_var('O(1D) + O2', units='eV/m3/s', shape3=[nAlts, 1, 1])
    call c%define_var('O(1D) + O', units='eV/m3/s', shape3=[nAlts, 1, 1])
    call c%define_var('O(1D) + e', units='eV/m3/s', shape3=[nAlts, 1, 1])
    call c%define_var('N(2D) + O2', units='eV/m3/s', shape3=[nAlts, 1, 1])
    call c%define_var('O+(2D) + e', units='eV/m3/s', shape3=[nAlts, 1, 1])
    call c%define_var('N(2D) + O', units='eV/m3/s', shape3=[nAlts, 1, 1])
    call c%define_var('N(2D) + e', units='eV/m3/s', shape3=[nAlts, 1, 1])
    call c%define_var('O+(2D) + N2', units='eV/m3/s', shape3=[nAlts, 1, 1])
    call c%define_var('O+(2P) + e', units='eV/m3/s', shape3=[nAlts, 1, 1])
    call c%define_var('O+(2P) + O', units='eV/m3/s', shape3=[nAlts, 1, 1])
    call c%define_var('O+(2P) + N2', units='eV/m3/s', shape3=[nAlts, 1, 1])
    call c%define_var('Chemical Heating', units='eV/m3/s', shape3=[nAlts, 1, 1])
  end subroutine define_schema_1dchm

  ! ---------------------------------------------------------------------------
  ! fill_1dchm: put 1DCHM variables at point (1,1).
  ! ---------------------------------------------------------------------------
  subroutine fill_1dchm(c, iBlock)
    use ModGITM
    use ModSources
    use ModConstants, only: Element_Charge
    type(OutputContainer), intent(inout) :: c
    integer, intent(in) :: iBlock
    real(output_kind) :: tmp_lon(nAlts), tmp_lat(nAlts), tmp_data(nAlts)
    integer :: iAlt, iiAlt, i

    call c%prepare(iBlock)
    if (.not. c%this_rank_writes) return

    tmp_lon = real(Longitude(1, 1), output_kind)
    tmp_lat = real(Latitude(1, 1), output_kind)
    call c%put('Longitude', tmp_lon)
    call c%put('Latitude', tmp_lat)

    do iAlt = 1, nAlts
      tmp_data(iAlt) = real(Altitude_GB(1, 1, iAlt, 1), output_kind)
    enddo
    call c%put('Altitude', tmp_data)

    do i = 1, nReactions
      do iAlt = 1, nAlts
        iiAlt = max(min(iAlt, nAlts), 1)
        tmp_data(iAlt) = real(per_dt(ChemicalHeatingSpecies(1, 1, iiAlt, i)/Element_Charge), output_kind)
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
      if (i == 14) call c%put('O+ + N2', tmp_data)
      if (i == 15) call c%put('O(1D) + N2', tmp_data)
      if (i == 16) call c%put('O(1D) + O2', tmp_data)
      if (i == 17) call c%put('O(1D) + O', tmp_data)
      if (i == 18) call c%put('O(1D) + e', tmp_data)
      if (i == 19) call c%put('N(2D) + O2', tmp_data)
      if (i == 20) call c%put('O+(2D) + e', tmp_data)
      if (i == 21) call c%put('N(2D) + O', tmp_data)
      if (i == 22) call c%put('N(2D) + e', tmp_data)
      if (i == 23) call c%put('O+(2D) + N2', tmp_data)
      if (i == 24) call c%put('O+(2P) + e', tmp_data)
      if (i == 25) call c%put('O+(2P) + O', tmp_data)
      if (i == 26) call c%put('O+(2P) + N2', tmp_data)
    enddo

    do iAlt = 1, nAlts
      iiAlt = max(min(iAlt, nAlts), 1)
      tmp_data(iAlt) = real(per_dt(ChemicalHeatingRate(1, 1, iiAlt)* &
                            cp(1, 1, iiAlt, 1)* &
                            Rho(1, 1, iiAlt, 1)* &
                            TempUnit(1, 1, iiAlt)/Element_Charge), output_kind)
    enddo
    call c%put('Chemical Heating', tmp_data)
  end subroutine fill_1dchm

  ! ---------------------------------------------------------------------------
  ! fill_1dnew: put 1DNEW interpolated vertical column variables.
  ! ---------------------------------------------------------------------------
  subroutine fill_1dnew(c, iBlock)
    use ModGITM
    use ModEUV, only: Sza
    use ModPlanet, only: nSpeciesTotal, nSpecies, nIons, cSpecies, cIons
    type(OutputContainer), intent(inout) :: c
    integer, intent(in) :: iBlock
    ! Route A1: probe location comes from module state (set by the driver).
    integer :: iiLon, iiLat
    real :: rLon, rLat
    real :: tmp_lon(nAlts), tmp_lat(nAlts), tmp_lt(nAlts), tmp_sza(nAlts), tmp_data(nAlts)
    real :: slice2d(0:nLons + 1, 0:nLats + 1)
    integer :: iAlt, i

    iiLon = iProbeLon; iiLat = iProbeLat
    rLon = rProbeLon; rLat = rProbeLat

    c%cType = '1DNEW'
    c%gridKind = GRID_GEO_2D

    call c%prepare(iBlock)
    if (.not. c%this_rank_writes) return

    tmp_lon = (rLon*Longitude(iiLon, iBlock) + (1.0 - rLon)*Longitude(iiLon + 1, iBlock))
    tmp_lat = (rLat*Latitude(iiLat, iBlock) + (1.0 - rLat)*Latitude(iiLat + 1, iBlock))
    tmp_lt = LocalTime(iiLon)
    tmp_sza = rLon*rLat*Sza(iiLon, iiLat, iBlock) + &
                   (1.0 - rLon)*rLat*Sza(iiLon + 1, iiLat, iBlock) + &
                   rLon*(1.0 - rLat)*Sza(iiLon, iiLat + 1, iBlock) + &
                   (1.0 - rLon)*(1.0 - rLat)*Sza(iiLon + 1, iiLat + 1, iBlock)

    call c%put_axis('Longitude', tmp_lon)
    call c%put_field('Local Time', tmp_lt, units='hr')
    call c%put_field('Latitude', tmp_lat)
    call c%put_field('Solar Zenith Angle', tmp_sza, units='rad')

    do iAlt = 1, nAlts
      tmp_data(iAlt) = Altitude_GB(iiLon, iiLat, iAlt, iBlock)
    enddo
    call c%put_field('Altitude', tmp_data)

    do iAlt = 1, nAlts
      slice2d = Rho(0:nLons + 1, 0:nLats + 1, iAlt, iBlock)
      tmp_data(iAlt) = inter2d(slice2d, iiLon, iiLat, rLon, rLat)
    enddo
    call c%put_field('Rho', tmp_data)

    do i = 1, nSpeciesTotal
      do iAlt = 1, nAlts
        slice2d = NDensityS(0:nLons + 1, 0:nLats + 1, iAlt, i, iBlock)
        tmp_data(iAlt) = inter2d(slice2d, iiLon, iiLat, rLon, rLat)
      enddo
      call c%put_field('['//trim(cSpecies(i))//']', tmp_data, units='/m3')
    enddo

    do iAlt = 1, nAlts
      slice2d = Temperature(0:nLons + 1, 0:nLats + 1, iAlt, iBlock)* &
                TempUnit(0:nLons + 1, 0:nLats + 1, iAlt)
      tmp_data(iAlt) = inter2d(slice2d, iiLon, iiLat, rLon, rLat)
    enddo
    call c%put_field('Temperature', tmp_data)

    do i = 1, 3
      do iAlt = 1, nAlts
        slice2d = Velocity(0:nLons + 1, 0:nLats + 1, iAlt, i, iBlock)
        tmp_data(iAlt) = inter2d(slice2d, iiLon, iiLat, rLon, rLat)
      enddo
      if (i == 1) call c%put_field('Vn (east)', tmp_data)
      if (i == 2) call c%put_field('Vn (north)', tmp_data)
      if (i == 3) call c%put_field('Vn (up)', tmp_data)
    enddo

    do i = 1, nSpecies
      do iAlt = 1, nAlts
        slice2d = VerticalVelocity(0:nLons + 1, 0:nLats + 1, iAlt, i, iBlock)
        tmp_data(iAlt) = inter2d(slice2d, iiLon, iiLat, rLon, rLat)
      enddo
      call c%put_field('Vn (up,'//trim(cSpecies(i))//')', tmp_data, units='m/s')
    enddo

    do i = 1, nIons
      do iAlt = 1, nAlts
        slice2d = IDensityS(0:nLons + 1, 0:nLats + 1, iAlt, i, iBlock)
        tmp_data(iAlt) = inter2d(slice2d, iiLon, iiLat, rLon, rLat)
      enddo
      call c%put_field('['//trim(cIons(i))//']', tmp_data, units='/m3')
    enddo

    do iAlt = 1, nAlts
      slice2d = eTemperature(0:nLons + 1, 0:nLats + 1, iAlt, iBlock)
      tmp_data(iAlt) = inter2d(slice2d, iiLon, iiLat, rLon, rLat)
    enddo
    call c%put_field('eTemperature', tmp_data)

    do iAlt = 1, nAlts
      slice2d = ITemperature(0:nLons + 1, 0:nLats + 1, iAlt, iBlock)
      tmp_data(iAlt) = inter2d(slice2d, iiLon, iiLat, rLon, rLat)
    enddo
    call c%put_field('iTemperature', tmp_data)

    do i = 1, 3
      do iAlt = 1, nAlts
        slice2d = IVelocity(0:nLons + 1, 0:nLats + 1, iAlt, i, iBlock)
        tmp_data(iAlt) = inter2d(slice2d, iiLon, iiLat, rLon, rLat)
      enddo
      if (i == 1) call c%put_field('Vi (east)', tmp_data)
      if (i == 2) call c%put_field('Vi (north)', tmp_data)
      if (i == 3) call c%put_field('Vi (up)', tmp_data)
    enddo

    do i = 1, nSpecies
      do iAlt = 1, nAlts
        slice2d = NDensityS(0:nLons + 1, 0:nLats + 1, iAlt, i, iBlock)/ &
                  NDensity(0:nLons + 1, 0:nLats + 1, iAlt, iBlock)
        tmp_data(iAlt) = inter2d(slice2d, iiLon, iiLat, rLon, rLat)
      enddo
      call c%put_field(trim(cSpecies(i))//' Mixing Ratio', tmp_data, units='')
    enddo
  end subroutine fill_1dnew

  ! ---------------------------------------------------------------------------
  ! define_schema_3dusr: define variables for the 3DUSR container.
  ! ---------------------------------------------------------------------------
  subroutine define_schema_3dusr(c)
    use ModSizeGitm, only: nLons, nLats, nAlts
    type(OutputContainer), intent(inout) :: c
    integer :: i

    c%cType = '3DUSR'
    c%gridKind = GRID_GEO_3D
    c%nGhostCells = 0

    call c%define_var('Longitude', units='degrees_east', shape3=[nLons, 1, 1], is_axis=.true.)
    call c%define_var('Latitude', units='degrees_north', shape3=[nLats, 1, 1], is_axis=.true.)
    call c%define_var('Altitude', units='m', shape3=[nLons, nLats, nAlts])
    do i = 1, nUsrVars3D
      call c%define_var(trim(UsrVars3D(i)%name), units=trim(UsrVars3D(i)%units), &
                        shape3=[nLons, nLats, nAlts], longName=trim(UsrVars3D(i)%longName))
    enddo
  end subroutine define_schema_3dusr

  ! ---------------------------------------------------------------------------
  ! fill_3dusr: put 3DUSR variables into the container.
  ! ---------------------------------------------------------------------------
  subroutine fill_3dusr(c, iBlock)
    use ModGITM
    use ModUserGITM, only: UserData3D
    type(OutputContainer), intent(inout) :: c
    integer, intent(in) :: iBlock
    integer :: i

    call c%prepare(iBlock)
    if (.not. c%this_rank_writes) return

    call c%put('Longitude', real(Longitude(1:nLons, iBlock), output_kind))
    call c%put('Latitude', real(Latitude(1:nLats, iBlock), output_kind))
    call c%put('Altitude', real(Altitude_GB(1:nLons, 1:nLats, 1:nAlts, iBlock), output_kind))
    do i = 1, nUsrVars3D
      call c%put(trim(UsrVars3D(i)%name), real(UserData3D(1:nLons, 1:nLats, 1:nAlts, i, iBlock), output_kind))
    enddo
  end subroutine fill_3dusr

  ! ---------------------------------------------------------------------------
  ! define_schema_2dusr: define variables for the 2DUSR container.
  ! ---------------------------------------------------------------------------
  subroutine define_schema_2dusr(c)
    use ModGITM, only: nLons, nLats
    type(OutputContainer), intent(inout) :: c
    integer :: i

    c%cType = '2DUSR'
    c%gridKind = GRID_GEO_2D
    c%nGhostCells = 0

    call c%define_var('Longitude', units='degrees_east', shape3=[nLons, 1, 1], is_axis=.true.)
    call c%define_var('Latitude', units='degrees_north', shape3=[nLats, 1, 1], is_axis=.true.)
    call c%define_var('Altitude', units='m', shape3=[nLons, nLats, 1])
    do i = 1, nUsrVars2D
      call c%define_var(trim(UsrVars2D(i)%name), units=trim(UsrVars2D(i)%units), &
                        shape3=[nLons, nLats, 1], longName=trim(UsrVars2D(i)%longName))
    enddo
  end subroutine define_schema_2dusr

  ! ---------------------------------------------------------------------------
  ! fill_2dusr: put 2DUSR variables into the container.
  ! ---------------------------------------------------------------------------
  subroutine fill_2dusr(c, iBlock)
    use ModGITM
    use ModUserGITM, only: UserData2D
    type(OutputContainer), intent(inout) :: c
    integer, intent(in) :: iBlock
    integer :: i

    call c%prepare(iBlock)
    if (.not. c%this_rank_writes) return

    call c%put('Longitude', real(Longitude(1:nLons, iBlock), output_kind))
    call c%put('Latitude', real(Latitude(1:nLats, iBlock), output_kind))
    call c%put('Altitude', real(Altitude_GB(1:nLons, 1:nLats, 1, iBlock), output_kind))
    do i = 1, nUsrVars2D
      call c%put(trim(UsrVars2D(i)%name), real(UserData2D(1:nLons, 1:nLats, 1, i, iBlock), output_kind))
    enddo
  end subroutine fill_2dusr

  ! ---------------------------------------------------------------------------
  ! define_schema_1dusr: define variables for the 1DUSR container.
  ! ---------------------------------------------------------------------------
  subroutine define_schema_1dusr(c)
    use ModGITM, only: nAlts
    type(OutputContainer), intent(inout) :: c
    integer :: i

    c%cType = '1DUSR'
    c%gridKind = GRID_GEO_2D
    c%nGhostCells = 0

    call c%define_var('Longitude', units='degrees_east', shape3=[nAlts, 1, 1], is_axis=.true.)
    call c%define_var('Latitude', units='degrees_north', shape3=[nAlts, 1, 1], is_axis=.true.)
    call c%define_var('Altitude', units='m', shape3=[nAlts, 1, 1])
    do i = 1, nUsrVars1D
      call c%define_var(trim(UsrVars1D(i)%name), units=trim(UsrVars1D(i)%units), &
                        shape3=[nAlts, 1, 1], longName=trim(UsrVars1D(i)%longName))
    enddo
  end subroutine define_schema_1dusr

  ! ---------------------------------------------------------------------------
  ! fill_1dusr: put 1DUSR variables into the container.
  ! ---------------------------------------------------------------------------
  subroutine fill_1dusr(c, iBlock)
    use ModGITM
    use ModUserGITM, only: UserData1D
    type(OutputContainer), intent(inout) :: c
    integer, intent(in) :: iBlock
    ! Route A1: probe location comes from module state (set by the driver).
    integer :: iiLon, iiLat
    real :: rLon, rLat
    real(output_kind) :: tmp(nAlts)
    integer :: i, iAlt

    iiLon = iProbeLon; iiLat = iProbeLat
    rLon = rProbeLon; rLat = rProbeLat

    call c%prepare(iBlock)
    if (.not. c%this_rank_writes) return

    tmp = real((rLon*Longitude(iiLon, iBlock) + (1.0 - rLon)*Longitude(iiLon + 1, iBlock)), output_kind)
    call c%put('Longitude', tmp)
    tmp = real((rLat*Latitude(iiLat, iBlock) + (1.0 - rLat)*Latitude(iiLat + 1, iBlock)), output_kind)
    call c%put('Latitude', tmp)
    do iAlt = 1, nAlts
      tmp(iAlt) = real(Altitude_GB(iiLon, iiLat, iAlt, iBlock), output_kind)
    enddo
    call c%put('Altitude', tmp)
    do i = 1, nUsrVars1D
      do iAlt = 1, nAlts
        tmp(iAlt) = real(UserData1D(iiLon, iiLat, iAlt, i), output_kind)
      enddo
      call c%put(trim(UsrVars1D(i)%name), tmp)
    enddo
  end subroutine fill_1dusr

  ! ---------------------------------------------------------------------------
  ! define_schema_0dusr: define variables for the 0DUSR container.
  ! ---------------------------------------------------------------------------
  subroutine define_schema_0dusr(c)
    type(OutputContainer), intent(inout) :: c
    integer :: i

    c%cType = '0DUSR'
    c%gridKind = GRID_GEO_2D
    c%nGhostCells = 0

    call c%define_var('Longitude', units='degrees_east', shape3=[1, 1, 1], is_axis=.true.)
    call c%define_var('Latitude', units='degrees_north', shape3=[1, 1, 1], is_axis=.true.)
    call c%define_var('Altitude', units='m', shape3=[1, 1, 1])
    do i = 1, nUsrVars0D
      call c%define_var(trim(UsrVars0D(i)%name), units=trim(UsrVars0D(i)%units), &
                        shape3=[1, 1, 1], longName=trim(UsrVars0D(i)%longName))
    enddo
  end subroutine define_schema_0dusr

  ! ---------------------------------------------------------------------------
  ! fill_0dusr: put 0DUSR variables into the container.
  ! ---------------------------------------------------------------------------
  subroutine fill_0dusr(c, iBlock)
    use ModGITM
    use ModUserGITM, only: UserData1D
    type(OutputContainer), intent(inout) :: c
    integer, intent(in) :: iBlock
    ! Route A1: probe location comes from module state (set by the driver).
    integer :: iiLon, iiLat, iiAlt
    real :: rLon, rLat, rAlt
    integer :: i

    iiLon = iProbeLon; iiLat = iProbeLat; iiAlt = iProbeAlt
    rLon = rProbeLon; rLat = rProbeLat; rAlt = rProbeAlt

    call c%prepare(iBlock)
    if (.not. c%this_rank_writes) return

    call c%put('Longitude', &
               real((rLon*Longitude(iiLon, iBlock) + (1.0 - rLon)*Longitude(iiLon + 1, iBlock)), &
                    output_kind))
    call c%put('Latitude', &
               real((rLat*Latitude(iiLat, iBlock) + (1.0 - rLat)*Latitude(iiLat + 1, iBlock)), &
                    output_kind))
    call c%put('Altitude', real(rAlt*Altitude_GB(iiLon, iiLat, iiAlt, iBlock) + &
                                (1.0 - rAlt)*Altitude_GB(iiLon + 1, iiLat + 1, iiAlt + 1, iBlock), output_kind))
    do i = 1, nUsrVars0D
      call c%put(trim(UsrVars0D(i)%name), real(UserData1D(iiLon, iiLat, iiAlt, i), output_kind))
    enddo
  end subroutine fill_0dusr

  ! ---------------------------------------------------------------------------
  ! per_dt: x/Dt, or 0 when Dt is 0.  The forced output before the first step
  ! runs with Dt = 0 and nothing accumulated yet, so the per-step rates are
  ! written as zeros there instead of 0/0 (which traps in -debug builds).
  ! For Dt > 0 this is exactly x/Dt, element by element.
  ! ---------------------------------------------------------------------------
  elemental real function per_dt(x)
    use ModGITM, only: Dt
    real, intent(in) :: x

    if (Dt > 0.0) then
      per_dt = x/Dt
    else
      per_dt = 0.0
    endif
  end function per_dt

  ! ---------------------------------------------------------------------------
  ! inter2d: Bilinear interpolation on a 2D grid slice.
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
  ! inter3d: Trilinear interpolation on a 3D grid.
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
