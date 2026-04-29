! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE
!
! ModOutputRegistry.f90
!
! Variable registry for GITM output types. Each output type registers its
! variables in a single place, eliminating the duplication between
! nvars_to_write calculations and header-writing code in the old system.
!
! Usage:
!   call init_output_registry()     ! once at startup
!   idx = find_output_type('3DALL') ! lookup by code
!   info = RegisteredTypes(idx)     ! get type descriptor
!   ! info%nVars, info%vars(:) now available for headers and gather

module ModOutputRegistry
  use ModInputs, only: nMaxOutputTypes

  implicit none

  integer, parameter :: MaxNameLen = 60
  integer, parameter :: MaxUnitsLen = 30
  integer, parameter :: MaxVarsPerType = 200

  type :: OutputVar
    character(len=MaxNameLen) :: name = ''
    character(len=MaxUnitsLen) :: units = ''
    character(len=MaxNameLen) :: longName = ''
    character(len=MaxNameLen) :: shortName = ''
  end type OutputVar

  type :: OutputTypeInfo
    character(len=5) :: code = ''
    integer :: nDims = 3
    integer :: nGhostCells = 2
    integer :: nVars = 0
    logical :: isUserType = .false.
    logical :: usesMagGrid = .false.     ! true if on magnetic (not geographic) coordinates
    logical :: isRegional = .false.      ! true if only blocks in specified region write
    logical :: needsInterp = .false.
    logical :: needsVtec = .false.
    type(OutputVar) :: vars(MaxVarsPerType)
  end type OutputTypeInfo

  type(OutputTypeInfo), target :: RegisteredTypes(nMaxOutputTypes)
  integer :: nRegisteredTypes = 0

contains

  ! -----------------------------------------------------------------
  ! Add a variable to an OutputTypeInfo descriptor
  ! -----------------------------------------------------------------
  !> Long and short names are optional and only written to netcdf's.
  subroutine add_var(info, varName, units, shortName, longName)
    type(OutputTypeInfo), intent(inout) :: info
    character(len=*), intent(in) :: varName
    character(len=*), intent(in) :: units
    character(len=*), intent(in), optional :: shortName, longName

    info%nVars = info%nVars + 1
    if (info%nVars > MaxVarsPerType) &
      call stop_gitm("ERROR: too many vars in output type "//info%code)
    info%vars(info%nVars)%name = varName
    if (present(shortName)) then
      info%vars(info%nVars)%shortName = shortName
    else
      info%vars(info%nVars)%shortName = varName
    endif
    if (present(longName)) then
      info%vars(info%nVars)%longName = longName
    else
      info%vars(info%nVars)%longName = varName
    endif
    info%vars(info%nVars)%units = units
  end subroutine add_var

  ! -----------------------------------------------------------------
  ! Find an output type by its 5-character code. Returns index into
  ! RegisteredTypes, or -1 if not found.
  ! -----------------------------------------------------------------
  integer function find_output_type(code)
    character(len=*), intent(in) :: code
    integer :: i

    do i = 1, nRegisteredTypes
      if (trim(RegisteredTypes(i)%code) == trim(code)) then
        find_output_type = i
        return
      endif
    enddo
    find_output_type = -1
  end function find_output_type

  ! -----------------------------------------------------------------
  ! Get a pointer to a new OutputTypeInfo slot for registration
  ! -----------------------------------------------------------------
  subroutine new_output_type(code, nDims, nGhostCells, info)
    character(len=*), intent(in) :: code
    integer, intent(in) :: nDims, nGhostCells
    type(OutputTypeInfo), pointer, intent(out) :: info

    nRegisteredTypes = nRegisteredTypes + 1
    info => RegisteredTypes(nRegisteredTypes)
    info%code = code
    info%nDims = nDims
    info%nGhostCells = nGhostCells
    info%nVars = 0
    info%isUserType = .false.
    info%usesMagGrid = .false.
    info%needsInterp = .false.
    info%needsVtec = .false.
  end subroutine new_output_type

  ! -----------------------------------------------------------------
  ! Geographic lon/lat/alt coordinates
  ! -----------------------------------------------------------------
  subroutine add_coord_vars(info)
    type(OutputTypeInfo), intent(inout) :: info
    call add_var(info, 'Longitude', 'rad')
    call add_var(info, 'Latitude', 'rad')
    call add_var(info, 'Altitude', 'm')
  end subroutine add_coord_vars

  ! =================================================================
  ! Master initialization: register all output types
  ! =================================================================
  subroutine init_output_registry()

    nRegisteredTypes = 0

    call register_3dall()
    call register_3dlst()
    call register_3dneu()
    call register_3dion()
    call register_3dthm()
    call register_3dchm()
    call register_3dusr()
    call register_3dglo()
    call register_3dmag()
    call register_3dhme()
    call register_3dmoh()
    call register_3dmov()
    call register_2dgel()
    call register_2dmel()
    call register_2dusr()
    call register_2dtec()
    call register_2danc()
    call register_2dhme()
    call register_1dall()
    call register_0dall()
    call register_1dglo()
    call register_1dthm()
    call register_1dnew()
    call register_1dchm()
    call register_1dcms()
    call register_1dusr()
    call register_0dusr()

  end subroutine init_output_registry

  ! =================================================================
  ! Individual registration subroutines
  ! =================================================================

  subroutine register_3dall()
    use ModPlanet, only: nSpeciesTotal, nSpecies, nIons, cSpecies, cIons
    type(OutputTypeInfo), pointer :: info
    integer :: i

    call new_output_type('3DALL', 3, 2, info)
    call add_coord_vars(info)
    call add_var(info, 'Rho', 'kg/m3')
    do i = 1, nSpeciesTotal
      call add_var(info, '['//trim(cSpecies(i))//']', '/m3')
    enddo
    call add_var(info, 'Temperature', 'K', 'T', 'neutral_temperature')
    call add_var(info, 'V!Dn!N (east)', 'm/s', 'Ve', 'eastward_neutral_velocity')
    call add_var(info, 'V!Dn!N (north)', 'm/s', 'Vn', 'northward_neutral_velocity')
    call add_var(info, 'V!Dn!N (up)', 'm/s', 'Vv', 'vertical_neutral_velocity')
    do i = 1, nSpecies
      call add_var(info, 'V!Dn!N (up,'//trim(cSpecies(i))//')', 'm/s', &
                   'Vv_'//trim(cSpecies(i)), 'vertical_'//trim(cSpecies(i))//'velocity')
    enddo
    do i = 1, nIons
      call add_var(info, '['//trim(cIons(i))//']', '/m3')
    enddo
    call add_var(info, 'eTemperature', 'K')
    call add_var(info, 'iTemperature', 'K')
    call add_var(info, 'V!Di!N (east)', 'm/s')
    call add_var(info, 'V!Di!N (north)', 'm/s')
    call add_var(info, 'V!Di!N (up)', 'm/s')
  end subroutine register_3dall

  subroutine register_3dlst()
    use ModPlanet, only: nSpeciesTotal, nIons, cSpecies, cIons
    use ModInputs, only: iRhoOutputList, iNeutralDensityOutputList, &
                         iNeutralWindOutputList, iIonDensityOutputList, &
                         iIonWindOutputList, iTemperatureOutputList
    type(OutputTypeInfo), pointer :: info
    integer :: i
    character(len=20) :: windNames(3)

    call new_output_type('3DLST', 3, 2, info)
    call add_coord_vars(info)

    if (iRhoOutputList) call add_var(info, 'Rho (kg/m3)', 'kg/m3')
    do i = 1, nSpeciesTotal
      if (iNeutralDensityOutputList(i)) &
        call add_var(info, '['//trim(cSpecies(i))//'] (/m3)', '/m3')
    enddo

    windNames = (/'Vn (east) (m/s) ', 'Vn (north) (m/s)', 'Vn (up) (m/s)   '/)
    do i = 1, 3
      if (iNeutralWindOutputList(i)) &
        call add_var(info, trim(windNames(i)), 'm/s')
    enddo

    do i = 1, nIons
      if (iIonDensityOutputList(i)) &
        call add_var(info, '['//trim(cIons(i))//'] (/m3)', '/m3')
    enddo

    windNames = (/'Vi (east) (m/s) ', 'Vi (north) (m/s)', 'Vi (up) (m/s)   '/)
    do i = 1, 3
      if (iIonWindOutputList(i)) &
        call add_var(info, trim(windNames(i)), 'm/s')
    enddo

    if (iTemperatureOutputList(1)) &
      call add_var(info, 'Neutral Temperature (K)', 'K')
    if (iTemperatureOutputList(2)) &
      call add_var(info, 'Ion Temperature (K)', 'K')
    if (iTemperatureOutputList(3)) &
      call add_var(info, 'Electron Temperature (K)', 'K')
  end subroutine register_3dlst

  subroutine register_3dneu()
    use ModPlanet, only: nSpeciesTotal, nSpecies, cSpecies
    type(OutputTypeInfo), pointer :: info
    integer :: i

    call new_output_type('3DNEU', 3, 2, info)
    call add_coord_vars(info)
    call add_var(info, 'Rho', 'kg/m3')
    do i = 1, nSpeciesTotal
      call add_var(info, '['//trim(cSpecies(i))//']', '/m3')
    enddo
    call add_var(info, 'Temperature', 'K')
    call add_var(info, 'V!Dn!N (east)', 'm/s')
    call add_var(info, 'V!Dn!N (north)', 'm/s')
    call add_var(info, 'V!Dn!N (up)', 'm/s')
    do i = 1, nSpecies
      call add_var(info, 'V!Dn!N (up,'//trim(cSpecies(i))//')', 'm/s')
    enddo
  end subroutine register_3dneu

  subroutine register_3dion()
    use ModPlanet, only: nIons, cIons
    type(OutputTypeInfo), pointer :: info
    integer :: i

    call new_output_type('3DION', 3, 2, info)
    call add_coord_vars(info)
    do i = 1, nIons
      call add_var(info, '['//trim(cIons(i))//']', '/m3')
    enddo
    call add_var(info, 'eTemperature', 'K')
    call add_var(info, 'iTemperature', 'K')
    call add_var(info, 'V!Di!N (east)', 'm/s')
    call add_var(info, 'V!Di!N (north)', 'm/s')
    call add_var(info, 'V!Di!N (up)', 'm/s')
    call add_var(info, 'Ed1', 'V/m')
    call add_var(info, 'Ed2', 'V/m')
    call add_var(info, 'Je1', 'A/m2')
    call add_var(info, 'Je2', 'A/m2')
    call add_var(info, 'Magnetic Latitude', 'deg')
    call add_var(info, 'Magnetic Longitude', 'deg')
    call add_var(info, 'B.F. East', 'T')
    call add_var(info, 'B.F. North', 'T')
    call add_var(info, 'B.F. Vertical', 'T')
    call add_var(info, 'B.F. Magnitude', 'T')
    call add_var(info, 'Potential', 'V')
    call add_var(info, 'E.F. East', 'V/m')
    call add_var(info, 'E.F. North', 'V/m')
    call add_var(info, 'E.F. Vertical', 'V/m')
    call add_var(info, 'E.F. Magnitude', 'V/m')
    call add_var(info, 'IN Collision Freq', 'Hz')
    call add_var(info, 'PressGrad (east)', 'm/s2')
    call add_var(info, 'PressGrad (north)', 'm/s2')
    call add_var(info, 'PressGrad (up)', 'm/s2')
  end subroutine register_3dion

  subroutine register_3dthm()
    type(OutputTypeInfo), pointer :: info

    call new_output_type('3DTHM', 3, 2, info)
    call add_coord_vars(info)
    call add_var(info, 'EUV Heating (K/s)', 'K/s')
    call add_var(info, 'Conduction (K/s)', 'K/s')
    call add_var(info, 'Molecular Conduction (K/s)', 'K/s')
    call add_var(info, 'Eddy Conduction (K/s)', 'K/s')
    call add_var(info, 'Eddy Adiabatic Conduction (K/s)', 'K/s')
    call add_var(info, 'Chemical Heating (K/s)', 'K/s')
    call add_var(info, 'Joule Heating (K/s)', 'K/s')
    call add_var(info, 'NO Cooling (K/s)', 'K/s')
    call add_var(info, 'O Cooling (K/s)', 'K/s')
    call add_var(info, 'Total Abs EUV', 'W/m2')
    call add_var(info, 'Cp', 'J/kg/K')
    call add_var(info, 'Rho', 'kg/m3')
    call add_var(info, 'E-Field Mag', 'V/m')
    call add_var(info, 'Sigma Ped', 'S/m')
    call add_var(info, 'Ionization Rate O_3P', '/s')
    call add_var(info, 'Ionization Rate O2', '/s')
    call add_var(info, 'Ionization Rate N2', '/s')
  end subroutine register_3dthm

  subroutine register_3dchm()
    type(OutputTypeInfo), pointer :: info

    call new_output_type('3DCHM', 3, 2, info)
    call add_coord_vars(info)
    call add_var(info, 'N!D2!U+!N + e', '/m3/s')
    call add_var(info, 'O!D2!U+!N + e', '/m3/s')
    call add_var(info, 'N!D2!U+!N + O', '/m3/s')
    call add_var(info, 'NO!U+!N + e', '/m3/s')
    call add_var(info, 'N!U+!N + O!D2!N', '/m3/s')
    call add_var(info, 'NO + N', '/m3/s')
    call add_var(info, 'O!U+!N + O!D2!N', '/m3/s')
    call add_var(info, 'N + O!D2!N', '/m3/s')
    call add_var(info, 'O!D2!U+!N + N', '/m3/s')
    call add_var(info, 'O!D2!U+!N + NO', '/m3/s')
    call add_var(info, 'O!D2!U+!N + N2', '/m3/s')
    call add_var(info, 'N!D2!U+!N + O!D2!N', '/m3/s')
    call add_var(info, 'N!U+!N + O', '/m3/s')
    call add_var(info, 'O!+!N + N!D2!N', '/m3/s')
    call add_var(info, 'O(1D) + N!D2!N', '/m3/s')
    call add_var(info, 'O(1D) + O!D2!N', '/m3/s')
    call add_var(info, 'O(1D) + O', '/m3/s')
    call add_var(info, 'O(1D) + e', '/m3/s')
    call add_var(info, 'N(2D) + O!D2!N', '/m3/s')
    call add_var(info, 'O!U+!N(2D)+e', '/m3/s')
    call add_var(info, 'N(2D) + O', '/m3/s')
    call add_var(info, 'N(2D) + e', '/m3/s')
    call add_var(info, 'O!U+!N(2D + N!D2!N', '/m3/s')
    call add_var(info, 'O!U+!N(2P) + e', '/m3/s')
    call add_var(info, 'O!U+!N(2P) + O', '/m3/s')
    call add_var(info, 'O!U+!N(2P) + N!D2!N', '/m3/s')
    call add_var(info, 'Chemical Heating Rate', 'K/s')
  end subroutine register_3dchm

  subroutine register_3dusr()
    type(OutputTypeInfo), pointer :: info
    call new_output_type('3DUSR', 3, 2, info)
    info%isUserType = .true.
  end subroutine register_3dusr

  subroutine register_3dglo()
    type(OutputTypeInfo), pointer :: info

    call new_output_type('3DGLO', 3, 2, info)
    call add_coord_vars(info)
    call add_var(info, '6300 A Emission', '')
    call add_var(info, 'PhotoElectronUp', '')
    call add_var(info, 'PhotoElectronDown', '')
  end subroutine register_3dglo

  subroutine register_3dmag()
    type(OutputTypeInfo), pointer :: info

    call new_output_type('3DMAG', 3, 2, info)
    call add_coord_vars(info)
    call add_var(info, 'Magnetic Latitude', 'deg')
    call add_var(info, 'Magnetic Longitude', 'deg')
    call add_var(info, 'B.F. East', 'T')
    call add_var(info, 'B.F. North', 'T')
    call add_var(info, 'B.F. Vertical', 'T')
    call add_var(info, 'B.F. Magnitude', 'T')
  end subroutine register_3dmag

  subroutine register_3dhme()
    use ModPlanet, only: nSpeciesTotal, nSpecies, nIons, cSpecies, cIons
    type(OutputTypeInfo), pointer :: info
    integer :: i

    call new_output_type('3DHME', 3, 0, info)
    call add_coord_vars(info)
    call add_var(info, 'Rho', 'kg/m3')
    do i = 1, nSpeciesTotal
      call add_var(info, '['//trim(cSpecies(i))//']', '/m3')
    enddo
    call add_var(info, 'Temperature', 'K')
    call add_var(info, 'V!Dn!N (east)', 'm/s')
    call add_var(info, 'V!Dn!N (north)', 'm/s')
    call add_var(info, 'V!Dn!N (up)', 'm/s')
    do i = 1, nSpecies
      call add_var(info, 'V!Dn!N (up,'//trim(cSpecies(i))//')', 'm/s')
    enddo
    do i = 1, nIons
      call add_var(info, '['//trim(cIons(i))//']', '/m3')
    enddo
    call add_var(info, 'eTemperature', 'K')
    call add_var(info, 'iTemperature', 'K')
    call add_var(info, 'V!Di!N (east)', 'm/s')
    call add_var(info, 'V!Di!N (north)', 'm/s')
    call add_var(info, 'V!Di!N (up)', 'm/s')
    call add_var(info, 'PhotoElectron Heating', 'K')
    call add_var(info, 'Joule Heating', 'K')
    call add_var(info, 'Specific Heat', 'J/kg/K')
    call add_var(info, 'Magnetic Latitude', 'deg')
    call add_var(info, 'Magnetic Longitude', 'deg')
    call add_var(info, 'B.F. East', 'T')
    call add_var(info, 'B.F. North', 'T')
    call add_var(info, 'B.F. Vertical', 'T')
    call add_var(info, 'B.F. Magnitude', 'T')
    call add_var(info, 'Potential', 'V')
    call add_var(info, 'PotentialY', 'V')
    call add_var(info, 'E.F. East', 'V/m')
    call add_var(info, 'E.F. North', 'V/m')
    call add_var(info, 'E.F. Vertical', 'V/m')
  end subroutine register_3dhme

  subroutine register_3dmoh()
    type(OutputTypeInfo), pointer :: info

    call new_output_type('3DMOH', 3, 2, info)
    call add_coord_vars(info)
    call add_var(info, 'Rho (kg/m3)', 'kg/m3')
    call add_var(info, 'Vn_east (m/s)', 'm/s')
    call add_var(info, 'Vn_north (m/s)', 'm/s')
    call add_var(info, 'Visc_Ve_rshear (m/s2)', 'm/s2')
    call add_var(info, 'Visc_Vn_rshear (m/s2)', 'm/s2')
    call add_var(info, 'IonDrag_east (m/s2)', 'm/s2')
    call add_var(info, 'IonDrag_north (m/s2)', 'm/s2')
    call add_var(info, 'Horiz_Adv_Ve (m/s2)', 'm/s2')
    call add_var(info, 'Horiz_Adv_Vn (m/s2)', 'm/s2')
    call add_var(info, 'PressGrad_east (m/s2)', 'm/s2')
    call add_var(info, 'PressGrad_north (m/s2)', 'm/s2')
    call add_var(info, 'Coriolis_east (m/s2)', 'm/s2')
    call add_var(info, 'Coriolis_north (m/s2)', 'm/s2')
    call add_var(info, 'Centrifugal_north (m/s2)', 'm/s2')
  end subroutine register_3dmoh

  subroutine register_3dmov()
    use ModPlanet, only: nSpecies, cSpecies
    type(OutputTypeInfo), pointer :: info
    integer :: i

    call new_output_type('3DMOV', 3, 2, info)
    call add_coord_vars(info)
    do i = 1, nSpecies
      call add_var(info, 'Vr_'//trim(cSpecies(i))//' (m/s)', 'm/s')
      call add_var(info, 'IonDrag_up_'//trim(cSpecies(i))//' (m/s2)', 'm/s2')
      call add_var(info, 'Horiz_Adv_Vr_'//trim(cSpecies(i))//' (m/s2)', 'm/s2')
    enddo
    call add_var(info, 'Coriolis_up (m/s2)', 'm/s2')
    call add_var(info, 'Centrif_up (m/s2)', 'm/s2')
    call add_var(info, 'EffectiveGravity (m/s2)', 'm/s2')
  end subroutine register_3dmov

  subroutine register_2dgel()
    type(OutputTypeInfo), pointer :: info

    call new_output_type('2DGEL', 2, 0, info)
    call add_coord_vars(info)
    call add_var(info, 'Potential', 'V')
    call add_var(info, 'Pedersen Conductance', 'S')
    call add_var(info, 'Hall Conductance', 'S')
    call add_var(info, 'Electron_Average_Energy_Diffuse', 'eV')
    call add_var(info, 'Electron_Energy_Flux_Diffuse', 'ergs/cm2/s')
    call add_var(info, 'Electron_Average_Energy_Wave', 'eV')
    call add_var(info, 'Electron_Energy_Flux_Wave', 'ergs/cm2/s')
    call add_var(info, 'Electron_Average_Energy_Mono', 'eV')
    call add_var(info, 'Electron_Energy_Flux_Mono', 'ergs/cm2/s')
    call add_var(info, 'Ion_Average_Energy', 'eV')
    call add_var(info, 'Ion_Energy_Flux', 'ergs/cm2/s')
    call add_var(info, 'DivJuAlt', 'A/m2')
    call add_var(info, 'Pedersen FL Conductance', 'S')
    call add_var(info, 'Hall FL Conductance', 'S')
    call add_var(info, 'DivJu FL', 'A/m2')
    call add_var(info, 'FL Length', 'm')
  end subroutine register_2dgel

  subroutine register_2dmel()
    type(OutputTypeInfo), pointer :: info

    call new_output_type('2DMEL', 2, 0, info)
    info%usesMagGrid = .true.
    call add_coord_vars(info)
    call add_var(info, 'MLT', 'hr')
    call add_var(info, 'GeoLat', 'deg')
    call add_var(info, 'GeoLon', 'deg')
    call add_var(info, 'Pedersen Conductance', 'S')
    call add_var(info, 'Hall Conductance', 'S')
    call add_var(info, 'DivJuAlt', 'A/m2')
    call add_var(info, 'Field Line Length', 'm')
    call add_var(info, 'Sigma PP', 'S')
    call add_var(info, 'Sigma LL', 'S')
    call add_var(info, 'Sigma H', 'S')
    call add_var(info, 'Sigma C', 'S')
    call add_var(info, 'Sigma PL', 'S')
    call add_var(info, 'Sigma LP', 'S')
    call add_var(info, 'K^D_{m\phi}', '')
    call add_var(info, 'K^D_{m\lamda}', '')
    call add_var(info, 'Solver A', '')
    call add_var(info, 'Solver B', '')
    call add_var(info, 'Solver C', '')
    call add_var(info, 'Solver D', '')
    call add_var(info, 'Solver E', '')
    call add_var(info, 'Solver S', '')
    call add_var(info, 'DynamoPotential', 'V')
    call add_var(info, 'Ed1new', 'V/m')
    call add_var(info, 'Ed2new', 'V/m')
    call add_var(info, 'Kphi', '')
    call add_var(info, 'Klamda', '')
  end subroutine register_2dmel

  subroutine register_2dusr()
    type(OutputTypeInfo), pointer :: info
    call new_output_type('2DUSR', 2, 0, info)
    info%isUserType = .true.
  end subroutine register_2dusr

  subroutine register_2dtec()
    type(OutputTypeInfo), pointer :: info

    call new_output_type('2DTEC', 2, 0, info)
    info%needsVtec = .true.
    call add_coord_vars(info)
    call add_var(info, 'Solar Zenith Angle', 'rad')
    call add_var(info, 'Vertical TEC', 'TECU')
  end subroutine register_2dtec

  subroutine register_2danc()
    type(OutputTypeInfo), pointer :: info

    call new_output_type('2DANC', 2, 0, info)
    info%needsVtec = .true.
    call add_coord_vars(info)
    call add_var(info, 'Local Time', 'hr')
    call add_var(info, 'Solar Zenith Angle', 'rad')
    call add_var(info, 'Vertical TEC', 'TECU')
    call add_var(info, 'AltIntJouleHeating (W/m2)', 'W/m2')
    call add_var(info, 'AltIntHeatingTransfer (W/m2)', 'W/m2')
    call add_var(info, 'AltIntEuvHeating (W/m2)', 'W/m2')
    call add_var(info, 'AltIntPhotoElectronHeating (W/m2)', 'W/m2')
    call add_var(info, 'AltIntChamicalHeating (W/m2)', 'W/m2')
    call add_var(info, 'AltIntRadCooling (W/m2)', 'W/m2')
    call add_var(info, 'AltIntCO2Cooling (W/m2)', 'W/m2')
    call add_var(info, 'AltIntNOCooling (W/m2)', 'W/m2')
    call add_var(info, 'AltIntOCooling (W/m2)', 'W/m2')
  end subroutine register_2danc

  subroutine register_2dhme()
    type(OutputTypeInfo), pointer :: info

    call new_output_type('2DHME', 2, 0, info)
    info%needsVtec = .true.
    info%isRegional = .true.
    call add_coord_vars(info)
    call add_var(info, 'Local Time', 'hr')
    call add_var(info, 'Vertical TEC', 'TECU')
  end subroutine register_2dhme

  subroutine register_1dall()
    use ModPlanet, only: nSpeciesTotal, nSpecies, nIons, cSpecies, cIons
    type(OutputTypeInfo), pointer :: info
    integer :: i

    call new_output_type('1DALL', 1, 0, info)
    info%needsInterp = .true.
    call add_coord_vars(info)
    call add_var(info, 'Rho', 'kg/m3')
    do i = 1, nSpeciesTotal
      call add_var(info, '['//trim(cSpecies(i))//']', '/m3')
    enddo
    call add_var(info, 'Temperature', 'K')
    call add_var(info, 'V!Dn!N (east)', 'm/s')
    call add_var(info, 'V!Dn!N (north)', 'm/s')
    call add_var(info, 'V!Dn!N (up)', 'm/s')
    do i = 1, nSpecies
      call add_var(info, 'V!Dn!N (up,'//trim(cSpecies(i))//')', 'm/s')
    enddo
    do i = 1, nIons
      call add_var(info, '['//trim(cIons(i))//']', '/m3')
    enddo
    call add_var(info, 'eTemperature', 'K')
    call add_var(info, 'iTemperature', 'K')
    call add_var(info, 'V!Di!N (east)', 'm/s')
    call add_var(info, 'V!Di!N (north)', 'm/s')
    call add_var(info, 'V!Di!N (up)', 'm/s')
  end subroutine register_1dall

  subroutine register_0dall()
    use ModPlanet, only: nSpeciesTotal, nSpecies, nIons, cSpecies, cIons
    type(OutputTypeInfo), pointer :: info
    integer :: i

    call new_output_type('0DALL', 0, 0, info)
    info%needsInterp = .true.
    call add_coord_vars(info)
    call add_var(info, 'Rho', 'kg/m3')
    do i = 1, nSpeciesTotal
      call add_var(info, '['//trim(cSpecies(i))//']', '/m3')
    enddo
    call add_var(info, 'Temperature', 'K')
    call add_var(info, 'V!Dn!N (east)', 'm/s')
    call add_var(info, 'V!Dn!N (north)', 'm/s')
    call add_var(info, 'V!Dn!N (up)', 'm/s')
    do i = 1, nSpecies
      call add_var(info, 'V!Dn!N (up,'//trim(cSpecies(i))//')', 'm/s')
    enddo
    do i = 1, nIons
      call add_var(info, '['//trim(cIons(i))//']', '/m3')
    enddo
    call add_var(info, 'eTemperature', 'K')
    call add_var(info, 'iTemperature', 'K')
    call add_var(info, 'V!Di!N (east)', 'm/s')
    call add_var(info, 'V!Di!N (north)', 'm/s')
    call add_var(info, 'V!Di!N (up)', 'm/s')
    ! 0DALL has extra vars: mixing ratios + heating
    do i = 1, nSpecies
      call add_var(info, trim(cSpecies(i))//' Mixing Ratio', '')
    enddo
    call add_var(info, 'RadCooling', 'K/s')
    call add_var(info, 'EuvHeating', 'K/s')
    call add_var(info, 'Conduction', 'K/s')
    call add_var(info, 'Heat Balance Total', 'K/s')
    call add_var(info, 'Heating Efficiency', '')
  end subroutine register_0dall

  subroutine register_1dglo()
    type(OutputTypeInfo), pointer :: info

    call new_output_type('1DGLO', 1, 0, info)
    call add_coord_vars(info)
    call add_var(info, '6300 A Emission', '')
    call add_var(info, 'PhotoElectronUp', '')
    call add_var(info, 'PhotoElectronDown', '')
  end subroutine register_1dglo

  subroutine register_1dthm()
    use ModPlanet, only: nSpeciesTotal, cSpecies
    type(OutputTypeInfo), pointer :: info
    integer :: i

    call new_output_type('1DTHM', 1, 0, info)
    call add_coord_vars(info)
    call add_var(info, 'EUV Heating (K/s)', 'K/s')
    call add_var(info, 'Conduction (K/s)', 'K/s')
    call add_var(info, 'Molecular Conduction (K/s)', 'K/s')
    call add_var(info, 'Eddy Conduction (K/s)', 'K/s')
    call add_var(info, 'Eddy Adiabatic Conduction (K/s)', 'K/s')
    call add_var(info, 'Chemical Heating (K/s)', 'K/s')
    call add_var(info, 'Joule Heating (K/s)', 'K/s')
    call add_var(info, 'NO Cooling (K/s)', 'K/s')
    call add_var(info, 'O Cooling (K/s)', 'K/s')
    call add_var(info, 'Total Abs EUV', 'W/m2')
    do i = 1, nSpeciesTotal
      call add_var(info, 'Production Rate '//trim(cSpecies(i)), '/m3/s')
    enddo
    do i = 1, nSpeciesTotal
      call add_var(info, 'Loss Rate '//trim(cSpecies(i)), '/m3/s')
    enddo
  end subroutine register_1dthm

  subroutine register_1dnew()
    use ModPlanet, only: nSpeciesTotal, nSpecies, nIons, cSpecies, cIons
    type(OutputTypeInfo), pointer :: info
    integer :: i

    call new_output_type('1DNEW', 1, 0, info)
    info%needsInterp = .true.
    ! 1DNEW has different coord layout: Lon, LocalTime, Lat, SZA, Alt
    call add_var(info, 'Longitude', 'rad')
    call add_var(info, 'Local Time', 'hr')
    call add_var(info, 'Latitude', 'rad')
    call add_var(info, 'Solar Zenith Angle', 'rad')
    call add_var(info, 'Altitude', 'm')
    call add_var(info, 'Rho', 'kg/m3')
    do i = 1, nSpeciesTotal
      call add_var(info, '['//trim(cSpecies(i))//']', '/m3')
    enddo
    call add_var(info, 'Temperature', 'K')
    call add_var(info, 'V!Dn!N (east)', 'm/s')
    call add_var(info, 'V!Dn!N (north)', 'm/s')
    call add_var(info, 'V!Dn!N (up)', 'm/s')
    do i = 1, nSpecies
      call add_var(info, 'V!Dn!N (up,'//trim(cSpecies(i))//')', 'm/s')
    enddo
    do i = 1, nIons
      call add_var(info, '['//trim(cIons(i))//']', '/m3')
    enddo
    call add_var(info, 'eTemperature', 'K')
    call add_var(info, 'iTemperature', 'K')
    call add_var(info, 'V!Di!N (east)', 'm/s')
    call add_var(info, 'V!Di!N (north)', 'm/s')
    call add_var(info, 'V!Di!N (up)', 'm/s')
    do i = 1, nSpecies
      call add_var(info, trim(cSpecies(i))//' Mixing Ratio', '')
    enddo
  end subroutine register_1dnew

  subroutine register_1dchm()
    type(OutputTypeInfo), pointer :: info

    call new_output_type('1DCHM', 1, 0, info)
    call add_coord_vars(info)
    call add_var(info, 'N!D2!U+!N + e', '/m3/s')
    call add_var(info, 'O!D2!U+!N + e', '/m3/s')
    call add_var(info, 'N!D2!U+!N + O', '/m3/s')
    call add_var(info, 'NO!U+!N + e', '/m3/s')
    call add_var(info, 'N!U+!N + O!D2!N', '/m3/s')
    call add_var(info, 'NO + N', '/m3/s')
    call add_var(info, 'O!U+!N + O!D2!N', '/m3/s')
    call add_var(info, 'N + O!D2!N', '/m3/s')
    call add_var(info, 'O!D2!U+!N + N', '/m3/s')
    call add_var(info, 'O!D2!U+!N + NO', '/m3/s')
    call add_var(info, 'O!D2!U+!N + N2', '/m3/s')
    call add_var(info, 'N!D2!U+!N + O!D2!N', '/m3/s')
    call add_var(info, 'N!U+!N + O', '/m3/s')
    call add_var(info, 'O!+!N + N!D2!N', '/m3/s')
    call add_var(info, 'O(1D) + N!D2!N', '/m3/s')
    call add_var(info, 'O(1D) + O!D2!N', '/m3/s')
    call add_var(info, 'O(1D) + O', '/m3/s')
    call add_var(info, 'O(1D) + e', '/m3/s')
    call add_var(info, 'N(2D) + O!D2!N', '/m3/s')
    call add_var(info, 'O!U+!N(2D)+e', '/m3/s')
    call add_var(info, 'N(2D) + O', '/m3/s')
    call add_var(info, 'N(2D) + e', '/m3/s')
    call add_var(info, 'O!U+!N(2D + N!D2!N', '/m3/s')
    call add_var(info, 'O!U+!N(2P) + e', '/m3/s')
    call add_var(info, 'O!U+!N(2P) + O', '/m3/s')
    call add_var(info, 'O!U+!N(2P) + N!D2!N', '/m3/s')
    call add_var(info, 'Chemical Heating Rate', 'K/s')
  end subroutine register_1dchm

  subroutine register_1dcms()
    ! 1DCMS is listed as valid but has no implementation in the dispatcher.
    ! Register as a placeholder.
    type(OutputTypeInfo), pointer :: info
    call new_output_type('1DCMS', 1, 0, info)
  end subroutine register_1dcms

  subroutine register_1dusr()
    type(OutputTypeInfo), pointer :: info
    call new_output_type('1DUSR', 1, 0, info)
    info%isUserType = .true.
    info%needsInterp = .true.
  end subroutine register_1dusr

  subroutine register_0dusr()
    type(OutputTypeInfo), pointer :: info
    call new_output_type('0DUSR', 0, 0, info)
    info%isUserType = .true.
    info%needsInterp = .true.
  end subroutine register_0dusr

end module ModOutputRegistry
