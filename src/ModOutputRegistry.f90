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
    call add_var(info, 'Rho', 'kg/m3', 'rho', 'Total neutral mass density')
    do i = 1, nSpeciesTotal
      call add_var(info, '['//trim(cSpecies(i))//']', '/m3', &
                   trim(cSpecies(i)), 'Number density of '//trim(cSpecies(i)))
    enddo
    call add_var(info, 'Temperature', 'K', 'Tn', 'Neutral temperature')
    call add_var(info, 'Vn (east)', 'm/s', 'Ve', 'Eastward neutral wind')
    call add_var(info, 'Vn (north)', 'm/s', 'Vn', 'Northward neutral wind')
    call add_var(info, 'Vn (up)', 'm/s', 'Vv', 'Vertical neutral wind')
    do i = 1, nSpecies
      call add_var(info, 'Vn (up,'//trim(cSpecies(i))//')', 'm/s', &
                   'Vv_'//trim(cSpecies(i)), 'Vertical neutral wind for '//trim(cSpecies(i)))
    enddo
    do i = 1, nIons
      call add_var(info, '['//trim(cIons(i))//']', '/m3', &
                   trim(cIons(i)), 'Number density of '//trim(cIons(i)))
    enddo
    call add_var(info, 'eTemperature', 'K', 'Te', 'Electron temperature')
    call add_var(info, 'iTemperature', 'K', 'Ti', 'Ion temperature')
    call add_var(info, 'Vi (east)', 'm/s', 'Vie', 'Eastward ion drift')
    call add_var(info, 'Vi (north)', 'm/s', 'Vin', 'Northward ion drift')
    call add_var(info, 'Vi (up)', 'm/s', 'Viv', 'Vertical ion drift')
  end subroutine register_3dall

  subroutine register_3dlst()
    use ModPlanet, only: nSpeciesTotal, nIons, cSpecies, cIons
    use ModInputs, only: iRhoOutputList, iNeutralDensityOutputList, &
                         iNeutralWindOutputList, iIonDensityOutputList, &
                         iIonWindOutputList, iTemperatureOutputList
    type(OutputTypeInfo), pointer :: info
    integer :: i

    call new_output_type('3DLST', 3, 2, info)
    call add_coord_vars(info)

    if (iRhoOutputList) &
      call add_var(info, 'Rho (kg/m3)', 'kg/m3', 'rho', 'Total neutral mass density')
    do i = 1, nSpeciesTotal
      if (iNeutralDensityOutputList(i)) &
        call add_var(info, '['//trim(cSpecies(i))//'] (/m3)', '/m3', &
                     trim(cSpecies(i)), 'Number density of '//trim(cSpecies(i)))
    enddo

    if (iNeutralWindOutputList(1)) &
      call add_var(info, 'Vn (east)', 'm/s', 'Ve', 'Eastward neutral wind')
    if (iNeutralWindOutputList(2)) &
      call add_var(info, 'Vn (north) (m/s)', 'm/s', 'Vn', 'Northward neutral wind')
    if (iNeutralWindOutputList(3)) &
      call add_var(info, 'Vn (up) (m/s)', 'm/s', 'Vv', 'Vertical neutral wind')

    do i = 1, nIons
      if (iIonDensityOutputList(i)) &
        call add_var(info, '['//trim(cIons(i))//'] (/m3)', '/m3', &
                     trim(cIons(i)), 'Number density of '//trim(cIons(i)))
    enddo

    if (iIonWindOutputList(1)) &
      call add_var(info, 'Vi (east) (m/s)', 'm/s', 'Vie', 'Eastward ion drift')
    if (iIonWindOutputList(2)) &
      call add_var(info, 'Vi (north) (m/s)', 'm/s', 'Vin', 'Northward ion drift')
    if (iIonWindOutputList(3)) &
      call add_var(info, 'Vi (up) (m/s)', 'm/s', 'Viv', 'Vertical ion drift')

    if (iTemperatureOutputList(1)) &
      call add_var(info, 'Neutral Temperature (K)', 'K', 'Tn', 'Neutral temperature')
    if (iTemperatureOutputList(2)) &
      call add_var(info, 'Ion Temperature (K)', 'K', 'Ti', 'Ion temperature')
    if (iTemperatureOutputList(3)) &
      call add_var(info, 'Electron Temperature (K)', 'K', 'Te', 'Electron temperature')
  end subroutine register_3dlst

  subroutine register_3dneu()
    use ModPlanet, only: nSpeciesTotal, nSpecies, cSpecies
    type(OutputTypeInfo), pointer :: info
    integer :: i

    call new_output_type('3DNEU', 3, 2, info)
    call add_coord_vars(info)
    call add_var(info, 'Rho', 'kg/m3', 'rho', 'Total neutral mass density')
    do i = 1, nSpeciesTotal
      call add_var(info, '['//trim(cSpecies(i))//']', '/m3', &
                   trim(cSpecies(i)), 'Number density of '//trim(cSpecies(i)))
    enddo
    call add_var(info, 'Temperature', 'K', 'Tn', 'Neutral temperature')
    call add_var(info, 'Vn (east)', 'm/s', 'Ve', 'Eastward neutral wind')
    call add_var(info, 'Vn (north)', 'm/s', 'Vn', 'Northward neutral wind')
    call add_var(info, 'Vn (up)', 'm/s', 'Vv', 'Vertical neutral wind')
    do i = 1, nSpecies
      call add_var(info, 'Vn (up,'//trim(cSpecies(i))//')', 'm/s', &
                   'Vv_'//trim(cSpecies(i)), 'Vertical neutral wind for '//trim(cSpecies(i)))
    enddo
  end subroutine register_3dneu

  subroutine register_3dion()
    use ModPlanet, only: nIons, cIons
    type(OutputTypeInfo), pointer :: info
    integer :: i

    call new_output_type('3DION', 3, 2, info)
    call add_coord_vars(info)
    do i = 1, nIons
      call add_var(info, '['//trim(cIons(i))//']', '/m3', &
                   trim(cIons(i)), 'Number density of '//trim(cIons(i)))
    enddo
    call add_var(info, 'eTemperature', 'K', 'Te', 'Electron temperature')
    call add_var(info, 'iTemperature', 'K', 'Ti', 'Ion temperature')
    call add_var(info, 'Vi (east)', 'm/s', 'Vie', 'Eastward ion drift')
    call add_var(info, 'Vi (north)', 'm/s', 'Vin', 'Northward ion drift')
    call add_var(info, 'Vi (up)', 'm/s', 'Viv', 'Vertical ion drift')
    call add_var(info, 'Ed1', 'V/m', 'Ed1', 'Electric field component 1 (magnetic east)')
    call add_var(info, 'Ed2', 'V/m', 'Ed2', 'Electric field component 2 (magnetic north)')
    call add_var(info, 'Je1', 'A/m2', 'Je1', 'Current density component 1 (magnetic east)')
    call add_var(info, 'Je2', 'A/m2', 'Je2', 'Current density component 2 (magnetic north)')
    call add_var(info, 'Magnetic Latitude', 'deg', 'MLat', 'Magnetic latitude')
    call add_var(info, 'Magnetic Longitude', 'deg', 'MLon', 'Magnetic longitude')
    call add_var(info, 'B.F. East', 'T', 'BF_e', 'Eastward magnetic field component')
    call add_var(info, 'B.F. North', 'T', 'BF_n', 'Northward magnetic field component')
    call add_var(info, 'B.F. Vertical', 'T', 'BF_v', 'Vertical magnetic field component')
    call add_var(info, 'B.F. Magnitude', 'T', 'BF_mag', 'Magnetic field magnitude')
    call add_var(info, 'Potential', 'V', 'pot', 'Electric potential')
    call add_var(info, 'E.F. East', 'V/m', 'EF_e', 'Eastward electric field component')
    call add_var(info, 'E.F. North', 'V/m', 'EF_n', 'Northward electric field component')
    call add_var(info, 'E.F. Vertical', 'V/m', 'EF_v', 'Vertical electric field component')
    call add_var(info, 'E.F. Magnitude', 'V/m', 'EF_mag', 'Electric field magnitude')
    call add_var(info, 'IN Collision Freq', 'Hz', 'nu_in', 'Ion-neutral collision frequency')
    call add_var(info, 'PressGrad (east)', 'm/s2', 'dP_east', 'Pressure gradient acceleration eastward')
    call add_var(info, 'PressGrad (north)', 'm/s2', 'dP_north', 'Pressure gradient acceleration northward')
    call add_var(info, 'PressGrad (up)', 'm/s2', 'dP_up', 'Pressure gradient acceleration upward')
  end subroutine register_3dion

  subroutine register_3dthm()
    type(OutputTypeInfo), pointer :: info

    call new_output_type('3DTHM', 3, 0, info)
    call add_coord_vars(info)
    call add_var(info, 'EUV Heating (K/s)', 'K/s', 'Q_euv', 'EUV heating rate')
    call add_var(info, 'Conduction (K/s)', 'K/s', 'Q_cond', 'Thermal conduction heating rate')
    call add_var(info, 'Molecular Conduction (K/s)', 'K/s', 'Q_mol_cond', 'Molecular conduction heating rate')
    call add_var(info, 'Eddy Conduction (K/s)', 'K/s', 'Q_eddy_cond', 'Eddy conduction heating rate')
    call add_var(info, 'Eddy Adiabatic Conduction (K/s)', 'K/s', 'Q_eddy_adiab', &
                 'Eddy adiabatic conduction heating rate')
    call add_var(info, 'Chemical Heating (K/s)', 'K/s', 'Q_chem', 'Chemical heating rate')
    call add_var(info, 'Joule Heating (K/s)', 'K/s', 'Q_joule', 'Joule heating rate')
    call add_var(info, 'NO Cooling (K/s)', 'K/s', 'Q_no_cool', 'NO infrared cooling rate')
    call add_var(info, 'O Cooling (K/s)', 'K/s', 'Q_o_cool', 'O infrared cooling rate')
    call add_var(info, 'Total Abs EUV', 'W/m2', 'EUV_abs', 'Total absorbed EUV flux')
    call add_var(info, 'Cp', 'J/kg/K', 'Cp', 'Specific heat at constant pressure')
    call add_var(info, 'Rho', 'kg/m3', 'rho', 'Total neutral mass density')
    call add_var(info, 'E-Field Mag', 'V/m', 'EF_mag', 'Electric field magnitude')
    call add_var(info, 'Sigma Ped', 'S/m', 'Sigma_P', 'Pedersen conductivity')
    call add_var(info, 'Ionization Rate O_3P', '/s', 'ion_rate_O3P', 'Ionization rate for O(3P)')
    call add_var(info, 'Ionization Rate O2', '/s', 'ion_rate_O2', 'Ionization rate for O2')
    call add_var(info, 'Ionization Rate N2', '/s', 'ion_rate_N2', 'Ionization rate for N2')
  end subroutine register_3dthm

  subroutine register_3dchm()
    type(OutputTypeInfo), pointer :: info

    call new_output_type('3DCHM', 3, 0, info)
    call add_coord_vars(info)
    call add_var(info, 'N2+ + e', '/m3/s', 'R_N2p_e', 'N2+ + e recombination rate')
    call add_var(info, 'O2+ + e', '/m3/s', 'R_O2p_e', 'O2+ + e recombination rate')
    call add_var(info, 'N2+ + O', '/m3/s', 'R_N2p_O', 'N2+ + O reaction rate')
    call add_var(info, 'NO+ + e', '/m3/s', 'R_NOp_e', 'NO+ + e recombination rate')
    call add_var(info, 'N+ + O2', '/m3/s', 'R_Np_O2', 'N+ + O2 reaction rate')
    call add_var(info, 'NO + N', '/m3/s', 'R_NO_N', 'NO + N reaction rate')
    call add_var(info, 'O+ + O2', '/m3/s', 'R_Op_O2', 'O+ + O2 reaction rate')
    call add_var(info, 'N + O2', '/m3/s', 'R_N_O2', 'N + O2 reaction rate')
    call add_var(info, 'O2+ + N', '/m3/s', 'R_O2p_N', 'O2+ + N reaction rate')
    call add_var(info, 'O2+ + NO', '/m3/s', 'R_O2p_NO', 'O2+ + NO reaction rate')
    call add_var(info, 'O2+ + N2', '/m3/s', 'R_O2p_N2', 'O2+ + N2 reaction rate')
    call add_var(info, 'N2+ + O2', '/m3/s', 'R_N2p_O2', 'N2+ + O2 reaction rate')
    call add_var(info, 'N+ + O', '/m3/s', 'R_Np_O', 'N+ + O reaction rate')
    call add_var(info, 'O+ + N2', '/m3/s', 'R_Op_N2', 'O+ + N2 reaction rate')
    call add_var(info, 'O(1D) + N2', '/m3/s', 'R_O1D_N2', 'O(1D) + N2 quench rate')
    call add_var(info, 'O(1D) + O2', '/m3/s', 'R_O1D_O2', 'O(1D) + O2 quench rate')
    call add_var(info, 'O(1D) + O', '/m3/s', 'R_O1D_O', 'O(1D) + O quench rate')
    call add_var(info, 'O(1D) + e', '/m3/s', 'R_O1D_e', 'O(1D) + e quench rate')
    call add_var(info, 'N(2D) + O2', '/m3/s', 'R_N2D_O2', 'N(2D) + O2 reaction rate')
    call add_var(info, 'O+(2D)+e', '/m3/s', 'R_O2Dp_e', 'O+(2D) + e recombination rate')
    call add_var(info, 'N(2D) + O', '/m3/s', 'R_N2D_O', 'N(2D) + O quench rate')
    call add_var(info, 'N(2D) + e', '/m3/s', 'R_N2D_e', 'N(2D) + e quench rate')
    call add_var(info, 'O+(2D + N2', '/m3/s', 'R_O2Dp_N2', 'O+(2D) + N2 reaction rate')
    call add_var(info, 'O+(2P) + e', '/m3/s', 'R_O2Pp_e', 'O+(2P) + e recombination rate')
    call add_var(info, 'O+(2P) + O', '/m3/s', 'R_O2Pp_O', 'O+(2P) + O reaction rate')
    call add_var(info, 'O+(2P) + N2', '/m3/s', 'R_O2Pp_N2', 'O+(2P) + N2 reaction rate')
    call add_var(info, 'Chemical Heating Rate', 'K/s', 'Q_chem_rate', 'Chemical heating rate')
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
    call add_var(info, '6300 A Emission', '', 'I_6300', '6300 Angstrom emission intensity')
    call add_var(info, 'PhotoElectronUp', '', 'PEflux_up', 'Upward photoelectron flux')
    call add_var(info, 'PhotoElectronDown', '', 'PEflux_dn', 'Downward photoelectron flux')
  end subroutine register_3dglo

  subroutine register_3dmag()
    type(OutputTypeInfo), pointer :: info

    call new_output_type('3DMAG', 3, 2, info)
    call add_coord_vars(info)
    call add_var(info, 'Magnetic Latitude', 'deg', 'MLat', 'Magnetic latitude')
    call add_var(info, 'Magnetic Longitude', 'deg', 'MLon', 'Magnetic longitude')
    call add_var(info, 'B.F. East', 'T', 'BF_e', 'Eastward magnetic field component')
    call add_var(info, 'B.F. North', 'T', 'BF_n', 'Northward magnetic field component')
    call add_var(info, 'B.F. Vertical', 'T', 'BF_v', 'Vertical magnetic field component')
    call add_var(info, 'B.F. Magnitude', 'T', 'BF_mag', 'Magnetic field magnitude')
  end subroutine register_3dmag

  subroutine register_3dhme()
    use ModPlanet, only: nSpeciesTotal, nSpecies, nIons, cSpecies, cIons
    type(OutputTypeInfo), pointer :: info
    integer :: i

    call new_output_type('3DHME', 3, 0, info)
    call add_coord_vars(info)
    call add_var(info, 'Rho', 'kg/m3', 'rho', 'Total neutral mass density')
    do i = 1, nSpeciesTotal
      call add_var(info, '['//trim(cSpecies(i))//']', '/m3', &
                   trim(cSpecies(i)), 'Number density of '//trim(cSpecies(i)))
    enddo
    call add_var(info, 'Temperature', 'K', 'Tn', 'Neutral temperature')
    call add_var(info, 'Vn (east)', 'm/s', 'Ve', 'Eastward neutral wind')
    call add_var(info, 'Vn (north)', 'm/s', 'Vn', 'Northward neutral wind')
    call add_var(info, 'Vn (up)', 'm/s', 'Vv', 'Vertical neutral wind')
    do i = 1, nSpecies
      call add_var(info, 'Vn (up,'//trim(cSpecies(i))//')', 'm/s', &
                   'Vv_'//trim(cSpecies(i)), 'Vertical neutral wind for '//trim(cSpecies(i)))
    enddo
    do i = 1, nIons
      call add_var(info, '['//trim(cIons(i))//']', '/m3', &
                   trim(cIons(i)), 'Number density of '//trim(cIons(i)))
    enddo
    call add_var(info, 'eTemperature', 'K', 'Te', 'Electron temperature')
    call add_var(info, 'iTemperature', 'K', 'Ti', 'Ion temperature')
    call add_var(info, 'Vi (east)', 'm/s', 'Vie', 'Eastward ion drift')
    call add_var(info, 'Vi (north)', 'm/s', 'Vin', 'Northward ion drift')
    call add_var(info, 'Vi (up)', 'm/s', 'Viv', 'Vertical ion drift')
    call add_var(info, 'PhotoElectron Heating', 'K', 'Q_pe', 'Photoelectron heating rate')
    call add_var(info, 'Joule Heating', 'K', 'Q_joule', 'Joule heating rate')
    call add_var(info, 'Specific Heat', 'J/kg/K', 'Cp', 'Specific heat at constant pressure')
    call add_var(info, 'Magnetic Latitude', 'deg', 'MLat', 'Magnetic latitude')
    call add_var(info, 'Magnetic Longitude', 'deg', 'MLon', 'Magnetic longitude')
    call add_var(info, 'B.F. East', 'T', 'BF_e', 'Eastward magnetic field component')
    call add_var(info, 'B.F. North', 'T', 'BF_n', 'Northward magnetic field component')
    call add_var(info, 'B.F. Vertical', 'T', 'BF_v', 'Vertical magnetic field component')
    call add_var(info, 'B.F. Magnitude', 'T', 'BF_mag', 'Magnetic field magnitude')
    call add_var(info, 'Potential', 'V', 'pot', 'Electric potential')
    call add_var(info, 'PotentialY', 'V', 'Phi_Y', 'Electric potential Y component')
    call add_var(info, 'E.F. East', 'V/m', 'EF_e', 'Eastward electric field component')
    call add_var(info, 'E.F. North', 'V/m', 'EF_n', 'Northward electric field component')
    call add_var(info, 'E.F. Vertical', 'V/m', 'EF_v', 'Vertical electric field component')
  end subroutine register_3dhme

  subroutine register_3dmoh()
    type(OutputTypeInfo), pointer :: info

    call new_output_type('3DMOH', 3, 0, info)
    call add_coord_vars(info)
    call add_var(info, 'Rho (kg/m3)', 'kg/m3', 'rho', 'Total neutral mass density')
    call add_var(info, 'Vn_east (m/s)', 'm/s', 'Ve', 'Eastward neutral wind')
    call add_var(info, 'Vn_north (m/s)', 'm/s', 'Vn', 'Northward neutral wind')
    call add_var(info, 'Visc_Ve_rshear (m/s2)', 'm/s2', 'Visc_east', &
                 'Viscous shear acceleration eastward')
    call add_var(info, 'Visc_Vn_rshear (m/s2)', 'm/s2', 'Visc_north', &
                 'Viscous shear acceleration northward')
    call add_var(info, 'IonDrag_east (m/s2)', 'm/s2', 'IonDrag_east', &
                 'Ion drag acceleration eastward')
    call add_var(info, 'IonDrag_north (m/s2)', 'm/s2', 'IonDrag_north', &
                 'Ion drag acceleration northward')
    call add_var(info, 'Horiz_Adv_Ve (m/s2)', 'm/s2', 'HorizAdv_east', &
                 'Horizontal advection acceleration eastward')
    call add_var(info, 'Horiz_Adv_Vn (m/s2)', 'm/s2', 'HorizAdv_north', &
                 'Horizontal advection acceleration northward')
    call add_var(info, 'PressGrad_east (m/s2)', 'm/s2', 'PressGrad_east', &
                 'Pressure gradient acceleration eastward')
    call add_var(info, 'PressGrad_north (m/s2)', 'm/s2', 'PressGrad_north', &
                 'Pressure gradient acceleration northward')
    call add_var(info, 'Coriolis_east (m/s2)', 'm/s2', 'Coriolis_east', &
                 'Coriolis acceleration eastward')
    call add_var(info, 'Coriolis_north (m/s2)', 'm/s2', 'Coriolis_north', &
                 'Coriolis acceleration northward')
    call add_var(info, 'Centrifugal_north (m/s2)', 'm/s2', 'Centrif_north', &
                 'Centrifugal acceleration northward')
  end subroutine register_3dmoh

  subroutine register_3dmov()
    use ModPlanet, only: nSpecies, cSpecies
    type(OutputTypeInfo), pointer :: info
    integer :: i

    call new_output_type('3DMOV', 3, 0, info)
    call add_coord_vars(info)
    do i = 1, nSpecies
      call add_var(info, 'Vr_'//trim(cSpecies(i))//' (m/s)', 'm/s', &
                   'Vv_'//trim(cSpecies(i)), 'Vertical neutral wind for '//trim(cSpecies(i)))
      call add_var(info, 'IonDrag_up_'//trim(cSpecies(i))//' (m/s2)', 'm/s2', &
                   'IonDrag_up_'//trim(cSpecies(i)), &
                   'Ion drag acceleration upward for '//trim(cSpecies(i)))
      call add_var(info, 'Horiz_Adv_Vr_'//trim(cSpecies(i))//' (m/s2)', 'm/s2', &
                   'HorizAdv_up_'//trim(cSpecies(i)), &
                   'Horizontal advection upward for '//trim(cSpecies(i)))
    enddo
    call add_var(info, 'Coriolis_up (m/s2)', 'm/s2', 'Coriolis_up', &
                 'Coriolis acceleration upward')
    call add_var(info, 'Centrif_up (m/s2)', 'm/s2', 'Centrif_up', &
                 'Centrifugal acceleration upward')
    call add_var(info, 'EffectiveGravity (m/s2)', 'm/s2', 'g_eff', &
                 'Effective gravitational acceleration')
  end subroutine register_3dmov

  subroutine register_2dgel()
    type(OutputTypeInfo), pointer :: info

    call new_output_type('2DGEL', 2, 0, info)
    call add_coord_vars(info)
    call add_var(info, 'Potential', 'V', 'pot', 'Electric potential')
    call add_var(info, 'Pedersen Conductance', 'S', 'PedCond', &
                 'Height-integrated Pedersen conductance')
    call add_var(info, 'Hall Conductance', 'S', 'HalCond', &
                 'Height-integrated Hall conductance')
    call add_var(info, 'Electron_Average_Energy_Diffuse', 'eV', 'AveE', &
                 'Diffuse electron precipitation average energy')
    call add_var(info, 'Electron_Energy_Flux_Diffuse', 'ergs/cm2/s', 'eFlux', &
                 'Diffuse electron precipitation energy flux')
    call add_var(info, 'Electron_Average_Energy_Wave', 'eV', 'AveE_W', &
                 'Wave-driven electron precipitation average energy')
    call add_var(info, 'Electron_Energy_Flux_Wave', 'ergs/cm2/s', 'eFlux_W', &
                 'Wave-driven electron precipitation energy flux')
    call add_var(info, 'Electron_Average_Energy_Mono', 'eV', 'AveE_M', &
                 'Monoenergetic electron precipitation average energy')
    call add_var(info, 'Electron_Energy_Flux_Mono', 'ergs/cm2/s', 'eFlux_M', &
                 'Monoenergetic electron precipitation energy flux')
    call add_var(info, 'Ion_Average_Energy', 'eV', 'AveE_I', 'Ion precipitation average energy')
    call add_var(info, 'Ion_Energy_Flux', 'ergs/cm2/s', 'eFlux_I', 'Ion precipitation energy flux')
    call add_var(info, 'DivJuAlt', 'A/m2', 'DivJuAlt', &
                 'Height-integrated divergence of horizontal current')
    call add_var(info, 'Pedersen FL Conductance', 'S', 'PedFLCond', &
                 'Field-line integrated Pedersen conductance')
    call add_var(info, 'Hall FL Conductance', 'S', 'HalFLCond', &
                 'Field-line integrated Hall conductance')
    call add_var(info, 'DivJu FL', 'A/m2', 'DivJuFL', 'Field-line divergence of current')
    call add_var(info, 'FL Length', 'm', 'FL_length', 'Field line length')
  end subroutine register_2dgel

  subroutine register_2dmel()
    type(OutputTypeInfo), pointer :: info

    call new_output_type('2DMEL', 2, 0, info)
    info%usesMagGrid = .true.
    ! We need mag coords, not geographic!
    call add_var(info, 'mlon', 'deg', 'mlon', 'Magnetic longitude')
    call add_var(info, 'mlat', 'deg', 'mlat', 'Magnetic latitude')

    call add_var(info, 'MLT', 'hr', 'MLT', 'Magnetic local time')
    call add_var(info, 'GeoLat', 'deg', 'GeoLat', 'Geographic latitude')
    call add_var(info, 'GeoLon', 'deg', 'GeoLon', 'Geographic longitude')
    call add_var(info, 'Pedersen Conductance', 'S', 'PedCond', &
                 'Height-integrated Pedersen conductance')
    call add_var(info, 'Hall Conductance', 'S', 'HalCond', &
                 'Height-integrated Hall conductance')
    call add_var(info, 'DivJuAlt', 'A/m2', 'DivJuAlt', &
                 'Height-integrated divergence of horizontal current')
    call add_var(info, 'Field Line Length', 'm', 'FL_length', 'Field line length')
    call add_var(info, 'Sigma PP', 'S', 'Sigma_PP', 'Sigma PP conductance component')
    call add_var(info, 'Sigma LL', 'S', 'Sigma_LL', 'Sigma LL conductance component')
    call add_var(info, 'Sigma H', 'S', 'Sigma_H', 'Hall conductance component')
    call add_var(info, 'Sigma C', 'S', 'Sigma_C', 'Cowling conductance component')
    call add_var(info, 'Sigma PL', 'S', 'Sigma_PL', 'Sigma PL conductance component')
    call add_var(info, 'Sigma LP', 'S', 'Sigma_LP', 'Sigma LP conductance component')
    call add_var(info, 'K^D_{m\phi}', '', 'K_mphi', 'Dynamo K_mphi coefficient')
    call add_var(info, 'K^D_{m\lamda}', '', 'K_mlambda', 'Dynamo K_mlambda coefficient')
    call add_var(info, 'Solver A', '', 'Solver_A', 'Dynamo solver coefficient A')
    call add_var(info, 'Solver B', '', 'Solver_B', 'Dynamo solver coefficient B')
    call add_var(info, 'Solver C', '', 'Solver_C', 'Dynamo solver coefficient C')
    call add_var(info, 'Solver D', '', 'Solver_D', 'Dynamo solver coefficient D')
    call add_var(info, 'Solver E', '', 'Solver_E', 'Dynamo solver coefficient E')
    call add_var(info, 'Solver S', '', 'Solver_S', 'Dynamo solver coefficient S')
    call add_var(info, 'DynamoPotential', 'V', 'Phi_dynamo', 'Dynamo electric potential')
    call add_var(info, 'Ed1new', 'V/m', 'Ed1new', 'Dynamo electric field component 1')
    call add_var(info, 'Ed2new', 'V/m', 'Ed2new', 'Dynamo electric field component 2')
    call add_var(info, 'Kphi', '', 'Kphi', 'Dynamo Kphi coefficient')
    call add_var(info, 'Klamda', '', 'Klamda', 'Dynamo Klambda coefficient')
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
    call add_var(info, 'Solar Zenith Angle', 'rad', 'SZA', 'Solar zenith angle')
    call add_var(info, 'Vertical TEC', 'TECU', 'TEC', 'Vertical total electron content')
  end subroutine register_2dtec

  subroutine register_2danc()
    type(OutputTypeInfo), pointer :: info

    call new_output_type('2DANC', 2, 0, info)
    info%needsVtec = .true.
    call add_coord_vars(info)
    call add_var(info, 'Local Time', 'hr', 'LT', 'Local solar time')
    call add_var(info, 'Solar Zenith Angle', 'rad', 'SZA', 'Solar zenith angle')
    call add_var(info, 'Vertical TEC', 'TECU', 'TEC', 'Vertical total electron content')
    call add_var(info, 'AltIntJouleHeating (W/m2)', 'W/m2', 'JH_int', &
                 'Altitude-integrated Joule heating')
    call add_var(info, 'AltIntHeatingTransfer (W/m2)', 'W/m2', 'HT_int', &
                 'Altitude-integrated heating transfer')
    call add_var(info, 'AltIntEuvHeating (W/m2)', 'W/m2', 'EUV_int', &
                 'Altitude-integrated EUV heating')
    call add_var(info, 'AltIntPhotoElectronHeating (W/m2)', 'W/m2', 'PE_int', &
                 'Altitude-integrated photoelectron heating')
    call add_var(info, 'AltIntChemicalHeating (W/m2)', 'W/m2', 'Chem_int', &
                 'Altitude-integrated chemical heating')
    call add_var(info, 'AltIntRadCooling (W/m2)', 'W/m2', 'RadCool_int', &
                 'Altitude-integrated radiative cooling')
    call add_var(info, 'AltIntCO2Cooling (W/m2)', 'W/m2', 'CO2Cool_int', &
                 'Altitude-integrated CO2 cooling')
    call add_var(info, 'AltIntNOCooling (W/m2)', 'W/m2', 'NOCool_int', &
                 'Altitude-integrated NO cooling')
    call add_var(info, 'AltIntOCooling (W/m2)', 'W/m2', 'OCool_int', &
                 'Altitude-integrated O cooling')
  end subroutine register_2danc

  subroutine register_2dhme()
    type(OutputTypeInfo), pointer :: info

    call new_output_type('2DHME', 2, 0, info)
    info%needsVtec = .true.
    info%isRegional = .true.
    call add_coord_vars(info)
    call add_var(info, 'Local Time', 'hr', 'LT', 'Local solar time')
    call add_var(info, 'Vertical TEC', 'TECU', 'TEC', 'Vertical total electron content')
  end subroutine register_2dhme

  subroutine register_1dall()
    use ModPlanet, only: nSpeciesTotal, nSpecies, nIons, cSpecies, cIons
    type(OutputTypeInfo), pointer :: info
    integer :: i

    call new_output_type('1DALL', 1, 0, info)
    info%needsInterp = .true.
    call add_coord_vars(info)
    call add_var(info, 'Rho', 'kg/m3', 'rho', 'Total neutral mass density')
    do i = 1, nSpeciesTotal
      call add_var(info, '['//trim(cSpecies(i))//']', '/m3', &
                   trim(cSpecies(i)), 'Number density of '//trim(cSpecies(i)))
    enddo
    call add_var(info, 'Temperature', 'K', 'Tn', 'Neutral temperature')
    call add_var(info, 'Vn (east)', 'm/s', 'Ve', 'Eastward neutral wind')
    call add_var(info, 'Vn (north)', 'm/s', 'Vn', 'Northward neutral wind')
    call add_var(info, 'Vn (up)', 'm/s', 'Vv', 'Vertical neutral wind')
    do i = 1, nSpecies
      call add_var(info, 'Vn (up,'//trim(cSpecies(i))//')', 'm/s', &
                   'Vv_'//trim(cSpecies(i)), 'Vertical neutral wind for '//trim(cSpecies(i)))
    enddo
    do i = 1, nIons
      call add_var(info, '['//trim(cIons(i))//']', '/m3', &
                   trim(cIons(i)), 'Number density of '//trim(cIons(i)))
    enddo
    call add_var(info, 'eTemperature', 'K', 'Te', 'Electron temperature')
    call add_var(info, 'iTemperature', 'K', 'Ti', 'Ion temperature')
    call add_var(info, 'Vi (east)', 'm/s', 'Vie', 'Eastward ion drift')
    call add_var(info, 'Vi (north)', 'm/s', 'Vin', 'Northward ion drift')
    call add_var(info, 'Vi (up)', 'm/s', 'Viv', 'Vertical ion drift')
  end subroutine register_1dall

  subroutine register_0dall()
    use ModPlanet, only: nSpeciesTotal, nSpecies, nIons, cSpecies, cIons
    type(OutputTypeInfo), pointer :: info
    integer :: i

    call new_output_type('0DALL', 0, 0, info)
    info%needsInterp = .true.
    call add_coord_vars(info)
    call add_var(info, 'Rho', 'kg/m3', 'rho', 'Total neutral mass density')
    do i = 1, nSpeciesTotal
      call add_var(info, '['//trim(cSpecies(i))//']', '/m3', &
                   trim(cSpecies(i)), 'Number density of '//trim(cSpecies(i)))
    enddo
    call add_var(info, 'Temperature', 'K', 'Tn', 'Neutral temperature')
    call add_var(info, 'Vn (east)', 'm/s', 'Ve', 'Eastward neutral wind')
    call add_var(info, 'Vn (north)', 'm/s', 'Vn', 'Northward neutral wind')
    call add_var(info, 'Vn (up)', 'm/s', 'Vv', 'Vertical neutral wind')
    do i = 1, nSpecies
      call add_var(info, 'Vn (up,'//trim(cSpecies(i))//')', 'm/s', &
                   'Vv_'//trim(cSpecies(i)), 'Vertical neutral wind for '//trim(cSpecies(i)))
    enddo
    do i = 1, nIons
      call add_var(info, '['//trim(cIons(i))//']', '/m3', &
                   trim(cIons(i)), 'Number density of '//trim(cIons(i)))
    enddo
    call add_var(info, 'eTemperature', 'K', 'Te', 'Electron temperature')
    call add_var(info, 'iTemperature', 'K', 'Ti', 'Ion temperature')
    call add_var(info, 'Vi (east)', 'm/s', 'Vie', 'Eastward ion drift')
    call add_var(info, 'Vi (north)', 'm/s', 'Vin', 'Northward ion drift')
    call add_var(info, 'Vi (up)', 'm/s', 'Viv', 'Vertical ion drift')
    ! 0DALL extra vars: mixing ratios + heating diagnostics
    do i = 1, nSpecies
      call add_var(info, trim(cSpecies(i))//' Mixing Ratio', '', &
                   'mr_'//trim(cSpecies(i)), 'Volume mixing ratio of '//trim(cSpecies(i)))
    enddo
    call add_var(info, 'RadCooling', 'K/s', 'Q_rad_cool', 'Radiative cooling rate')
    call add_var(info, 'EuvHeating', 'K/s', 'Q_euv', 'EUV heating rate')
    call add_var(info, 'Conduction', 'K/s', 'Q_cond', 'Thermal conduction heating rate')
    call add_var(info, 'Heat Balance Total', 'K/s', 'Q_total', 'Total heat balance')
    call add_var(info, 'Heating Efficiency', '', 'HeatEff', 'Heating efficiency')
  end subroutine register_0dall

  subroutine register_1dglo()
    type(OutputTypeInfo), pointer :: info

    call new_output_type('1DGLO', 1, 0, info)
    call add_coord_vars(info)
    call add_var(info, '6300 A Emission', '', 'I_6300', '6300 Angstrom emission intensity')
    call add_var(info, 'PhotoElectronUp', '', 'PEflux_up', 'Upward photoelectron flux')
    call add_var(info, 'PhotoElectronDown', '', 'PEflux_dn', 'Downward photoelectron flux')
  end subroutine register_1dglo

  subroutine register_1dthm()
    use ModPlanet, only: nSpeciesTotal, cSpecies
    type(OutputTypeInfo), pointer :: info
    integer :: i

    call new_output_type('1DTHM', 1, 0, info)
    call add_coord_vars(info)
    call add_var(info, 'EUV Heating (K/s)', 'K/s', 'Q_euv', 'EUV heating rate')
    call add_var(info, 'Conduction (K/s)', 'K/s', 'Q_cond', 'Thermal conduction heating rate')
    call add_var(info, 'Molecular Conduction (K/s)', 'K/s', 'Q_mol_cond', &
                 'Molecular conduction heating rate')
    call add_var(info, 'Eddy Conduction (K/s)', 'K/s', 'Q_eddy_cond', 'Eddy conduction heating rate')
    call add_var(info, 'Eddy Adiabatic Conduction (K/s)', 'K/s', 'Q_eddy_adiab', &
                 'Eddy adiabatic conduction heating rate')
    call add_var(info, 'Chemical Heating (K/s)', 'K/s', 'Q_chem', 'Chemical heating rate')
    call add_var(info, 'Joule Heating (K/s)', 'K/s', 'Q_joule', 'Joule heating rate')
    call add_var(info, 'NO Cooling (K/s)', 'K/s', 'Q_no_cool', 'NO infrared cooling rate')
    call add_var(info, 'O Cooling (K/s)', 'K/s', 'Q_o_cool', 'O infrared cooling rate')
    call add_var(info, 'Total Abs EUV', 'W/m2', 'EUV_abs', 'Total absorbed EUV flux')
    do i = 1, nSpeciesTotal
      call add_var(info, 'Production Rate '//trim(cSpecies(i)), '/m3/s', &
                   'P_'//trim(cSpecies(i)), 'Production rate of '//trim(cSpecies(i)))
    enddo
    do i = 1, nSpeciesTotal
      call add_var(info, 'Loss Rate '//trim(cSpecies(i)), '/m3/s', &
                   'L_'//trim(cSpecies(i)), 'Loss rate of '//trim(cSpecies(i)))
    enddo
  end subroutine register_1dthm

  subroutine register_1dnew()
    use ModPlanet, only: nSpeciesTotal, nSpecies, nIons, cSpecies, cIons
    type(OutputTypeInfo), pointer :: info
    integer :: i

    call new_output_type('1DNEW', 1, 0, info)
    info%needsInterp = .true.
    ! 1DNEW has different coord layout: Lon, LocalTime, Lat, SZA, Alt
    call add_var(info, 'Longitude', 'rad', 'Lon', 'Geographic longitude')
    call add_var(info, 'Local Time', 'hr', 'LT', 'Local solar time')
    call add_var(info, 'Latitude', 'rad', 'Lat', 'Geographic latitude')
    call add_var(info, 'Solar Zenith Angle', 'rad', 'SZA', 'Solar zenith angle')
    call add_var(info, 'Altitude', 'm', 'Alt', 'Altitude above sea level')
    call add_var(info, 'Rho', 'kg/m3', 'rho', 'Total neutral mass density')
    do i = 1, nSpeciesTotal
      call add_var(info, '['//trim(cSpecies(i))//']', '/m3', &
                   trim(cSpecies(i)), 'Number density of '//trim(cSpecies(i)))
    enddo
    call add_var(info, 'Temperature', 'K', 'Tn', 'Neutral temperature')
    call add_var(info, 'Vn (east)', 'm/s', 'Ve', 'Eastward neutral wind')
    call add_var(info, 'Vn (north)', 'm/s', 'Vn', 'Northward neutral wind')
    call add_var(info, 'Vn (up)', 'm/s', 'Vv', 'Vertical neutral wind')
    do i = 1, nSpecies
      call add_var(info, 'Vn (up,'//trim(cSpecies(i))//')', 'm/s', &
                   'Vv_'//trim(cSpecies(i)), 'Vertical neutral wind for '//trim(cSpecies(i)))
    enddo
    do i = 1, nIons
      call add_var(info, '['//trim(cIons(i))//']', '/m3', &
                   trim(cIons(i)), 'Number density of '//trim(cIons(i)))
    enddo
    call add_var(info, 'eTemperature', 'K', 'Te', 'Electron temperature')
    call add_var(info, 'iTemperature', 'K', 'Ti', 'Ion temperature')
    call add_var(info, 'Vi (east)', 'm/s', 'Vie', 'Eastward ion drift')
    call add_var(info, 'Vi (north)', 'm/s', 'Vin', 'Northward ion drift')
    call add_var(info, 'Vi (up)', 'm/s', 'Viv', 'Vertical ion drift')
    do i = 1, nSpecies
      call add_var(info, trim(cSpecies(i))//' Mixing Ratio', '', &
                   'mr_'//trim(cSpecies(i)), 'Volume mixing ratio of '//trim(cSpecies(i)))
    enddo
  end subroutine register_1dnew

  subroutine register_1dchm()
    type(OutputTypeInfo), pointer :: info

    call new_output_type('1DCHM', 1, 0, info)
    call add_coord_vars(info)
    call add_var(info, 'N2+ + e', '/m3/s', 'R_N2p_e', 'N2+ + e recombination rate')
    call add_var(info, 'O2+ + e', '/m3/s', 'R_O2p_e', 'O2+ + e recombination rate')
    call add_var(info, 'N2+ + O', '/m3/s', 'R_N2p_O', 'N2+ + O reaction rate')
    call add_var(info, 'NO+ + e', '/m3/s', 'R_NOp_e', 'NO+ + e recombination rate')
    call add_var(info, 'N+ + O2', '/m3/s', 'R_Np_O2', 'N+ + O2 reaction rate')
    call add_var(info, 'NO + N', '/m3/s', 'R_NO_N', 'NO + N reaction rate')
    call add_var(info, 'O+ + O2', '/m3/s', 'R_Op_O2', 'O+ + O2 reaction rate')
    call add_var(info, 'N + O2', '/m3/s', 'R_N_O2', 'N + O2 reaction rate')
    call add_var(info, 'O2+ + N', '/m3/s', 'R_O2p_N', 'O2+ + N reaction rate')
    call add_var(info, 'O2+ + NO', '/m3/s', 'R_O2p_NO', 'O2+ + NO reaction rate')
    call add_var(info, 'O2+ + N2', '/m3/s', 'R_O2p_N2', 'O2+ + N2 reaction rate')
    call add_var(info, 'N2+ + O2', '/m3/s', 'R_N2p_O2', 'N2+ + O2 reaction rate')
    call add_var(info, 'N+ + O', '/m3/s', 'R_Np_O', 'N+ + O reaction rate')
    call add_var(info, 'O!+ + N2', '/m3/s', 'R_Op_N2', 'O+ + N2 reaction rate')
    call add_var(info, 'O(1D) + N2', '/m3/s', 'R_O1D_N2', 'O(1D) + N2 quench rate')
    call add_var(info, 'O(1D) + O2', '/m3/s', 'R_O1D_O2', 'O(1D) + O2 quench rate')
    call add_var(info, 'O(1D) + O', '/m3/s', 'R_O1D_O', 'O(1D) + O quench rate')
    call add_var(info, 'O(1D) + e', '/m3/s', 'R_O1D_e', 'O(1D) + e quench rate')
    call add_var(info, 'N(2D) + O2', '/m3/s', 'R_N2D_O2', 'N(2D) + O2 reaction rate')
    call add_var(info, 'O+(2D)+e', '/m3/s', 'R_O2Dp_e', 'O+(2D) + e recombination rate')
    call add_var(info, 'N(2D) + O', '/m3/s', 'R_N2D_O', 'N(2D) + O quench rate')
    call add_var(info, 'N(2D) + e', '/m3/s', 'R_N2D_e', 'N(2D) + e quench rate')
    call add_var(info, 'O+(2D + N2', '/m3/s', 'R_O2Dp_N2', 'O+(2D) + N2 reaction rate')
    call add_var(info, 'O+(2P) + e', '/m3/s', 'R_O2Pp_e', 'O+(2P) + e recombination rate')
    call add_var(info, 'O+(2P) + O', '/m3/s', 'R_O2Pp_O', 'O+(2P) + O reaction rate')
    call add_var(info, 'O+(2P) + N2', '/m3/s', 'R_O2Pp_N2', 'O+(2P) + N2 reaction rate')
    call add_var(info, 'Chemical Heating Rate', 'K/s', 'Q_chem_rate', 'Chemical heating rate')
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
