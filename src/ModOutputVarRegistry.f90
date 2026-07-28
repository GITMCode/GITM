! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE
!
! ModOutputVarRegistry.f90
!
! Central metadata registry for container outputs: one canonical (units,
! longName) pair per physical variable, keyed by the variable's short name.
! Before this, each define_schema_<type> re-typed its own units/longName
! strings, so the same quantity could be described differently in two types
! ('Rho' had two longNames).  Naming a variable once here avoids that drift.
!
! put_field/put_axis consult lookup_var_meta on first touch only, so it costs
! nothing per timestep.  Computed names built in loops (e.g.
! '['//trim(cSpecies(i))//']') are not in the table; those call sites pass
! units=/longName= to put_field instead.  A name that is neither registered nor
! given an explicit units= is an error, so nothing reaches a file with empty
! metadata.  The table covers the fixed variables of the converted one-shot
! types and grows as more are converted; entries are shared across types by name.

module ModOutputVarRegistry

  implicit none

  private
  public :: lookup_var_meta

  ! Component lengths match ModOutputContainer's Variable type.
  integer, parameter :: RegNameLen = 64
  integer, parameter :: RegUnitsLen = 32
  integer, parameter :: RegLongNameLen = 128

  type :: VarMeta
    character(len=RegNameLen)     :: name = ''
    character(len=RegUnitsLen)    :: units = ''
    character(len=RegLongNameLen) :: longName = ''
  end type VarMeta

  ! One row per PHYSICAL variable, not per schema appearance.  Order is
  ! irrelevant; lookup is by name.
  type(VarMeta), parameter :: &
    VarTable(*) = [ &
    ! -- coordinates (shared by every geographic type) ----------------------
    VarMeta('Longitude', 'degrees_east', 'Geographic longitude'), &
    VarMeta('Latitude', 'degrees_north', 'Geographic latitude'), &
    VarMeta('Altitude', 'm', 'Altitude above surface'), &
    ! -- neutral state ------------------------------------------------------
    VarMeta('Rho', 'kg/m3', 'Total neutral mass density'), &
    VarMeta('Temperature', 'K', 'Neutral temperature'), &
    VarMeta('Vn (east)', 'm/s', 'Eastward neutral wind'), &
    VarMeta('Vn (north)', 'm/s', 'Northward neutral wind'), &
    VarMeta('Vn (up)', 'm/s', 'Vertical neutral wind'), &
    ! -- ionosphere ---------------------------------------------------------
    VarMeta('eTemperature', 'K', 'Electron temperature'), &
    VarMeta('iTemperature', 'K', 'Ion temperature'), &
    VarMeta('Vi (east)', 'm/s', 'Eastward ion drift'), &
    VarMeta('Vi (north)', 'm/s', 'Northward ion drift'), &
    VarMeta('Vi (up)', 'm/s', 'Vertical ion drift'), &
    ! -- 2D ancillary -------------------------------------------------------
    VarMeta('LT', 'hr', 'Local Time'), &
    VarMeta('SZA', 'rad', 'Solar zenith angle'), &
    VarMeta('TEC', 'TECU', 'Total electron content'), &
    VarMeta('JH_int', 'W/m2', 'Alt-integrated Joule heating'), &
    VarMeta('HT_int', 'W/m2', 'Alt-integrated Heat transfer'), &
    VarMeta('EUV_int', 'W/m2', 'Alt-integrated EUV heating'), &
    VarMeta('PE_int', 'W/m2', 'Alt-integrated PE heating'), &
    VarMeta('Chem_int', 'W/m2', 'Alt-integrated Chem heating'), &
    VarMeta('RadCool_int', 'W/m2', 'Alt-integrated Rad cooling'), &
    VarMeta('CO2Cool_int', 'W/m2', 'Alt-integrated CO2 cooling'), &
    VarMeta('NOCool_int', 'W/m2', 'Alt-integrated NO cooling'), &
    VarMeta('OCool_int', 'W/m2', 'Alt-integrated O cooling'), &
    ! -- electrodynamics (introduced by 3DION) ------------------------------
    VarMeta('Ed1', 'V/m', 'Electric field component 1 (magnetic east)'), &
    VarMeta('Ed2', 'V/m', 'Electric field component 2 (magnetic north)'), &
    VarMeta('Je1', 'A/m2', 'Current density component 1 (magnetic east)'), &
    VarMeta('Je2', 'A/m2', 'Current density component 2 (magnetic north)'), &
    VarMeta('Magnetic Latitude', 'deg', 'Magnetic latitude'), &
    VarMeta('Magnetic Longitude', 'deg', 'Magnetic longitude'), &
    VarMeta('B.F. East', 'T', 'Eastward magnetic field'), &
    VarMeta('B.F. North', 'T', 'Northward magnetic field'), &
    VarMeta('B.F. Vertical', 'T', 'Vertical magnetic field'), &
    VarMeta('B.F. Magnitude', 'T', 'Magnetic field magnitude'), &
    VarMeta('Potential', 'V', 'Electric potential'), &
    VarMeta('E.F. East', 'V/m', 'Eastward electric field'), &
    VarMeta('E.F. North', 'V/m', 'Northward electric field'), &
    VarMeta('E.F. Vertical', 'V/m', 'Vertical electric field'), &
    VarMeta('E.F. Magnitude', 'V/m', 'Electric field magnitude'), &
    VarMeta('IN Collision Freq', 'Hz', 'Ion-neutral collision frequency'), &
    VarMeta('PressGrad (east)', 'm/s2', 'Neutral Pressure gradient acceleration eastward'), &
    VarMeta('PressGrad (north)', 'm/s2', 'Neutral Pressure gradient acceleration northward'), &
    VarMeta('IonPressGrad (east)', 'm/s2', 'Ion Pressure gradient acceleration eastward'), &
    VarMeta('IonPressGrad (north)', 'm/s2', 'Ion Pressure gradient acceleration northward'), &
    VarMeta('IonPressGrad (up)', 'm/s2', 'Ion Pressure gradient acceleration upward'), &
    ! -- 2D electrodynamics (introduced by 2DGEL) ---------------------------
    VarMeta('PedCond', 'S', 'Height-integrated Pedersen conductance'), &
    VarMeta('HalCond', 'S', 'Height-integrated Hall conductance'), &
    VarMeta('AveE', 'keV', 'Diffuse electron precipitation average energy'), &
    VarMeta('eFlux', 'ergs/cm2/s', 'Diffuse electron precipitation energy flux'), &
    VarMeta('AveE_W', 'keV', 'Wave-driven electron precipitation average energy'), &
    VarMeta('eFlux_W', 'ergs/cm2/s', 'Wave-driven electron precipitation energy flux'), &
    VarMeta('AveE_M', 'keV', 'Monoenergetic electron precipitation average energy'), &
    VarMeta('eFlux_M', 'ergs/cm2/s', 'Monoenergetic electron precipitation energy flux'), &
    VarMeta('AveE_I', 'keV', 'Ion precipitation average energy'), &
    VarMeta('eFlux_I', 'ergs/cm2/s', 'Ion precipitation energy flux'), &
    VarMeta('DivJuAlt', 'A/m2', 'Height-integrated divergence of horizontal current'), &
    VarMeta('PedFLCond', 'S', 'Field-line integrated Pedersen conductance'), &
    VarMeta('HalFLCond', 'S', 'Field-line integrated Hall conductance'), &
    VarMeta('DivJuFL', 'A/m2', 'Field-line divergence of current'), &
    VarMeta('FL_length', 'm', 'Field line length'), &
    ! -- thermal budget (introduced by 3DTHM) -------------------------------
    VarMeta('EUV Heating Rate', 'K/s', 'EUV heating rate'), &
    VarMeta('Conduction', 'K/s', 'Thermal conduction heating rate'), &
    VarMeta('Molecular Conduction', 'K/s', 'Molecular conduction heating rate'), &
    VarMeta('Chemical Heating Rate', 'K/s', 'Chemical heating rate'), &
    VarMeta('Joule Heating Rate', 'K/s', 'Joule heating rate'), &
    VarMeta('NO Cooling Rate', 'K/s', 'NO infrared cooling rate'), &
    VarMeta('O Cooling Rate', 'K/s', 'O infrared cooling rate'), &
    VarMeta('Total Abs EUV', 'W/m2', 'Total absorbed EUV flux'), &
    VarMeta('Cp', 'J/kg/K', 'Specific heat at constant pressure'), &
    VarMeta('Sigma Ped', 'S/m', 'Pedersen conductivity'), &
    VarMeta('Ionization Rate O_3P', '/s', 'Ionization rate for O(3P)'), &
    VarMeta('Ionization Rate O2', '/s', 'Ionization rate for O2'), &
    VarMeta('Ionization Rate N2', '/s', 'Ionization rate for N2'), &
    ! -- chemistry (introduced by 3DCHM) -------------------------------------
    ! Same energy source as 3DTHM's 'Chemical Heating Rate' but volumetric
    ! (eV/m3/s) instead of a temperature tendency (K/s).
    VarMeta('Chemical Heating', 'eV/m3/s', 'Volumetric chemical heating rate'), &
    ! -- horizontal momentum budget (introduced by 3DMOH) --------------------
    VarMeta('Visc_Ve_rshear', 'm/s2', 'Viscous east shear accel'), &
    VarMeta('Visc_Vn_rshear', 'm/s2', 'Viscous north shear accel'), &
    VarMeta('IonDrag_east', 'm/s2', 'Ion drag east accel'), &
    VarMeta('IonDrag_north', 'm/s2', 'Ion drag north accel'), &
    VarMeta('Horiz_Adv_Ve', 'm/s2', 'Horiz advection east accel'), &
    VarMeta('Horiz_Adv_Vn', 'm/s2', 'Horiz advection north accel'), &
    VarMeta('Coriolis_east', 'm/s2', 'Coriolis east accel'), &
    VarMeta('Coriolis_north', 'm/s2', 'Coriolis north accel'), &
    VarMeta('Centrifugal_north', 'm/s2', 'Centrifugal north accel'), &
    ! -- vertical momentum budget (introduced by 3DMOV) -----------------------
    VarMeta('Coriolis_up', 'm/s2', 'Coriolis upward accel'), &
    VarMeta('Centrif_up', 'm/s2', 'Centrifugal upward accel'), &
    VarMeta('EffectiveGravity', 'm/s2', 'Effective gravity accel'), &
    ! -- glow placeholders (introduced by 3DGLO) -------------------------------
    VarMeta('6300 A Emission', '', '6300 Angstrom emission'), &
    VarMeta('PhotoElectronUp', '', 'Photoelectron flux up'), &
    VarMeta('PhotoElectronDown', '', 'Photoelectron flux down')]

contains

  ! ---------------------------------------------------------------------------
  ! lookup_var_meta: canonical metadata for name.
  !
  ! Returns .true. and fills units/longName when name is registered. Returns
  ! .false. otherwise (computed species/ion names) — units/longName are
  ! intent(out), so on a miss they come back undefined, not preserved;
  ! callers must check the result before using them.
  ! ---------------------------------------------------------------------------
  logical function lookup_var_meta(name, units, longName)
    character(len=*), intent(in)  :: name
    character(len=*), intent(out) :: units
    character(len=*), intent(out) :: longName
    integer :: i

    lookup_var_meta = .false.
    do i = 1, size(VarTable)
      if (trim(VarTable(i)%name) == trim(name)) then
        units = trim(VarTable(i)%units)
        longName = trim(VarTable(i)%longName)
        lookup_var_meta = .true.
        return
      endif
    enddo
  end function lookup_var_meta

end module ModOutputVarRegistry
