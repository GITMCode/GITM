! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE
!
! ModOutputVarRegistry.f90
!
! Central metadata registry for container outputs: one canonical (units,
! longName) pair per physical variable, keyed by the variable's short name.
!
! Why this exists: before the registry, every define_schema_<type> re-typed the
! units and longName strings for the variables it happened to contain, so the
! same physical quantity could be described differently in two output types
! ('Rho' had two longNames; 'Altitude' still does).  Naming a variable once,
! here, makes that drift structurally impossible.
!
! Contract: put_field/put_axis (ModOutputContainer) consult lookup_var_meta on
! FIRST TOUCH only, so this costs nothing per timestep.  Names that cannot live
! in a static table -- the computed species/ion names built in loops, e.g.
! '['//trim(cSpecies(i))//']' -- are not registered; those call sites pass
! units=/longName= override arguments to put_field instead.  A name that is
! neither registered nor given an explicit units= is an error, so a variable can
! never reach a file with silently empty metadata.
!
! Scope note: this table currently covers the fixed (non-computed) variables of
! the output types converted to the one-shot put API.  It grows as more types
! are converted; entries are shared across types by name, which is the point.

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
  ! irrelevant -- lookup is by name.
  !
  ! 'Altitude' carries 3DALL's longName.  2DANC currently writes the same
  ! variable as plain 'Altitude'; that call site overrides longName= to keep
  ! output byte-identical, and the override goes away when the drifted
  ! longNames are canonicalized as a deliberate, one-time output change.
  type(VarMeta), parameter :: VarTable(25) = [ &
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
                              VarMeta('OCool_int', 'W/m2', 'Alt-integrated O cooling')]

contains

  ! ---------------------------------------------------------------------------
  ! lookup_var_meta — canonical metadata for name.
  !
  ! Returns .true. and fills units/longName when name is registered; returns
  ! .false. and leaves them untouched otherwise (computed species/ion names).
  ! Callers supply their own metadata in that case.
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
