subroutine run_potential_model(ie, potential)

  class(ieModel) :: ie
  real, dimension(ie%neednMlts, &
                  ie%neednLats), intent(out) :: potential
  real :: currentTilt = rBadValue, lastTilt = rBadValue
  real :: potVal

  integer :: iMLT, iLat
  !-----------------------------------------------------------------------------
  if (ie%iDebugLevel > 2) &
          write(*, *) "=> Getting Potential"
  potential = 0.0

  if (ie%isCoupleInitialized) then
    call ie%get_ie_to_ua('pot', potential)
  else
    if (ie%iDebugLevel > 0) &
            write(*, *) "No potential to get before coupling has been "// &
                        "performed with IE for the first time."
  end if

end subroutine run_potential_model
!============================================================================
subroutine run_aurora_model_electron_diffuse(ie, eflux, avee)
  class(ieModel) :: ie
  real, dimension(ie%neednMlts, &
          ie%neednLats), intent(out) :: eFlux
  real, dimension(ie%neednMlts, &
          ie%neednLats), intent(out) :: AveE
  integer :: iError = 0
  !-----------------------------------------------------------------------------
  if (ie%isCoupleInitialized) then
    call ie%get_ie_to_ua('def', eflux)
    call ie%get_ie_to_ua('dae', avee)
  else
    if (ie%iDebugLevel > 0) &
            write(*, *) "No diffuse auroras to get before coupling has been "//&
                        "performed with IE for the first time."
  endif
end subroutine run_aurora_model_electron_diffuse
!============================================================================
subroutine run_aurora_model_electron_mono(ie, eflux, avee)
  ! At the moment, this only works for ovation & AMIE Files...
  class(ieModel) :: ie
  real, dimension(ie%neednMlts, &
          ie%neednLats), intent(out) :: eFlux
  real, dimension(ie%neednMlts, &
          ie%neednLats), intent(out) :: AveE
  integer :: iError = 0
  !-----------------------------------------------------------------------------
  if (ie%isCoupleInitialized) then
    call ie%get_ie_to_ua('mef', eflux)
    call ie%get_ie_to_ua('mae', avee)
  else
    if (ie%iDebugLevel > 0) &
            write(*, *) "No mono auroras to get before coupling has been "// &
                        "performed with IE for the first time."
  endif
end subroutine run_aurora_model_electron_mono
!============================================================================
subroutine run_aurora_model_electron_wave(ie, eflux, avee)
  ! At the moment, this only works for ovation & AMIE files...
  class(ieModel) :: ie
  real, dimension(ie%neednMlts, &
          ie%neednLats), intent(out) :: eFlux
  real, dimension(ie%neednMlts, &
          ie%neednLats), intent(out) :: AveE
  integer :: iError = 0
  !-----------------------------------------------------------------------------
  if (ie%isCoupleInitialized) then
    call ie%get_ie_to_ua('wef', eflux)
    call ie%get_ie_to_ua('wae', avee)
  else
    if (ie%iDebugLevel > 0) &
            write(*, *) "No wave auroras to get before coupling has been "// &
                        "performed with IE for the first time."
  endif
end subroutine run_aurora_model_electron_wave
!============================================================================
subroutine run_aurora_model_ion_diffuse(ie, eflux, avee)
  ! At the moment, this only works for ovation and AMIE files...
  class(ieModel) :: ie
  real, dimension(ie%neednMlts, &
          ie%neednLats), intent(out) :: eFlux
  real, dimension(ie%neednMlts, &
          ie%neednLats), intent(out) :: AveE
  integer :: iError = 0
  !-----------------------------------------------------------------------------
  if (ie%isCoupleInitialized) then
    call ie%get_ie_to_ua('ief', eflux)
    call ie%get_ie_to_ua('iae', avee)
  else
    if (ie%iDebugLevel > 0) &
            write(*, *) "No ion auroras to get before coupling has been "// &
                        "performed with IE for the first time."
  endif
end subroutine run_aurora_model_ion_diffuse
!============================================================================
!! Set model name and interpret. It needs to be swmf though...
subroutine set_efield_model(this, efield_model)
  class(ieModel) :: this
  character(len=*), intent(in) :: efield_model
  ! --------------------------------------------------------------------------
  if (this%iDebugLevel > 0) &
          write(*, *) "=> Setting efield model to : SWMF"
  this%iEfield_ = 1
  if (this%iDebugLevel > 0) &
          write(*, *) "=> That is model : ", this%iEfield_
end subroutine set_efield_model
!============================================================================
subroutine set_aurora_model(this, aurora_model)
  class(ieModel) :: this
  character(len=*), intent(in) :: aurora_model
  ! --------------------------------------------------------------------------
  if (this%iDebugLevel > 0) &
          write(*, *) "=> Setting aurora model to : SWMF"
  this%iAurora_ = 1
  if (this%iDebugLevel > 0) &
          write(*, *) "=> That is model : ", this%iAurora_
end subroutine set_aurora_model
!============================================================================
integer function efield_interpret_name(efieldString)
  use ModErrors
  implicit none
  character(len=*), intent(in) :: efieldString
  character(len=len(efieldString)) :: efieldLower
  character(len=*), parameter :: NameSub="efield_interpret_name"
  ! --------------------------------------------------------------------------
  efield_interpret_name = -1

  efieldLower = efieldString
  call lower_case(efieldLower)
  ! why even do this when it has to be the same no matter what
  if (trim(efieldLower) == "SWMF" .or. (trim(efieldLower) == "swmf")) then
      efield_interpret_name = 1
  else
      call CON_stop(NameSub//" GITM can only use the SWMF IE potential "//&
                             "when a component of the SWMF.")
  endif
  return
end function efield_interpret_name
!============================================================================
integer function aurora_interpret_name(auroraString)
  use ModErrors
  implicit none
  character(len=*), intent(in) :: auroraString
  character(len=len(auroraString)) :: auroraLower
  character(len=*), parameter :: NameSub="aurora_interpret_name"
  ! --------------------------------------------------------------------------
  aurora_interpret_name = -1

  auroraLower = auroraString
  call lower_case(auroraLower)
  ! why even do this when it has to be the same no matter what
  if (trim(auroraLower) == "SWMF" .or. (trim(auroraLower) == "swmf")) then
    aurora_interpret_name = 1
  else
    call CON_stop(NameSub//" GITM can only use the SWMF IE auroras when"//&
            " a component of the SWMF.")
  endif
  return
end function aurora_interpret_name
!============================================================================
subroutine initialize(this)
  use ModErrors
  class(ieModel) :: this
  character(len=iCharLenIE_) :: modelDirTotal
  character(len=iCharLenIE_) :: inFileNameTotal
  character(len=iCharLenIE_) :: name
  integer :: UnitTmp_ = 76
  integer :: iError = 0
  if (this%iDebugLevel >= 0) &
          write(*, *) "> Initializing IE Library"

  if (this%iDebugLevel > 1) &
          write(*, *) "==> Model data directory : ", trim(this%modelDir)

end subroutine initialize
