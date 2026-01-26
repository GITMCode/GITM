submodule (ModIEGITM) ModIEDummy
  implicit none

contains

  !============================================================================
  !! Set model name and interpret. It needs to be swmf though...
  module subroutine set_efield_model(this, efield_model)
    class(ieModel) :: this
    character(len=*), intent(in) :: efield_model
    if (this%iDebugLevel > 0) &
      write(*, *) "=> Setting efield model to : ", trim(efield_model)
    this%iEfield_ = efield_interpret_name(efield_model)
    if (this%iDebugLevel > 0) &
      write(*, *) "=> That is model : ", this%iEfield_
  end subroutine set_efield_model

  module subroutine set_aurora_model(this, aurora_model)
    class(ieModel) :: this
    character(len=*), intent(in) :: aurora_model
    if (this%iDebugLevel > 0) &
      write(*, *) "=> Setting aurora model to : ", trim(aurora_model)
    this%iEfield_ = aurora_interpret_name(aurora_model)
    if (this%iDebugLevel > 0) &
      write(*, *) "=> That is model : ", this%iAurora_
  end subroutine set_aurora_model

  integer function efield_interpret_name(efieldString)
    use ModErrors
    implicit none
    character(len=*), intent(in) :: efieldString
    character(len=len(efieldString)) :: efieldLower
    return
  end function efield_interpret_name
  ! ------------------------------------------------------------
  integer function aurora_interpret_name(auroraString)
    use ModErrors
    implicit none
    character(len=*), intent(in) :: auroraString
    character(len=len(auroraString)) :: auroraLower
    return
  end function aurora_interpret_name
  ! ------------------------------------------------------------
  subroutine run_aurora_model(ie)
    class(ieModel) :: ie
    integer :: iError = 0
    return
  end subroutine run_aurora_model
  ! ------------------------------------------------------------
  subroutine run_aurora_model_electron_diffuse(ie, eflux, avee)
    class(ieModel) :: ie
    real, dimension(ie%neednMlts, &
            ie%neednLats), intent(out) :: eFlux
    real, dimension(ie%neednMlts, &
            ie%neednLats), intent(out) :: AveE
    integer :: iError = 0
    return
  end subroutine run_aurora_model_electron_diffuse
  ! ------------------------------------------------------------
  subroutine run_aurora_model_electron_mono(ie, eflux, avee)
    ! At the moment, this only works for ovation & AMIE Files...
    class(ieModel) :: ie
    real, dimension(ie%neednMlts, &
            ie%neednLats), intent(out) :: eFlux
    real, dimension(ie%neednMlts, &
            ie%neednLats), intent(out) :: AveE
    integer :: iError = 0
    return
  end subroutine run_aurora_model_electron_mono
  ! ------------------------------------------------------------
  subroutine run_aurora_model_electron_wave(ie, eflux, avee)
    ! At the moment, this only works for ovation & AMIE files...
    class(ieModel) :: ie
    real, dimension(ie%neednMlts, &
            ie%neednLats), intent(out) :: eFlux
    real, dimension(ie%neednMlts, &
            ie%neednLats), intent(out) :: AveE
    integer :: iError = 0
    return
  end subroutine run_aurora_model_electron_wave
  ! ------------------------------------------------------------
  subroutine run_aurora_model_ion_diffuse(ie, eflux, avee)
    ! At the moment, this only works for ovation and AMIE files...
    class(ieModel) :: ie
    real, dimension(ie%neednMlts, &
            ie%neednLats), intent(out) :: eFlux
    real, dimension(ie%neednMlts, &
            ie%neednLats), intent(out) :: AveE
    integer :: iError = 0
    return
  end subroutine run_aurora_model_ion_diffuse

  ! ------------------------------------------------------------
  ! These functions are for getting information that was
  ! derived when the auroral model was run
  ! ------------------------------------------------------------

  ! Supply the polar cap to the user:
  ! polarcap = 1 inside the polar cap, and 0 elsewhere.
  subroutine get_polarcap_results(ie, polarcap)
    class(ieModel) :: ie
    real, dimension(ie%neednMlts, &
            ie%neednLats), intent(out) :: polarcap
    return
  end subroutine get_polarcap_results
  ! ------------------------------------------------------------
  subroutine set_bz(this, value)
    class(ieModel), intent(inout) :: this
    real, intent(in) :: value
  end subroutine set_bz
  ! ------------------------------------------------------------
  subroutine set_by(this, value)
    class(ieModel), intent(inout) :: this
    real, intent(in) :: value
  end subroutine set_by
  ! ------------------------------------------------------------
  subroutine set_swv(this, value)
    class(ieModel), intent(inout) :: this
    real, intent(in) :: value
  end subroutine set_swv
  ! ------------------------------------------------------------
  subroutine set_swn(this, value)
    class(ieModel), intent(inout) :: this
    real, intent(in) :: value
  end subroutine set_swn
  ! ------------------------------------------------------------
  subroutine set_hp(this, value)
    class(ieModel), intent(inout) :: this
    real, intent(in) :: value
  end subroutine set_hp
  ! ------------------------------------------------------------
  subroutine set_hpn(this, value)
    class(ieModel), intent(inout) :: this
    real, intent(in) :: value
  end subroutine set_hpn
  ! ------------------------------------------------------------
  subroutine set_hps(this, value)
    class(ieModel), intent(inout) :: this
    real, intent(in) :: value
  end subroutine set_hps
  ! ------------------------------------------------------------
  subroutine set_hp_from_ae(this, value)
    class(ieModel), intent(inout) :: this
    real, intent(in) :: value
  end subroutine set_hp_from_ae
  ! ------------------------------------------------------------
  subroutine set_au(this, value)
    class(ieModel), intent(inout) :: this
    real, intent(in) :: value
  end subroutine set_au
  ! ------------------------------------------------------------
  subroutine set_al(this, value)
    class(ieModel), intent(inout) :: this
    real, intent(in) :: value
  end subroutine set_al
  ! ------------------------------------------------------------
  subroutine set_ae(this, value)
    class(ieModel), intent(inout) :: this
    real, intent(in) :: value
  end subroutine set_ae
  ! ------------------------------------------------------------
  subroutine set_useAeForHp(this)
    class(ieModel), intent(inout) :: this
  end subroutine set_useAeForHp
  ! ------------------------------------------------------------
  subroutine unset_useAeForHp(this)
    class(ieModel), intent(inout) :: this
  end subroutine unset_useAeForHp
  ! ------------------------------------------------------------
  subroutine set_kp(this, value)
    class(ieModel), intent(inout) :: this
    real, intent(in) :: value
  end subroutine set_kp
  ! ------------------------------------------------------------
  module subroutine set_nMlts(this, iValue)
    class(ieModel) :: this
    integer, intent(in) :: iValue
  end subroutine set_nMlts
  ! ------------------------------------------------------------
  module subroutine set_nLats(this, iValue)
    class(ieModel) :: this
    integer, intent(in) :: iValue
  end subroutine set_nLats
  ! ------------------------------------------------------------
  module subroutine set_grid(this, MltsIn, LatsIn)
    !      use ModAMIE_Interface, only: set_interpolation_indices
    class(ieModel) :: this
    real, dimension(this%neednMlts, this%neednLats), intent(in) :: MltsIn
    real, dimension(this%neednMlts, this%neednLats), intent(in) :: LatsIn
  end subroutine set_grid
  ! ------------------------------------------------------------

end submodule ModIEDummy
