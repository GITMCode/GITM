submodule (ModIEGITM) ModIECoupling
    implicit none
    save

contains
  !============================================================================
  module subroutine run_potential_model(ie, potential)

    class(ieModel) :: ie
    real, dimension(ie%neednMlts, &
                    ie%neednLats), intent(out) :: potential
    real :: currentTilt = rBadValue, lastTilt = rBadValue
    real :: potVal

    integer :: iMLT, iLat
    potential = 0.0

    ! call put_from_buffer

    return
  end subroutine run_potential_model
  !============================================================================
!  subroutine initialize()
!  end subroutine initialize
  !============================================================================
end submodule ModIECoupling