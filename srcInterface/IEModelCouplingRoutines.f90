
  subroutine run_potential_model(ie, potential)

    use UA_wrapper, only: UA_put_from_ie
    class(ieModel) :: ie
    real, dimension(ie%neednMlts, &
                    ie%neednLats), intent(out) :: potential
    real :: currentTilt = rBadValue, lastTilt = rBadValue
    real :: potVal

    integer :: iMLT, iLat
    potential = 0.0

    call UA_put_from_ie()

    return
  end subroutine run_potential_model

  subroutine initialize()