submodule (ModIEGITM) ModGridRoutines
    implicit none

    contains
    ! ------------------------------------------------------------
    ! Set the number of mlts in the grid to get
    subroutine set_nMlts(this, iValue)
        class(ieModel) :: this
        integer, intent(in) :: iValue
        if (this%iDebugLevel > 3) &
            write(*, *) "=> Setting neednMlts to : ", iValue
        this%neednMlts = iValue
    end subroutine set_nMlts

    ! ------------------------------------------------------------
    ! Set the number of lats in the grid to get
    subroutine set_nLats(this, iValue)
            class(ieModel) :: this
            integer, intent(in) :: iValue
            if (this%iDebugLevel > 3) &
            write(*, *) "=> Setting neednMlats to : ", iValue
    this%neednLats = iValue
    end subroutine set_nLats

    ! ------------------------------------------------------------
    ! set grid
    subroutine set_grid(this, MltsIn, LatsIn)
            class(ieModel) :: this
    real, dimension(this%neednMlts, this%neednLats), intent(in) :: MltsIn
    real, dimension(this%neednMlts, this%neednLats), intent(in) :: LatsIn
    if (this%iDebugLevel > 2) &
    write(*, *) "=> Setting Grid"
    call this%mlts(MltsIn)
            call this%lats(LatsIn)

    if (this%iEfield_ == iAmiePot_ .or. &
    this%iAurora_ == iAmieAur_) then
!    call set_interpolation_indices( &
!    this%neednMlts, this%neednLats, this%needMlts, this%needLats)
            endif
    end subroutine set_grid

            ! ------------------------------------------------------------
            ! set mlts
    subroutine set_mlts(this, MltsIn)
            class(ieModel), intent(inout) :: this
        real, dimension(this%neednMlts, this%neednLats), intent(in) :: MltsIn
        integer :: iError = 0, iMlt, iLat
        if (this%iDebugLevel > 2) &
        write(*, *) "=> Setting Mlts"
        if (allocated(this%needMlts)) deallocate(this%needMlts)
            allocate(this%needMlts(this%neednMlts, this%neednLats), stat=iError)
            if (iError /= 0) then
            call set_error("Error allocating Mlts!")
            return
            else
            do iMlt = 1, this%neednMlts
            do iLat = 1, this%neednLats
    this%needMlts(iMlt, iLat) = mod((MltsIn(iMlt, iLat) + 24.0), 24.0)
            enddo
            enddo
            endif
    end subroutine set_mlts

            ! ------------------------------------------------------------
            ! set lats
    subroutine set_lats(this, LatsIn)
            class(ieModel) :: this
            real, dimension(this%neednMlts, this%neednLats), intent(in) :: LatsIn
            integer :: iError = 0, iMlt, iLat
            if (this%iDebugLevel > 2) &
            write(*, *) "=> Setting Lats"
            if (allocated(this%needLats)) deallocate(this%needLats)
            allocate(this%needLats(this%neednMlts, this%neednLats), stat=iError)
                    if (iError /= 0) then
            call set_error("Error allocating Lats!")
            return
            else
            do iMlt = 1, this%neednMlts
    do iLat = 1, this%neednLats
    this%needLats(iMlt, iLat) = LatsIn(iMlt, iLat)
    enddo
    enddo
    endif
    end subroutine set_lats

end submodule ModGridRoutines