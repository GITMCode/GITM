!============================================================================
subroutine set_nMlts(this, iValue)
    class(ieModel) :: this
    integer, intent(in) :: iValue
    if (this%iDebugLevel > 3) &
            write(*, *) "=> Setting neednMlts to : ", iValue
    this%neednMlts = iValue
end subroutine set_nMlts
!============================================================================
subroutine set_nLats(this, iValue)
    class(ieModel) :: this
    integer, intent(in) :: iValue
    if (this%iDebugLevel > 3) &
            write(*, *) "=> Setting neednMlats to : ", iValue
    this%neednLats = iValue
end subroutine set_nLats
!============================================================================
subroutine set_grid(this, MltsIn, LatsIn)
    class(ieModel) :: this
    real, dimension(this%neednMlts, this%neednLats), intent(in) :: MltsIn
    real, dimension(this%neednMlts, this%neednLats), intent(in) :: LatsIn
    if (this%iDebugLevel > 2) &
            write(*, *) "=> Setting Grid"
    call this%mlts(MltsIn)
    call this%lats(LatsIn)
end subroutine set_grid
!============================================================================
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
                this%needMlts(iMlt, iLat) = mod((MltsIn(iMlt, iLat) + 24.0), &
                        24.0)
            enddo
        enddo
    endif
end subroutine set_mlts
!============================================================================
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
!============================================================================
subroutine set_ie_nMlts(this, iValue)
    class(ieModel) :: this
    integer, intent(in) :: iValue
    if (this%iDebugLevel > 3) &
            write(*, *) "=> Setting havenMlts to : ", iValue
    this%havenMlts = iValue
end subroutine set_ie_nMlts
!============================================================================
subroutine set_ie_nLats(this, iValue)
    class(ieModel) :: this
    integer, intent(in) :: iValue
    if (this%iDebugLevel > 3) &
            write(*, *) "=> Setting havenMlats to : ", iValue
    this%havenLats = iValue
end subroutine set_ie_nLats
!============================================================================
subroutine set_ie_nBlks(this, iValue)
    class(ieModel) :: this
    integer, intent(in) :: iValue
    if (this%iDebugLevel > 3) &
            write(*, *) "=> Setting havenMlats to : ", iValue
    this%havenBlks = iValue
end subroutine set_ie_nBlks
!============================================================================
subroutine set_ie_mlts(this, MltsIn)
    class(ieModel), intent(inout) :: this
    real, dimension(this%havenMlts - 1), intent(in) :: MltsIn
    integer :: iError = 0, iLat

    if(allocated(this%haveMlts)) deallocate(this%haveMlts)
    allocate(this%haveMlts(this%havenMlts, this%havenLats, this%havenBlks))

    ! Set MLT for both hemispheres:
    do iLat=1,this%havenLats
        this%haveMlts(1,iLat,1) = mod(MltsIn(this%havenMLTs-2) - 180.0, 360.0)
        this%HaveMlts(2:this%havenMLTs,iLat,1) = mod(MltsIn + 180.0, 360.0)
        this%HaveMlts(1,iLat,2) = mod(MltsIn(this%havenMLTs-2) - 180.0, 360.0)
        this%HaveMlts(2:this%havenMLTs,iLat,2) = mod(MltsIn + 180.0, 360.0)
    enddo

end subroutine set_ie_mlts
!============================================================================
subroutine set_ie_lats(this, LatsIn)
    class(ieModel) :: this
    real, dimension(this%havenLats), intent(in) :: LatsIn
    integer :: iError = 0, iMlt, iLat, ii

    if(allocated(this%haveLats)) deallocate(this%haveLats)
    allocate(this%haveLats(this%havenMlts, this%havenLats, this%havenBlks))

    do iMlt=1,this%havenMLTs
        ! Northern hemisphere lats:
        this%haveLats(iMlt,:,1) = LatsIn
        ! Southern hemisphere lats (flipped):
        do iLat=1, this%havenLats
            ii = this%havenLats - iLat + 1
            this%haveLats(iMlt,iLat,2) = 0 - LatsIn(ii)
        enddo
    end do

end subroutine set_ie_lats
!============================================================================
subroutine initializeCouple(this, UseDiffuse, UseMono, UseWave, UseIon)
    class(ieModel) :: this
    logical, intent(in) :: UseDiffuse, UseMono, UseWave, UseIon
    ! Resizes all input arrays from IE to be correct
    ! Field-aligned Currents:
    if(allocated(this%haveFac)) deallocate(this%haveFac)
    allocate(this%haveFac(this%havenMlts, this%havenLats, this%havenBlks))
    this%haveFac = 0
    ! Potentials:
    if(allocated(this%havePotential)) deallocate(this%havePotential)
    allocate(this%havePotential(this%havenMlts, this%havenLats, this%havenBlks))
    this%havePotential = 0
    ! Electron diffuse:
    if(allocated(this%haveDiffuseEeFlux)) deallocate(this%haveDiffuseEeFlux)
    if(allocated(this%haveDiffuseEAveE)) deallocate(this%haveDiffuseEAveE)
    if(UseDiffuse) then
        allocate(this%haveDiffuseEeFlux(this%havenMlts, &
                this%havenLats, this%havenBlks))
        allocate(this%haveDiffuseEAveE(this%havenMlts, &
                this%havenLats, this%havenBlks))
        this%haveDiffuseEeFlux = 0
        this%haveDiffuseEAveE = 0
    end if
    ! Ion diffuse:
    if(allocated(this%haveDiffuseIeFlux)) deallocate(this%haveDiffuseIeFlux)
    if(allocated(this%haveDiffuseIAveE)) deallocate(this%haveDiffuseIAveE)
    if(UseIon) then
        allocate(this%haveDiffuseIeFlux(this%havenMlts, &
                this%havenLats, this%havenBlks))
        allocate(this%haveDiffuseIAveE(this%havenMlts, &
                this%havenLats, this%havenBlks))
        this%haveDiffuseIeFlux = 0
        this%haveDiffuseIAveE = 0
    end if
    ! Discrete or Monoenergetic:
    if(allocated(this%haveMonoEeFlux)) deallocate(this%haveMonoEeFlux)
    if(allocated(this%haveMonoEAveE)) deallocate(this%haveMonoEAveE)
    if(UseMono) then
        allocate(this%haveMonoEeFlux(this%havenMlts, &
                this%havenLats, this%havenBlks))
        allocate(this%haveMonoEAveE(this%havenMlts, &
                this%havenLats, this%havenBlks))
        this%haveMonoEeFlux = 0
        this%haveMonoEAveE = 0
    end if
    ! Broadband or Wave-drive:
    if(allocated(this%haveWaveEeFlux)) deallocate(this%haveWaveEeFlux)
    if(allocated(this%haveWaveEAveE)) deallocate(this%haveWaveEAveE)
    if(UseWave) then
        allocate(this%haveWaveEeFlux(this%havenMlts, &
                this%havenLats, this%havenBlks))
        allocate(this%haveWaveEAveE(this%havenMlts, &
                this%havenLats, this%havenBlks))
        this%haveWaveEeFlux = 0
        this%haveWaveEAveE = 0
    end if
    ! Is Polar Cap (1 if is polar cap, 0 otherwise):
    if(allocated(this%havePolarCap)) deallocate(this%havePolarCap)
    allocate(this%havePolarCap(this%havenMlts, this%havenLats, this%havenBlks))
    this%havePolarCap = 0

    this%isCoupleInitalized = .true.

end subroutine initializeCouple