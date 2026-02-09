! --------------------------------------------------------------------
!\
! This routine finds a point on in the spherical file system, given
! a Theta, Phi:
! LocIn(1) = Phi
! LocIn(2) = Theta
! It returns a 5 element array:
! LocOut(1) = Index of Longitude
! LocOut(2) = Index of Latitude
! LocOut(4) = Multiplication factor for Longitude
! LocOut(5) = Multiplication factor for Latitude
!/
! --------------------------------------------------------------------

subroutine find_ua_point(this, LocIn, LocOut)

    use ModErrors

    implicit none

    class(ieModel) :: this
    real, dimension(2), intent(in)  :: LocIn
    real, dimension(4), intent(out) :: LocOut
    real :: MLTIn, LatIn, MLTUp, MLTDown
    integer :: LatIndex, MltIndex, iPt=0
    real :: MLTMin=0, MLTMax=24.0, MLThalf=12.0, LatMin=-90.0, LatMax=90.0

    logical :: IsFound

    character(len=*), parameter :: NameSub = "find_point"

    ! this was taken from AMIE, not sure will be correct for IE
    ! Might need to convert MLTs to Lons or vice versa
    LocOut = -1.0

    !\
    ! Check to see if the point is even on the grid.
    !/

    MLTIn = mod(LocIn(1) + MLTMax, MLTMax)

    LatIn = LocIn(2)
    if (LatIn > LatMax) then
        LatIn = 2 * LatMax - LatIn
        MLTIn = mod(MLTIn + MLThalf, MLTMax)
    endif
    if (LatIn < LatMin) then
        LatIn = 2 * LatMin - LatIn
        MLTIn = mod(MLTIn + MLThalf, MLTMax)
    endif

    if (MLTIn > MLTMax .or. MLTIn < MLTMin .or. &
        LatIn > LatMax .or. LatIn < LatMin) then
            ! should update to be var based
        call set_error("Input lat / mlt is outside of -90-90 and 0-24 range! " &
                // NameSub)
        return
    endif

    ! only works for regular grid, should eventually check in case IE uses
    ! something else

    do iPt=1, this%havenMlts
        if ((this%haveMLTs(iPt+1, 2) > MLTIn) .and. (this%haveMLTs(iPt, 2) <= MLTIn)) &
        MltIndex=iPt
    enddo

    do iPt=1, this%havenLats
        if ((this%haveLats(2, iPt+1) < LatIn) .and. (this%haveLats(2, iPt) >= LatIn)) &
        LatIndex=iPt
    enddo
    
    ! check my work to be sure!
    MLTUp = this%haveMLTs(MltIndex+1, LatIndex)
    MLTDown = this%haveMLTs(MltIndex, LatIndex)
 
    ! store needed values
    LocOut(1) = MltIndex
    LocOut(2) = LatIndex
    
    LocOut(3) = (MLTIn - MLTDown)/(MLTUp - MLTDown)
    ! To keep the MLT and LAT ratios uses consistent, we need to calculate
    ! the ratios backward, since lat shrinks with increasing index:
    LocOut(4) = (this%haveLats(MltIndex,LatIndex) - LatIn)/ &
            (this%haveLats(MltIndex,LatIndex) - &
            this%haveLats(MltIndex,LatIndex+1))

    if (this%iDebugLevel > 4) &
        write(*, *) 'file point!', LatIn, MltIn, LocOut
    if (this%iDebugLevel > 7) then
        if (LatIn > this%haveLats(MltIndex,LatIndex)) &
        write(*,*) "LatIn, LatUp", LatIn, this%haveLats(MltIndex,LatIndex)
        if (LatIn <= this%haveLats(MltIndex,LatIndex+1)) &
        write(*,*) "LatIn, LatDown", LatIn, this%haveLats(MltIndex,LatIndex+1)
        if (MLTIn >= MLTUp) write(*,*) "MLTIn, MLTUp", MltIn, MltUP
        if (MltIn < MltDown) write(*,*) "MltIn, MltDown", MltIn, MltDown
    endif

end subroutine find_ua_point
!==============================================================================
subroutine set_ie_ua_interpolation_indices(this, mltsIn, latsIn)

    class(ieModel) :: this
    real, intent(in) :: mltsIn(this%neednMlts, this%neednLats), &
            latsIn(this%neednMlts, this%neednLats)
    real, dimension(2) :: mlt_and_lat
    real, dimension(4) :: interpolation_info
    integer :: iError, iLat, iMlt, im, il

    if (this%iDebugLevel > 2) &
            write(*, *) "=> Getting IE->UA interpolation indices", &
                    this%neednMlts, this%neednLats
    if (allocated(this%IeUaInterpolationIndices)) then
        deallocate(this%IeUaInterpolationIndices)
        deallocate(this%IeUaInterpolationRatios)
    endif
    allocate(this%IeUaInterpolationIndices(this%neednMlts, this%neednLats, 2), &
            stat=iError)
    if (iError /= 0) then
        call set_error("Error allocating IeUaInterpolationIndices!")
        return
    endif
    allocate(this%IeUaInterpolationRatios(this%neednMlts, this%neednLats, 2), &
            stat=iError)
    if (iError /= 0) then
        call set_error("Error allocating IeUaInterpolationRatios!")
        return
    endif

    this%IeUaInterpolationIndices = -1

    do iMlt = 1, this%neednMlts
        do iLat = 1, this%neednLats
            mlt_and_lat(1) = mltsIn(iMlt, iLat)
            mlt_and_lat(2) = latsIn(iMlt, iLat)
            call this%find_ua_point(mlt_and_lat, interpolation_info)
            if (iError == 0) then
                this%IeUaInterpolationIndices(iMlt, iLat, 1:2) = &
                        interpolation_info(1:2)
                this%IeUaInterpolationRatios(iMlt, iLat, 1:2) = &
                        interpolation_info(3:4)

                iM = this%IeUaInterpolationIndices(iMLT, iLat, 1)
                iL = this%IeUaInterpolationIndices(iMLT, iLat, 2)   
            else
                this%IeUaInterpolationIndices(iMlt, iLat, 1:2) = -1
            endif
        enddo
    enddo

end subroutine set_ie_ua_interpolation_indices
!==============================================================================
subroutine get_ie_values_for_ua(this, iVarToGetIn, valueOut)

    class(ieModel) :: this
    character(len=*), intent(in) :: iVarToGetIn
    real, intent(inout) ::  valueOut(this%neednMlts, this%neednLats)
    integer :: iMlt, iLat, iM, iL, iFile
    real :: dM, dL
    logical :: isNorth
    real :: current_var(this%havenMlts, this%havenLats)
    character(len=*), parameter :: NameSub = "get_ie_values_for_ua"

    valueOut = 0.0

    do iMlt = 1, this%neednMlts; do iLat = 1, this%neednLats

        iM = this%IeUaInterpolationIndices(iMLT, iLat, 1)
        iL = this%IeUaInterpolationIndices(iMLT, iLat, 2)
        dM = this%IeUaInterpolationRatios(iMLT, iLat, 1)
        dL = this%IeUaInterpolationRatios(iMLT, iLat, 2)
        
        if(iM < 0) cycle

        ! I don't think this is necessary when coupling from IE
        ! I hope I'm not wrong
        !! We need to figure out which block we should use, since
        !! the find indices could only be called once.  Ugh.
        !! At this point, assume all files / blocks have the same grid
        !if (allFiles(iB)%nTimesInMemory == 0) then
        !    isNorth = allFiles(iB)%isNorth
        !    ! We have nothing in memory for this file!
        !    iB = 0
        !    do iFile = 1, nFiles
        !        if ((allFiles(iFile)%nTimesInMemory > 0) .and. &
        !                (allFiles(iFile)%isNorth .eqv. isNorth)) &
        !                iB = iFile
        !    enddo
        !endif

        ! I don't know enough fortran to know if we can do this better
        ! Please tell me we can do better, this *is* a cry for help
        select case (iVarToGetIn)
            case("pot")
                current_var = this%havePotential
            case("def")
                current_var = this%haveDiffuseEeFlux
            case("dae")
                current_var = this%haveDiffuseEaveE
            case("mef")
                current_var = this%haveMonoEeFlux
            case("mae")
                current_var = this%haveMonoEaveE
            case("wef")
                current_var = this%haveWaveEeFlux
            case("wae")
                current_var = this%haveWaveEaveE
            case("ief")
                current_var = this%haveDiffuseIeFlux
            case("iae")
                current_var = this%haveDiffuseIaveE
            case default
                call CON_stop(NameSub//": "//iVarToGetIn//&
                        "is not a valid variable to get")

        end select
        ValueOut(iMLT, iLat) = &
                    (1.0 - dM)*(1.0 - dL)*current_var(iL, iM) + &
                    (1.0 - dM)*(dL)*current_var(iL + 1, iM) + &
                    (dM)*(dL)*current_var(iL + 1, iM + 1) + &
                    (dM)*(1.0 - dL)*current_var(iL, iM + 1)
    enddo; enddo
    return
end subroutine get_ie_values_for_ua
