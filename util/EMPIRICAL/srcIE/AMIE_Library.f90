!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!----------------------------------------------------------------------

subroutine AMIE_SetFileName(cFileNameIn)
  use ModAMIE_Interface
  implicit none
  character (len=iCharLenIE_), intent(in) :: cFileNameIn
  AMIE_FileName = cFileNameIn
end subroutine AMIE_SetFileName

!----------------------------------------------------------------------

subroutine AMIE_GetFileName(cFileNameOut)
  use ModAMIE_Interface
  implicit none
  character (len=iCharLenIE_), intent(out) :: cFileNameOut
  cFileNameOut = AMIE_FileName
end subroutine AMIE_GetFileName

!----------------------------------------------------------------------

subroutine AMIE_GetnTimes(nTimesOut)
  use ModAMIE_Interface
  implicit none
  integer, intent(out) :: nTimesOut
  nTimesOut = AMIE_nTimes
end subroutine AMIE_GetnTimes

!----------------------------------------------------------------------

subroutine AMIE_GetnMLTs(nMLTsOut)
  use ModAMIE_Interface
  implicit none
  integer, intent(out) :: nMLTsOut
  nMLTsOut = AMIE_nMLTs
end subroutine AMIE_GetnMLTs

!----------------------------------------------------------------------

subroutine AMIE_GetnLats(nLatsOut)
  use ModAMIE_Interface
  implicit none
  integer, intent(out) :: nLatsOut
  nLatsOut = AMIE_nLats
end subroutine AMIE_GetnLats

!----------------------------------------------------------------------

subroutine AMIE_GetLats(LatsOut)

  use ModAMIE_Interface
  implicit none

  real, dimension(AMIE_nMlts, AMIE_nLats, AMIE_nBLKs), intent(out) :: LatsOut
  integer :: i,j

  do i = 1, AMIE_nMLTs
     do j = AMIE_nLats, 1, -1
        LatsOut(i, j, AMIE_North_) = AMIE_Lats(AMIE_nLats - j + 1)
     enddo
     LatsOut(i, 1:AMIE_nLats, AMIE_South_) = -AMIE_Lats(1:AMIE_nLats)
  enddo

end subroutine AMIE_GetLats

!----------------------------------------------------------------------

subroutine AMIE_GetMLTs(MLTsOut)

  use ModAMIE_Interface
  implicit none

  real, dimension(AMIE_nMlts, AMIE_nLats, AMIE_nBLKs), intent(out) :: MLTsOut
  integer :: j

  do j = 1, AMIE_nLats
     MLTsOut(1:AMIE_nMLTs,j,1) = AMIE_MLTs(1:AMIE_nMLTs)
     MLTsOut(1:AMIE_nMLTs,j,2) = AMIE_MLTs(1:AMIE_nMLTs)
  enddo

end subroutine AMIE_GetMLTs

!----------------------------------------------------------------------
! If there is an error, report it and crash
!----------------------------------------------------------------------

subroutine report_error_and_crash(iError, StringSource)

  use ModErrors
  implicit none

  integer, intent(in) :: iError
  character (len=*) :: StringSource

  if (iError /= 0) then
     write(*,*) "Error in AMIE_Library!"
     write(*,*) "Source : ", StringSource
     write(*,*) cErrorCodes(iError)
     isOk = .false.
  endif

end subroutine report_error_and_crash

!----------------------------------------------------------------------
! New method!
! East/West potential(y)
!----------------------------------------------------------------------

subroutine AMIE_GetPotentialY_New(TimeIn, Method, PotentialYOut, iError)

  use ModAMIE_Interface
  implicit none

  real*8, intent(in) :: TimeIn
  integer, intent(in) :: Method
  real, dimension(AMIE_nMLTs, AMIE_nLats, AMIE_nBLKs), intent(out) :: PotentialYOut
  integer, intent(out) :: iError

  call AMIE_GetValue_New2(iPotentialY_, TimeIn, Method, iError)
  call report_error_and_crash(iError, "AMIE_GetPotentialY_New")
  PotentialYOut = AMIE_Interpolated

end subroutine AMIE_GetPotentialY_New

!----------------------------------------------------------------------
! New method!
! Overall potential
!----------------------------------------------------------------------

subroutine AMIE_GetPotential_New(TimeIn, Method, PotentialOut, iError)

  use ModAMIE_Interface
  implicit none

  real*8, intent(in) :: TimeIn
  integer, intent(in) :: Method
  real, dimension(AMIE_nMLTs, AMIE_nLats, AMIE_nBLKs), intent(out) :: PotentialOut
  integer, intent(out) :: iError

  call AMIE_GetValue_New2(iPotential_, TimeIn, Method, iError)
  call report_error_and_crash(iError, "AMIE_GetPotential_New2")
  PotentialOut = AMIE_Interpolated

end subroutine AMIE_GetPotential_New

!----------------------------------------------------------------------
! New method!
! Electron Energy Flux
!----------------------------------------------------------------------

subroutine AMIE_GetEFlux_New(TimeIn, Method, EFluxOut, iError)

  use ModAMIE_Interface
  implicit none

  real*8, intent(in) :: TimeIn
  integer, intent(in) :: Method
  real, dimension(AMIE_nMLTs, AMIE_nLats, AMIE_nBLKs), intent(out) :: EFluxOut
  integer, intent(out) :: iError

  call AMIE_GetValue_New2(iEle_diff_eflux_, TimeIn, Method, iError)
  call report_error_and_crash(iError, "AMIE_GetEFlux_New")
  EFluxOut = AMIE_Interpolated

end subroutine AMIE_GetEFlux_New

!----------------------------------------------------------------------
! New method!
! Electron Average Energy
!----------------------------------------------------------------------

subroutine AMIE_GetAveE_New(TimeIn, Method, AveEOut, iError)

  use ModAMIE_Interface
  implicit none

  real*8, intent(in) :: TimeIn
  integer, intent(in) :: Method
  real, dimension(AMIE_nMLTs, AMIE_nLats, AMIE_nBLKs), intent(out) :: AveEOut
  integer, intent(out) :: iError

  call AMIE_GetValue_New2(iEle_diff_avee_, TimeIn, Method, iError)
  call report_error_and_crash(iError, "AMIE_GetAveE_New")
  AveEOut = AMIE_Interpolated

end subroutine AMIE_GetAveE_New

!----------------------------------------------------------------------
! Ion Energy Flux
!----------------------------------------------------------------------

subroutine AMIE_GetIonEFlux_New(TimeIn, Method, IonEFluxOut, iError)

  use ModAMIE_Interface
  implicit none

  real*8, intent(in) :: TimeIn
  integer, intent(in) :: Method
  real, dimension(AMIE_nMLTs, AMIE_nLats, AMIE_nBLKs), intent(out) :: IonEFluxOut
  integer, intent(out) :: iError

  if (useIons) then
     call AMIE_GetValue_New2(iIon_diff_eFlux_, TimeIn, Method, iError)
     call report_error_and_crash(iError, "AMIE_GetIonEFlux_New")
     IonEFluxOut = AMIE_Interpolated
  else
     IonEFluxOut = rDummyeFlux
  endif

end subroutine AMIE_GetIonEFlux_New

!----------------------------------------------------------------------
! Ion Average Energy
!----------------------------------------------------------------------

subroutine AMIE_GetIonAveE_New(TimeIn, Method, IonAveEOut, iError)

  use ModAMIE_Interface
  implicit none

  real*8, intent(in) :: TimeIn
  integer, intent(in) :: Method
  real, dimension(AMIE_nMLTs, AMIE_nLats, AMIE_nBLKs), intent(out) :: IonAveEOut
  integer, intent(out) :: iError

  if (useIons) then
     call AMIE_GetValue_New2(iIon_diff_avee_, TimeIn, Method, iError)
     call report_error_and_crash(iError, "AMIE_GetIonAveE_New")
     IonAveEOut = AMIE_Interpolated
  else
     IonAveEOut = rDummyIonAveE
  endif

end subroutine AMIE_GetIonAveE_New

!----------------------------------------------------------------------
! Electron Wave Energy Flux
!----------------------------------------------------------------------

subroutine AMIE_GetEleWaveEflux_New(TimeIn, Method, EleWaveEfluxOut, iError)

  use ModAMIE_Interface
  implicit none

  real*8, intent(in) :: TimeIn
  integer, intent(in) :: Method
  real, dimension(AMIE_nMLTs, AMIE_nLats, AMIE_nBLKs), intent(out) :: EleWaveEFluxOut
  integer, intent(out) :: iError

  if (useIons) then
     call AMIE_GetValue_New2(iEle_wave_eflux_, TimeIn, Method, iError)
     call report_error_and_crash(iError, "AMIE_GetEleWaveEflux_New")
     EleWaveEfluxOut = AMIE_Interpolated
  else
     EleWaveEfluxOut = rDummyeFlux
  endif

end subroutine AMIE_GetEleWaveEflux_New

!----------------------------------------------------------------------
! Electron Wave Average Energy
!----------------------------------------------------------------------

subroutine AMIE_GetEleWaveAveE_New(TimeIn, Method, EleWaveAveEOut, iError)

  use ModAMIE_Interface
  implicit none

  real*8, intent(in) :: TimeIn
  integer, intent(in) :: Method
  real, dimension(AMIE_nMLTs, AMIE_nLats, AMIE_nBLKs), intent(out) :: EleWaveAveEOut
  integer, intent(out) :: iError

  if (useIons) then
     call AMIE_GetValue_New2(iEle_wave_avee_, TimeIn, Method, iError)
     call report_error_and_crash(iError, "AMIE_GetEleWaveAveE_New")
     EleWaveAveEOut = AMIE_Interpolated
  else
     EleWaveAveEOut = rDummyAveE
  endif

end subroutine AMIE_GetEleWaveAveE_New

!----------------------------------------------------------------------
! Electron Mono Energy Flux
!----------------------------------------------------------------------

subroutine AMIE_GetEleMonoEflux_New(TimeIn, Method, EleMonoEfluxOut, iError)

  use ModAMIE_Interface
  implicit none

  real*8, intent(in) :: TimeIn
  integer, intent(in) :: Method
  real, dimension(AMIE_nMLTs, AMIE_nLats, AMIE_nBLKs), intent(out) :: EleMonoEFluxOut
  integer, intent(out) :: iError

  if (useIons) then
     call AMIE_GetValue_New2(iEle_mono_eflux_, TimeIn, Method, iError)
     call report_error_and_crash(iError, "AMIE_GetEleMonoEflux_New")
     EleMonoEfluxOut = AMIE_Interpolated
  else
     EleMonoEfluxOut = rDummyeFlux
  endif

end subroutine AMIE_GetEleMonoEflux_New

!----------------------------------------------------------------------
! Electron Mono Average Energy
!----------------------------------------------------------------------

subroutine AMIE_GetEleMonoAveE_New(TimeIn, Method, EleMonoAveEOut, iError)

  use ModAMIE_Interface
  implicit none

  real*8, intent(in) :: TimeIn
  integer, intent(in) :: Method
  real, dimension(AMIE_nMLTs, AMIE_nLats, AMIE_nBLKs), intent(out) :: EleMonoAveEOut
  integer, intent(out) :: iError

  if (useIons) then
     call AMIE_GetValue_New2(iEle_mono_avee_, TimeIn, Method, iError)
     call report_error_and_crash(iError, "AMIE_GetEleMonoAveE_New")
     EleMonoAveEOut = AMIE_Interpolated
  else
     EleMonoAveEOut = rDummyAveE
  endif

end subroutine AMIE_GetEleMonoAveE_New

!----------------------------------------------------------------------
! Get General Values
!   - This has no output (except error)
!   - it sets the variable AMIE_Interpolated 
!----------------------------------------------------------------------

subroutine AMIE_GetValue_New2(iVarIn, TimeIn, Method, iError)

  use ModErrors
  use ModAMIE_Interface
  implicit none

  real*8, intent(in) :: TimeIn
  integer, intent(in) :: iVarIn
  integer, intent(in) :: Method
  integer, intent(out) :: iError

  integer :: iTime, i, j, iLat, iBLK, iL
  logical :: IsDone
  real*8  :: dT, VerySmall = 1.0e-6

  iError = 0

  do iBLK = AMIE_South_, AMIE_North_

     IsDone = .false.
     iTime = 1

     do while (.not. IsDone)
        if (TimeIn - AMIE_Time(iTime,iBLK) < VerySmall) IsDone = .true.
        if ((iTime == AMIE_nTimes) .and. (.not.IsDone)) then
           iTime = iTime + 1 
           IsDone = .true.
        endif
        iTime = iTime + 1
     enddo

     if (iTime <= AMIE_nTimes+1) then

        iTime = iTime - 1

        if (iTime == 1) then

           ! If we are before the start time, allow users to extrapolate
           ! up to 5 dT.

           dT = AMIE_Time(2,iBLK) - AMIE_Time(1,iBLK)
           if (TimeIn + 5*dt < AMIE_Time(1,iBLK)) then
              AMIE_Interpolated = -1.0e32
              iError = ecBeforeStartTime_
              return
           endif
        endif
     else
        dT = AMIE_Time(2,iBLK) - AMIE_Time(1,iBLK)

        ! If we are after the end time, allow users to extrapolate
        ! up to 5 dT.

        if (TimeIn - 5*dt < AMIE_Time(AMIE_nTimes,iBLK)) then
           iTime = AMIE_nTimes
        else
           AMIE_Interpolated = -1.0e32
           iError = ecAfterEndTime_
           return
        endif
     endif

     if ( Method == AMIE_After_ .or. &
          Method == AMIE_Closest_) then

        if (Method == AMIE_Closest_) then
           if (iTime > 1) then
              if (abs(TimeIn-AMIE_Time(iTime,iBLK)) > &
                   abs(TimeIn-AMIE_Time(iTime-1,iBLK))) &
                   iTime = iTime - 1
           endif
        endif

        if (iBLK == AMIE_South_) then
           AMIE_Interpolated(1:AMIE_nMLTs, 1:AMIE_nLats,iBLK) =  &
                AMIE_Storage(iVarIn, 1:AMIE_nMLTs, 1:AMIE_nLats,iTime,iBLK)
        else
           ! Reverse the North block of AMIE data for now...
           do iLat = AMIE_nLats,1,-1
              AMIE_Interpolated(1:AMIE_nMLTs, iLat,iBLK) =  &
                   AMIE_Storage(iVarIn, 1:AMIE_nMLTs, AMIE_nLats - iLat + 1,iTime,iBLK)
           enddo
        endif

     endif

     if (Method == AMIE_Interpolate_) then
        ! This will do extrapolation if it is before the first time
        if (iTime == 1) iTime = iTime + 1
        ! dT is the percentage of the way away from the current point
        dT = (AMIE_Time(iTime,iBLK) - TimeIn) / &
             (AMIE_Time(iTime,iBLK) - AMIE_Time(iTime-1,iBLK))

        ! Use 1-dT for the selected point, since dt = 0 if you are exactly
        ! on the selected point
        if (iBLK == AMIE_South_) then
           AMIE_Interpolated(1:AMIE_nMLTs, 1:AMIE_nLats,iBLK) =  &
                (1.0 - dt)*AMIE_Storage(iVarIn, 1:AMIE_nMLTs, 1:AMIE_nLats,iTime,iBLK)+&
                        dt*AMIE_Storage(iVarIn, 1:AMIE_nMLTs, 1:AMIE_nLats,iTime-1,iBLK)
        else
           ! Reverse the 2nd block of AMIE data for now...
           do iLat = AMIE_nLats,1,-1
              iL = AMIE_nLats-iLat+1
              AMIE_Interpolated(1:AMIE_nMLTs, iLat,iBLK) =  &
                   (1.0 - dt)*AMIE_Storage(iVarIn, 1:AMIE_nMLTs, iL, iTime, iBLK) + &
                           dt*AMIE_Storage(iVarIn, 1:AMIE_nMLTs, iL, iTime-1, iBLK)
           enddo
        endif
     endif

  enddo

end subroutine AMIE_GetValue_New2

!----------------------------------------------------------------------
! Get All Values:
! - potential
! - e- average energy
! - e- total energy flux
! - ion average energy
! - ion total energy flux
!----------------------------------------------------------------------

subroutine get_AMIE_values(rtime)

  use ModEIE_Interface
  implicit none

  real*8, intent(in) :: rtime
  integer :: iError

  call AMIE_GetPotential_New(rtime, EIE_Interpolate_, EIEr3_HavePotential, iError)
  call AMIE_GetAveE_New(rtime, EIE_Closest_, EIEr3_HaveAveE, iError)
  call AMIE_GetEFlux_New(rtime, EIE_Closest_, EIEr3_HaveEFlux, iError)
  if (useIons) then
     call AMIE_GetIonAveE_New(rtime, EIE_Closest_, EIEr3_HaveIonAveE, iError)
     call AMIE_GetIonEFlux_New(rtime, EIE_Closest_, EIEr3_HaveIonEFlux, iError)
  else
     EIEr3_HaveIonAveE = EIE_fill_IonAveE
     EIEr3_HaveIonEFlux = EIE_fill_eFlux
  endif
  if (useWave) then
     call AMIE_GetEleWaveAveE_New(rtime, EIE_Closest_, EIEr3_HaveWaveAveE, iError)
     call AMIE_GetEleWaveEFlux_New(rtime, EIE_Closest_, EIEr3_HaveWaveEFlux, iError)
  else
     EIEr3_HaveWaveAveE = EIE_fill_AveE
     EIEr3_HaveWaveEFlux = EIE_fill_eFlux
  endif
  if (useMono) then
     call AMIE_GetEleMonoAveE_New(rtime, EIE_Closest_, EIEr3_HaveMonoAveE, iError)
     call AMIE_GetEleMonoEFlux_New(rtime, EIE_Closest_, EIEr3_HaveMonoEFlux, iError)
  else
     EIEr3_HaveMonoAveE = EIE_fill_AveE
     EIEr3_HaveMonoEFlux = EIE_fill_eFlux
  endif
  
end subroutine get_AMIE_values

!----------------------------------------------------------------------
! Get the secondary potential (this is the potential for the E/W direction)
!----------------------------------------------------------------------

subroutine get_AMIE_PotentialY(rtime)

  use ModEIE_Interface
  implicit none

  real*8, intent(in) :: rtime
  integer :: iError

  call AMIE_GetPotentialY_New(rtime, EIE_Interpolate_, EIEr3_HavePotential, iError)

end subroutine get_AMIE_PotentialY
