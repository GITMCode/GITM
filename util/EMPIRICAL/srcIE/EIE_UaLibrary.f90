!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf

!----------------------------------------------------------------------

subroutine UA_SetnMLTs(nMLTsIn)
  use ModEIE_Interface
  implicit none
  integer, intent(in) :: nMLTsIn
  UAi_NeednMLTs = nMLTsIn
end subroutine UA_SetnMLTs

!----------------------------------------------------------------------

subroutine UA_SetnLats(nLatsIn)
  use ModEIE_Interface
  implicit none
  integer, intent(in) :: nLatsIn
  UAi_NeednLats = nLatsIn
end subroutine UA_SetnLats

!----------------------------------------------------------------------

!subroutine UA_SetTime(TimeIn)
!!  use ModKind
!  use ModEIE_Interface
!  implicit none
!!  Integer, parameter :: Real8_ = selected_real_kind(14,200)
!  real (kind=Real8_), intent(in) :: TimeIn
!  UAd_NeedTime = TimeIn
!end subroutine UA_SetTime

!----------------------------------------------------------------------

subroutine UA_SetNorth
!  use ModKind
  use ModEIE_Interface
  implicit none
  UAl_IsNorth = .true.
end subroutine UA_SetNorth

!----------------------------------------------------------------------

subroutine UA_SetSouth
!  use ModKind
  use ModEIE_Interface
  implicit none
  UAl_IsNorth = .false.
end subroutine UA_SetSouth

!----------------------------------------------------------------------

subroutine UA_SetMLTs(MLTsIn, iError)

  use ModErrors
  use ModEIE_Interface

  implicit none

  integer, intent(out) :: iError
  real, dimension(UAi_NeednMLTs, UAi_NeednLats), intent(in) :: MLTsIn

  integer :: i, j

  iError = 0
  if (allocated(UAr2_NeedMLTs)) deallocate(UAr2_NeedMLTs)
  allocate(UAr2_NeedMLTs(UAi_NeednMLTs, UAi_NeednLats), stat=iError)
  if (iError /= 0) then
    iError = ecAllocationError_
    return
  else
    do i = 1, UAi_NeednMLTs
      do j = 1, UAi_NeednLats
        UAr2_NeedMLTs(i, j) = mod((MLTsIn(i, j) + 24.0), 24.0)
      end do
    end do
  end if

end subroutine UA_SetMLTs

!----------------------------------------------------------------------

subroutine UA_SetLats(LatsIn, iError)

  use ModErrors
  use ModEIE_Interface

  implicit none

  integer, intent(out) :: iError
  real, dimension(UAi_NeednMLTs, UAi_NeednLats), intent(in) :: LatsIn

  iError = 0
  if (allocated(UAr2_NeedLats)) deallocate(UAr2_NeedLats)
  allocate(UAr2_NeedLats(UAi_NeednMLTs, UAi_NeednLats), stat=iError)
  if (iError /= 0) then
    iError = ecAllocationError_
    return
  else
    UAr2_NeedLats(1:UAi_NeednMLTs, 1:UAi_NeednLats) = &
      LatsIn(1:UAi_NeednMLTs, 1:UAi_NeednLats)
  end if

end subroutine UA_SetLats

!----------------------------------------------------------------------

subroutine UA_SetGrid(MLTsIn, LatsIn, iError)

  use ModErrors
  use ModEIE_Interface

  implicit none

  integer, intent(out) :: iError
  real, dimension(UAi_NeednMLTs, UAi_NeednLats), intent(in) :: MLTsIn, LatsIn
  real, dimension(2) :: UAr1_Location
  real, dimension(5) :: EIEr1_Location

  integer :: i, j

  iError = 0

  call UA_SetMLTs(MLTsIn, iError)
  if (iError /= 0) return

  call UA_SetLats(LatsIn, iError)
  if (iError /= 0) return

  if (UAl_UseGridBasedEIE) then

    if (allocated(UAi3_InterpolationIndices)) &
      deallocate(UAi3_InterpolationIndices)
    allocate(UAi3_InterpolationIndices(UAi_NeednMLTs, UAi_NeednLats, 3), &
             stat=iError)
    if (iError /= 0) then
      iError = ecAllocationError_
      return
    end if

    if (allocated(UAr3_InterpolationRatios)) &
      deallocate(UAr3_InterpolationRatios)
    allocate(UAr3_InterpolationRatios(UAi_NeednMLTs, UAi_NeednLats, 2), &
             stat=iError)
    if (iError /= 0) then
      iError = ecAllocationError_
      return
    end if

    do i = 1, UAi_NeednLats
      do j = 1, UAi_NeednMLTs

        UAr1_Location(1) = mod((UAr2_NeedMLTs(j, i) + 24.0), 24.0)
        UAr1_Location(2) = UAr2_NeedLats(j, i)

        call EIE_FindPoint(UAr1_Location, EIEr1_Location, iError)

        if (iError == 0) then
          UAi3_InterpolationIndices(j, i, 1:3) = EIEr1_Location(1:3)
          UAr3_InterpolationRatios(j, i, 1:2) = EIEr1_Location(4:5)
        else
          UAi3_InterpolationIndices(j, i, 1) = -1
        end if

      end do
    end do

  end if

end subroutine UA_SetGrid

!----------------------------------------------------------------------
! General get value code - just do linear interpolation with the value
!----------------------------------------------------------------------

subroutine UA_GetValue(ValueIn, ValueOut, Filler, iError)

  use ModEIE_Interface

  implicit none

  integer, intent(out) :: iError
  real, intent(in)     :: Filler
  real, dimension(EIEi_HavenMlts, EIEi_HavenLats, EIEi_HavenBLKs), intent(in) :: &
    ValueIn
  real, dimension(UAi_NeednMLTs, UAi_NeednLats), intent(out) :: ValueOut

  integer :: iMLT, iLat, iB, iM, iL
  real    :: dM, dL

  iError = 0

  do iMLT = 1, UAi_NeednMLTs
    do iLat = 1, UAi_NeednLats

      iB = UAi3_InterpolationIndices(iMLT, iLat, 1)
      iM = UAi3_InterpolationIndices(iMLT, iLat, 2)
      iL = UAi3_InterpolationIndices(iMLT, iLat, 3)

      if (iB < 0 .or. iM < 0 .or. iL < 0) then
        ValueOut(iMLT, iLat) = Filler
      else

        dM = UAr3_InterpolationRatios(iMLT, iLat, 1)
        dL = UAr3_InterpolationRatios(iMLT, iLat, 2)

        ValueOut(iMLT, iLat) = &
          (1.0 - dM)*(1.0 - dL)*ValueIn(iM, IL, iB) + &
          (1.0 - dM)*(dL)*ValueIn(iM, IL + 1, iB) + &
          (dM)*(dL)*ValueIn(iM + 1, IL + 1, iB) + &
          (dM)*(1.0 - dL)*ValueIn(iM + 1, IL, iB)

      end if

    end do
  end do

end subroutine UA_GetValue

!----------------------------------------------------------------------
! Potential
!----------------------------------------------------------------------

subroutine UA_GetPotential(PotentialOut, iError)

  use ModEIE_Interface
  use ModErrors

  implicit none

  integer, intent(out) :: iError
  real, dimension(EIEi_HavenMlts, EIEi_HavenLats, EIEi_HavenBLKs) :: EIE_Value
  real, dimension(UAi_NeednMLTs, UAi_NeednLats)               :: ValueOut
  real, dimension(UAi_NeednMLTs, UAi_NeednLats), intent(out)  :: PotentialOut
  real :: Filler = 0.0
  integer :: iMLT, iLat

  iError = 0

  if (UAl_UseGridBasedEIE) then

    EIE_Value = EIEr3_HavePotential

    call UA_GetValue(EIE_Value, ValueOut, Filler, iError)

    if (iError /= 0) then
      write(*, *) "Error in routine UA_GetPotential:"
      write(*, *) cErrorCodes(iError)
      stop
    else
      PotentialOut = ValueOut
    end if

  else

    call UA_GetNonGridBasedPotential(ValueOut, iError)

    ! In Weimer, there are sometimes 0 potentials.  I don't know why.
    if (maxval(abs(UAr2_NeedLats)) > 80.0) then

      do iMlt = 2, UAi_NeednMLTs - 1
        do iLat = 2, UAi_NeednLats - 1
          if (ValueOut(iMlt, iLat) == 0.0) then
            if (ValueOut(iMlt, iLat - 1)*ValueOut(iMlt, iLat + 1) > 0) then
              ! this is a "hole" in the potential, since the potential is
              ! the same sign both above and below the current point.
              ValueOut(iMlt, iLat) = ( &
                                     ValueOut(iMlt - 1, iLat) + &
                                     ValueOut(iMlt, iLat - 1) + &
                                     ValueOut(iMlt + 1, iLat) + &
                                     ValueOut(iMlt, iLat + 1))/4.0
            end if
          end if
        end do
      end do

    end if

    if (iError == 0) then
      PotentialOut = ValueOut
    else
      write(*, *) 'error in UA_GetPotential (nongrid): ', cErrorCodes(iError)
      stop
      PotentialOut = -1.0e32
    end if

  end if

end subroutine UA_GetPotential

!----------------------------------------------------------------------
! Potential - non-grid-based
!----------------------------------------------------------------------

subroutine UA_GetNonGridBasedPotential(PotentialOut, iError)

  use EIE_ModWeimer, only: WeiEpot96, get_tilt
  use w05sc, only: setmodel, epotval

  use ModEIE_Interface
  use ModErrors
  use ModTimeConvert, ONLY: time_real_to_int

  implicit none

  integer, intent(out)                                      :: iError
  real, dimension(UAi_NeednMLTs, UAi_NeednLats), intent(out) :: PotentialOut

  integer, dimension(7) :: itime
  integer :: iYear, iMonth, iDay, iHour, iMinute, iSecond

  integer :: iChange, iHemisphere

  integer :: iMLT, iLat

  real    :: Potential, ETheta, EPhi, MLT, Lat, tilt, hr

  iError = 0

  PotentialOut = 0.0

  if (index(EIE_NameOfEFieldModel, 'zero') > 0) then
    return
  end if

  if (IOd_NeedTime < 0.0) then
    iError = ecTimeNotSet_
    return
  end if

  if (index(EIE_NameOfEFieldModel, 'hpi') <= 0) then
    if (IOr_NeedIMFBz < -1000.0) then
      iError = ecIMFBzNotSet_
      return
    end if
    if (IOr_NeedIMFBy < -1000.0) then
      iError = ecIMFByNotSet_
      return
    end if
  end if

  if (index(EIE_NameOfEFieldModel, 'weimer05') > 0) then
    if (IOr_NeedSWV < -1000.0) then
      iError = ecSWVNotSet_
      return
    end if
    if (IOr_NeedSWN < -1000.0) then
      iError = ecSWNNotSet_
      return
    end if
    iChange = 1
  end if

  if (index(EIE_NameOfEFieldModel, 'weimer01') > 0) then
    if (IOr_NeedSWV < -1000.0) then
      iError = ecSWVNotSet_
      return
    end if
    iChange = 1
  end if

  if (index(EIE_NameOfEFieldModel, 'weimer96') > 0) then
    if (IOr_NeedSWV < -1000.0) then
      iError = ecSWVNotSet_
      return
    end if
    iChange = 1
  end if

  if (index(EIE_NameOfEFieldModel, 'millstone') > 0) then
    if (IOr_NeedHPI < -1000.0) then
      iError = ecHPINotSet_
      return
    end if
  end if

  if (index(EIE_NameOfEFieldModel, 'hmr89') > 0) then
    if (IOr_NeedKP < -1000.0) then
      iError = ecKPNotSet_
      return
    end if
    iChange = 1
  end if

  call time_real_to_int(IOd_NeedTime, itime)

  iYear = itime(1)
  iMonth = itime(2)
  iDay = itime(3)
  iHour = itime(4)
  iMinute = itime(5)
  iSecond = itime(6)

  Lat = 0.0

  do iMLT = 1, UAi_NeednMLTs
    do iLat = 1, UAi_NeednLats

      if (iLat == 1) then
        iChange = 1
      else
        if (UAr2_NeedLats(iMLT, iLat)*UAr2_NeedLats(iMLT, iLat - 1) <= 0.0) &
          iChange = 1
      end if

      MLT = UAr2_NeedMLTs(iMLT, iLat)
      Lat = UAr2_NeedLats(iMLT, iLat)

      if (UAl_IsNorth) then
        if (Lat > 0.0) then
          iHemisphere = 1
        else
          Lat = abs(Lat)
          iHemisphere = -1
        end if
      else
        iHemisphere = -1
        if (Lat < 0.0) Lat = abs(Lat)
      end if

      if (index(EIE_NameOfEFieldModel, 'weimer05') > 0) then
        if (abs(lat) >= 45.0) then
          hr = float(iHour) + float(iMinute)/60.0
          if (IsFixedTilt) then
            tilt = 0.0
          else
            tilt = iHemisphere*get_tilt(iYear, iMonth, iDay, hr)
          end if
          call setmodel(IOr_NeedIMFBy, IOr_NeedIMFBz, tilt, &
                        IOr_NeedSWV, IOr_NeedSWN, 'epot')
          call epotval(lat, mlt, 0.0, Potential)
          Potential = Potential*1000.0
          iChange = 0
        else
          Potential = 0.0
          ETheta = 0.0
          EPhi = 0.0
        end if
      end if

      if (index(EIE_NameOfEFieldModel, 'weimer96') > 0) then
        if (abs(lat) >= 45.0) then
          !call timing_start("weimer")
          call WEIEPOT96(iYear, iMonth, iDay, iHour, iMinute, IOr_NeedSWV, &
                         IOr_NeedIMFBy, IOr_NeedIMFBz, iHemisphere, iChange, &
                         Lat, MLT, ETheta, EPhi, Potential)
          !call timing_stop("weimer")
          Potential = Potential*1000.0
          iChange = 0
        else
          Potential = 0.0
          ETheta = 0.0
          EPhi = 0.0
        end if
      end if

      if (index(EIE_NameOfEFieldModel, 'millstone_hpi') > 0) then
        call MHEMODL(Lat, mlt, IOr_NeedHPI, IOr_NeedIMFBy, &
                     IOr_NeedIMFBz, 1, ETheta, EPhi, Potential)
        Potential = Potential*1000.0
      end if

      if (index(EIE_NameOfEFieldModel, 'millstone_imf') > 0) then
        call MHEMODL(lat, mlt, IOr_NeedHPI, IOr_NeedIMFBy, &
                     IOr_NeedIMFBz, 2, ETheta, EPhi, Potential)
        Potential = Potential*1000.0
      end if

      if (index(EIE_NameOfEFieldModel, 'hmr89') > 0) then
        call HMREPOT(lat, mlt, IOr_NeedIMFBy, IOr_NeedIMFBz, &
                     IOr_NeedKp, iChange, ETheta, EPhi, Potential)
        Potential = Potential*1000.0
        iChange = 0
      end if

      if (index(EIE_NameOfEFieldModel, 'izmem') > 0) then
        call IZEPOT(iMonth, lat, mlt, iHemisphere*IOr_NeedIMFBy, &
                    IOr_NeedIMFBz, ETheta, EPhi, Potential)
        Potential = Potential*1000.0
      end if

      PotentialOut(iMLT, iLat) = Potential

    end do
  end do

end subroutine UA_GetNonGridBasedPotential

!----------------------------------------------------------------------
! Electron Average Energy
!----------------------------------------------------------------------

subroutine UA_GetAveE(AveEOut, iError)

  use ModEIE_Interface
  use ModErrors

  implicit none

  integer, intent(out) :: iError
  real, dimension(EIEi_HavenMlts, EIEi_HavenLats, EIEi_HavenBLKs) :: EIE_Value
  real, dimension(UAi_NeednMLTs, UAi_NeednLats)               :: ValueOut
  real, dimension(UAi_NeednMLTs, UAi_NeednLats), intent(out)  :: AveEOut

  iError = 0

  if (UAl_UseGridBasedEIE) then

    EIE_Value = EIEr3_HaveAveE
    call UA_GetValue(EIE_Value, ValueOut, 0.1, iError)

    if (iError /= 0) then
      write(*, *) "Error in routine UA_GetAveE:"
      write(*, *) cErrorCodes(iError)
    else
      AveEOut = ValueOut
    end if

  else

    call UA_GetNonGridBasedAveE(ValueOut, iError)

    if (iError == 0) then
      AveEOut = ValueOut
    else
      AveEOut = -1.0e32
    end if

  end if

end subroutine UA_GetAveE

!----------------------------------------------------------------------
! Electron Average Energy - non-grid-based
!----------------------------------------------------------------------

subroutine UA_GetNonGridBasedAveE(AveEOut, iError)

  use ModEIE_Interface
  use ModErrors

  implicit none

  integer, intent(out) :: iError
  real, dimension(UAi_NeednMLTs, UAi_NeednLats) :: AveEOut

  integer :: iMLT, iLat
  real :: Lat, MLT, hpi
  real :: ped, hal, avkev, eflx

  iError = 0
  if (IOr_NeedHPINorm < -1000.0) then
    iError = ecHPINotSet_
    return
  end if

  do iMLT = 1, UAi_NeednMLTs
    do iLat = 1, UAi_NeednLats

      MLT = UAr2_NeedMLTs(iMLT, iLat)
      Lat = UAr2_NeedLats(iMLT, iLat)
      hpi = IOr_NeedHPINorm

      call get_auroral_conductance(Lat, MLT, hpi, &
                                   ped, hal, avkev, eflx)

      AveEOut(iMLT, iLat) = avkev

    end do
  end do

end subroutine UA_GetNonGridBasedAveE

!----------------------------------------------------------------------
! Electron Energy Flux
!----------------------------------------------------------------------

subroutine UA_GetEFlux(EFluxOut, iError)

  use ModEIE_Interface
  use ModErrors

  implicit none

  integer, intent(out) :: iError
  real, dimension(EIEi_HavenMlts, EIEi_HavenLats, EIEi_HavenBLKs) :: EIE_Value
  real, dimension(UAi_NeednMLTs, UAi_NeednLats)               :: ValueOut
  real, dimension(UAi_NeednMLTs, UAi_NeednLats), intent(out)  :: EFluxOut

  iError = 0

  if (UAl_UseGridBasedEIE) then

    EIE_Value = EIEr3_HaveEFlux
    call UA_GetValue(EIE_Value, ValueOut, 1.0e-10, iError)

    if (iError /= 0) then
      write(*, *) "Error in routine UA_GetEFlux:"
      write(*, *) cErrorCodes(iError)
      stop
    else
      EFluxOut = ValueOut
    end if

  else

    call UA_GetNonGridBasedEFlux(ValueOut, iError)

    if (iError == 0) then
      EFluxOut = ValueOut
    else
      EFluxOut = -1.0e32
    end if

  end if

end subroutine UA_GetEFlux

!----------------------------------------------------------------------
! Electron Energy Flux - non-grid-based
!----------------------------------------------------------------------

subroutine UA_GetNonGridBasedEFlux(EFluxOut, iError)

  use ModEIE_Interface
  use ModErrors

  implicit none

  integer, intent(out) :: iError
  real, dimension(UAi_NeednMLTs, UAi_NeednLats) :: EFluxOut

  integer :: iMLT, iLat
  real :: Lat, MLT, hpi
  real :: ped, hal, avkev, eflx

  iError = 0

  if (IOr_NeedHPINorm < -1000.0) then
    iError = ecHPINotSet_
    return
  end if

  do iMLT = 1, UAi_NeednMLTs
    do iLat = 1, UAi_NeednLats

      MLT = UAr2_NeedMLTs(iMLT, iLat)
      Lat = UAr2_NeedLats(iMLT, iLat)

      call get_auroral_conductance(Lat, MLT, IOr_NeedHPINorm, &
                                   ped, hal, avkev, eflx)

      EFluxOut(iMLT, iLat) = eflx

    end do
  end do

end subroutine UA_GetNonGridBasedEFlux

!----------------------------------------------------------------------
! Ions Average Energy
!----------------------------------------------------------------------

subroutine UA_GetIonAveE(IonAveEOut, iError)

  use ModEIE_Interface
  use ModErrors

  implicit none

  integer, intent(out) :: iError
  real, dimension(EIEi_HavenMlts, EIEi_HavenLats, EIEi_HavenBLKs) :: EIE_Value
  real, dimension(UAi_NeednMLTs, UAi_NeednLats)               :: ValueOut
  real, dimension(UAi_NeednMLTs, UAi_NeednLats), intent(out)  :: IonAveEOut

  iError = 0

  if (UAl_UseGridBasedEIE) then

    EIE_Value = EIEr3_HaveIonAveE
    call UA_GetValue(EIE_Value, ValueOut, 0.1, iError)

    if (iError /= 0) then
      write(*, *) "Error in routine UA_GetAveE:"
      write(*, *) cErrorCodes(iError)
    else
      IonAveEOut = ValueOut
    end if

  else

    call UA_GetNonGridBasedAveE(ValueOut, iError)

    if (iError == 0) then
      IonAveEOut = ValueOut
    else
      IonAveEOut = 0.0
    end if

  end if

end subroutine UA_GetIonAveE

!----------------------------------------------------------------------
! Ion Average Energy - non-grid-based
!----------------------------------------------------------------------

subroutine UA_GetNonGridBasedIonAveE(IonAveEOut, iError)

  use ModEIE_Interface
  use ModErrors

  implicit none

  integer, intent(out) :: iError
  real, dimension(UAi_NeednMLTs, UAi_NeednLats) :: IonAveEOut

  write(*, *) "UA_GetNonGridBasedIonAveE doen't work yet!"
  IonAveEOut = 0.0

end subroutine UA_GetNonGridBasedIonAveE

!----------------------------------------------------------------------
! Ion Energy Flux
!----------------------------------------------------------------------

subroutine UA_GetIonEFlux(IonEFluxOut, iError)

  use ModEIE_Interface
  use ModErrors

  implicit none

  integer, intent(out) :: iError
  real, dimension(EIEi_HavenMlts, EIEi_HavenLats, EIEi_HavenBLKs) :: EIE_Value
  real, dimension(UAi_NeednMLTs, UAi_NeednLats)               :: ValueOut
  real, dimension(UAi_NeednMLTs, UAi_NeednLats), intent(out)  :: IonEFluxOut

  iError = 0

  if (UAl_UseGridBasedEIE) then

    EIE_Value = EIEr3_HaveIonEFlux
    call UA_GetValue(EIE_Value, ValueOut, 1.0e-10, iError)

    if (iError /= 0) then
      write(*, *) "Error in routine UA_GetEFlux:"
      write(*, *) cErrorCodes(iError)
      stop
    else
      IonEFluxOut = ValueOut
    end if

  else

    call UA_GetNonGridBasedEFlux(ValueOut, iError)

    if (iError == 0) then
      IonEFluxOut = ValueOut
    else
      IonEFluxOut = -1.0e32
    end if

  end if

end subroutine UA_GetIonEFlux

!----------------------------------------------------------------------
! Ion Energy Flux - non-grid-based
!----------------------------------------------------------------------

subroutine UA_GetNonGridBasedIonEFlux(IonEFluxOut, iError)

  use ModEIE_Interface
  use ModErrors

  implicit none

  integer, intent(out) :: iError
  real, dimension(UAi_NeednMLTs, UAi_NeednLats) :: IonEFluxOut

  write(*, *) "UA_GetNonGridBasedIonEFlux doesn't work yet!"
  IonEFluxOut = 0.0

end subroutine UA_GetNonGridBasedIonEFlux

