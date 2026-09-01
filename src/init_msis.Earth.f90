! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

!-----------------------------------------------------------------------------
! $Id: init_msis.Earth.f90,v 1.40 2017/08/09 15:18:05 ridley Exp $
! Author: Aaron Ridley, UMichigan
!
! Modified:
!           AGB Oct 2013 - Adapted to use specified realistic F10.7 value when
!                          driving F10.7 in RCMR
!           Asad Feb 2013 - Adapted to use F10.7 = 150 sfu when using RCMR
!
! Comments: Routines to initialize the GITM thermosphere with MSIS
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! This is the driver code for MSIS, since there are now two version
! that can be called. Since NRL has hard coded the reals to be real(4)
! for MSIS2.1, we have to translate variables.  This does not have to
! be done for MSIS00.
! ------------------------------------------------------------------------------

subroutine call_msis(lonDeg, latDeg, altKm, f107, f107a, densities10, temp)

  use ModTime
  use EUA_ModMsis00, only: meters, gtd7
  use ModMsis21, only: gtd8d
  use ModInputs, only: useMsis21
  use ModConstants, only: Boltzmanns_Constant, AMU
  use ModIndicesInterfaces, only: get_HPI

  implicit none

  real, intent(in) :: lonDeg, latDeg, altKm, f107, f107a
  real, intent(out) :: densities10(10)
  real, intent(out) :: temp

  ! MSIS-2.1 hard-codes the size of reals, and needs the following:
  integer :: iyd
  real(4) :: sec
  real(4) :: alt
  real(4) :: glat
  real(4) :: glong
  real(4) :: stl
  real(4) :: f107a_4
  real(4) :: f107_4
  real(4) :: ap_4(7)
  integer :: mass
  real(4) :: d(10), t(2)

  real, dimension(1:2) :: msis_temp
  real, dimension(1:9) :: msis_dens9
  real, dimension(7) :: AP

  real :: Lst
  real :: ffactor, no, h, hp
  integer :: iError

  LST = mod(utime/3600.0 + LonDeg/15.0, 24.0)
  AP = 10

  ! We don't often have Ap, but have hemispheric power. So, use that:
  call get_HPI(CurrentTime, HP, iError)
  if (iError > 0) hp = 40.0
  Ap = min(200., max(-40.72 + 1.3*HP, 10.))

  if (useMsis21) then
    iyd = iJulianDay
    sec = utime
    alt = altKm
    glat = latDeg
    glong = lonDeg
    stl = LST
    f107a_4 = f107a
    f107_4 = f107
    ap_4 = AP
    ! mass is not used, but passed anyways
    mass = -1
    call gtd8d( &
      iyd, sec, alt, glat, glong, stl, f107a_4, f107_4, ap_4, mass, &
      d, t)
    temp = t(2)
    ! Convert to /m3
    ! 10th density is NO now!
    densities10 = d*1e6
  else
    CALL GTD7(iJulianDay, utime, AltKm, LatDeg, LonDeg, LST, &
              f107a, f107, AP, 48, msis_dens9, msis_temp)
    temp = msis_temp(2)
    densities10(1:9) = msis_dens9
    ffactor = 6.36*log(f107) - 13.8
    no = (ffactor*1.0e13 + 8.0e13)*1.24 ! 12.4 ! 12.4 is roughly exp
    ! This is obviously an approximation:
    h = Boltzmanns_Constant*msis_temp(2)/ &
        (9.5*28.0*AMU)/1000.0
    densities10(10) = no*exp(-(altKm - 100.0)/h)
  endif

end subroutine call_msis

subroutine get_msis_temperature(lon, lat, alt, t, h)

  use ModTime
  use ModInputs
  use ModPlanet
  use ModGITM
  use ModRCMR, only: RCMRFlag, RCMROutType
  use ModIndicesInterfaces

  use EUA_ModMsis00, only: meters, gtd7
  use ModMsis21, only: gtd8d

  implicit none

  real, intent(in) :: lon, lat, alt
  real, intent(out) :: t, h

  real, dimension(1:2) :: msis_temp
  real, dimension(1:9) :: msis_dens
  real :: msis_dens10(10)
  real :: msis_temp1

  real :: LonDeg, LatDeg, AltKm, LST
  real, dimension(7)    :: AP
  real :: nO, nO2, nN2, m, r, g

  integer :: iError

  !-------------------------------------------------------

  ap = 10.0

  call meters(.true.)

  LonDeg = lon*180.0/pi
  LatDeg = lat*180.0/pi
  AltKm = alt/1000.0
  LST = mod(utime/3600.0 + LonDeg/15.0, 24.0)
  iError = 0

  call get_f107(CurrentTime, f107, iError)

  if (iError /= 0) then
    write(*, *) "Error in getting F107 value.  Is this set?"
    write(*, *) "Code : ", iError
    call stop_gitm("Stopping in euv_ionization_heat")
  endif

  call get_f107a(CurrentTime, f107a, iError)
  if (iError /= 0) then
    write(*, *) "Error in getting F107a value.  Is this set?"
    write(*, *) "Code : ", iError
    call stop_gitm("Stopping in euv_ionization_heat")
  endif

  if (RCMRFlag .and. RCMROutType == "F107") then
    call call_msis(lonDeg, latDeg, altKm, f107_msis, f107a_msis, &
                   msis_dens10, msis_temp1)
    msis_dens = msis_dens10(1:9)
    msis_temp = msis_temp1
  else
    call call_msis(lonDeg, latDeg, altKm, f107, f107a, msis_dens10, msis_temp1)
    msis_dens = msis_dens10(1:9)
    msis_temp = msis_temp1
  endif

  t = msis_temp(2)
  nO = msis_dens(2)
  nN2 = msis_dens(3)
  nO2 = msis_dens(4)

  m = (nO*mass(iO_3P_) + nO2*mass(iO2_) + nN2*mass(iN2_))/ &
      (nO + nO2 + nN2)

  r = RBody + alt
  g = Gravitational_Constant*(RBody/r)**2
  h = Boltzmanns_Constant*t/(m*g)

end subroutine get_msis_temperature

!--------------------------------------------------------------
! Initialize MSIS itself:
!--------------------------------------------------------------

subroutine initialize_msis_routines

  use ModInputs, only: UseMsis, &
                       UseMSISDiurnal, UseMSISSemidiurnal, UseMSISTerdiurnal, useMsis21, sw_msis
  use EUA_ModMsis00, ONLY: meters, tselec
  use msis_init, only: msisinit

  implicit none

  real*4 :: sw_msis4x25(25)
  logical, save :: isFirstTIme = .true.

  ! We only need to initialize MSIS once!
  if (.not. isFirstTime) then
    return
  else
    isFirstTime = .false.
  endif

  ! We want units of /m3 and not /cm3

  call meters(.true.)

  if (UseMsis) then
    sw_msis = 1
    if (.not. UseMSISDiurnal) then
      CALL report("...Using MSIS without diurnal tidal variations...", 0)
      sw_msis(7) = 0
    endif
    if (.not. UseMSISSemidiurnal) then
      CALL report("...Using MSIS without semi-diurnal tidal variations...", 0)
      sw_msis(8) = 0
    endif
    if (.not. UseMSISTerdiurnal) then
      CALL report("...Using MSIS without terdiurnal tidal variations...", 0)
      sw_msis(14) = 0
    endif
  else
    sw_msis = 0
    sw_msis(1) = 1
    sw_msis(9) = 1
  endif

  sw_msis(9) = 0
  sw_msis(2) = 0

  if (useMsis21) then
    sw_msis(2) = 1
    sw_msis4x25 = sw_msis
    call msisinit(parmpath='UA/DataIn/LowerBCs/', switch_legacy=sw_msis4x25)
  else
    call tselec(sw_msis)
  endif

end subroutine initialize_msis_routines

!--------------------------------------------------------------
!
!--------------------------------------------------------------

subroutine init_msis

  use ModGITM
  use ModInputs
  use ModConstants
  use ModPlanet
  use ModTime
  use ModIndicesInterfaces

  use EUA_ModMsis00, ONLY: meters, gtd7, tselec
  use msis_init, only: msisinit

  implicit none

  ! msis variables

  real, dimension(1:2) :: msis_temp
  real, dimension(1:9) :: msis_dens
  real :: msis_dens10(10)
  real :: msis_temp1

  integer :: iBlock, iAlt, iLat, iLon, iSpecies, iyd, iError
  real :: geo_lat, geo_lon, geo_alt, geo_lst, m, k, ut
  real :: ffactor, h, no
  real, dimension(7)  :: ap = 10.0

  real*4 :: hwm_utime, hwm_alt, hwm_lat, hwm_lon, hwm_lst
  real*4 :: hwm_f107a, hwm_f107, hwm_ap(2), qw(2)

  character(250) :: path = './DataIn/LowerBCs/'

  call report("init_msis", 0)

  iError = 0
  call get_f107(CurrentTime, f107, iError)
  if (iError /= 0) then
    write(*, *) "Error in getting F107 value (init_msis).  Is this set?"
    write(*, *) "Code : ", iError
    call stop_gitm("Stopping in advance")
  endif

  call get_f107a(CurrentTime, f107a, iError)
  if (iError /= 0) then
    write(*, *) "Error in getting F107a value (init_msis).  Is this set?"
    write(*, *) "Code : ", iError
    call stop_gitm("Stopping in advance")
  endif

  call initialize_msis_routines

  !--------------------------------------------------------------------------
  !
  !  From the msis90 library:
  !
  !     OUTPUT:
  !        D(1) - HE NUMBER DENSITY(CM-3)
  !        D(2) - O NUMBER DENSITY(CM-3)
  !        D(3) - N2 NUMBER DENSITY(CM-3)
  !        D(4) - O2 NUMBER DENSITY(CM-3)
  !        D(5) - AR NUMBER DENSITY(CM-3)
  !        D(6) - TOTAL MASS DENSITY(GM/CM3)
  !        D(7) - H NUMBER DENSITY(CM-3)
  !        D(8) - N NUMBER DENSITY(CM-3)
  !        T(1) - EXOSPHERIC TEMPERATURE
  !        T(2) - TEMPERATURE AT ALT
  !
  !      TO GET OUTPUT IN M-3 and KG/M3:   CALL METERS(.TRUE.)
  !
  !      O, H, and N set to zero below 72.5 km
  !      Exospheric temperature set to average for altitudes below 120 km.
  !
  !--------------------------------------------------------------------------

  if (DoRestart) return

  !           The following is for test and special purposes:
  !            TO TURN ON AND OFF PARTICULAR VARIATIONS CALL TSELEC(SW)
  !               WHERE SW IS A 25 ELEMENT ARRAY CONTAINING 0. FOR OFF, 1.
  !               FOR ON, OR 2. FOR MAIN EFFECTS OFF BUT CROSS TERMS ON
  !               FOR THE FOLLOWING VARIATIONS
  !               1 - F10.7 EFFECT ON MEAN  2 - TIME INDEPENDENT
  !               3 - SYMMETRICAL ANNUAL    4 - SYMMETRICAL SEMIANNUAL
  !               5 - ASYMMETRICAL ANNUAL   6 - ASYMMETRICAL SEMIANNUAL
  !               7 - DIURNAL               8 - SEMIDIURNAL
  !               9 - DAILY AP             10 - ALL UT/LONG EFFECTS
  !              11 - LONGITUDINAL         12 - UT AND MIXED UT/LONG
  !              13 - MIXED AP/UT/LONG     14 - TERDIURNAL
  !              15 - DEPARTURES FROM DIFFUSIVE EQUILIBRIUM
  !              16 - ALL TINF VAR         17 - ALL TLB VAR
  !              18 - ALL TN1 VAR           19 - ALL S VAR
  !              20 - ALL TN2 VAR           21 - ALL NLB VAR
  !              22 - ALL TN3 VAR           23 - TURBO SCALE HEIGHT VAR

  ! Initialize data

  iyd = iTimeArray(1)*1000 + iJulianDay

  do iBlock = 1, nBlocks
    do iAlt = -1, nAlts + 2
      do iLon = -1, nLons + 2
        do iLat = -1, nLats + 2

          geo_lon = mod(Longitude(iLon, iBlock)*180.0/pi + 360.0, 360.0)

          geo_lat = Latitude(iLat, iBlock)*180.0/pi
          if (geo_lat < -90.0) then
            geo_lat = -180.0 - geo_lat
            geo_lon = mod(geo_lon + 180.0, 360.0)
          endif
          if (geo_lat > 90.0) then
            geo_lat = 180.0 - geo_lat
            geo_lon = mod(geo_lon + 180.0, 360.0)
          endif

          geo_alt = Altitude_GB(iLon, iLat, iAlt, iBlock)/1000.0
          geo_lst = mod(utime/3600.0 + geo_lon/15.0, 24.0)

          !
          ! Call MSIS (results will be im mks units)
          !

          !CALL GTD7(iJulianDay, utime, geo_alt, geo_lat, geo_lon, geo_lst, &
          !          F107A, F107, AP, 48, msis_dens, msis_temp)
          call call_msis(geo_lon, geo_lat, geo_alt, f107, f107a, &
                         msis_dens10, msis_temp1)
          msis_dens = msis_dens10(1:9)
          msis_temp = msis_temp1

          ! Initialize densities to zero in case msis does not set it
          NDensityS(iLon, iLat, iAlt, :, iBlock) = 1.0

          NDensityS(iLon, iLat, iAlt, iHe_, iBlock) = &
            max(msis_dens10(1), 100.0)
          NDensityS(iLon, iLat, iAlt, iO_3P_, iBlock) = &
            max(msis_dens10(2), 100.0)
          NDensityS(iLon, iLat, iAlt, iN2_, iBlock) = &
            max(msis_dens10(3), 100.0)
          NDensityS(iLon, iLat, iAlt, iO2_, iBlock) = &
            max(msis_dens10(4), 100.0)
          NDensityS(iLon, iLat, iAlt, iH_, iBlock) = &
            max(msis_dens10(7), 100.0)
          NDensityS(iLon, iLat, iAlt, iN_4S_, iBlock) = &
            max(msis_dens10(8), 100.0)
          NDensityS(iLon, iLat, iAlt, iNO_, iBlock) = &
            max(msis_dens10(10), 100.0)

          NDensityS(iLon, iLat, iAlt, iN_2P_, iBlock) = &
            NDensityS(iLon, iLat, iAlt, iN_4S_, iBlock)/10000.0
          NDensityS(iLon, iLat, iAlt, iN_2D_, iBlock) = &
            NDensityS(iLon, iLat, iAlt, iN_4S_, iBlock)/100.0
          NDensityS(iLon, iLat, iAlt, iO_1D_, iBlock) = &
            NDensityS(iLon, iLat, iAlt, iO_3P_, iBlock)*0.0 + 1

          MeanMajorMass(iLon, iLat, iAlt) = 0

          do iSpecies = 1, nSpecies
            MeanMajorMass(iLon, iLat, iAlt) = &
              MeanMajorMass(iLon, iLat, iAlt) + &
              Mass(iSpecies)* &
              NDensityS(iLon, iLat, iAlt, iSpecies, iBlock)/ &
              sum(NDensityS(iLon, iLat, iAlt, 1:nSpecies, iBlock))
          enddo

          TempUnit(iLon, iLat, iAlt) = &
            MeanMajorMass(iLon, iLat, iAlt)/Boltzmanns_Constant

          Temperature(iLon, iLat, iAlt, iBlock) = &
            msis_temp(2)/TempUnit(iLon, iLat, iAlt)

          Rho(iLon, iLat, iAlt, iBlock) = msis_dens(6)

          NDensity(iLon, iLat, iAlt, iBlock) = &
            sum(NDensityS(iLon, iLat, iAlt, 1:nSpecies, iBlock))

          LogNS(iLon, iLat, iAlt, 1:nSpecies, iBlock) = &
            log(NDensityS(iLon, iLat, iAlt, 1:nSpecies, iBlock))

          hwm_utime = utime
          hwm_alt = geo_alt
          hwm_lat = geo_lat
          hwm_lon = geo_lon
          hwm_lst = geo_lst
          hwm_f107a = f107a
          hwm_f107 = f107
          hwm_ap(1) = -1.0
          hwm_ap(2) = 4.0

          if (UseMSISDiurnal .and. &
              UseMSISSemidiurnal .and. &
              UseMSISTerdiurnal) then

            call hwm14(iyd, hwm_utime, hwm_alt, hwm_lat, hwm_lon, hwm_lst, &
                       hwm_f107a, hwm_f107, hwm_ap, path, qw)

            ! qw is north&east
            Velocity(iLon, iLat, iAlt, iEast_, iBlock) = qw(2)
            Velocity(iLon, iLat, iAlt, iNorth_, iBlock) = qw(1)

          else

            Velocity(iLon, iLat, iAlt, iEast_, iBlock) = 0.0
            Velocity(iLon, iLat, iAlt, iNorth_, iBlock) = 0.0

          endif

        enddo
      enddo
    enddo

    Rho(:, :, :, iBlock) = 0.0
    NDensity(:, :, :, iBlock) = 0.0

    do iSpecies = 1, nSpecies

      NDensity(:, :, :, iBlock) = NDensity(:, :, :, iBlock) + &
                                  NDensityS(:, :, :, iSpecies, iBlock)

      Rho(:, :, :, iBlock) = Rho(:, :, :, iBlock) + &
                             Mass(iSpecies)*NDensityS(:, :, :, iSpecies, iBlock)

    enddo

    call calc_co2(iBlock)

  enddo

end subroutine init_msis

!--------------------------------------------------------------
!
!--------------------------------------------------------------

subroutine msis_bcs(iJulianDay, UTime, Alt, LatIn, LonIn, Lst, &
                    F107A, F107, AP, LogNS, Temp, LogRho, v)

  use ModTime, only: iTimeArray
  use ModPlanet
  use ModInputs, only: &
    UseHmeTides, UseFileTides, &
    UseOBCExperiment, sw_msis, UseMSIS21, co2ppm
  use EUA_ModMsis00, ONLY: gtd7, tselec

  implicit none

  real, dimension(25) :: sw_tmp = 1.0

  integer, intent(in) :: iJulianDay
  real, intent(in) :: uTime, Alt, LatIn, LonIn, LST, f107a, f107
  real, intent(in):: ap
  real, intent(out) :: LogNS(nSpecies), Temp, LogRho, v(2)

  real :: lat, lon

  real :: msis_temp(2)
  real :: msis_dens(9), oMSIS, oCurrentSeason, oOffsetSeason
  real :: msis_dens10(10)
  real :: msis_temp1
  real :: AP_I(7), ffactor, no
  integer :: iyd, iJulianDayOffset

  real*4 :: hwm_utime, hwm_alt, hwm_lat, hwm_lon, hwm_lst
  real*4 :: hwm_f107a, hwm_f107, hwm_ap(2), qw(2)

  real :: base, season, vari
  
  character(250) :: path = './DataIn/LowerBCs/'

  lat = LatIn
  lon = mod(LonIn + 360.0, 360.0)

  if (lat > 90) then
    lat = 180 - lat
    lon = mod(lon + 180.0, 360.0)
  endif
  if (lat < -90) then
    lat = -180 - lat
    lon = mod(lon + 180.0, 360.0)
  endif

  !----------------------------------------------------------------------------
  AP_I = AP

  call call_msis(lon, lat, alt, f107, f107a, msis_dens10, msis_temp1)
  msis_dens = msis_dens10(1:9)
  msis_temp = msis_temp1

  LogNS(iO_3P_) = alog(max(msis_dens10(2), 1.0))
  LogNS(iO2_) = alog(max(msis_dens10(4), 1.0))
  LogNS(iN2_) = alog(max(msis_dens10(3), 1.0))
  LogNS(iN_4S_) = alog(max(msis_dens10(8), 1.0))
  LogNS(iHe_) = alog(max(msis_dens10(1), 1.0))
  LogNS(iNO_) = alog(max(msis_dens10(10), 1.0))
  logNS(iCO2_) = alog(CO2ppm*1e-6/(1.0 - CO2ppm*1e-6)* &
                      (msis_dens(1) + msis_dens(2) + msis_dens(3)))

  Temp = msis_temp(2)
  LogRho = alog(msis_dens(6))

  iyd = iTimeArray(1)*1000 + iJulianDay
  hwm_utime = utime
  hwm_alt = alt
  hwm_lat = lat
  hwm_lon = lon
  hwm_lst = lst
  hwm_f107a = f107a
  hwm_f107 = f107
  hwm_ap(1) = -1.0
  hwm_ap(2) = -1.0

  if (UseHmeTides .or. UseFileTides) then
    ! set_vertical_bcs replaces these with TidesEast/TidesNorth.
    V(1) = 0.0
    V(2) = 0.0
  else
    call hwm14(iyd, hwm_utime, hwm_alt, hwm_lat, hwm_lon, hwm_lst, &
               hwm_f107a, hwm_f107, hwm_ap, path, qw)
    ! qw is north&east
    V(1) = qw(2)
    V(2) = qw(1)
  endif

  ! Do some O experimentation (ONLY for MSIS00!):

  if (UseOBCExperiment .and. .not. UseMSIS21) then

    oMSIS = msis_dens10(2)

    sw_tmp = sw_msis

    ! Experimentation:
    ! This will make MSIS only have latitudinal variation:
    sw_tmp(7) = 0
    sw_tmp(8) = 0
    sw_tmp(11) = 0
    sw_tmp(14) = 0

    call tselec(sw_tmp)

    CALL GTD7(iJulianDay, uTime, Alt, Lat, Lon, LST, &
              F107A, F107, AP_I, 48, msis_dens, msis_temp)

    oCurrentSeason = msis_dens(2)

    ! Offset by 6 months:
    iJulianDayOffset = mod(iJulianDay + 182, 365)

    CALL GTD7(iJulianDayOffset, uTime, Alt, Lat, Lon, LST, &
              F107A, F107, AP_I, 48, msis_dens, msis_temp)

    oOffsetSeason = msis_dens(2)

    ! Put MSIS back to normal
    call tselec(sw_msis)

    LogNS(iO_3P_) = alog(max(oMSIS - oCurrentSeason + oOffsetSeason, 1.0))

  endif

  if (UseOBCExperiment .and. UseMSIS21) then
     oMSIS = exp(LogNS(iO_3P_))
     ! base should be around 1:
     base = 1.2
     ! Add more O to the summer hemisphere:
     season = sin((iJulianDay - 90) * 2 * PI / 365)
     vari = 0.2 * sin(Lat * Pi / 180) * season
     LogNS(iO_3P_) = alog( oMSIS * (base + vari))
  endif
  
end subroutine msis_bcs

subroutine calc_co2(iBlock)

  use ModPlanet
  use ModGITM
  use ModInputs, only: CO2ppm

  implicit none

  integer, intent(in) :: iBlock
  integer :: iAlt
  real, dimension(-1:nLons + 2, -1:nLats + 2, -1:nAlts + 2) :: &
    Have, Hn2, Ho, r, Hco2, Hco2t

  Have = -Boltzmanns_Constant* &
         Temperature(:, :, :, iBlock)*TempUnit/( &
         MeanMajorMass*Gravity_GB(:, :, :, iBlock))
  Hn2 = -Boltzmanns_Constant* &
        Temperature(:, :, :, iBlock)*TempUnit/( &
        Mass(iN2_)*Gravity_GB(:, :, :, iBlock))
  Ho = -Boltzmanns_Constant* &
       Temperature(:, :, :, iBlock)*TempUnit/( &
       Mass(iO_3P_)*Gravity_GB(:, :, :, iBlock))
  Hco2t = -Boltzmanns_Constant* &
          Temperature(:, :, :, iBlock)*TempUnit/( &
          Mass(iCO2_)*Gravity_GB(:, :, :, iBlock))

  ! This calculates the ratio between the current average scale height
  ! and the Oxygen scale height.  If the scale height is the oxygen
  ! scale height, then the atmosphere is in molecular diffusion.  If it
  ! is far away from the oxygen scale height (closer to the N2 scale height)
  ! the atmosphere is well mixed and we should use the Have.

  r = (Ho - Have)/(Ho - Hn2)
  where (r > 1.0) r = 1.0
  where (r < 0.9) r = 0.0

  Hco2 = (1.0 - r)*Hco2t + r*Have

  ! Start at the bottom of the model:

  do iAlt = -1, 0
    NDensityS(:, :, iAlt, iCO2_, iBlock) = &
      CO2ppm*1e-6/(1.0 - CO2ppm*1e-6)*NDensity(:, :, iAlt, iBlock)
  enddo

  ! Then do hydrostatic to the top using the inferred scale height.

  do iAlt = 1, nAlts + 2
    NDensityS(:, :, iAlt, iCO2_, iBlock) = &
      NDensityS(:, :, iAlt - 1, iCO2_, iBlock)* &
      Temperature(:, :, iAlt - 1, iBlock)/Temperature(:, :, iAlt, iBlock)* &
      exp(-dAlt_GB(:, :, iAlt, iBlock)/Hco2(:, :, iAlt))
  enddo

end subroutine calc_co2
