! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

!--------------------------------------------------------------
!  Corrections by S. W. Bougher (9/28/07)
!  -- atomic oxygen replaces O2 in mean mass and scale height
!  -- artificially set all ion densities to 1.0 m-3 (not 1.0E+06 m-3)
!  Corrections by S. W. Bougher (11/07/12)
!  -- artificially set all ion densities to 1.0 m-3 (not 1.0E+04 m-3)
!  Corrections by S. W. Bougher (11/28/12)
!  -- artificially set O2+ and CO2+ ion densities to 1.0 m-3 (LBC at 80 km)
!--------------------------------------------------------------

subroutine get_msis_temperature(lon, lat, alt, t, h)

  use ModTime
  use ModInputs
  use ModPlanet
  use ModGITM

  implicit none

  real, intent(in) :: lon, lat, alt
  real, intent(out) :: t, h

! real :: nCO2, nO2, nN2, nCO, m, r, g
! real :: nCO2, nO,  nN2, nCO, m, r, g
  real :: nCO2, nOX, nN2, nCO, m, r, g
  integer :: i

  i = 1
! do while (alt >= newalt(i))
  do while ((alt >= newalt(i)) .and. (i <= nAlts + 2))
    i = i + 1
  enddo
  i = i - 1

  t = InTemp(i)
  nCO2 = InNDensityS(i, iCO2_)
! nO2 = InNDensityS(i,iO2_)
! nO   = InNDensityS(i,iO_)
  nOX = InNDensityS(i, iO_)
  nCO = InNDensityS(i, iCO_)
  nN2 = InNDensityS(i, iN2_)

! m = (nCO2 * mass(iCO2_) + &
!      nO2 * mass(iO2_) + &
!      nN2 * mass(iN2_) + &
!      nCO * mass(iCO_)) / (nCO2 + nO2 + nN2 + nCO)
! m = (nCO2 * mass(iCO2_) + &
!      nO   * mass(iO_) + &
!      nN2 * mass(iN2_) + &
!      nCO * mass(iCO_)) / (nCO2 + nO  + nN2 + nCO)

  m = (nCO2*mass(iCO2_) + &
       nOX*mass(iO_) + &
       nN2*mass(iN2_) + &
       nCO*mass(iCO_))/(nCO2 + nOX + nN2 + nCO)

  r = RBody + alt
!  g = Gravitational_Constant * (RBody/r) ** 2
  g = Gravitational_Constant

  h = Boltzmanns_Constant*t/(m*g)

end subroutine get_msis_temperature

subroutine init_msis

  use ModPlanet
  use ModGITM
  use ModEUV
  use ModInputs

  implicit none

  integer, parameter:: ninitialAlts = 25
  integer :: iBlock
  integer :: iiLon, iiLat, iiAlt, iminiono, ialtlow(1), TimeArray(7), ilatlow(1)
  integer :: iLon, iLat, iAlt, iSpecies, iIon, iError, jlat, klon, iline

  logical :: Done = .False., NotStarted = .True.
  character(len=iCharLen_) :: cLine
  real :: inDensities(9), altlow, althigh, latlow, lathigh
  real :: ralt, invAltDiff, altFind, altdiff, LogElectronDensity, dalt(nspeciestotal), alttemp(nInAlts)
  real, dimension(nInitialAlts) :: tempalt, LogInitialDensity, InitialEDensity, InitialAlt

  SurfaceAlbedo(:, :, :) = 0.0

  do iblock = 1, nblocks
    do ilon = 1, nlons
      do ilat = 1, nlats
        jlat = nint((latitude(ilat, iblock)*180.0/pi + 93.75)/7.5)
        klon = nint((longitude(ilon, iblock)*180.0/pi + 5.0)/10.0)
        SurfaceAlbedo(ilon, ilat, iblock) = dummyalbedo(jlat, klon)
        tinertia(ilon, ilat, iblock) = dummyti(jlat, klon)
      enddo
    enddo
  enddo

  if (useDustDistribution) call read_dust

  if (DoRestart) return

  do iBlock = 1, nBlocks
    !write(*,*) "init_msis"
    !write(*,*) '==> Now Initializing Mars Background Composition', iBlock

    dSurfaceTemp(:, :, :) = 0.0
    dSubSurfaceTemp(:, :, :) = 0.0
    SurfaceTemp(:, :, :) = 170.0
    SubsurfaceTemp(:, :, :) = 180.0
    !\
    ! Initializes the Planet with the same Chemistry as Above
    !/

    initialEDensity = 0.0
    initialAlt = 0.0

    open(iInputUnit_, file='DataIn/MarsInitialIonosphere.txt')
    do while (.not. Done)
      read(iInputUnit_, *) cLine
      if (cline .eq. '#START') Done = .True.
    enddo

    Done = .False.
    ialt = 1

    do while (.not. Done)
      read(iInputUnit_, *, iostat=iError) inDensities
      if (iError .ne. 0) then
        Done = .True.
      else
        initialAlt(ialt) = inDensities(1)
        !Convert to m^-3
        initialEDensity(ialt) = inDensities(9)*1.0e6
      endif

      iAlt = iAlt + 1

    enddo
    close(iInputUnit_)

    LogInitialDensity = log10(InitialEDensity)

    do iLat = -1, nLats + 2
      iiLat = min(max(iLat, 1), nLats)
      do iLon = -1, nLons + 2
        iiLon = min(max(iLon, 1), nLons)

        iMiniono = iAltMinIono(iiLon, iiLat, iBlock)

        do iAlt = iMinIono - 2, nalts + 2
          iialt = max(-1, ialt)
          tempalt = initialAlt
          altFind = altitude_GB(iLon, iLat, iiAlt, iBlock)/1000.0

          where (altfind - tempalt .lt. 0) tempalt = -1.0e9

          ialtlow = maxloc(tempalt)

          if (ialtlow(1) .eq. ninitialalts) ialtlow(1) = ialtlow(1) - 1

          altlow = initialAlt(ialtlow(1))
          althigh = initialAlt(ialtlow(1) + 1)

          invaltdiff = 1/(althigh - altlow)

          if (altFind .lt. initialAlt(1)) then

            ralt = (altlow - altFind)*(LogInitialDensity(ialtlow(1) + 1) - LogInitialDensity(ialtlow(1)))* &
                   invAltDiff
            LogElectronDensity = LogInitialDensity(ialtlow(1)) - ralt

          else

            if (altFind .gt. initialAlt(ninitialalts)) then

              ralt = (altFind - althigh)* &
                     (LogInitialDensity(ialtlow(1) + 1) - LogInitialDensity(ialtlow(1)))*invAltDiff
              LogElectronDensity = LogInitialDensity(ialtlow(1) + 1) + ralt

            else
              ralt = (Althigh - altFind)*(LogInitialDensity(ialtlow(1) + 1) - &
                                          LogInitialDensity(ialtlow(1)))*invAltDiff
              LogElectronDensity = LogInitialDensity(ialtlow(1) + 1) - ralt

            endif
          endif
          IDensityS(iLon, iLat, iialt, iE_, iBlock) = 10**LogElectronDensity

        enddo

        where (IDensityS(iLon, iLat, :, iE_, iBlock) .lt. 1.0) IDensityS(iLon, iLat, :, iE_, iBlock) = 1.0

        do iAlt = -1, nalts + 2
          alttemp = newalt

          altFind = altitude_GB(iLon, iLat, iAlt, iBlock)/1000.0
          where (altfind - alttemp .lt. 0) alttemp = -1.0e9

          ialtlow = maxloc(alttemp)

          if (ialtlow(1) .eq. ninalts) ialtlow(1) = ialtlow(1) - 1

          altlow = newalt(ialtlow(1))
          althigh = newalt(ialtlow(1) + 1)

          invaltdiff = 1/(althigh - altlow)

          if (altFind .lt. newalt(1)) then

            dalt = (altlow - altFind)*(InNDensityS(ialtlow(1) + 1, :) - inNDensityS(ialtlow(1), :))* &
                   invAltDiff
            NDensityS(iLon, iLat, ialt, :, iBlock) = inndensitys(ialtlow(1), :) - dalt
            ralt = (altlow - altFind)*(InTemp(ialtlow(1) + 1) - inTemp(ialtlow(1)))* &
                   invAltDiff
            Temperature(iLon, iLat, iAlt, iBlock) = InTemp(ialtlow(1)) - ralt
          else

            if (altFind .gt. newalt(ninalts)) then

              dalt = (altFind - althigh)* &
                     (inndensitys(ialtlow(1) + 1, :) - inndensitys(ialtlow(1), :))*invAltDiff
              NDensityS(iLon, iLat, ialt, :, iBlock) = Inndensitys(ialtlow(1) + 1, :) + dalt
              ralt = (altFind - althigh)*(InTemp(ialtlow(1) + 1) - inTemp(ialtlow(1)))* &
                     invAltDiff
              Temperature(iLon, iLat, iAlt, iBlock) = InTemp(ialtlow(1)) + ralt

            else
              dalt = (Althigh - altFind)*(inNDensitys(ialtlow(1) + 1, :) - &
                                          Inndensitys(ialtlow(1), :))*invAltDiff
              NDensityS(iLon, iLat, ialt, :, iBlock) = inNDensityS(ialtlow(1) + 1, :) - dalt
              ralt = (althigh - altFind)*(InTemp(ialtlow(1) + 1) - inTemp(ialtlow(1)))* &
                     invAltDiff
              Temperature(iLon, iLat, iAlt, iBlock) = InTemp(ialtlow(1)) - ralt
            endif
          endif

        enddo

        !           NDensityS(iLon,iLat,-1:nAlts + 2,1:nSpeciesTotal,iBlock) =  1.0

        !           do iSpecies = 1, nSpeciesTotal
        !
        !              ! If Mars_input.f90 then *1.0e+06 Converts from cm^-3 to m^-3
        !              ! Do not multiply by 1.0e+06 for MarsAtmosphere.txt - the
        !              ! densities in the file have already been done so
        !              NDensityS(iLon,iLat,-1:nAlts+2,iSpecies,iBlock) =  &
        !                   InNDensityS(-1:nAlts+2,iSpecies)!*(1.0e+06)
        !
        !           enddo ! end inner ispecies loop
!           NDensityS(iLon,iLat,-1:nAlts+2,iN2_,iBlock) = &
!                NDensityS(iLon,iLat,-1:nAlts+2,iN2_,iBlock) *0.1
!           NDensityS(iLon,iLat,-1:nAlts+2,iO_,iBlock) = &
!                NDensityS(iLon,iLat,-1:nAlts+2,iO_,iBlock) *0.1
!           NDensityS(iLon,iLat,-1:nAlts+2,iAr_,iBlock) = &
!                NDensityS(iLon,iLat,-1:nAlts+2,iAr_,iBlock) *0.1
!           NDensityS(iLon,iLat,-1:nAlts+2,iO2_,iBlock) = &
!                NDensityS(iLon,iLat,-1:nAlts+2,iO2_,iBlock) *0.1
!           NDensityS(iLon,iLat,-1:nAlts+2,iNO_,iBlock) = &
!                NDensityS(iLon,iLat,-1:nAlts+2,iNO_,iBlock) *0.1

        !\
        ! These first few are from the input file read in earlier
        !

        !\
        ! This just arbitrarily sets the ionospheric densities to 1.0
        ! and adds them up to get the electron density
        !

        !     IDensityS(iLon,iLat,-1:nAlts + 2,ie_,iBlock) = 0.0

        !do iIon = 1, nIons-1
        !
           !! IDensityS(iLon,iLat,-1:nAlts + 2,iIon,iBlock) =  1.0e6
        !  IDensityS(iLon,iLat,-1:nAlts + 2,iIon,iBlock) =  1.0
        !
        !   IDensityS(iLon,iLat,-1:nAlts + 2,ie_,iBlock) = &
        !        IDensityS(iLon,iLat,-1:nAlts + 2,ie_,iBlock) + &
        !        IDensityS(iLon,iLat,-1:nAlts + 2,iIon,iBlock)
        !
        !enddo ! end inner ispecies loop

        !
        ! End ion loop
        !/

      enddo! end iLon loop
    enddo ! end iLat loop

!  Initialization of Major Ions for Ion Calculation
    IDensityS(:, :, :, iO2P_, iBlock) = 0.9*IDensityS(:, :, :, iE_, iBlock)
    IDensityS(:, :, :, iCO2P_, iBlock) = 0.1*IDensityS(:, :, :, iE_, iBlock)
!  LBC for all  Ions for Ion Calculation
!    IDensityS(:,:,:,iO2P_,iBlock) = 1.0e0
!    IDensityS(:,:,:,iCO2P_,iBlock) = 1.0e0
    IDensityS(:, :, :, iOP_, iBlock) = 1.0e0
    IDensityS(:, :, :, iN2P_, iBlock) = 1.0e0
    IDensityS(:, :, :, iNOP_, iBlock) = 1.0e0

    !BP (initialize all to 1)
    IDensityS(:, :, :, :, :) = 1.0e0

    !write(*,*) '============> init_msis.Mars.f90 Major Diagnostics:  Begin'
!     Temperature(:,:,:,iBlock) = 175.
    !\
    ! Altitude Ghost Cells

    Temperature(:, :, -1, iBlock) = Temperature(:, :, 1, iBlock)
    Temperature(:, :, 0, iBlock) = Temperature(:, :, 1, iBlock)

    Temperature(:, :, nAlts + 1, iBlock) = Temperature(:, :, nAlts, iBlock)
    Temperature(:, :, nAlts + 2, iBlock) = Temperature(:, :, nAlts, iBlock)

    !\
    ! Longitude Ghost Cells

    Temperature(-1, :, :, iBlock) = Temperature(1, :, :, iBlock)
    Temperature(0, :, :, iBlock) = Temperature(1, :, :, iBlock)

    Temperature(nLons + 1, :, :, iBlock) = Temperature(nLons, :, :, iBlock)
    Temperature(nLons + 2, :, :, iBlock) = Temperature(nLons, :, :, iBlock)

    !\
    ! Latitude Ghost Cells

    Temperature(:, -1, :, iBlock) = Temperature(:, 1, :, iBlock)
    Temperature(:, 0, :, iBlock) = Temperature(:, 1, :, iBlock)

    Temperature(:, nLats + 1, :, iBlock) = Temperature(:, nLats, :, iBlock)
    Temperature(:, nLats + 2, :, iBlock) = Temperature(:, nLats, :, iBlock)

    !\
    ! Calculating MeanMajorMass -----------------------------
    !/

    !\
    ! Initialize MeanMajorMass to 0.0
    !/
    NDensityS = exp(nDensityS)
    MeanMajorMass(-1:nLons + 2, -1:nLats + 2, -1:nAlts + 2) = 0.0
    MeanIonMass(-1:nLons + 2, -1:nLats + 2, -1:nAlts + 2) = 0.0

    ! Calculate MeanMajorMass -----------------------------
    ! Calculate TempUnit -----------------------------

    do iLat = -1, nLats + 2
      do iLon = -1, nLons + 2
        do iAlt = -1, nAlts + 2

          NDensity(iLon, iLat, iAlt, iBlock) = 0.0

          do iSpecies = 1, nSpeciesTotal
            NDensity(iLon, iLat, iAlt, iBlock) = &
              NDensity(iLon, iLat, iAlt, iBlock) + &
              NDensityS(iLon, iLat, iAlt, iSpecies, iBlock)
          enddo

          do iSpecies = 1, nSpeciesTotal
            MeanMajorMass(iLon, iLat, iAlt) = &
              MeanMajorMass(iLon, iLat, iAlt) + &
              Mass(iSpecies)*NDensityS(iLon, iLat, iAlt, iSpecies, iBlock)/ &
              NDensity(iLon, iLat, iAlt, iBlock)
          enddo

          do iIon = 1, nIons - 1
            MeanIonMass(iLon, iLat, iAlt) = &
              MeanIonMass(iLon, iLat, iAlt) + &
              MassI(iIon)*IDensityS(iLon, iLat, iAlt, iIon, iBlock)/ &
              IDensityS(iLon, iLat, iAlt, ie_, iBlock)
          enddo

        enddo
      enddo
    enddo

    TempUnit(-1:nLons + 2, -1:nLats + 2, -1:nAlts + 2) = &
      MeanMajorMass(-1:nLons + 2, -1:nLats + 2, -1:nAlts + 2)/ &
      Boltzmanns_Constant

! More recent data

    NDensityS(-1:nLons + 2, -1:nLats + 2, -1:nAlts + 2, iHe_, iBlock) = &
      (2.0e-06)*NDensity(-1:nLons + 2, -1:nLats + 2, -1:nAlts + 2, iBlock)

! Previous EUV He observations
!     NDensityS(-1:nLons+2,-1:nLats+2,-1:nAlts+2,iHe_,iBlock) = &
!        (0.7e-06)*NDensity(-1:nLons+2,-1:nLats+2,-1:nAlts+2,iBlock)
    !\
    ! Initialize Rho to 0.0
    !/

    Rho(-1:nLons + 2, -1:nLats + 2, -1:nAlts + 2, iBlock) = 0.0

    Temperature(-1:nLons + 2, -1:nLats + 2, -1:nAlts + 2, iBlock) = &
      Temperature(-1:nLons + 2, -1:nLats + 2, -1:nAlts + 2, iBlock)/ &
      TempUnit(-1:nLons + 2, -1:nLats + 2, -1:nAlts + 2)

    Rho(-1:nLons + 2, -1:nLats + 2, -1:nAlts + 2, iBlock) = &
      MeanMajorMass(-1:nLons + 2, -1:nLats + 2, -1:nAlts + 2)* &
      NDensity(-1:nLons + 2, -1:nLats + 2, -1:nAlts + 2, iBlock)

    !write(*,*) '==> Now Completing Mars Background Composition: END', iBlock
    !BP: checking out the initialization of the ions
    call calc_electron_temperature(iBlock)
    !if (iProc == 0) then
    !  write(*,*) "Ion Densities in init_msis.Venus.f90:"
    !  write(*,*) IDensityS(:,:,:,ie_,iBlock)
    !endif

  enddo

end subroutine init_msis

subroutine msis_bcs(iJulianDay, UTime, Alt, Lat, Lon, Lst, &
                    F107A, F107, AP, LogNS, Temp, LogRho)

  write(*, *) "You can not use MSIS with any planet except Earth!!!"
  write(*, *) "If you ARE running Earth, then make the code again, using"
  write(*, *) "configure Earth ; make"
  call stop_gitm("I can not continue...")

end subroutine msis_bcs

subroutine read_dust
  use ModPlanet
  use ModInputs
  use ModGITM

  implicit none

  integer :: iBlock
  integer :: ialtlow(1), TimeArray(7), ilatlow(1)
  integer :: iLon, iLat, iAlt, iError, jlat, klon, iline

  logical :: Done = .False., NotStarted = .True.
  character(len=iCharLen_) :: cLine
  real :: altlow, althigh, latlow, lathigh
  real :: ralt, invAltDiff, altFind, altdiff, dalt(nspeciestotal), alttemp(nInAlts)

  call readDustHeader
  do iBlock = 1, nBlocks
    call setTau(iBlock)
  enddo

  call cleanDust

end subroutine read_dust

subroutine readDustHeader
  use ModPlanet
  use ModInputs, only: cDustFile, iCharLen_, DustFileType

  character(len=iCharLen_) :: cLine
  logical :: notstarted = .True.

  open(unit=iInputUnit_, file=cDustFile)

  do while (notstarted)
    read(iInputUnit_, *, iostat=iError) cLine

    if (iError .ne. 0) then
      write(*, *) "Error reading Dust file"
      write(*, *) "Is the header missing?"
      call stop_GITM("In init_msis_Mars")
    endif
    if (cline(1:7) .eq. '#HEADER') notstarted = .False.

  enddo

  if (DustFileType .eq. "FullHorizontal") then
    read(iInputUnit_, *, iostat=iError) nDustLats
    read(iInputUnit_, *, iostat=iError) nDustLons

    allocate(DustLatitude(nDustLats))
    allocate(DustLongitude(nDustLons))

    read(iInputUnit_, *, iostat=iError) DustLatitude
    read(iInputUnit_, *, iostat=iError) DustLongitude

  else if (DustFileType .eq. "MCSVertical") then
    read(iInputUnit_, *, iostat=iError) nDustLats
    read(iInputUnit_, *, iostat=iError) nDustTimes
    read(iInputUnit_, *, iostat=iError) nDustAlts

    nDustLons = 1

    allocate(DustLatitude(nDustLats))
  endif

  close(iInputUnit_)

end subroutine readDustHeader

subroutine setTau(iBlock)

  use ModInputs
  use ModPlanet
  use ModGITM, only: Latitude, Longitude, iproc

  integer, intent(IN) :: iBlock
  real :: TempDust(ndustlats, ndustlons), tempconrath(ndustlats), Temp(nDustLats*nDustLons), MCSTemp(5)
  real :: rlat, invLatDiff, LatFind, Latdiff, Dust, templat(ndustlats), templon(ndustlons), conrath, invlondiff
  real :: lathigh, latlow, lonhigh, lonlow, V11, V12, V21, V22
  character(len=iCharLen_) :: cLine
  logical :: notstarted

  real, dimension(nDustTimes, nDustLats, nDustAlts) :: CumulativeTau, DustMixingRatio

  integer :: TimeArray(7), i, ilat, ilon, ilatlow(1), ilonlow(1), ialt, itime

  iline = 1

  open(unit=iInputUnit_, file=cDustFile)

  notstarted = .True.

  do while (notstarted)
    read(iInputUnit_, *, iostat=iError) cLine
    if (cline(1:6) .eq. '#START') notstarted = .False.
  enddo

  if (DustFileType .eq. "FullHorizontal") then

    !The flies are ordered by lats first then lons.  I.e.,
    do while (iError .eq. 0)
      read(iInputUnit_, *, iostat=iError) TimeArray(1:6), Temp
      if (iproc .eq. 0) then
        write(*, *) TimeArray(1:6)
      endif
      i = 1
      do iLat = 1, ndustlats
        do iLon = 1, ndustlons
          TempDust(iLat, iLon) = Temp(i)
          i = i + 1
        enddo
      enddo
      TimeArray(7) = 0
      call time_int_to_real(TimeArray, TimeDust(iLine))

      do iLat = 1, nLats
        do iLon = 1, nLons

          latFind = Latitude(ilat, iBlock)*180/pi
          lonFind = Longitude(ilon, iBlock)*180/pi

          templat = DustLatitude
          templon = DustLongitude

          where (LatFind - tempLat .lt. -0.00001) tempLat = -1.0e9
          ilatlow = maxloc(tempLat)

          where (LonFind - tempLon < -0.00001) tempLon = -1.0e9
          ilonlow = maxloc(tempLon)

          if (ilatlow(1) .eq. nDustLats) ilatlow = ilatlow - 1
          if (ilonlow(1) .eq. nDustLons) ilonlow = ilonlow - 1

          Latlow = DustLatitude(ilatlow(1))
          LatHigh = DustLatitude(ilatlow(1) + 1)
          Lonlow = DustLongitude(ilonlow(1))
          LonHigh = DustLongitude(ilonlow(1) + 1)

          invLatdiff = 1/(Lathigh - Latlow)
          invLondiff = 1/(Lonhigh - Lonlow)

          if (LatFind .lt. DustLatitude(1) .or. LonFind .lt. DustLongitude(1) &
              .or. LatFind .gt. DustLatitude(nDustLats) .or. LonFind .ge. DustLongitude(nDustLons)) then
            write(*, *) 'Dust grid does not cover GITM grid'
            write(*, *) 'Stopping...'
            call stop_gitm('Stopping in init_msis.Mars')
          endif

          V11 = TempDust(ilatlow(1), ilonlow(1))
          V12 = TempDust(ilatlow(1) + 1, ilonlow(1))
          V21 = TempDust(ilatlow(1), ilonlow(1) + 1)
          V22 = TempDust(ilatlow(1) + 1, ilonlow(1) + 1)

          rlat1 = (Lathigh - LatFind)*V11*invLatDiff + &
                  (LatFind - LatLow)*V21*invLatDiff

          rlat2 = (Lathigh - LatFind)*V12*invLatDiff + &
                  (LatFind - LatLow)*V22*invLatDiff

          Dust = (Lonhigh - LonFind)*invLonDiff*rlat1 + (LonFind - LonLow)*invLonDiff*rlat2
!           if (iproc .eq. 1 .and. ilat .eq. 1 .and. ilon .eq. 5 .and. iline .eq. 54) then
!              write(*,*) iproc,ilat,ilon,iline,dust,lonhigh,lonfind,lonlow,lathigh,latfind,latlow
!              write(*,*) V11, V12, V21,V22
!              write(*,*) "lats: ",ilatlow(1)+1, ilonlow(1),ndustlats,ndustlons
!              stop
!           endif
          HorizontalDustProfile(iline, iLat, iLon, iblock) = Dust

        enddo
      enddo
      iline = iline + 1

    enddo
    nDustTimes = iline - 1

    nConrathTimes = nDustTimes
    TimeConrath = TimeDust
    HorizontalConrathProfile = 0.03

  else if (DustFileType .eq. "MCSVertical") then
    !Read in data
    do iTime = 1, nDustTimes
      do iLat = 1, nDustLats
        do iAlt = 1, nDustAlts
          read(iInputUnit_, *, iostat=iError) TimeArray(1:6), MCSTemp
          if (iError .ne. 0) then
            write(*, *) "Error reading dustfile"
            call stop_gitm('Stopping in init_msis.Mars')
          endif
          DustPressureLevel(iAlt) = MCSTemp(3)
          CumulativeTau(iTime, iLat, iAlt) = MCSTemp(4)
          DustMixingRatio(iTime, iLat, iAlt) = MCSTemp(5)

        enddo
        DustLatitude(iLat) = MCSTemp(1)
      enddo

      TimeArray(7) = 0
      call time_int_to_real(TimeArray, rTime)
      TimeDust(iTime) = rtime
    enddo

    do iLat = 1, nLats
      latFind = Latitude(ilat, iBlock)*180/pi
      templat = DustLatitude

      where (LatFind - tempLat .lt. -0.00001) tempLat = -1.0e9
      ilatlow = maxloc(tempLat)

      if (ilatlow(1) .eq. nDustLats) ilatlow = ilatlow - 1

      Latlow = DustLatitude(ilatlow(1))
      LatHigh = DustLatitude(ilatlow(1) + 1)
      invLatdiff = 1/(Lathigh - Latlow)

      if (LatFind .lt. DustLatitude(1)) then
        CumulativeTauProfile(1:nDustTimes, ilat, 1:nDustAlts, iblock) = &
          CumulativeTau(1:nDustTimes, 1, 1:nDustAlts)

        DustMixingRatioProfile(1:nDustTimes, ilat, 1:nDustAlts, iblock) = &
          DustMixingRatio(1:nDustTimes, 1, 1:nDustAlts)

      else if (LatFind .gt. DustLatitude(nDustLats)) then

        CumulativeTauProfile(1:nDustTimes, ilat, 1:nDustAlts, iblock) = &
          CumulativeTau(1:nDustTimes, nDustLats, 1:nDustAlts)

        DustMixingRatioProfile(1:nDustTimes, ilat, 1:nDustAlts, iblock) = &
          DustMixingRatio(1:nDustTimes, nDustLats, 1:nDustAlts)

      else
        !Multiply by .999 to handle case when grids match up.
        !Rounding error can result in negative numbers

        CumulativeTauProfile(1:nDustTimes, ilat, 1:nDustAlts, iblock) = &
          CumulativeTau(1:nDustTimes, ilatlow(1) + 1, 1:nDustAlts) - &
          (Lathigh - LatFind)*0.999*invLatDiff* &
          (CumulativeTau(1:nDustTimes, ilatlow(1) + 1, 1:nDustAlts) - &
           CumulativeTau(1:nDustTimes, ilatlow(1), 1:nDustAlts))

        DustMixingRatioProfile(1:nDustTimes, ilat, 1:nDustAlts, iblock) = &
          DustMixingRatio(1:nDustTimes, ilatlow(1) + 1, 1:nDustAlts) - &
          (Lathigh - LatFind)*0.999*invLatDiff* &
          (DustMixingRatio(1:nDustTimes, ilatlow(1) + 1, 1:nDustAlts) - &
           DustMixingRatio(1:nDustTimes, ilatlow(1), 1:nDustAlts))
      endif

    enddo
  endif

  close(iInputUnit_)

  !Horizontal Conrath Parameter Distribution

  !     TimeArray(1:6) = temp(1:6)
  !     TimeArray(7) = 0
  !     call time_int_to_real(TimeArray,TimeConrath(iLine))
  !     TempConrath = Temp(7:ndustlats+6)
  !
  !     do iLat = 1, nLats
  !        latFind = Latitude(ilat,iBlock)*180/pi
  !        templat = DustLatitude
  !        where(LatFind - tempLat .lt. -0.00001) tempLat = -1.0e9
  !
  !        ilatlow =  maxloc(tempLat)
  !
  !        if (ilatlow(1) .eq. nDustLats) ialtlow = ialtlow - 1
  !
  !        Latlow = DustLatitude(ilatlow(1))
  !        LatHigh = DustLatitude(ilatlow(1)+1)
  !
  !        invLatdiff = 1/(Lathigh - Latlow)
  !
  !        if (LatFind .lt. DustLatitude(1)) then
  !
  !           rlat = (latlow-latFind)*(TempConrath(ilatlow(1) + 1)-TempConrath(ilatlow(1))) * &
  !                invLatDiff
  !           Conrath = TempConrath(ilat) - rlat
  !
  !        else
  !
  !           if (LatFind .ge. DustLatitude(nDustLats)) then
  !
  !              rlat = (LatFind-Lathigh)* &
  !                   (TempConrath(ilatlow(1) + 1) - TempConrath(ilatlow(1))) * invLatDiff
  !
  !              Conrath = TempConrath(ilatlow(1) + 1) + rlat
  !
  !           else
  !              rlat = (Lathigh - LatFind)*(TempConrath(ilatlow(1) + 1) - &
  !                   TempConrath(ilatlow(1))) * invLatDiff
  !              Conrath = TempConrath(ilatlow(1) + 1) - rlat
  !
  !           endif
  !        endif
  !        HorizontalConrathProfile(iline,iLat,iblock) = Conrath
  !
  !     enddo
  !
  !     iline = iline + 1

end subroutine setTau

subroutine cleanDust
  use ModPlanet, only: DustLatitude, DustLongitude

  if (allocated(DustLatitude)) then
    deallocate(DustLatitude)
  endif
  if (allocated(DustLongitude)) then
    deallocate(DustLongitude)
  endif

end subroutine cleanDust
