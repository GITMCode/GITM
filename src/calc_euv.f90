! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

subroutine euv_ionization_heat(iBlock)

  use ModGITM
  use ModEUV
  use ModPlanet
  use ModConstants
  use ModInputs
  use ModSources
  use ModTime, only: tSimulation, CurrentTime
  use ModIndicesInterfaces

  implicit none

  integer, intent(in) :: iBlock

  integer :: iAlt, iWave, iSpecies, iNeutral, iIon, iError, iLon, iLat
  integer :: iRxn
  real, dimension(nLons, nLats) :: Tau, Intensity

  logical :: IsFirstTime(nBlocksMax) = .true.

  real :: NeutralDensity(nLons, nLats, nSpeciesTotal)
  real :: ChapmanLittle(nLons, nLats, nSpecies)
  real :: EHeat(nLons, nLats)
  real :: nEuvHeating(nLons, nLats, nAlts), neEuvHeating(nLons, nLats, nAlts)

  if (IsFirstTime(iBlock)) then

    IsFirstTime(iBlock) = .false.

    ! This transfers the specific photo absorption and ionization cross
    ! sections into general variables, so we can use loops...

    call fill_photo

  else
    if (floor((tSimulation - dT)/dTEUV) == &
        floor(tSimulation/dTEUV)) return
  endif

  call report("euv_ionization_heat", 2)
  call start_timing("euv_ionization_heat")
  call chapman_integrals(iBlock)

  EuvIonRate = 0.0
  EuvHeating(:, :, :, iBlock) = 0.0
  PhotoElectronHeating(:, :, :, iBlock) = 0.0
  eEuvHeating(:, :, :, iBlock) = 0.0
  EuvIonRateS(:, :, :, :, iBlock) = 0.0
  EuvDissRateS(:, :, :, :, iBlock) = 0.0
  nEuvHeating(:, :, :) = 0.0
  neEuvHeating(:, :, :) = 0.0

  do iAlt = 1, nAlts

    NeutralDensity = NDensityS(1:nLons, 1:nLats, iAlt, 1:nSpeciesTotal, iBlock)
    ChapmanLittle = Chapman(:, :, iAlt, 1:nSpecies, iBlock)
    EHeat = 0.0

    do iWave = 1, Num_WaveLengths_High

      Tau = 0.0
      do iSpecies = 1, nSpecies
        Tau = Tau + &
              photoabs(iWave, iSpecies)*ChapmanLittle(:, :, iSpecies)
      enddo

      Intensity = Flux_of_EUV(iWave)*exp(-1.0*Tau)

      do iIon = 1, nIons - 1
        iNeutral = PhotoIonFrom(iIon)
        ! If we have an ionization cross section, use it with the
        ! appropriate neutral density, otherwise, there is no ionization
        if (iNeutral > 0) then
          EuvIonRateS(:, :, iAlt, iIon, iBlock) = &
            EuvIonRateS(:, :, iAlt, iIon, iBlock) + &
            Intensity*PhotoIon(iWave, iIon)* &
            NeutralDensity(:, :, iNeutral)* &
            (1.0 + PhotoElecIon(iWave, iIon))
        else
          EuvIonRateS(:, :, iAlt, iIon, iBlock) = 0.0
        endif
      enddo

      do iSpecies = 1, nSpecies
        EuvDissRateS(:, :, iAlt, iSpecies, iBlock) = &
          EuvDissRateS(:, :, iAlt, iSpecies, iBlock) + &
          Intensity*PhotoDis(iWave, iSpecies)* &
          NeutralDensity(:, :, iSpecies)* &
          (1.0 + PhotoElecDiss(iWave, iSpecies))
      enddo

      do iSpecies = 1, nSpecies
        EHeat = EHeat + &
                Intensity*PhotonEnergy(iWave)* &
                photoabs(iWave, iSpecies)*NeutralDensity(:, :, iSpecies)
      enddo

    enddo

    EuvHeating(:, :, iAlt, iBlock) = EHeat*HeatingEfficiency_CB(:, :, iAlt, iBlock)
    eEuvHeating(:, :, iAlt, iBlock) = EHeat*eHeatingEfficiency_CB(:, :, iAlt, iBlock)

    ! Why this? Is this for Mars or Venus or something?
    !do ilon = 1, nlons
    !   do ilat =1 ,nlats
    !      if (Altitude_GB(iLon,iLat,iAlt,iBlock) .lt. 80000.0) then
    !         EUVHeating(iLon,iLat,iAlt,iBlock) =0.0
    !         eEUVHeating(iLon,iLat,iAlt,iBlock) =0.0
    !      endif
    !   enddo
    !enddo

  enddo

  if (IncludeEclipse) call calc_eclipse_effects
  if (IsEarth) call night_euv_ionization

  EuvIonRateS = EuvIonRateS + nEuvIonRateS
  EuvHeating(:, :, :, iBlock) = EuvHeating(:, :, 1:nAlts, iBlock) + nEuvHeating
  eEuvHeating(:, :, :, iBlock) = eEuvHeating(:, :, 1:nAlts, iBlock) + neEuvHeating

  !\
  ! Zero out EuvHeating if specified not to use it.
  !/

  if (UseSolarHeating) then

    EuvHeating2d = 0.0

    do iAlt = 1, nAlts

      EuvHeating2d(1:nLons, 1:nLats) = &
        EuvHeating2d(1:nLons, 1:nLats) + &
        EuvHeating(:, :, iAlt, iBlock)* &
        dAlt_GB(1:nLons, 1:nLats, iAlt, iBlock)

      EuvHeating(:, :, iAlt, iBlock) = EuvHeating(:, :, iAlt, iBlock)/ &
                                       Rho(1:nLons, 1:nLats, iAlt, iBlock)/ &
                                       cp(1:nLons, 1:nLats, iAlt, iBlock)/ &
                                       TempUnit(1:nLons, 1:nLats, iAlt)

      EuvTotal(:, :, iAlt, iBlock) = EuvHeating(:, :, iAlt, iBlock)* &
                                     TempUnit(1:nLons, 1:nLats, iAlt)/ &
                                     HeatingEfficiency_CB(:, :, iAlt, iBlock)

    enddo

  else
    EuvHeating = 0.0
  endif

  call end_timing("euv_ionization_heat")

contains

  !\
  ! -------------------------------------------------------------
  !/

  subroutine calc_eclipse_effects

    real :: x, y, z, xp, yp, zp, distance, PercentDone
    real :: yPos, zPos, e
    real :: factor(nLons, nLats, nAlts)

    Factor = 1.0

    if (CurrentTime > EclipseStartTime .and. &
        CurrentTime < EclipseEndTime) then

      PercentDone = (CurrentTime - EclipseStartTime)/ &
                    (EclipseEndTime - EclipseStartTime)

      yPos = EclipseStartY + (EclipseEndY - EclipseStartY)*PercentDone
      zPos = EclipseStartZ + (EclipseEndZ - EclipseStartZ)*PercentDone

      ! Need to rotate lat/lon/alt into solar coordinate X,Y,Z system

      !   x = r * cos(localtime-!pi) * cos(lats*!dtor)
      !   y = r * sin(localtime-!pi) * cos(lats*!dtor)
      !   z = r * sin(lats*!dtor)
      !
      !   xp = x * cos(tilt) - z * sin(tilt)
      !   zp = x * sin(tilt) + z * cos(tilt)

      do iLon = 1, nLons
        do iLat = 1, nLats
          do iAlt = 1, nAlts
            ! Calculate X,Y,Z for points on Earth (in lat/localtime/alt)
            x = &
              RadialDistance_GB(iLon, iLat, iAlt, iBlock)* &
              cos(localtime(iLon)*pi/12.0 - pi)* &
              cos(Latitude(iLat, iBlock))
            y = &
              RadialDistance_GB(iLon, iLat, iAlt, iBlock)* &
              sin(localtime(iLon)*pi/12.0 - pi)* &
              cos(Latitude(iLat, iBlock))
            z = &
              RadialDistance_GB(iLon, iLat, iAlt, iBlock)* &
              sin(Latitude(iLat, iBlock))

            ! rotate around y-axis for the Earth Tilt:
            xp = x*cos(-SunDeclination) - z*sin(-SunDeclination)
            yp = y
            zp = x*sin(-SunDeclination) + z*cos(-SunDeclination)

            distance = sqrt((yPos - yp)**2 + (zPos - zp)**2)

            ! If it is close to the eclipse, do calculations

            if (Distance < EclipseMaxDistance + EclipseExpWidth*2) then

              ! Assume a linearly decreasing effect:
              e = EclipsePeak - EclipsePeak*distance/EclipseMaxDistance
              if (e < 0) e = 0

              ! Towards the edge, assume an exponential:
              e = e + EclipseExpAmp*exp(-((Distance - EclipseMaxDistance)/EclipseExpWidth)**2)
              Factor(iLon, iLat, iAlt) = 1.0 - e

              ! Change the ionization rate and heating rate by the factor:
              do iIon = 1, nIons - 1
                EuvIonRateS(iLon, iLat, iAlt, iIon, iBlock) = &
                  EuvIonRateS(iLon, iLat, iAlt, iIon, iBlock)* &
                  Factor(iLon, iLat, iAlt)
              enddo
              EuvHeating(iLon, iLat, iAlt, iBlock) = &
                EuvHeating(iLon, iLat, iAlt, iBlock)* &
                Factor(iLon, iLat, iAlt)

              N2PERateS(iLon, iLat, iAlt, :, iBlock) = &
                N2PERateS(iLon, iLat, iAlt, :, iBlock)*Factor(iLon, iLat, iAlt)
              O2PERateS(iLon, iLat, iAlt, :, iBlock) = &
                O2PERateS(iLon, iLat, iAlt, :, iBlock)*Factor(iLon, iLat, iAlt)
              OPERateS(iLon, iLat, iAlt, :, iBlock) = &
                OPERateS(iLon, iLat, iAlt, :, iBlock)*Factor(iLon, iLat, iAlt)
              CH4PERateS(iLon, iLat, iAlt, :, iBlock) = &
                CH4PERateS(iLon, iLat, iAlt, :, iBlock)*Factor(iLon, iLat, iAlt)

            endif

          enddo
        enddo
      enddo

    endif

  end subroutine calc_eclipse_effects

  !\
  ! -------------------------------------------------------------
  !/

  subroutine night_euv_ionization

    ! Nighttime EUV Ionization

    real :: night_col(nLons, nLats, nAlts, nSpecies)
    real :: dAlts(nLons, nLats)
    real, dimension(nLons, nLats) :: nTau, nIntensity
    real :: nNeutralDensity(nLons, nLats, nSpeciesTotal)
    real :: tempflux(nLons, nLats), SZALocal(nLons, nLats)
    real :: ScaleHeightS(nLons, nLats)
    real :: scatterfac = 0.01
    integer :: sinexpo = 15

    real :: nEHeat(nLons, nLats)
    real :: nPhotonEnergy(Num_NightWaveLens)
    real :: tempPhotonEnergy = 0.0

    nEuvIonRates = 0.0
    night_col = 0

    SZALocal = SZA(1:nLons, 1:nLats, iBlock)

    do iWave = 1, Num_NightWaveLens
      tempflux = 0.0
      select case (iWave)
      case (1)
        where (SZALocal .GE. Pi/2.) tempflux = -0.573563*SZALocal**3 &
                                               + 4.66159*SZALocal**2 - 12.6936*SZALocal + 11.66

        where (SZALocal .LT. Pi/2.) tempflux = (-0.573563*(Pi/2.)**3 &
                                                + 4.66159*(Pi/2.)**2 - 12.6936*(Pi/2.) + 11.66) &
                                               *sin(SZALocal)**sinexpo

        tempflux = tempflux*Flux_of_EUV(18)*scatterfac

        where (tempflux .LE. 0.) tempflux = 0.0

        tempPhotonEnergy = PhotonEnergy(18)

      case (2)
        where (SZALocal .GE. Pi/2.) tempflux = -4.39459*SZALocal**3 &
                                               + 28.8112*SZALocal**2 - 62.9068*SZALocal + 45.7575

        where (SZALocal .LT. Pi/2.) tempflux = (-4.39459*(Pi/2.)**3 &
                                                + 28.8112*(Pi/2.)**2 - 62.9068*(Pi/2.) + 45.7575) &
                                               *sin(SZALocal)**sinexpo

        tempflux = tempflux*Flux_of_EUV(35)*scatterfac

        where (tempflux .LE. 0.0) tempflux = 0.0

        tempPhotonEnergy = PhotonEnergy(35)

      case (4)
        where (SZALocal .GE. Pi/2.) tempflux = -0.356061*SZALocal**3 &
                                               + 3.0883*SZALocal**2 - 8.98066*SZALocal + 8.86673

        where (SZALocal .LT. Pi/2.) tempflux = (-0.356061*(Pi/2.)**3 &
                                                + 3.0883*(Pi/2.)**2 - 8.98066*(Pi/2.) + 8.86673) &
                                               *sin(SZALocal)**sinexpo

        tempflux = tempflux*Flux_of_EUV(12)*scatterfac

        where (tempflux .LE. 0.0) tempflux = 0.0

        tempPhotonEnergy = PhotonEnergy(12)

      case default
        tempflux = 0.0
        tempPhotonEnergy = 0.0
      end select

      nighteuvflux(iWave, :, :, iBlock) = tempflux
      nPhotonEnergy(iWave) = tempPhotonEnergy

    enddo

    !! Calculate Tau
    do iAlt = nAlts, 1, -1
      nNeutralDensity = NDensityS(1:nLons, 1:nLats, iAlt, 1:nSpeciesTotal, iBlock)
      dAlts = dAlt_GB(1:nLons, 1:nLats, iAlt, iBlock)

      do iSpecies = 1, nSpecies

        if (iAlt .LT. nAlts) then
          night_col(:, :, iAlt, iSpecies) = night_col(:, :, iAlt + 1, iSpecies) &
                                            + nNeutralDensity(:, :, iSpecies)*dAlts
        else
          ScaleHeightS = &
            Temperature(1:nLons, 1:nLats, iAlt, iBlock) &
            *TempUnit(1:nLons, 1:nLats, iAlt)*Boltzmanns_Constant &
            /(-Gravity_GB(1:nLons, 1:nLats, iAlt, iBlock)*Mass(iSpecies))
          night_col(:, :, iAlt, iSpecies) = nNeutralDensity(:, :, iSpecies) &
                                            *ScaleHeightS
        endif

      enddo
    enddo

    do iAlt = 1, nAlts
      nEHeat = 0.0
      nNeutralDensity = NDensityS( &
                        1:nLons, 1:nLats, iAlt, 1:nSpeciesTotal, iBlock)
      do iWave = 1, Num_NightWaveLens
        nTau = 0.0

        do iSpecies = 1, nSpecies
          nTau = nTau + &
                 night_photoabs(iWave, iSpecies)*night_col(:, :, iAlt, iSpecies)
        enddo

        nIntensity = nighteuvflux(iWave, :, :, iBlock)*exp(-1.0*nTau)

        do iIon = 1, nIons - 1
          iNeutral = PhotoIonFrom(iIon)
          ! If we have an ionization cross section, use it with the
          ! appropriate neutral density, otherwise, there is no ionization
          if (iNeutral > 0) then
            nEuvIonRateS(:, :, iAlt, iIon, iBlock) = &
              nEuvIonRateS(:, :, iAlt, iIon, iBlock) + &
              nIntensity*night_photoion(iWave, iIon)* &
              nNeutralDensity(:, :, iNeutral)
          endif
        enddo

        do iSpecies = 1, nSpecies
          nEHeat = nEHeat + &
                   nIntensity* &
                   nPhotonEnergy(iWave)* &
                   night_photoabs(iWave, iSpecies)* &
                   nNeutralDensity(:, :, iSpecies)

        enddo

      enddo

      nEuvHeating(:, :, iAlt) = nEHeat*HeatingEfficiency_CB(:, :, iAlt, iBlock)
      neEuvHeating(:, :, iAlt) = nEHeat*eHeatingEfficiency_CB(:, :, iAlt, iBlock)
      do ilon = 1, nlons
        do ilat = 1, nlats
          if (Altitude_GB(iLon, iLat, iAlt, iBlock) .lt. 80000.0) then
            nEUVHeating(iLon, iLat, iAlt) = 0.0
            neEUVHeating(iLon, iLat, iAlt) = 0.0
          endif
        enddo
      enddo

    enddo

  end subroutine night_euv_ionization

end subroutine euv_ionization_heat

!-------------------------------------------------------------------
! Subroutine for initializing photoabsorption, photoionization,
! and branching ratio quantities.
!-------------------------------------------------------------------

subroutine init_euv

  use ModEUV
  use ModTime
  use ModInputs

  implicit none

  integer :: N, NN, iError
  real :: wavelength_ave

  call report("init_euv", 2)

  do N = 1, Num_WaveLengths_High
    ! Calculate the energy in the bin:
    waveLengthMid(n) = (WAVEL(N) + WAVES(N))/2.0
    PhotonEnergy(N) = &
      Planck_Constant*Speed_Light/ &
      (waveLengthMid(n)*1.0e-10)
  enddo

  EUVEFF = 0.05

  DO N = 1, Num_WaveLengths_Low

    NN = N + 15    ! 16:52

    BranchingRatio_OPlus2P(N) = 0.
    IF (N .GT. 14) then
      BranchingRatio_OPlus2P(N) = &
        1.-BranchingRatio_OPlus2D(N) - BranchingRatio_OPlus4S(N)
    endif

    PhotoAbs_O2(NN) = Photoabsorption_O2(N)
    PhotoAbs_O(NN) = Photoabsorption_O(N)
    PhotoAbs_N2(NN) = Photoabsorption_N2(N)
    PhotoAbs_CO2(NN) = Photoabsorption_CO2(N)
    PhotoAbs_CO(NN) = Photoabsorption_CO(N)

    PhotoIon_O2(NN) = Photoionization_O2(N)
    PhotoIon_CO2(NN) = Photoionization_CO2(N)
    PhotoIon_CO(NN) = Photoionization_CO(N)
    PhotoIon_OPlus4S(NN) = Photoionization_O(N)*BranchingRatio_OPlus4S(N)
    PhotoIon_N2(NN) = Photoionization_N2(N)
    PhotoIon_N(NN) = Photoionization_N(N)
    PhotoIon_OPlus2D(NN) = Photoionization_O(N)*BranchingRatio_OPlus2D(N)
    PhotoIon_OPlus2P(NN) = Photoionization_O(N)*BranchingRatio_OPlus2P(N)

    BranchingRatio_N2(NN) = BranchingRatio_N2_to_NPlus(N)
    BranchingRatio_O2(NN) = BranchingRatio_O2_to_OPlus(N)

  enddo

  ! Convert everything from /cm2 to /m2

  PhotoAbs_O2 = PhotoAbs_O2/10000.0
  PhotoAbs_O = PhotoAbs_O/10000.0
  PhotoAbs_N2 = PhotoAbs_N2/10000.0
  PhotoAbs_CH4 = PhotoAbs_CH4/10000.0
  PhotoAbs_CO2 = PhotoAbs_CO2/10000.0
  PhotoAbs_CO = PhotoAbs_CO/10000.0
  PhotoAbs_He = PhotoAbs_He/10000.0
  PhotoAbs_H2 = PhotoAbs_H2/10000.0
  PhotoAbs_NO = PhotoAbs_NO/10000.0

  PhotoIon_O2 = PhotoIon_O2/10000.0
  PhotoIon_CO2 = PhotoIon_CO2/10000.0
  PhotoIon_CO = PhotoIon_CO/10000.0
  PhotoIon_OPlus4S = PhotoIon_OPlus4S/10000.0
  PhotoIon_N2 = PhotoIon_N2/10000.0
  PhotoIon_N = PhotoIon_N/10000.0
  PhotoIon_OPlus2D = PhotoIon_OPlus2D/10000.0
  PhotoIon_OPlus2P = PhotoIon_OPlus2P/10000.0
  PhotoIon_H2 = PhotoIon_H2/10000.0
  PhotoIon_NO = PhotoIon_NO/10000.0

  if (UseEUVData) call read_fism_file(iError, CurrentTime, EndTime)

  if (UseRidleyEUV) then
    if (iDebugLevel >= 0) then
      write(*,*) '> Setting EUV Model to Ridley Version ', RidleyEUVModelNumber
    endif
    if (RidleyEUVModelNumber == 0) then
      call read_euv_waves(iError)
    else
      call read_euv2026_model(iError)
    endif
  endif

end subroutine init_euv

