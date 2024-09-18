! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

!------------------------------------------------------------------------------
! Calculate the vertical advection terms for all lat / lon
!------------------------------------------------------------------------------

subroutine advance_vertical_all

  use ModInputs
  use ModGitm

  implicit none

  integer :: iBlock, iLon, iLat

  call report("advance_vertical_all", 2)
  call start_timing("vertical_all")

  do iBlock = 1, nBlocks

    call calc_rates(iBlock)
    call calc_viscosity(iBlock)

    do iLon = 1, nLons; do iLat = 1, nLats
        call advance_vertical(iLon, iLat, iBlock)
      end do; end do

  end do

  if (DoCheckForNans) then
    call check_for_nans_ions("After Vertical")
    call check_for_nans_neutrals("After Vertical")
    call check_for_nans_temps("After Vertical")
  end if

  call end_timing("vertical_all")

end subroutine advance_vertical_all

!------------------------------------------------------------------------------
! Calculate the vertical advection terms for an individual lat / lon
!------------------------------------------------------------------------------

subroutine advance_vertical(iLon, iLat, iBlock)
  use ModEUV, only: SZA
  use ModGITM
  use ModPlanet, only: nSpecies, nIonsAdvect
  use ModConstants, only: pi
  use ModSources, only: EUVHeating, KappaEddyDiffusion
  use ModInputs
  use ModVertical, ONLY: &
    LogRho, &
    cMax1 => cMax, &
    LogNS1 => LogNS, &
    Heating, &
    Vel_GD, &
    Temp, iLon1D, iLat1D, iBlock1D, &
    Centrifugal, Coriolis, &
    LogINS, &
    IVel, Lat, Lon, &
    VertVel, &
    MeanMajorMass_1d, &
    gamma_1d, &
    EddyCoef_1d, &
    ViscCoef_1d, &
    ViscCoefS_1d, &
    KappaTemp_1d, &
    Gravity_G, Altitude_G, dAlt_C, InvRadialDistance_C, dAlt_F, InvDAlt_F, &
    Cv_1D, dAltdLon_1D, dAltdLat_1D, SZAVertical, MLatVertical

  implicit none

  integer, intent(in) :: iLon, iLat, iBlock

  integer :: iIon, iSpecies, iAlt, iDim
  real    :: KappaTemp1(0:nAlts + 1)

  iLon1D = iLon
  iLat1D = iLat
  iBlock1D = iBlock

  SZAVertical = SZA(iLon, iLat, iBlock)
  MLatVertical = MLatitude(iLon, iLat, nAlts, iBlock)

  KappaTemp1 = KappaTemp(iLon, iLat, :, iBlock)
  dAltdLon_1D = dAltdLon_CB(iLon, iLat, 1, iBlock)
  dAltdLat_1D = dAltdLat_CB(iLon, iLat, 1, iBlock)
  EddyCoef_1d(0:nAlts + 1) = KappaEddyDiffusion(iLon, iLat, 0:nAlts + 1, iBlock)
  ViscCoef_1d(0:nAlts + 1) = ViscCoef(iLon, iLat, 0:nAlts + 1)
  ViscCoefS_1d(0:nAlts + 1, 1:nSpecies) = ViscCoefS(iLon, iLat, 0:nAlts + 1, 1:nSpecies)
  KappaTemp_1d(0:nAlts + 1) = KappaTemp(iLon, iLat, 0:nAlts + 1, iBlock)
  Cv_1D(-1:nAlts + 2) = cp(iLon, iLat, -1:nAlts + 2, iBlock)

  if (minval(NDensityS(iLon, iLat, :, 1:nSpecies, iBlock)) <= 0.0) then
    write(*, *) "negative density found!"
    write(*, *) NDensityS(iLon, iLat, 1, 1:nSpecies, iBlock)
    call stop_gitm("Can't Continue")
  end if

  Heating = EuvHeating(iLon, iLat, :, iBlock)
  Centrifugal = (CosLatitude(iLat, iBlock)*OmegaBodyInput)**2
  Coriolis = 2*CosLatitude(iLat, iBlock)*OmegaBodyInput
  LogRho = log(Rho(iLon, iLat, :, iBlock))
  do iDim = 1, 3
    Vel_GD(:, iDim) = Velocity(iLon, iLat, :, iDim, iBlock)
  end do

  !!!! CHANGE !!!!

  Temp = Temperature(iLon, iLat, :, iBlock)*TempUnit(iLon, iLat, :)
  do iSpecies = 1, nSpecies
    LogNS1(:, iSpecies) = log(NDensityS(iLon, iLat, :, iSpecies, iBlock))
    VertVel(:, iSpecies) = VerticalVelocity(iLon, iLat, :, iSpecies, iBlock)
  end do

  cMax1 = cMax_GDB(iLon, iLat, :, iUp_, iBlock)

  do iDim = 1, 3
    IVel(:, iDim) = IVelocity(iLon, iLat, :, iDim, iBlock)
  end do

  do iSpecies = 1, nIons - 1 !Advect
    if (UseImprovedIonAdVection) then
      LogINS(:, iSpecies) = IDensityS(iLon, iLat, :, iSpecies, iBlock)
    else
      LogINS(:, iSpecies) = log(IDensityS(iLon, iLat, :, iSpecies, iBlock))
    end if
  end do

  MeanMajorMass_1d = MeanMajorMass(iLon, iLat, :)
  gamma_1d = gamma(ilon, ilat, :, iBlock)

!!!!  LogINS  = IDensityS(iLon,iLat,:,1:nIonsAdvect,iBlock)

  Lat = Latitude(iLat, iBlock)*180.0/pi
  Lon = Longitude(iLon, iBlock)*180.0/pi

  ! Cell centered variables
  Gravity_G = Gravity_GB(iLon, iLat, :, iBlock)
  Altitude_G = Altitude_GB(iLon, iLat, :, iBlock)
  dAlt_C = dAlt_GB(iLon, iLat, -1:nAlts + 2, iBlock)
  InvRadialDistance_C = InvRadialDistance_GB(iLon, iLat, -1:nAlts + 2, iBlock)

  ! Face centered variables
  ! This is the distance between cell centers.
  ! Note that face(i) is between cells i and i-1 (like in BATSRUS)
  dAlt_F(0:nAlts + 2) = Altitude_G(0:nAlts + 2) - Altitude_G(-1:nAlts + 1)
  dAlt_F(-1) = dAlt_F(0)
  InvDAlt_F = 1.0/dAlt_F
  SpeciesDensityOld(iLon, iLat, :, 1:nSpeciesTotal, iBlock) = &
    NDensityS(iLon, iLat, :, :, iBlock)
  SpeciesDensityOld(iLon, iLat, :, nSpeciesTotal + 1:nSpeciesAll, iBlock) = &
    IDensityS(iLon, iLat, :, 1:nIons - 1, iBlock)

  if (UseAUSMSolver) then
    call advance_vertical_1d_ausm
  else
    call advance_vertical_1d_rusanov
  end if

  Rho(iLon, iLat, :, iBlock) = exp(LogRho)

  do iDim = 1, 3
    Velocity(iLon, iLat, :, iDim, iBlock) = Vel_GD(:, iDim)
  end do

  Temperature(iLon, iLat, :, iBlock) = Temp/TempUnit(iLon, iLat, :)

  do iSpecies = 1, nSpecies
    LogNS(iLon, iLat, :, iSpecies, iBlock) = LogNS1(:, iSpecies)
    VerticalVelocity(iLon, iLat, :, iSpecies, iBlock) = VertVel(:, iSpecies)
  end do

  do iSpecies = nSpecies + 1, nSpecies
    LogNS(iLon, iLat, :, iSpecies, iBlock) = LogNS1(:, iSpecies)
  end do

  nDensity(iLon, iLat, :, iBlock) = 0.0
  do iSpecies = 1, nSpecies
    nDensityS(iLon, iLat, :, iSpecies, iBlock) = exp(LogNS1(:, iSpecies))
    nDensity(iLon, iLat, :, iBlock) = nDensity(iLon, iLat, :, iBlock) + &
                                      nDensityS(iLon, iLat, :, iSpecies, iBlock)
  end do

  if (UseIonAdvection) then

    do iIon = 1, nIons - 1 !Advect
      if (UseImprovedIonAdvection) then
        IDensityS(iLon, iLat, :, iIon, iBlock) = LogINS(:, iIon)
        ! Put in a floor on ion densities....
        do iAlt = 1, nAlts
          if (IDensityS(iLon, iLat, iAlt, iIon, iBlock) < 1) then
            IDensityS(iLon, iLat, iAlt, iIon, iBlock) = 1.0
          end if
        end do
      else
        IDensityS(iLon, iLat, :, iIon, iBlock) = exp(LogINS(:, iIon))
      end if
    end do

    !\
    ! New Electron Density
    !/
    IDensityS(iLon, iLat, :, ie_, iBlock) = 0.0
    do iIon = 1, nIons - 1
      IDensityS(iLon, iLat, :, ie_, iBlock) = &
        IDensityS(iLon, iLat, :, ie_, iBlock) + &
        IDensityS(iLon, iLat, :, iIon, iBlock)
    end do
  end if

  SpeciesDensity(iLon, iLat, :, 1:nSpeciesTotal, iBlock) = &
    NDensityS(iLon, iLat, :, :, iBlock)
  SpeciesDensity(iLon, iLat, :, nSpeciesTotal + 1:nSpeciesAll, iBlock) = &
    IDensityS(iLon, iLat, :, 1:nIons - 1, iBlock)

end subroutine advance_vertical
