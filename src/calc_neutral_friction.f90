! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

subroutine calc_neutral_friction(DtIn, oVel, EddyCoef_1d, NDensity_1d, NDensityS_1d, &
                                 GradLogCon, Temp)

  use ModGITM
  use ModSources
  use ModPlanet, only: Diff0, DiffExp, IsEarth
  use ModInputs, only: UseNeutralFriction, DoCheckForNans
  use ieee_arithmetic

  implicit none

  real, intent(in) :: DtIn
  real, intent(inout) :: oVel(1:nAlts, 1:nSpecies)
  real, intent(in) :: EddyCoef_1d(1:nAlts)
  real, intent(in) :: NDensity_1d(1:nAlts)
  real, intent(in) :: NDensityS_1d(1:nAlts, 1:nSpecies)
  real, intent(in) :: GradLogCon(1:nAlts, 1:nSpecies)
  real, intent(in) :: Temp(1:nAlts)

  integer :: iSpecies, jSpecies
  real :: CoefMatrix(nSpecies, nSpecies), kTOverM
  real :: Matrix(nSpecies, nSpecies)
  real :: Vel(nSpecies), Parity
  integer :: iPivot(nSpecies)

  real :: TempDij
  real :: InvDij(nSpecies)
  real :: Dij(nSpecies)
  real :: denscale
  real :: mscale
  real :: mms
  real :: mmwos(nSpecies)

  real :: EddyContribution(nSpecies)

  integer :: iAlt

  call report("calc_neutral_friction", 4)

  if (.not. UseNeutralFriction) return

  EddyContribution(1:nSpecies) = 0.0

  do iAlt = 1, nAlts

    Vel = oVel(iAlt, 1:nSpecies)
    CoefMatrix = 0.0

    mms = 0.0
    mmwos = 0.0
    InvDij = 0.0

    EddyContribution(1:nSpecies) = 0.0

    do iSpecies = 1, nSpecies

      InvDij(iSpecies) = 0.0
      kTOverM = Boltzmanns_Constant*Temp(iAlt)/Mass(iSpecies)
      denscale = 1.0/NDensity_1d(iAlt)
      do jSpecies = 1, nSpecies
        if (jSpecies == iSpecies) cycle

        if (DoCheckForNans) then
           if (ieee_is_nan(Temp(iAlt))) &
                write(*, *) "Friction : Temp is nan", iAlt
           if (ieee_is_nan(NDensity_1d(iAlt))) &
                write(*, *) "Friction : NDen is nan", iAlt
           if (ieee_is_nan(NDensityS_1d(iAlt, jSpecies))) &
                write(*, *) "Friction : NDenS is nan", iAlt, jSpecies
        endif

        ! TempDij are the Dij binary coefficients
        ! Based upon the formulation by Banks and Kokarts.
        ! These coefficients demand that
        ! (1) NDensity be in cm^-3 (hence the 1.0e-06) factor below
        ! (2) Additionally, the Dij's are in cm^2/s, thus the 1.0e-04 factor
        TempDij = (1.0e-04)* &              ! Scales the Dij from cm^2/s -> m^2/s
                  (Diff0(iSpecies, jSpecies)*(Temp(iAlt)**DiffExp(iSpecies, jSpecies)))/ &
                  (NDensity_1d(iAlt)*(1.0e-06)) ! Converts to #/cm^-3

        ! Add the eddy coefficient:
        TempDij = TempDij + EddyCoef_1d(iAlt)
        CoefMatrix(iSpecies, jSpecies) = &
          - DtIn * kTOverM * denscale* NDensityS_1d(iAlt, jSpecies) / TempDij

        InvDij(iSpecies) = InvDij(iSpecies) + &
                           denscale*NDensityS_1d(iAlt, jSpecies)/ &
                           (TempDij)

      enddo  ! End DO over jSpecies

      EddyContribution(iSpecies) = &
        -1.0*EddyCoef_1d(iAlt)*GradLogCon(iAlt, iSpecies)

    enddo  !End DO Over iSpecies

    Matrix = CoefMatrix
    do iSpecies = 1, nSpecies
      Matrix(iSpecies, iSpecies) = &
        1.0 - sum(CoefMatrix(iSpecies, :))
    enddo
    call ludcmp(Matrix, nSpecies, nSpecies, iPivot, Parity)
    call lubksb(Matrix, nSpecies, nSpecies, iPivot, Vel)

    oVel(iAlt, 1:nSpecies) = Vel(1:nSpecies) + &
                             EddyContribution(1:nSpecies)
  enddo

end subroutine calc_neutral_friction
