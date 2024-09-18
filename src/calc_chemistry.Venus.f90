! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

subroutine calc_chemistry(iBlock)

!  No override of CO2 densities; explicit calculation

  use ModGITM
  use ModInputs, only: iDebugLevel, UseIonChemistry, UseNeutralChemistry, f107, f107a, &
                       useImplicitChemistry
  use ModConstants
  use ModSources
  use ModChemistry
  use ModGITMImplicit
  use ModTime, only: iStep
  use ModPlanet, only: ialtminiono
  use ModUserGITM

  implicit none

  integer, intent(in) :: iBlock

  real  :: IonSources(nIons), IonLosses(nIons), nSources(nSpeciesTotal), ISources(nIons)
  real  :: NeutralSources(nSpeciesTotal), NeutralLosses(nSpeciesTotal)
  real  :: ChemicalHeatingSub, Emission(nEmissions)

  integer :: iLon, iLat, iAlt, niters, iIon, iNeutral, ivar, nImplicitSourceCalls, iminiono

  real :: dttotal, dtsub, dtMin
  real :: tli(nIons), tsi(nIons), tln(nSpeciesTotal), tsn(nSpeciesTotal)

  real :: EmissionTotals(nEmissions), F

  real :: ionso, ionlo, neuso, neulo

  logical :: doImplicit = .true.
  !---------------------------------------------------------------------------
  UseNeutralConstituent = .true.
  UseIonConstituent = .true.

  nImplicitSourceCalls = nSpeciesAll*2 + 1

  call report("Chemistry", 2)
  call start_timing("calc_chemistry")
  if (iDebugLevel > 3) then
    do iIon = 1, nIons
      write(*, *) "====> start calc_chemistry: Max Ion Density: ", iIon, &
        maxval(IDensityS(1:nLons, 1:nLats, (nAlts*4)/5, iIon, iBlock))
    end do
  end if

  ChemicalHeatingRate = 0.0
  ChemicalHeatingRateIon = 0.0
  ChemicalHeatingRateEle = 0.0
  ChemicalHeatingSpecies = 0.0

  !   if (istep .lt. 10000) then
  !      useimplicitchemistry = .false.
  !   else
  !      useimplicitchemistry = .true.
  !   endif
  !  neutralsourcestotal = 0.0
  !  neutrallossestotal = 0.0

  do iLon = 1, nLons
    do iLat = 1, nLats

      iMinIono = max(iAltMinIono(iLon, iLat, iBlock), 1.0)

      do ialt = iminiono, nAlts

        DtSub = Dt

        Ions = IDensityS(iLon, iLat, iAlt, 1:nIons, iBlock)
        Neutrals = NDensityS(iLon, iLat, iAlt, :, iBlock)

        call calc_reaction_rates(iLon, iLat, iAlt, iBlock)
        call calc_chemical_sources(iLon, iLat, iAlt, iBlock, &
                                   IonSources, IonLosses, &
                                   NeutralSources, NeutralLosses, &
                                   ChemicalHeatingSub, Emission)

        if (.not. UseIonChemistry) then
          IonSources = 0.0
          IonLosses = 0.0
        else
          do iIon = 1, nIons - 1
            if (.not. UseIonConstituent(iIon)) then
              IonSources(iIon) = 0.0
              IonLosses(iIon) = 0.0
            end if
          end do
        end if

        if (.not. UseNeutralChemistry) then
          NeutralSources = 0.0
          NeutralLosses = 0.0
        else
          do iNeutral = 1, nSpeciesTotal
            if (.not. UseNeutralConstituent(iNeutral)) then
              NeutralSources(iNeutral) = 0.0
              NeutralLosses(iNeutral) = 0.0
            end if
          end do
        end if

        ! Take Implicit time step
        Ions(ie_) = 0.0
        do iIon = 1, nIons - 1
          ionso = IonSources(iIon)
          ionlo = IonLosses(iIon)/(Ions(iIon) + 1.0e-6)
          Ions(iIon) = (Ions(iIon) + ionso*DtSub)/ &
                       (1 + DtSub*ionlo)
          ! sum for e-
          Ions(ie_) = Ions(ie_) + Ions(iIon)
        end do

        do iNeutral = 1, nSpeciesTotal

          neuso = NeutralSources(iNeutral)
          neulo = NeutralLosses(iNeutral)/(Neutrals(iNeutral) + 1.0e-6)

          Neutrals(iNeutral) = (Neutrals(iNeutral) + neuso*DtSub)/ &
                               (1 + DtSub*neulo)

        end do

        NDensityS(iLon, iLat, iAlt, :, iBlock) = Neutrals
        IDensityS(iLon, iLat, iAlt, :, iBlock) = Ions

        userdata1d(1, 1, ialt, 25) = 0
        userdata1d(1, 1, ialt, 25) = (NeutralSources(1) - NeutralLosses(1))*DtSub

        ChemicalHeatingRate(iLon, iLat, iAlt) = &
          ChemicalHeatingRate(iLon, iLat, iAlt) + &
          ChemicalHeatingSub*DtSub

        ChemicalHeatingSpecies(iLon, iLat, iAlt, :) = &
          ChemicalHeatingSpecies(iLon, iLat, iAlt, :) + &
          ChemicalHeatingS*DtSub

        EmissionTotals = EmissionTotals + Emission*DtSub

      end do
    end do
  end do

  ChemicalHeatingRate(:, :, :) = &
    ChemicalHeatingRate(:, :, :)*Element_Charge/ &
    TempUnit(1:nLons, 1:nLats, 1:nAlts)/cp(1:nLons, 1:nLats, 1:nAlts, iBlock)/ &
    rho(1:nLons, 1:nLats, 1:nAlts, iBlock)

  ChemicalHeatingRateIon(:, :, :) = &
    ChemicalHeatingRateIon(:, :, :)*Element_Charge

  ChemicalHeatingSpecies = ChemicalHeatingSpecies*Element_Charge

  if (iDebugLevel > 3) then
    do iIon = 1, nIons
      write(*, *) "====> calc_chemistry: Max Ion Density: ", iIon, &
        maxval(IDensityS(1:nLons, 1:nLats, (nAlts*4)/5, iIon, iBlock))
    end do
  end if
  call end_timing("calc_chemistry")

end subroutine calc_chemistry

