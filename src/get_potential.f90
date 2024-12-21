! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

subroutine init_get_potential

  use ModGITM
  use ModTime
  use ModInputs
  use ModIE
  use ModErrors
  use ModElectrodynamics, only: IEModel_

  implicit none

  integer :: iError = 0

  call report("init_get_potential", 2)

  allocate(IEModel_)
  IEModel_ = ieModel()

  call IEModel_%verbose(iDebugLevel)

  
  call IEModel_ % efield_model(cPotentialModel)
  call IEModel_ % aurora_model(cAuroralModel)

  ! Most likely do not need to change, use the default.
  call IEModel_%model_dir("extIE/")

  ! If we are using AMIE files, set north and south files:
  call IEModel_%filename_north(cAMIEFileNorth)
  call IEModel_%filename_south(cAMIEFileSouth)

  ! Initialize the IE library after setting it up:
  call IEModel_%init()

  ! Load in indices for the 0th time-step
  call set_ie_indices(IEModel_, CurrentTime)

  ! Check that correct indices are present:
  call IEModel_%check_indices

  call report("ieModel indices were checked", 5)

  ! If there were errors initializing, stop here
  if (.not. isOk) then
    call set_error("Failed to initialize ieModel! Exiting!")
    call report_errors
    call stop_gitm("Failed to initialize ieModel in get_potential. Check indices, probably.")
  endif

  ! Initialize the grid:
  call IEModel_%nMlts(nLons + 4)
  call IEModel_%nLats(nLats + 4)

  ! if (UseBarriers) call MPI_BARRIER(iCommGITM, iError)
  ! Initialize the IE library after setting it up:
  call IEModel_%init()
  ! then also regional amie... maybe there's a way to do it cleanly

  call report("done with init_get_potential", 2)

  if (UseBarriers) call MPI_BARRIER(iCommGITM, iError)
  if (iError /= 0) call set_error("MPI Barrier falied in init_get_potential")

end subroutine init_get_potential

!--------------------------------------------------------------------
! get_potential
!--------------------------------------------------------------------

subroutine get_potential(iBlock)

  use ModGITM
  use ModTime
  use ModInputs
  use ModUserGITM
  use ModErrors
  use ModElectrodynamics, only: IEModel_
  use ModMpi

  implicit none

  integer, intent(in) :: iBlock

  integer :: iError, iLat, iLon, iAlt, iPot, nPot, iDir, nDir = 1
  logical :: IsFirstTime = .true.
  logical :: IsFirstPotential(nBlocksMax) = .true.
  logical :: IsFirstAurora(nBlocksMax) = .true.
  real    :: mP, dis, TempWeight
  real    ::  LocalSumDiffPot, MeanDiffPot, LatBoundNow

  real, dimension(-1:nLons + 2, -1:nLats + 2) :: TempPotential2d
  real, dimension(-1:nLons + 2, -1:nLats + 2, 2) :: TempPotential, AMIEPotential
  real, dimension(-1:nLons + 2, -1:nLats + 2) :: Grid, dynamo, SubMLats, SubMLons
  real, dimension(-1:nLons + 2, -1:nLats + 2) :: lats, mlts, EFlux
  real :: by, bz, CuspLat, CuspMlt

  logical :: UAl_UseGridBasedEIE

  call start_timing("get_potential")
  call report("get_potential", 2)

  iError = 0

  if (index(cPlanet, "Earth") == 0) then

    potential = 0.0
    ElectronAverageEnergy = 0.1
    ElectronEnergyFlux = 0.0001
    return

  endif

  ! Grid might be reset by calc_electrodynamics
  call IEModel_%nMlts(nLons + 4)
  call IEModel_%nLats(nLats + 4)

  call IEModel_%time_real(tSimulation)

  if (floor((tSimulation - dt)/DtPotential) /= &
      floor((tsimulation)/DtPotential) .or. IsFirstPotential(iBlock)) then
    
    call report("Getting Potential", 1)
    call set_ie_indices(IEModel_, CurrentTime)

    Potential(:, :, :, iBlock) = 0.0

    do iAlt = -1, nAlts + 2

      call iemodel_%grid(MLT(-1:nLons + 2, -1:nLats + 2, iAlt), &
                           MLatitude(-1:nLons + 2, -1:nLats + 2, iAlt, iBlock))

      if (.not. isOk) then
        call set_error("Error in routine get_potential (SetGrid):")
        call report_errors
        call stop_gitm("Stopping in get_potential")
      endif

      if (iDebugLevel > 1 .and. iAlt == 1) &
        write(*, *) "==> Getting IE potential"

      TempPotential = 0.0
      TempPotential2d = 0.0

      call iemodel_%get_potential(TempPotential2d)
      TempPotential(:, :, 1) = TempPotential2d

      if (.not. isOk) then
        call set_error("Error in routine get_potential (getting potential):")
        call report_errors
        call stop_gitm("Stopping in get_potential")
      endif

      nDir = 1

      if (UseDynamo .and. .not. Is1D) then
        dynamo = 0.0
        call get_dynamo_potential( &
          MLongitude(-1:nLons + 2, -1:nLats + 2, iAlt, iBlock), &
          MLatitude(-1:nLons + 2, -1:nLats + 2, iAlt, iBlock), dynamo)

        ! Set latitude boundary between region of high lat convection
        ! and region of neutral wind dyanmo based on if SWMF potential
        ! is being used:
        if (IsFramework) then
          LatBoundNow = 45.
        else
          LatBoundNow = DynamoHighLatBoundary
        endif

        do iDir = 1, nDir
          do iLon = -1, nLons + 2
            do iLat = -1, nLats + 2
              if (abs(MLatitude(iLon, iLat, iAlt, iBlock)) < LatBoundNow) then
                dis = (LatBoundNow - &
                       abs(MLatitude(iLon, iLat, iAlt, iBlock)))/20.0
                if (dis > 1.0) then
                  TempPotential(iLon, iLat, iDir) = dynamo(iLon, iLat)
                else
                  TempPotential(iLon, iLat, iDir) = &
                    (1.0 - dis)*TempPotential(iLon, iLat, iDir) + &
                    dis*dynamo(iLon, iLat)
                endif
              endif
            enddo
          enddo
        enddo

      endif

      Potential(:, :, iAlt, iBlock) = TempPotential(:, :, 1)
      if (UseTwoAMIEPotentials) then
        PotentialY(:, :, iAlt, iBlock) = TempPotential(:, :, 2)
      else
        PotentialY(:, :, iAlt, iBlock) = Potential(:, :, iAlt, iBlock)
      endif

      !----------------------------------------------
      ! Another example of user output

      if (iAlt == 1) then

        UserData2d(1:nLons, 1:nLats, 1, 1, iBlock) = &
          TempPotential(1:nLons, 1:nLats, 1)/1000.0
      endif

    enddo

    IsFirstPotential(iBlock) = .false.

  endif

  ! -----------------------------------------------------
  ! Now get the aurora.
  ! This assumes that the field lines are basically
  ! vertical starting at the top of the model.
  ! -----------------------------------------------------

  if (floor((tSimulation - dt)/DtAurora) /= &
      floor((tsimulation)/DtAurora) .or. IsFirstAurora(iBlock)) then

    call report("Getting Aurora", 1)

    iAlt = nAlts + 1

    call ieModel_%get_aurora(ElectronEnergyFlux, ElectronAverageEnergy)

    if (.not. isOk) then
      call set_error("Error in routine get_potential (getting aurora):")
      call report_errors
      call stop_gitm("Stopping in get_potential")
    endif

    ! Sometimes, in AMIE, things get messed up in the
    ! Average energy, so go through and fix some of these.
    do iLat = -1, nLats + 2
      do iLon = -1, nLons + 2
        if (ElectronAverageEnergy(iLon, iLat) < 0.0) then
          ElectronAverageEnergy(iLon, iLat) = 0.1
          write(*, *) "ave e i,j Negative : ", iLon, iLat, &
            ElectronAverageEnergy(iLon, iLat)
        endif
        if (ElectronAverageEnergy(iLon, iLat) > 100.0) then
          write(*, *) "ave e i,j Positive : ", iLon, iLat, &
            ElectronAverageEnergy(iLon, iLat)
          ElectronAverageEnergy(iLon, iLat) = 0.1
        endif
      enddo
    enddo

    ! -----------------------------------------------------
    ! Get Ion Precipitation if desired
    ! -----------------------------------------------------

    if (UseIonPrecipitation) then
      call stop_gitm("I don't know how to do that yet, sorry!")
      !   call UA_GetIonAveE(IonAverageEnergy, iError)
      !   if (iError /= 0) then
      !     write(*, *) "Error in get_potential (UA_GetAveE):"
      !     write(*, *) iError
      !     IonAverageEnergy = 1.0
      !   endif

      !   call UA_GetIonEFlux(IonEnergyFlux, iError)
      !   if (iError /= 0) then
      !     write(*, *) "Error in get_potential (UA_GetEFlux):"
      !     write(*, *) iError
      !     IonEnergyFlux = 0.1
      !   endif

      ! endif

    endif

    if (UseCusp) then

      lats = abs(MLatitude(-1:nLons + 2, -1:nLats + 2, iAlt, iBlock))

      if (maxval(lats) > 50) then

        mlts = mod(MLT(-1:nLons + 2, -1:nLats + 2, iAlt) + 24.0, 24.0)

        by = iemodel_%needImfBy
        bz = iemodel_%needImfBz

        ! If we are in the southern hemisphere, reverse by:
        if (lats(nLons/2, nLats/2) < 0.0) by = -by

        if (bz > 0) then
          ! Newell et al., 1988:
          CuspLat = 77.2 + 0.11*bz
          ! Asai et al., Earth Planets Space, 2005:
          CuspMlt = 11.755 + 0.169*by
        else
          ! Asai et al., Earth Planets Space, 2005:
          CuspMlt = 11.949 + 0.0826*by
          ! Zhang et al., JGR, 2005:
          if (Bz > -10) then
            CuspLat = 77.2 + 1.1*bz
          else
            CuspLat = 21.7*exp(0.1*bz) + 58.2
          endif
        endif

        EFlux = CuspEFlux* &
                exp(-abs(lats - CuspLat)/CuspLatHalfWidth)* &
                exp(-abs(mlts - CuspMlt)/CuspMltHalfWidth)

        do iLat = -1, nLats + 2
          do iLon = -1, nLons + 2
            if (EFlux(iLon, iLat) > 0.1) then
              ElectronEnergyFlux(iLon, iLat) = EFlux(iLon, iLat)
              ElectronAverageEnergy(iLon, iLat) = CuspAveE
            endif
          enddo
        enddo

      endif

    endif

    if (iDebugLevel > 2) &
      write(*, *) "==> Max, electron_ave_ene : ", &
      maxval(ElectronAverageEnergy), &
      maxval(ElectronEnergyFlux)

    IsFirstAurora(iBlock) = .false.

  endif

  if (iDebugLevel > 1) &
    write(*, *) "==> Min, Max, CPC Potential : ", &
    minval(Potential(:, :, :, iBlock))/1000.0, &
    maxval(Potential(:, :, :, iBlock))/1000.0, &
    (maxval(Potential(:, :, :, iBlock)) - minval(Potential(:, :, :, iBlock)))/1000.0

  call end_timing("get_potential")

end subroutine get_potential

subroutine get_dynamo_potential(lons, lats, pot)

  use ModGITM
  use ModInputs, only: iDebugLevel, DynamoHighLatBoundary
  use ModElectrodynamics

  implicit none

  real, dimension(-1:nLons + 2, -1:nLats + 2), intent(in) :: lons
  real, dimension(-1:nLons + 2, -1:nLats + 2), intent(in) :: lats
  real, dimension(-1:nLons + 2, -1:nLats + 2), intent(out) :: pot

  integer :: iLon, iLat, iL, iM
  real    :: dM, dL, LatIn, LonIn, iError

  logical :: IsFound

  pot = 0.0

  IsFound = .false.

  do iLon = -1, nLons + 2
    do iLat = -1, nLats + 2

      LatIn = lats(iLon, iLat)
      LonIn = mod(lons(iLon, iLat) + 360.0, 360.0)

      IsFound = .false.

      if (abs(LatIn) <= MagLatMC(nMagLons, nMagLats)) then

        iM = 1
        do while (iM < nMagLons + 1)
          iL = 1
          do while (iL < nMagLats)

            !\
            ! Check to see if the point is within the current cell
            !/

            if (LatIn < MagLatMC(iM, iL + 1) .and. &
                LatIn >= MagLatMC(iM, iL) .and. &
                LonIn < MagLonMC(iM + 1, iL) .and. &
                LonIn >= MagLonMC(iM, iL)) then

              dM = (LonIn - MagLonMC(iM, iL))/ &
                   (MagLonMC(iM + 1, iL) - MagLonMC(iM, iL))

              dL = (LatIn - MagLatMC(iM, iL))/ &
                   (MagLatMC(iM, iL + 1) - MagLatMC(iM, iL))

              pot(iLon, iLat) = &
                (1.0 - dM)*(1.0 - dL)*DynamoPotentialMC(iM, iL) + &
                (1.0 - dM)*(dL)*DynamoPotentialMC(iM, iL + 1) + &
                (dM)*(dL)*DynamoPotentialMC(iM + 1, iL + 1) + &
                (dM)*(1.0 - dL)*DynamoPotentialMC(iM + 1, iL)

              iL = nMagLats
              iM = nMagLons

              IsFound = .true.

            endif

            iL = iL + 1

          enddo

          iM = iM + 1

        enddo

        ! Check for near pole
        if (.not. IsFound) then

          if (LatIn > 88.0) then
            IsFound = .true.
            pot(iLon, iLat) = sum(DynamoPotentialMC(:, nMagLats))/(nMagLons + 1)
          endif

          if (LatIn < -88.0) then
            IsFound = .true.
            pot(iLon, iLat) = sum(DynamoPotentialMC(:, 0))/(nMagLons + 1)
          endif

          write(*, *) "Inside the low latitude, but can't find the point!"
          write(*, *) LatIn, LonIn

        endif

      endif

      if (.not. IsFound) then
        if (abs(LatIn) < MagLatMC(nMagLons, nMagLats)) &
          write(*, *) "=====> Could not find point : ", &
          LatIn, LonIn, DynamoHighLatBoundary, &
          MagLatMC(nMagLons, nMagLats)
        pot(iLon, iLat) = 0.0
      endif

    enddo

  enddo

end subroutine get_dynamo_potential

! Subroutine to automatically set indices for the given IEModel
! Will need updateing when HPn/s are introduced
subroutine set_ie_indices(IEModel_, TimeIn)

  use ModKind
  use ModIndicesInterfaces
  use ModIE, only: ieModel
  ! use ModElectrodynamics, only: ieModel_
  use ModErrors
  use ModTime, only: EndTime !could pull current time too, but better to be explicit
  use ModInputs, only: TimeDelayHighLat

  implicit none

  type(ieModel), intent(inout) :: IEModel_
  real, intent(in) :: TimeIn

  integer :: iError = 0
  real    :: val


  call IEModel_%time_real(TimeIn)

  if (IEModel_%doReadMHD) then

    call read_MHDIMF_Indices_new( &
      iError, &
      TimeIn + TimeDelayHighLat, &
      EndTime + TimeDelayHighLat)

    if (iError /= 0) call set_error("Issue reading IMF file in get_potential")

    call get_IMF_Bz(TimeIn + TimeDelayHighLat, val, iError)
    if (val < -50.0) val = -50.0
    if (val > 50.0) val = 50.0
    call IEModel_%imfBz(val)

    call get_IMF_By(TimeIn, val, iError)
    if (val < -50.0) val = -50.0
    if (val > 50.0) val = 50.0
    call IEModel_%imfBy(val)

    call get_SW_V(TimeIn, val, iError)
    if (val < -1800.0) val = -1800.0
    if (val > 1800.0) val = 1800.0
    call IEModel_%swV(val)

    call get_SW_N(TimeIn, val, iError)
    if (val > 80) val = 80
    call IEModel_%swN(val)

    if (iError /= 0 .or. .not. isOk) then
      call set_error("IMF values could not be set!")
      return
    endif
  endif

  if (IEModel_%doReadSME) then
    ! read_sme works differently, it reads the entire file, always. So does not
    ! need to be called again.

    call get_AU(TimeIn, val, iError)
    call IEModel_%au(val)
    call get_AL(TimeIn, val, iError)
    call IEModel_%al(val)

    if (iError /= 0 .or. .not. isOk) then
      call set_error("In set_ie_indices SME values could not be set!")
      return
    endif
  endif

  if (IEModel_%doReadHPI) then

    if (ieModel_%useAeForHp) then
      call read_NOAAHPI_Indices_new(iError, &
      TimeIn + TimeDelayHighLat, &
                                    EndTime + TimeDelayHighLat)

      if (iError /= 0) call set_error("HPI values could not be read.")

    endif

    call get_HPI(TimeIn, val, iError)
    call IEModel_%hp(val)
    call IEModel_%hpN(val)
    call IEModel_%hpS(val)

    if (iError /= 0 .or. .not. isOk) then
      call set_error("HPI values could not be set!")
      return
    endif
  endif

  if (IEModel_%doReadKP) then
    ! I don't think there's a way to set KP unless you use a single value. No time:
    call get_KP(val, iError)
    call IEModel_%kp(val)

    if (iError /= 0 .or. .not. isOk) then
       call set_error("KP values could not be set!")
       return
    endif
  endif

end subroutine set_ie_indices
