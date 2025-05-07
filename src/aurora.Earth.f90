! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

subroutine aurora(iBlock)

  use ModGITM
  use ModSources
  use ModTime, only: tSimulation, CurrentTime
  use ModInputs
  use ModConstants
  use ModUserGITM
  use ModErrors
  use ModIE
  use ModMpi
  use ModIndicesInterfaces
  use ModElectrodynamics, only: IEModel_

  implicit none

  integer, intent(in) :: iBlock

  real :: ion_av_kev, ion_eflx_ergs, ion_eflux, ion_avee
  real :: factor, p, Q0
  integer :: i, j, k, n, iError, iED, iErr, iEnergy
  logical :: IsDone, IsTop, HasSomeAurora!, UseMono, UseWave
  real, dimension(ED_N_Energies) :: &
    e_diffuse_ED_flux, i_diffuse_ED_flux, mono_ED_flux, wave_ED_flux, &
    ED_Flux ! for temp values

  real :: hpi, hpi_NH, hpi_SH

  real, dimension(nLons, nLats, nAlts) :: temp, AuroralBulkIonRate, &
                                          IonPrecipitationBulkIonRate, IonPrecipitationHeatingRate

  logical :: IsFirstTime(nBlocksMax) = .true.

  real :: f1, f2, f3, f4, f5, power
  real :: de1, de2, de3, de4, de5, detotal, h

  real :: LocalVar, HPn, HPs, avepower, ratio, ratio_NH, ratio_SH

  if (UseFangEnergyDeposition .and. IsFirstTime(iBlock)) then
    call initialize_fang_arrays
  else
    if (floor((tSimulation - dT)/dTAurora) == &
        floor(tSimulation/dTAurora)) return
  endif

  AuroralBulkIonRate = 0.0
  ! AuroralHeatingRate(:, :, :, iBlock) = 0.0
  AuroralIonRateS = 0.0
  HPn = 0
  HPs = 0

  call report("Aurora", 1)
  call start_timing("Aurora")

  if (iBlock == 1) then
    HemisphericPowerNorth = 0.0
    HemisphericPowerSouth = 0.0
  endif

  if (NormalizeAuroraToHP) then

    ! Let's scale our hemispheric power so it is roughly the same as what
    ! is measured.

    ! Calculate HP from e- flux
    call calculate_HP(iBlock, HPn, HPs)

    ! Collect all of the powers by summing them together

    LocalVar = HemisphericPowerNorth/1.0e9
    call MPI_REDUCE(LocalVar, HPn, 1, MPI_REAL, MPI_SUM, &
                    0, iCommGITM, iError)

    LocalVar = HemisphericPowerSouth/1.0e9
    call MPI_REDUCE(LocalVar, HPs, 1, MPI_REAL, MPI_SUM, &
                    0, iCommGITM, iError)

    ! Average north and south together
    avepower = (HPn + HPs)/2.0

    ! If we are only have one hemisphere or the other, assign to avepower
    if (HPs < 0.1*HPn) avepower = HPn
    if (HPn < 0.1*HPs) avepower = HPs

    call MPI_Bcast(avepower, 1, MPI_Real, 0, iCommGITM, ierror)
    call MPI_Bcast(Hps, 1, MPI_Real, 0, iCommGITM, ierror)
    call MPI_Bcast(HPn, 1, MPI_Real, 0, iCommGITM, ierror)

    if (DoSeparateHPI) then
      call get_hpi_n(CurrentTime, hpi_NH, iError)
      call get_hpi_s(CurrentTime, hpi_SH, iError)

      if (Hps == 0) call set_error("HPs is zero!")
      if (Hpn == 0) call set_error("HPn is zero!")

      ratio_NH = Hpi_NH/HPn
      ratio_SH = Hpi_SH/HPs

    else

      call get_hpi(CurrentTime, Hpi, iError)
      if (avepower == 0) then
        call set_error("AvePower is zero!")
        avepower = 1.0
      endif
      ratio_NH = Hpi/avepower
      ratio_SH = Hpi/avepower

    endif

    if (iDebugLevel >= 0) then
      if ((iDebugLevel == 0) .and. IsFirstTime(iBlock) .and. .not. doSeparateHPI) then
        write(*, *) '---------------------------------------------------'
        write(*, *) 'Using auroral normalizing ratios!!! '
        write(*, *) 'no longer reporting!'
        write(*, *) '---------------------------------------------------'
      elseif ((iDebugLevel == 0) .and. IsFirstTime(iBlock) .and. doSeparateHPI) then
        write(*, *) '---------------------------------------------------'
        write(*, *) 'Using HPI from each hemisphere to normalize aurora!!'
        write(*, *) 'no longer reporting!'
        write(*, *) '---------------------------------------------------'
      endif
      if (iDebugLevel >= 2) then
        if (doSeparateHPI) then
          write(*, *) 'Auroral normalizing ratios: '
          write(*, *) 'Hpi(NH)  HPI(SH)  HPI(NH-modeled)  HPI(SH-modeled)  ratio_n  ratio_s'
          write(*, *) Hpi_NH, Hpi_SH, HPn, HPs, ratio_NH, ratio_sh
        else
          write(*, *) 'Auroral normalizing ratio: ', Hpi_NH, ratio_SH, ratio_NH
        endif
      endif

    endif
    do i = 1, nLats
      do j = 1, nLons
        if (ElectronEnergyFluxDiffuse(j, i) > 0.1) then
          if (latitude(i, iBlock) < 0.0) then
            ElectronEnergyFluxDiffuse(j, i) = ElectronEnergyFluxDiffuse(j, i)*ratio_sh
          else
            ElectronEnergyFluxDiffuse(j, i) = ElectronEnergyFluxDiffuse(j, i)*ratio_NH
          endif
        endif
      enddo
    enddo

  endif

  ! Reset the hemispheric power

  if (iBlock == 1) then
    HemisphericPowerNorth = 0.0
    HemisphericPowerSouth = 0.0
  endif

  if (iProc == 0 .and. AveEFactor /= 1.0 .and. IsFirstTime(iBlock)) then
    write(*, *) "Auroral Experiments!!!!"
    write(*, *) "AveEFactor : ", AveEFactor
  endif
  if (iProc == 0 .and. IsKappaAurora .and. IsFirstTime(iBlock)) then
    write(*, *) "Auroral Experiments!!!!"
    write(*, *) "kappa : ", AuroraKappa
  endif

  do i = 1, nLats
    do j = 1, nLons

      UserData2d(j, i, 1, 2:nUserOutputs, iBlock) = 0.0

      if (UseIonAurora) then
        ion_eflx_ergs = IonEnergyFlux(j, i)
        ion_av_kev = IonAverageEnergy(j, i)
      else
        ion_eflx_ergs = 0.001
        ion_av_kev = 10.0
      endif

      e_diffuse_ED_flux = 0.0
      i_diffuse_ED_flux = 0.0
      wave_ED_flux = 0.0
      mono_ED_flux = 0.0

      HasSomeAurora = .false.

      ! For diffuse auroral models (default)
      if (ElectronEnergyFluxDiffuse(j, i) > 0.1 &
          .and. ElectronAverageEnergyDiffuse(j, i) > 0.1 &
          .and. UseDiffuseAurora) then
        call do_diffuse_aurora(ElectronEnergyFluxDiffuse(j, i), &
                               ElectronAverageEnergyDiffuse(j, i)*AveEFactor, &
                               e_diffuse_ED_flux)
        HasSomeAurora = .true.
      endif
      if (ion_eflx_ergs > 0.1) then
        call do_diffuse_aurora(ion_av_kev, ion_eflx_ergs, i_diffuse_ED_flux)
        HasSomeAurora = .true.
      endif

      ! Monoenergetic aurora
      if (UseMonoAurora .and. &
          ElectronAverageEnergyMono(j, i) > 1.0e4 .and. &
          ElectronEnergyFluxMono(j, i) > 0.1) then

        call do_mono_aurora(ElectronEnergyFluxMono(j, i), &
                            ElectronAverageEnergyMono(j, i), &
                            mono_ED_flux)
        HasSomeAurora = .true.

      endif

      ! Wave (broadband) aurora
      if (UseWaveAurora .and. &
          ElectronEnergyFluxWave(j, i) > 0.1) then
        call do_wave_aurora(ElectronEnergyFluxWave(j, i), &
                            ElectronAverageEnergyWave(j, i), &
                            wave_ED_flux)
        HasSomeAurora = .true.
      endif

      ED_EnergyFlux = e_diffuse_ED_flux + i_diffuse_ED_flux + wave_ED_flux + mono_ED_flux

      if (HasSomeAurora) then

        call calc_fang_rates(j, i, iBlock, AuroralBulkIonRate)

      endif

    enddo
  enddo

  ! From Rees's book:

  temp = 0.92*NDensityS(1:nLons, 1:nLats, 1:nAlts, iN2_, iBlock) + &
         1.00*NDensityS(1:nLons, 1:nLats, 1:nAlts, iO2_, iBlock) + &
         0.56*NDensityS(1:nLons, 1:nLats, 1:nAlts, iO_3P_, iBlock)

  AuroralIonRateS(:, :, :, iO_3P_, iBlock) = &
    0.56*AuroralBulkIonRate* &
    NDensityS(1:nLons, 1:nLats, 1:nAlts, iO_3P_, iBlock)/temp
  AuroralIonRateS(:, :, :, iO2_, iBlock) = &
    1.00*AuroralBulkIonRate* &
    NDensityS(1:nLons, 1:nLats, 1:nAlts, iO2_, iBlock)/temp
  AuroralIonRateS(:, :, :, iN2_, iBlock) = &
    0.92*AuroralBulkIonRate* &
    NDensityS(1:nLons, 1:nLats, 1:nAlts, iN2_, iBlock)/temp

  ! if (UseAuroralHeating) then
  !   AuroralHeating = AuroralHeatingRate(:, :, :, iBlock)/ &
  !                    TempUnit(1:nLons, 1:nLats, 1:nAlts)/cp(:, :, 1:nAlts, iBlock)/ &
  !                    rho(1:nLons, 1:nLats, 1:nAlts, iBlock)
  ! else
  ! AuroralHeating = 0.0
  ! endif
  ! write(*,*) AuroralHeating

  IsFirstTime(iBlock) = .false.

  call end_timing("Aurora")

contains

  ! --------------------------
  ! Diffuse Aurora can be represented by kappa or maxwellian
  ! - This is for Newell diffuse, or all other auroral models!
  ! --------------------------
  subroutine do_diffuse_aurora(eflux_ergs, av_kev, diff_ED_Energy_Flux)
    real, intent(in) :: av_kev, eflux_ergs
    real, intent(inout), dimension(:) :: diff_ED_Energy_Flux
    real :: avee, eflux

    avee = av_kev*1000.0        ! keV -> eV
    eflux = eflux_ergs*6.242e11  ! ergs/cm2/s -> eV/cm2/s

    ! 100 * 100 is for (eV/cm2/s -> J/m2/s)
    power = (eflux)*Element_Charge*100.0*100.0* &
            dLatDist_FB(j, i, nAlts, iBlock)* &
            dLonDist_FB(j, i, nAlts, iBlock)

    if ((latitude(i, iBlock) < 0.0)) then
      HemisphericPowerSouth = HemisphericPowerSouth + power
    else
      HemisphericPowerNorth = HemisphericPowerNorth + power
    endif

    if (IsKappaAurora) then
      call calc_kappa(AuroraKappa, eflux, avee, diff_ED_Energy_Flux)

    else
      ! This calls the Maxwellian from Fang et al. [2010]
      call calc_maxwellian(eflux, avee, diff_ED_Energy_Flux)
    endif

  end subroutine do_diffuse_aurora

  subroutine calc_kappa(kappaVal, ieflux, iavee, kap_ED_EnergyFlux)
    ! This is a Kappa Function from Fang et al. [2010]:

    real, intent(in) :: kappaVal
    real, intent(in) :: ieflux, iavee
    real, intent(inout), dimension(:) :: kap_ED_EnergyFlux

    do n = 1, ED_N_Energies
      ED_Flux(n) = ieflux/2/(iavee/2)**3* & ! a=Q0/2/E0**3
                   (AuroraKappa - 1)*(AuroraKappa - 2)/ &
                   (AuroraKappa**2)* &
                   ed_energies(n)* &
                   (1 + ed_energies(n)/(AuroraKappa*(iavee/2)))** &
                   (-AuroraKappa - 1)

      kap_ED_EnergyFlux(n) = &
        ED_flux(n)* &
        ED_Energies(n)* &
        ED_delta_energy(n)
    enddo

  end subroutine calc_kappa

  subroutine calc_maxwellian(total_flux, avg_e, max_ed_flux)
    ! From (Fang et al., [2010]), a
    ! Maxwellian is defined as:
    ! DifferentialNumberFlux = Q0/2/E0**3 * E * exp(-E/E0),
    ! where:
    ! Q0 = Total Energy Flux
    ! E0 = Characteristic Energy (0.5*avee)
    ! E = mid-point of energy bin
    !
    ! Calculate as a * E * exp(-E/E0)

    real, intent(in) :: total_flux, avg_e
    real, intent(out), dimension(:) :: max_ed_flux

    real :: a, E0

    E0 = (0.5*avg_e)
    a = (total_flux)/2/E0**3

    do n = 1, ED_N_Energies
      max_ed_flux(n) = a*ED_Energies(n)*exp(-ED_Energies(n)/E0)
      ! Maxwellian_Energy-Deposition_flux
      max_ed_flux(n) = &
        max_ed_flux(n)* &
        ED_Energies(n)* &
        ED_delta_energy(n)
    enddo

  end subroutine calc_maxwellian

  ! --------------------------
  ! Monoenergetic aurora
  ! --------------------------
  subroutine do_mono_aurora(eflux_ergs, av_kev, ED_MonoEnergyFlux)

    real, intent(in) :: eflux_ergs, av_kev
    real, intent(inout), dimension(:) :: ED_MonoEnergyFlux
    real :: avee, eflux

    avee = av_kev*1000.0        ! keV -> eV
    eflux = eflux_ergs*6.242e11  ! ergs/cm2/s -> eV/cm2/s

    power = eflux*Element_Charge*100.0*100.0* &    ! (eV/cm2/s -> J/m2/s)
            dLatDist_FB(j, i, nAlts, iBlock)* &
            dLonDist_FB(j, i, nAlts, iBlock)

    if (latitude(i, iBlock) < 0.0) then
      HemisphericPowerSouth = HemisphericPowerSouth + power
    else
      HemisphericPowerNorth = HemisphericPowerNorth + power
    endif

    call calc_mono(eflux, av_kev, ED_MonoEnergyFlux)

  end subroutine do_mono_aurora

  subroutine calc_mono(ENumberFluxMono, mono_av_kev, ED_MonoEnergyFlux)
    real, intent(in) :: ENumberFluxMono, mono_av_kev
    real, intent(inout), dimension(:) :: ED_MonoEnergyFlux
    ! Mono-Energetic goes into one bin only!
    do n = 2, ED_N_Energies - 1
      if (mono_av_kev < ED_energies(n - 1) .and. mono_av_kev >= ED_energies(n)) then
        ED_flux(n) = ED_Flux(n) + &
                     ENumberFluxMono/ &
                     (ED_Energies(n - 1) - ED_Energies(n))
        ED_MonoEnergyFlux(n) = &
          ED_flux(n)* &
          ED_Energies(n)* &
          ED_delta_energy(n)
      endif
    enddo
  end subroutine calc_mono

  ! --------------------------
  ! Wave aurora
  ! --------------------------
  subroutine do_wave_aurora(eflux_ergs, av_kev, wave_EDEnergyFlux)
    real, intent(in) :: eflux_ergs, av_kev
    real, intent(inout), dimension(:) :: wave_EDEnergyFlux
    real :: avee, eflux

    avee = av_kev*1000.0        ! keV -> eV
    eflux = eflux_ergs*6.242e11  ! ergs/cm2/s -> eV/cm2/s

    if (latitude(i, iBlock) < 0.0) then
      HemisphericPowerSouth = HemisphericPowerSouth + power
    else
      HemisphericPowerNorth = HemisphericPowerNorth + power
    endif

    ! #TODO: eventually use a gaussian.

    ! Leaving old code in just in case...
    ! Using the same maxwellian as the diffuse aurora.
    if (eflux > 0.1 .and. avee > 0.1) &
    call calc_maxwellian(eflux, avee, wave_EDEnergyFlux)

    ! k = 0
    ! do n = 3, ED_N_Energies - 3
    !   if (wave_av_kev < ED_energies(n - 1) .and. wave_av_kev >= ED_energies(n)) then
    !     k = n
    !   endif
    ! enddo
    ! if (k > 3) then
    !   f1 = 1.0
    !   f2 = 1.2
    !   f3 = 1.3
    !   f4 = f2
    !   f5 = f1
    !   de1 = ED_energies(k - 3) - ED_energies(k - 2)
    !   de2 = ED_energies(k - 2) - ED_energies(k - 1)
    !   de3 = ED_energies(k - 1) - ED_energies(k)
    !   de4 = ED_energies(k) - ED_energies(k + 1)
    !   de5 = ED_energies(k + 1) - ED_energies(k + 2)
    !   detotal = (f1 + f2 + f3 + f4 + f5)

    !   ED_flux(k - 2) = f1*ElectronNumberFluxWave(j, i)/detotal/de1
    !   ED_flux(k - 1) = f2*ElectronNumberFluxWave(j, i)/detotal/de2
    !   ED_flux(k) = f3*ElectronNumberFluxWave(j, i)/detotal/de3
    !   ED_flux(k + 1) = f4*ElectronNumberFluxWave(j, i)/detotal/de4
    !   ED_flux(k + 2) = f5*ElectronNumberFluxWave(j, i)/detotal/de5

    !   do n = k - 2, n + 2
    !     wave_EDEnergyFlux(n) = ED_flux(n)*ED_Energies(n)*ED_delta_energy(n)
    !   enddo
    ! endif
  end subroutine do_wave_aurora

end subroutine aurora

! --------------------------------
! Initialize the variables necessary to run Fang ED code
! --------------------------------

subroutine initialize_fang_arrays
! Allocate variables necessary for Fang Energy Deposition

  use ModInputs
  use ModSources

  allocate(Fang_Ci(ED_N_Energies, 8), stat=iErr)
  allocate(Fang_y(ED_N_Energies, nAlts), stat=iErr)
  allocate(Fang_f(ED_N_Energies, nAlts), stat=iErr)
  allocate(Fang_Pij(8, 4), stat=iErr)

  ! Ions
  if (UseIonAurora) then
    allocate(Fang_Ion_Ci(ED_N_Energies, 12), stat=iErr)
    allocate(Fang_Ion_y(ED_N_Energies, nAlts), stat=iErr)
    allocate(Fang_Ion_f(ED_N_Energies, nAlts), stat=iErr)
    allocate(Fang_Ion_Pij(12, 4), stat=iErr)
  endif

  if (iErr /= 0) then
    call stop_gitm("Error allocating Fang arrays in aurora")
  endif

  !electrons
  Fang_Pij(1, :) = (/1.25E+00, 1.45903, -2.42E-01, 5.95E-02/)
  Fang_Pij(2, :) = (/2.24E+00, -4.23E-07, 1.36E-02, 2.53E-03/)
  Fang_Pij(3, :) = (/1.42E+00, 1.45E-01, 1.70E-02, 6.40E-04/)
  Fang_Pij(4, :) = (/0.248775, -1.51E-01, 6.31E-09, 1.24E-03/)
  Fang_Pij(5, :) = (/-0.465119, -1.05E-01, -8.96E-02, 1.22E-02/)
  Fang_Pij(6, :) = (/3.86E-01, 1.75E-03, -7.43E-04, 4.61E-04/)
  Fang_Pij(7, :) = (/-6.45E-01, 8.50E-04, -4.29E-02, -2.99E-03/)
  Fang_Pij(8, :) = (/9.49E-01, 1.97E-01, -2.51E-03, -2.07E-03/)

  ! This is from:
  ! Fang, X., D. Lummerzheim, and C. H. Jackman (2013),
  !           Proton impact ionization and a fast calculation method,
  !           J. Geophys. Res. Space Physics, 118, 5369–5378,
  !           doi:10.1002/jgra.50484:

  !ions
  if (UseIonAurora) then
    Fang_Ion_Pij(1, :) = (/2.55050E+00, 2.69476e-01, -2.58425E-01, 4.43190E-02/)
    Fang_Ion_Pij(2, :) = (/6.39287E-01, -1.85817e-01, -3.15636E-02, 1.01370E-02/)
    Fang_Ion_Pij(3, :) = (/1.63996E+00, 2.43580e-01, 4.29873E-02, 3.77803E-02/)
    Fang_Ion_Pij(4, :) = (/-2.13479E-01, 1.42464e-01, 1.55840E-02, 1.97407E-03/)
    Fang_Ion_Pij(5, :) = (/-1.65764E-01, 3.39654e-01, -9.87971E-03, 4.02411E-03/)
    Fang_Ion_Pij(6, :) = (/-3.59358E-02, 2.50330e-02, -3.29365E-02, 5.08057E-03/)
    Fang_Ion_Pij(7, :) = (/-6.26528E-01, 1.46865e+00, 2.51853E-01, -4.57132E-02/)
    Fang_Ion_Pij(8, :) = (/1.01384E+00, 5.94301e-02, -3.27839E-02, 3.42688E-03/)
    Fang_Ion_Pij(9, :) = (/-1.29454E-06, -1.43623e-01, 2.82583E-01, 8.29809E-02/)
    Fang_Ion_Pij(10, :) = (/-1.18622E-01, 1.79191e-01, 6.49171E-02, -3.99715E-03/)
    Fang_Ion_Pij(11, :) = (/2.94890E+00, -5.75821e-01, 2.48563E-02, 8.31078E-02/)
    Fang_Ion_Pij(12, :) = (/-1.89515E-01, 3.53452e-02, 7.77964E-02, -4.06034E-03/)
  endif

  ! Electrons
  do iEnergy = 1, ED_N_Energies
    do i = 1, 8
      Fang_Ci(iEnergy, i) = 0.0
      do j = 0, 3
        Fang_Ci(iEnergy, i) = Fang_Ci(iEnergy, i) + &
                              Fang_Pij(i, j + 1)*log(ED_Energies(iEnergy)/1000.0)**j
      enddo
    enddo
  enddo
  Fang_Ci = exp(Fang_Ci)

  ! Ions
  if (UseIonAurora) then
    do iEnergy = 1, ED_N_Energies
      do i = 1, 12
        Fang_Ion_Ci(iEnergy, i) = 0.0
        do j = 0, 3
          Fang_Ion_Ci(iEnergy, i) = Fang_Ion_Ci(iEnergy, i) + &
                                    Fang_Ion_Pij(i, j + 1)*log(ED_Energies(iEnergy)/1000.0)**j
        enddo
      enddo
    enddo
    Fang_Ion_Ci = exp(Fang_Ion_Ci)
  endif

end subroutine initialize_fang_arrays

! -----------------------------------------
! Calculate energy deposition rates
! -----------------------------------------

subroutine calc_fang_rates(j, i, iBlock, AuroralBulkIonRate)
  ! Given (j,i), on an iBlock, calculate AuroralBulkIonRate

  use ModInputs
  use ModSources
  use ModGITM

  integer, intent(in) :: i, j, iBlock
  real, dimension(nLons, nLats, nAlts), intent(out) :: AuroralBulkIonRate
  real :: BulkScaleHeight1d(nAlts)
  real :: fac(nAlts)
  real :: Ci(8) ! e-
  real :: Ion_Ci(12) ! ions

  real :: Fang_de = 0.035

  BulkScaleHeight1d = &
    Temperature(j, i, 1:nAlts, iBlock) &
    *TempUnit(j, i, 1:nAlts)*Boltzmanns_Constant &
    /(-Gravity_GB(j, i, 1:nAlts, iBlock)* &
      MeanMajorMass(j, i, 1:nAlts))*100.0 ! Convert to cm

  do iEnergy = 1, ED_N_Energies

    ! /10.0 in this statement is for kg/m2 to g/cm2
    ! /1000.0 is conversion from eV to keV
    ! Fang doesn't include the dip angle, be we do.
    Fang_y(iEnergy, :) = 2.0/(ED_Energies(iEnergy)/1000.0)* &
                         (ColumnIntegralRho(j, i, 1:nAlts)/10.0/6e-6)**0.7
    !sinDipAngle(j,i,1:nAlts,iBlock) / 6e-6) ** 0.7

    Ci = Fang_Ci(iEnergy, :)
    Fang_f(iEnergy, :) = &
      Ci(1)*Fang_y(iEnergy, :)**Ci(2)* &
      exp(-Ci(3)*Fang_y(iEnergy, :)**Ci(4)) + &
      Ci(5)*Fang_y(iEnergy, :)**Ci(6)* &
      exp(-Ci(7)*Fang_y(iEnergy, :)**Ci(8))

    ! Energy flux is in eV/cm2/s and Fang needs keV/cm2/s:
    fac = ED_energyflux(iEnergy)/1000.0/ &
          Fang_de/ &
          BulkScaleHeight1d

    ! I think that the 1e6 is cm3 to m3
    AuroralBulkIonRate(j, i, 1:nAlts) = &
      AuroralBulkIonRate(j, i, 1:nAlts) + 1e6*Fang_f(iEnergy, :)*fac

    if (UseIonAurora) then

      ! /10.0 in this statement is for kg/m2 to g/cm2
      ! /1000.0 is conversion from eV to keV
      Fang_Ion_y(iEnergy, :) = 7.5/(ED_Energies(iEnergy)/1000.0)* &
                               (ColumnIntegralRho(j, i, 1:nAlts)/10.0/1e-4)**0.9

      Ion_Ci = Fang_Ion_Ci(iEnergy, :)

      Fang_Ion_f(iEnergy, :) = &
        Ion_Ci(1)*Fang_Ion_y(iEnergy, :)**Ion_Ci(2)* &
        exp(-Ion_Ci(3)*Fang_Ion_y(iEnergy, :)**Ion_Ci(4)) + &
        Ion_Ci(5)*Fang_Ion_y(iEnergy, :)**Ion_Ci(6)* &
        exp(-Ion_Ci(7)*Fang_Ion_y(iEnergy, :)**Ion_Ci(8)) + &
        Ion_Ci(9)*Fang_Ion_y(iEnergy, :)**Ion_Ci(10)* &
        exp(-Ion_Ci(11)*Fang_Ion_y(iEnergy, :)**Ion_Ci(12))

      fac = ED_ion_energyflux(iEnergy)/1000.0/ &
            Fang_de/ &
            BulkScaleHeight1d

      AuroralBulkIonRate(j, i, 1:nAlts) = &
        AuroralBulkIonRate(j, i, 1:nAlts) + 1e6*Fang_Ion_f(iEnergy, :)*fac

    endif

  enddo

  ! AuroralHeatingRate(j, i, 1:nAlts, iBlock) = 0.0 ! ?? what is this for? move back to aurora?

end subroutine calc_fang_rates

subroutine calculate_HP(iBlock, HPn, HPs)
  ! Calculate hemispheric power in n/s (NPn, HPs)
  ! Scaling is done in aurora

  use ModGITM
  use ModUserGITM
  use ModSources
  use ModIndicesInterfaces
  use ModInputs
  use ModMpi

  real :: LocalVar
  real, intent(out) :: HPn, HPs
  integer, intent(in) :: iBlock

  do i = 1, nLats
    do j = 1, nLons

      eflx_ergs = ElectronEnergyFluxDiffuse(j, i) !/ (1.0e-7 * 100.0 * 100.0)

      if (eflx_ergs > 0.1) then
        eflux = eflx_ergs*6.242e11  ! ergs/cm2/s -> eV/cm2/s

        !(eV/cm2/s -> J/m2/s)
        power = eflux*Element_Charge*100.0*100.0* &
                dLatDist_FB(j, i, nAlts, iBlock)* &
                dLonDist_FB(j, i, nAlts, iBlock)

        if (latitude(i, iBlock) < 0.0) then
          HemisphericPowerSouth = HemisphericPowerSouth + power
        else
          HemisphericPowerNorth = HemisphericPowerNorth + power
        endif

      endif

    enddo
  enddo

end subroutine calculate_HP

