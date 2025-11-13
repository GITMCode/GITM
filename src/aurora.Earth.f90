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

  real :: Q0
  integer :: i, j, n, iError, iErr, iEnergy
  logical :: HasSomeAurora
  real, dimension(ED_N_Energies) :: &
    e_diffuse_ED_flux, i_diffuse_ED_flux, mono_ED_flux, wave_ED_flux, &
    ED_Flux ! for temp values

  real :: power

  real, dimension(nLons, nLats, nAlts) :: temp, AuroralBulkIonRate

  logical :: IsFirstTime(nBlocksMax) = .true.

  if (IsFirstTime(iBlock)) then
    call initialize_fang_arrays
  else
    if (floor((tSimulation - dT)/dTAurora) == &
        floor(tSimulation/dTAurora)) return
  endif

  AuroralBulkIonRate = 0.0
  AuroralIonRateS = 0.0

  call report("Aurora", 1)
  call start_timing("Aurora")

  if (iBlock == 1) then
    HemisphericPowerNorth = 0.0
    HemisphericPowerSouth = 0.0
  endif

  if (NormalizeAuroraToHP) then

    ! Let's scale our hemispheric power so it is roughly the same as what
    ! is measured.

    call NormalizeDiffuseAuroraToHP

  endif

  ! Reset the hemispheric power

  if (iBlock == 1) then
    HemisphericPowerNorth = 0.0
    HemisphericPowerSouth = 0.0
    HemisphericPowerNorth_diffuse = 0.0
    HemisphericPowerSouth_diffuse = 0.0
    HemisphericPowerNorth_wave = 0.0
    HemisphericPowerSouth_wave = 0.0
    HemisphericPowerSouth_mono = 0.0
    HemisphericPowerNorth_mono = 0.0
    HemisphericPowerNorth_ion = 0.0
    HemisphericPowerSouth_ion = 0.0
  endif

  if (iProc == 0 .and. AveEFactor /= 1.0 .and. IsFirstTime(iBlock)) then
    write(*, *) "Reminder / Warning!!!"
    write(*, *) " --> Auroral Experiment, AveEFactor is set"
    write(*, *) " --> AveEFactor : ", AveEFactor
    write(*, *) " --> This factor is currently applied in get_potential!"
  endif
  if (iProc == 0 .and. IsKappaAurora .and. IsFirstTime(iBlock)) then
    write(*, *) "Reminder / Warning!!!"
    write(*, *) " --> Auroral Experiments, AuroraKappa is set"
    write(*, *) " --> kappa : ", AuroraKappa
  endif

  do i = 1, nLats
    do j = 1, nLons

      UserData2d(j, i, 1, 2:nUserOutputs, iBlock) = 0.0

      e_diffuse_ED_flux = 0.0
      i_diffuse_ED_flux = 0.0
      wave_ED_flux = 0.0
      mono_ED_flux = 0.0

      HasSomeAurora = .false.

      ! For diffuse auroral models (default)
      if (ElectronEnergyFluxDiffuse(j, i) > 0.001 &
          .and. ElectronAverageEnergyDiffuse(j, i) > 0.01 &
          .and. ElectronAverageEnergyDiffuse(j, i) < MaxAveEAurora &
          .and. UseDiffuseAurora &
          ) then
        call do_diffuse_aurora(ElectronEnergyFluxDiffuse(j, i), &
                               ElectronAverageEnergyDiffuse(j, i), &
                               e_diffuse_ED_flux)
        HasSomeAurora = .true.
      endif

      if (IonEnergyFlux(j, i) > 0.001 &
          .and. IonAverageEnergy(j, i) > 0.1 &
          .and. IonAverageEnergy(j, i) < MaxAveEAurora &
          .and. UseIonAurora &
          ) then
        call do_diffuse_aurora(IonEnergyFlux(j, i), &
                               IonAverageEnergy(j, i), &
                               i_diffuse_ED_flux)
        HasSomeAurora = .true.
      endif

      ! Monoenergetic aurora
      if (UseMonoAurora &
          .and. ElectronAverageEnergyMono(j, i) > 0.1 &
          .and. ElectronEnergyFluxMono(j, i) > 0.1 &
          .and. ElectronAverageEnergyMono(j, i) < MaxAveEAurora &
          ) then
        if (HasSomeAurora .or. AllowAurWODiffuse) &
          call do_mono_aurora(ElectronEnergyFluxMono(j, i), &
                              ElectronAverageEnergyMono(j, i), &
                              mono_ED_flux)

      endif

      ! Wave (broadband) aurora
      if (UseWaveAurora &
          .and. ElectronEnergyFluxWave(j, i) > 0.1 &
          .and. ElectronAverageEnergyWave(j, i) > 0.1 &
          .and. ElectronAverageEnergyWave(j, i) < MaxAveEAurora &
          ) then
        if (HasSomeAurora .or. AllowAurWODiffuse) &
          call do_wave_aurora(ElectronEnergyFluxWave(j, i), &
                              ElectronAverageEnergyWave(j, i), &
                              wave_ED_flux)
      endif

      do n = 1, ED_N_Energies
        ED_EnergyFlux(n) = e_diffuse_ED_flux(n) + wave_ED_flux(n) + mono_ED_flux(n)
        ED_Ion_EnergyFlux(n) = i_diffuse_ED_flux(n) ! for consistency
      enddo

      if (HasSomeAurora) &
        call calc_fang_rates(j, i, iBlock, AuroralBulkIonRate)

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

  IsFirstTime(iBlock) = .false.

  HemisphericPowerNorth = HemisphericPowerNorth_diffuse + HemisphericPowerNorth_ion &
                          + HemisphericPowerNorth_mono + HemisphericPowerNorth_wave
  HemisphericPowerSouth = HemisphericPowerSouth_diffuse + HemisphericPowerSouth_ion &
                          + HemisphericPowerSouth_mono + HemisphericPowerSouth_wave

  call end_timing("Aurora")

contains

  ! --------------------------
  ! Diffuse Aurora can be represented by kappa or maxwellian
  ! - This is for Newell diffuse, or all other auroral models
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
      HemisphericPowerSouth_diffuse = HemisphericPowerSouth_diffuse + power
    else
      HemisphericPowerNorth_diffuse = HemisphericPowerNorth_diffuse + power
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

  subroutine calc_maxwellian(eflux, avg_e, max_ed_flux)
    ! From (Fang et al., [2010]), a
    ! Maxwellian is defined as:
    ! DifferentialNumberFlux = Q0/2/E0**3 * E * exp(-E/E0),
    ! where:
    ! Q0 = Total Energy Flux
    ! E0 = Characteristic Energy (0.5*avee)
    ! E = mid-point of energy bin
    !
    ! Calculate as a * E * exp(-E/E0)

    real, intent(in) :: eflux, avg_e
    real, intent(inout), dimension(:) :: max_ed_flux

    real :: a, E0

    E0 = (0.5*avg_e)
    a = (eflux)/2/E0**3

    do n = 1, ED_N_Energies
      max_ed_flux(n) = a*ED_Energies(n)*exp(-ED_Energies(n)/E0)
      ! Maxwellian_Energy-Deposition_flux
      max_ed_flux(n) = &
        max_ed_flux(n)* &
        ED_Energies(n)* &
        ED_delta_energy(n)
    enddo

  end subroutine calc_maxwellian

  ! Gaussian distribution (in log coordinates)
  ! sig controls the width (std dev)
  ! Is not normalized to eflux! Done for mono/wave differently
  subroutine calc_gaussian(avee, gaus_ed_flux, sig)
    real, intent(in) :: avee, sig
    real, intent(inout), dimension(:) :: gaus_ed_flux

    do n = 1, ED_N_Energies
      gaus_ed_flux(n) = exp(-(log10(ED_Energies(n)) - log10(avee))**2 &
                            /(2.0*sig*sig))
    enddo

  end subroutine calc_gaussian

  ! --------------------------
  ! Monoenergetic aurora
  ! --------------------------
  subroutine do_mono_aurora(eflux_ergs, av_kev, Mono_ED_EnergyFlux)

    real, intent(in) :: eflux_ergs, av_kev
    real, intent(inout), dimension(:) :: Mono_ED_EnergyFlux
    real :: avee, eflux, sum_gaus

    avee = av_kev*1000.0        ! keV -> eV
    eflux = eflux_ergs*6.242e11  ! ergs/cm2/s -> eV/cm2/s

    power = eflux*Element_Charge*100.0*100.0* &    ! (eV/cm2/s -> J/m2/s)
            dLatDist_FB(j, i, nAlts, iBlock)* &
            dLonDist_FB(j, i, nAlts, iBlock)

    if (latitude(i, iBlock) < 0.0) then
      HemisphericPowerSouth_mono = HemisphericPowerSouth_mono + power
    else
      HemisphericPowerNorth_mono = HemisphericPowerNorth_mono + power
    endif

    ! Mono is treated as a (skinny) gaussian, width of 0.1 keV
    if (avee > 0.0) &
      call calc_gaussian(avee, Mono_ED_EnergyFlux, 0.1)

    ! Normalize to E-Flux
    sum_gaus = sum(Mono_ED_EnergyFlux)
    do n = 1, ED_N_Energies
      Mono_ED_EnergyFlux(n) = Mono_ED_EnergyFlux(n)*eflux/sum_gaus
    enddo

  end subroutine do_mono_aurora

  ! --------------------------
  ! Wave aurora
  ! --------------------------
  subroutine do_wave_aurora(eflux_ergs, av_kev, wave_EDEnergyFlux)
    real, intent(in) :: eflux_ergs, av_kev
    real, intent(inout), dimension(:) :: wave_EDEnergyFlux
    real, dimension(ED_N_Energies) :: tmp_wavflux
    real :: avee, eflux, sum_gaus
    integer :: nearBin

    avee = av_kev*1000.0        ! keV -> eV
    eflux = eflux_ergs*6.242e11  ! ergs/cm2/s -> eV/cm2/s

    power = eflux*Element_Charge*100.0*100.0* &    ! (eV/cm2/s -> J/m2/s)
            dLatDist_FB(j, i, nAlts, iBlock)* &
            dLonDist_FB(j, i, nAlts, iBlock)

    if (latitude(i, iBlock) < 0.0) then
      HemisphericPowerSouth_wave = HemisphericPowerSouth_wave + power
    else
      HemisphericPowerNorth_wave = HemisphericPowerNorth_wave + power
    endif

    ! Wave is a gaussian, centered at aveE, width of 0.4 keV
    if (avee > 0.1) &
      call calc_gaussian(avee, tmp_wavflux, 0.4)

    ! But only in a few bins!
    ! Flux only exists in bin holding aveE & 3 bins to either side
    do iEnergy = 1, ED_N_Energies - 1
      if (ED_Energies(iEnergy) < avee .and. ED_Energies(iEnergy + 1) > avee) then
        ! center bin
        wave_EDEnergyFlux(iEnergy) = tmp_wavflux(iEnergy)
        do nearBin = 1, 3 ! 3 on either side
          if (iEnergy + nearBin < ED_N_Energies) &
            wave_EDEnergyFlux(iEnergy + nearBin) = tmp_wavflux(iEnergy + nearBin)
          if (iEnergy - nearBin > 0) &
            wave_EDEnergyFlux(iEnergy - nearBin) = tmp_wavflux(iEnergy - nearBin)
        enddo
      endif
    enddo

    ! Then normalize
    sum_gaus = sum(wave_EDEnergyFlux)
    ! write(*,*) sum_gaus,
    do n = 1, ED_N_Energies
      wave_EDEnergyFlux(n) = wave_EDEnergyFlux(n)*eflux/sum_gaus
    enddo
  end subroutine do_wave_aurora

  subroutine NormalizeDiffuseAuroraToHP

    real :: hpi_NH, hpi_SH, Hpi
    real :: LocalVar, HPn, HPs, avepower, ratio, ratio_NH, ratio_SH

    HPn = 0
    HPs = 0
    ! Calculate HP from (diffuse) e- flux
    call calculate_HP(HPn, HPs)

    ! Collect all of the powers by summing them together
    LocalVar = HPn/1.0e9
    call MPI_REDUCE(LocalVar, HPn, 1, MPI_REAL, MPI_SUM, &
                    0, iCommGITM, iError)

    LocalVar = HPs/1.0e9
    call MPI_REDUCE(LocalVar, HPs, 1, MPI_REAL, MPI_SUM, &
                    0, iCommGITM, iError)

    ! Average north and south together
    avepower = (HPn + HPs)/2.0

    if (avepower == 0) then
      call set_error("AvePower is zero!")
      avepower = 1.0
    endif

    ! If we are only have one hemisphere or the other, assign to avepower
    if (HPs < 0.1*HPn) then
      avepower = HPn
      HPs = HPn
    endif
    if (HPn < 0.1*HPs) then
      avepower = HPs
      HPn = HPs
    endif

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
      ratio_NH = Hpi/HPn
      ratio_SH = Hpi/HPs

    endif

    if (iDebugLevel >= 0) then
      if ((iDebugLevel == 0) .and. IsFirstTime(iBlock) .and. .not. doSeparateHPI) then
        write(*, *) '---------------------------------------------------'
        write(*, *) 'Using auroral normalizing ratios!!! '
        write(*, *) '  -> forcing the HP to be the same in both hemispheres!'
        write(*, *) '  -> no longer reporting!'
        write(*, *) '---------------------------------------------------'
      elseif ((iDebugLevel == 0) .and. IsFirstTime(iBlock) .and. doSeparateHPI) then
        write(*, *) '---------------------------------------------------'
        write(*, *) 'Using HPI from each hemisphere to normalize aurora!!'
        write(*, *) '  -> Normalizing to mean, then calculating ratios seperately'
        write(*, *) '  -> no longer reporting!'
        write(*, *) '---------------------------------------------------'
      endif
      if (iDebugLevel >= 1) then
        if (doSeparateHPI) then
          write(*, *) 'Auroral normalizing ratios: '
          write(*, *) 'HP(N)  HP(S)  HP(N-model)  HP(S-model)  ratio_n  ratio_s'
          write(*, *) Hpi_NH, Hpi_SH, HPn, HPs, ratio_NH, ratio_sh
        else
          write(*, *) 'Auroral normalizing ratio (hp, hp model averaged, ratios (S/N)): ', &
            Hpi, avepower, ratio_SH, ratio_NH
        endif
      endif

    endif
    do i = 1, nLats
      do j = 1, nLons
        if (ElectronEnergyFluxDiffuse(j, i) > 0.001) then
          if (latitude(i, iBlock) < 0.0) then
            ElectronEnergyFluxDiffuse(j, i) = ElectronEnergyFluxDiffuse(j, i)*ratio_sh
          else
            ElectronEnergyFluxDiffuse(j, i) = ElectronEnergyFluxDiffuse(j, i)*ratio_NH
          endif
        endif
      enddo
    enddo
  end subroutine NormalizeDiffuseAuroraToHP

  subroutine calculate_HP(HPn, HPs)
    ! Calculate hemispheric power in n/s (NPn, HPs)
    ! Scaling is done in NormalizeDiffuseAuroraToHP

    real, intent(out) :: HPn, HPs
    real :: eflx_ergs, eflux

    HPn = 0.0
    HPs = 0.0

    do i = 1, nLats
      do j = 1, nLons

        eflx_ergs = ElectronEnergyFluxDiffuse(j, i) !/ (1.0e-7 * 100.0 * 100.0)

        if (eflx_ergs > 0.001) then
          eflux = eflx_ergs*6.242e11  ! ergs/cm2/s -> eV/cm2/s

          !(eV/cm2/s -> J/m2/s)
          power = eflux*Element_Charge*100.0*100.0* &
                  dLatDist_FB(j, i, nAlts, iBlock)* &
                  dLonDist_FB(j, i, nAlts, iBlock)

          if (latitude(i, iBlock) < 0.0) then
            HPs = HPs + power
          else
            HPn = HPn + power
          endif

        endif

      enddo
    enddo

  end subroutine calculate_HP

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

  call init_energy_deposition()

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
  real, dimension(nLons, nLats, nAlts), intent(inout) :: AuroralBulkIonRate
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

      fac = ED_Ion_EnergyFlux(iEnergy)/1000.0/ &
            Fang_de/ &
            BulkScaleHeight1d

      AuroralBulkIonRate(j, i, 1:nAlts) = &
        AuroralBulkIonRate(j, i, 1:nAlts) + 1e6*Fang_Ion_f(iEnergy, :)*fac

    endif

  enddo

end subroutine calc_fang_rates

! ================================== !
! create the energy deposition stuff !
! ================================== !

subroutine init_energy_deposition()
  use ModSources

  ! temporary
  real, dimension(ED_N_Energies):: energy_edges
  real ed_max_energy, ed_min_energy, de, logmin

  allocate(ED_Energies(ED_N_Energies), stat=ierr)
  allocate(ED_EnergyFlux(ED_N_Energies), stat=ierr)
  allocate(ED_Ion_EnergyFlux(ED_N_Energies), stat=ierr)
  allocate(ED_delta_energy(ED_N_Energies), stat=ierr)

  ! min/max energy bins
  ED_Max_Energy = 5.0e5
  ED_Min_Energy = 10.0

  ! log space the energy bins
  de = (log10(ED_Max_Energy) - log10(ED_Min_Energy))/(ED_N_Energies - 1)
  logmin = log10(ED_Min_Energy)

  do i = 1, ED_N_Energies
    ED_Energies(i) = 10.0**(logmin + (i - 1)*de)
  enddo

  ! Calculate delta btwn bins
  ED_delta_energy(1) = ED_Energies(1) - ED_Energies(2)
  energy_edges(1) = ED_Energies(1) + ED_delta_energy(1)/2.0
  iEnergy = 1
  do iEnergy = 2, ED_N_Energies
    energy_edges(iEnergy) = &
      (ED_Energies(iEnergy - 1) + ED_Energies(iEnergy))/2
    ED_delta_energy(iEnergy - 1) = &
      abs(energy_edges(iEnergy - 1) - energy_edges(iEnergy))
  enddo

  iEnergy = ED_N_Energies
  ED_delta_energy(ED_N_Energies) = &
    ED_Energies(ED_N_Energies) - ED_Energies(ED_N_Energies - 1)

end subroutine init_energy_deposition

