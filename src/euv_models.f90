! Copyright 2026, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

!-------------------------------------------------------------------
! Subroutine for calculating scaled solar flux.
!-------------------------------------------------------------------

subroutine calc_scaled_euv

  use ModEUV
  use ModInputs
  use ModTime
  use ModIndicesInterfaces
  use ModGITM, only: dt
  use ModB107, only: b107_convert

  implicit none

  integer :: iError, N, NN
  ! Solar index for the empirical models: F10.7, or B10.7 with #B107
  real :: f107e, f107ae, f107, f107a
  real :: EuvScale

  real, dimension(Num_WaveLengths_High) :: FISM_Flux
  real, dimension(Num_WaveLengths_Low) :: EUVAC_Flux
  real, dimension(Num_WaveLengths_High) :: Tobiska_Flux
  real, dimension(Num_WaveLengths_High) :: Empirical_Flux

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

  ! The empirical models (EUVAC, Tobiska) use f107e/f107ae: B10.7
  call b107_convert(f107, f107a, f107e, f107ae)

  ! Defaults, and explainations:
  !   UseEUVData = FISM. False by default
  !   UseRidleyEUV is False by default
  !
  !   UseEUVAC is False by default
  !   UseTobiska is False by default

  FISM_Flux = 0.0
  Empirical_Flux = 0.0
  Tobiska_Flux = 0.0
  EUVAC_Flux = 0.0

  ! --------------------------------------------------------
  ! These are Empirical Models:
  if (UseRidleyEUV) call run_ridley_euv(Empirical_Flux, f107, f107a)
  if (UseEUVAC .or. UseTobiska) then
     ! The old way of doing this is to run both, then blend them
     call run_euvac(EUVAC_Flux, f107e, f107ae)
     call run_tobiska_euv(Tobiska_Flux, f107e, f107ae)

     ! EUVAC is stored in EUVAC_Flux
     if (UseEUVACBlend .or. UseTobiska) &
          Empirical_Flux = Tobiska_Flux
     do N = 1, Num_WaveLengths_Low
        NN = N + 15
        if (UseTobiska) then
           Empirical_Flux(NN) = 0.5*(EUVAC_Flux(N) + Tobiska_Flux(NN))
        else
           Empirical_Flux(NN) = EUVAC_Flux(N)
        endif
     enddo
  endif

  ! --------------------------------------------------------
  ! This is for FISM:
  if (UseEUVData) then 
    call get_FISM_Fluxes(FISM_Flux)

    ! #EUVSCALE: flat F10.7a-dependent scaling of the spectrum. It was fit
    ! against FISM, so it multiplies only the FISM term
    ! If you want PURE FISM, then set:
    !   #EUVSCALE
    !   EuvScaleBase = 1.0
    !   EuvScaleSlope = 0.0
    !   EuvScaleF107aRef = 100.0 (doesn't matter)
    !
    !   and
    !
    !   EUV_Ratio_Empirical = 0.0 (should be default)
    EuvScale = max(0.0, EuvScaleBase + EuvScaleSlope*(f107a - EuvScaleF107aRef))

    ! Blend Data and Empirical (set in set_inputs) and
    ! take into account the sun distance to the planet.  Applied after
    ! Flux_of_EUV is assigned, so repeated calls cannot compound it:
    Flux_of_EUV = &
      (EUV_Ratio_Empirical * Empirical_Flux + &
      (1.0 - EUV_Ratio_Empirical) * FISM_Flux*EuvScale)/ &
      (SunPlanetDistance**2)
  else
    Flux_of_EUV = &
      Empirical_Flux / (SunPlanetDistance**2)
  endif

  if (iDebugLevel >= 2) then
     !do N = 1, Num_WaveLengths_High
     n = 1
     write(*, *) 'flux of euv in bin : ', n, Flux_of_euv(n)
     !enddo
  endif

  ! --------------------------------------------------------
  ! Calculate total integrated energy for log file:

  TotalIntegratedEuvEnergy = 0.0
  do N = 1, Num_WaveLengths_High
    TotalIntegratedEuvEnergy = &
      TotalIntegratedEuvEnergy + &
      Flux_of_EUV(N)*(Planck_Constant*Speed_Light)/ &
      (wavelengthmid(N)*1.0e-10)
  enddo

end subroutine calc_scaled_euv

!-------------------------------------------------------------------
! Set EUV Flags
!-------------------------------------------------------------------

subroutine set_euv_flags(iError)

  use ModInputs
  implicit none

  logical :: isFound

  integer, intent(out) :: iError

  iError = 0
  isFound = .false.

  ! If this is called, we need to turn everything off, and allow
  ! the user to set a value
  UseEUVData = .false.
  UseEUVAC = .false.
  UseTobiska = .false.
  UseAboveHigh = .false.
  UseBelowLow = .false.
  UseRidleyEUV = .false.
  RidleyEUVModelNumber = 2

  ! Turn everything to lower case:
  call lower_case(cEUVModel)

  if (cEUVModel(1:4) == 'fism') then
     call report(' -> Setting EUV Model to FISM', 0)
    UseEUVData = .true.
    isFound = .true.
  endif

  if (cEUVModel(1:6) == 'ridley') then
    UseRidleyEUV = .true.
    if (cEUVModel(7:7) == '0') then
      RidleyEUVModelNumber = 0
      isFound = .true.
      call report(' -> Setting EUV Model to Ridley0', 0)
    endif
    if (cEUVModel(7:7) == '1') then
      RidleyEUVModelNumber = 1
      isFound = .true.
      call report(' -> Setting EUV Model to Ridley1', 0)
    endif
    if (cEUVModel(7:7) == '2') then
      RidleyEUVModelNumber = 2
      isFound = .true.
      call report(' -> Setting EUV Model to Ridley2', 0)
    endif
  endif
  
  if (cEUVModel(1:5) == 'euvac') then
    UseEUVAC = .true.
    isFound = .true.
    call report(' -> Setting EUV Model to EUVAC', 0)
  endif

  ! Augmentation to above
  if (cEUVModel(1:10) == 'euvacblend') then
    UseEUVACBlend = .true.
    call report(' -> Augmenting EUVAC Model', 0)
  endif

  if (cEUVModel(1:5) == 'tobiska') then
    UseTobiska = .true.
    isFound = .true.
    call report(' -> Setting EUV Model to Tobiska', 0)
  endif

  if (.not. isFound) then
    if (iDebugLevel >= 0) then
      write(*,*) 'EUV Models Available: '
      write(*,*) '  - EUVAC: the OG, unaltered'
      write(*,*) '  - EUVACBlend: the OG with some additions from Tobiska'
      write(*,*) '  - Tobiska: a model made by Kent Tobiska?'
      write(*,*) '  - Ridley0: the first Ridley model that is unverified'
      write(*,*) '  - Ridley1: Uses F107-FISM relationship, F107 & F107A have power'
      write(*,*) '  - Ridley2: Uses F107-FISM relationship, has F107-F107A'
    endif
    iError = 1
  else
    iError = 0
  endif

end subroutine set_euv_flags


!-------------------------------------------------------------------
! Subroutine for calculating the Tobiska-EUV flux in a vacuum.
!-------------------------------------------------------------------

subroutine run_tobiska_euv(tobiska_flux, f107, f107a)

  use ModEUV
  use ModInputs

  implicit none

  real, dimension(1:Num_WaveLengths_High), intent(out) :: tobiska_flux
  real, intent(in) :: f107, f107a

  real :: f107_Ratio, r1, r2, hlybr, fexvir, hlya, heiew
  real :: heimod, hlymod, xuvf, xuvfac

  integer :: i, ii, n

  ! This stuff is for Tobiska:
  hlybr = 0.
  fexvir = 0.
  hlya = 3.E+11 + 0.4E+10*(f107 - 70.)
  heiew = 0.

  ! I think that this is the Hinteregger stuff, which goes from
  ! 2A - 1750A, so it is used to fill in the Above and Below spectrum
  f107_ratio = (f107 - 68.0)/(243.0 - 68.0)
  do N = 1, Num_WaveLengths_High
    Tobiska_Flux(N) = RFLUX(N) + (XFLUX(N) - RFLUX(N))*f107_Ratio
  enddo

  xuvfac = 4.0 - f107_ratio
  if (xuvfac < 1.0) xuvfac = 1.0

  ! This is a total hack, just comparing FISM to these fluxes:
  tobiska_flux(Num_WaveLengths_High - 4) = &
    tobiska_flux(Num_WaveLengths_High - 4)*3.0
  tobiska_flux(Num_WaveLengths_High - 3) = &
    tobiska_flux(Num_WaveLengths_High - 3)*3.0
  tobiska_flux(Num_WaveLengths_High - 2) = &
    tobiska_flux(Num_WaveLengths_High - 2)*5.0
  tobiska_flux(Num_WaveLengths_High - 1) = &
    tobiska_flux(Num_WaveLengths_High - 1)*100.0
  tobiska_flux(Num_WaveLengths_High) = &
    tobiska_flux(Num_WaveLengths_High)*800.0

  !   UseAboveHigh is False by default
  !   UseBelowLow is False by default
  if (.not. UseAboveHigh) tobiska_flux(56:Num_WaveLengths_High) = 0
  if (.not. UseBelowLow) tobiska_flux(1:Num_WaveLengths_Low) = 0

  if (HLYA > 0.001) then
    hlymod = hlya
  else
    if (heiew > 0.001) then
      hlymod = heiew*3.77847e9 + 8.40317e10
    else
      hlymod = 8.70e8*f107 + 1.90e11
    endif
  endif

  if (heiew > 0.001) then
    heimod = heiew*3.77847e9 + 8.40317e10
  else
    heimod = hlymod
  endif

  ! hlymod is SME Lyman-alpha
  ! heimod is He I 10,830A, scaled to Lyman-alpha
  do N = 16, 55
    Tobiska_Flux(N) = &
      TCHR0(N) + &
      TCHR1(N)*hlymod + &
      TCHR2(N)*heimod + &
      TCOR0(N) + &
      TCOR1(N)*f107 + &
      TCOR2(N)*f107a
  enddo

  !
  ! Substitute in H Lyman-alpha and XUVFAC if provided:
  !

  if (hlya > 0.001) Tobiska_Flux(12) = hlya/1.E9
  if (xuvfac > 0.001) THEN
    xuvf = xuvfac
  else
    xuvf = 1.0
  endif

  !
  ! Convert from gigaphotons to photons, cm^-2 to m^-2, etc.:
  !

  do N = 1, Num_WaveLengths_High
    IF (Tobiska_Flux(N) < 0.0) Tobiska_Flux(N) = 0.0
    ! I don't know why we scale this...
    IF ((WAVEL(N) < 251.0) .AND. (WAVES(N) > 15.0)) then
      Tobiska_Flux(N) = Tobiska_Flux(N)*xuvf
    endif
    ! Convert to photons/m^2/s
    Tobiska_Flux(N) = Tobiska_Flux(N)*1.E9*10000.0
  enddo

end subroutine run_tobiska_euv

!-------------------------------------------------------------------
! Subroutine for calculating the Ridley-EUV flux in a vacuum.
!-------------------------------------------------------------------

subroutine run_ridley_euv(ridley_flux, f107, f107a)

  use ModEUV
  use ModInputs, only: iDebugLevel

  implicit none

  integer :: i, ii
  real, dimension(1:Num_WaveLengths_High), intent(out) :: ridley_flux
  ! Daily and 81-day solar index to use 
  real, intent(in) :: f107, f107a
  real :: f107e, f107ae, flux, wavelength
  integer :: n

  ! This model was derived WITHOUT the B107 stuff, so let's not
  ! use that stuff. But, the model subtracted off a low F107
  ! value to allow the exponential to be stronger:
  f107e = f107 - RidleyF107Subtract
  f107ae = f107a - RidleyF107Subtract
  if (iDebugLevel > 2) &
    write(*, *) 'Running Ridley EUV model: ',f107e, f107ae, RidleyF107Subtract
  do N = 1, Num_WaveLengths_High
    flux = &
      RidleySlopes(1, N)*(f107e**RidleyPowers(1, N)) + &
      RidleySlopes(2, N)*(f107ae**RidleyPowers(2, N)) + &
      RidleySlopes(3, N)*(f107ae - f107e) + &
      RidleyIntercepts(N)
    if (flux < 0) flux = 0.0
    ridley_Flux(N) = &
      flux * wavelengthmid(n) * 1.0e-10/ &
      (Planck_Constant*Speed_Light)
  enddo

end subroutine run_ridley_euv

!-------------------------------------------------------------------
! Subroutine for calculating the EUV flux in a vacuum.
!-------------------------------------------------------------------

subroutine run_euvac(euvac_flux, f107e, f107ae)

  use ModEUV
  use ModInputs

  implicit none

  integer :: i, ii
  real    :: flxfac
  real, dimension(1:Num_WaveLengths_Low), intent(out) :: euvac_flux
  ! Daily and 81-day solar index to use (F10.7, or B10.7 with #B107)
  real, intent(in) :: f107e, f107ae

  !:::::::::::::::::::::::::::::::: EUVAC :::::::::::::::::::::::
  !------ This EUV flux model uses the F74113 solar reference spectrum and
  !------ ratios determined from Hinteregger's SERF1 model. It uses the daily
  !------ F10.7 flux (F107) and the 81 day mean (F107A) as a proxy for
  !------ scaling
  !------ The fluxes are returned in EUVFLX and correspond to the 37
  !------ wavelength
  !------ bins of Torr et al. [1979] Geophys. Res. Lett. p771.
  !------ See Richards et al. [1994] J. Geophys. Res. p8981 for details.
  !
  !...... F107   = input daily 10.7 cm flux index. (e.g. 74)
  !...... F107A  = input 81 day ave. of daily F10.7 centered on current day
  !...... EUVFLX = output array for EUV flux in units of photons/cm2/sec.
  !
  !
  !----- loop through the wavelengths calculating the scaling factors and
  !----- the resulting solar flux.
  !----- The scaling factors are restricted to be greater than 0.8
  !

  do i = 1, Num_waveLengths_Low

    ! This has to be backwards, due to the way that the wavelengths are defined:
    ii = Num_waveLengths_Low - i + 1
    FLXFAC = (1.0 + AFAC(II)*(0.5*(f107e + f107ae) - 80.0))
    IF (FLXFAC .LT. 0.8) FLXFAC = 0.8
    EUVAC_Flux(i) = F74113(II)*FLXFAC*1.0E9*10000.

  enddo

end subroutine run_euvac

!-------------------------------------------------------------------
! Get FISM Fluxes
!-------------------------------------------------------------------

subroutine get_FISM_fluxes(FISM_Flux)

  use ModEUV
  use ModInputs
  use ModTime
  use ModIndicesInterfaces
  use ModGITM, only: dt

  implicit none

  real, dimension(Num_WaveLengths_High), intent(out) :: FISM_Flux

  !DAVES:
  real :: SeeTime(nSeeTimes), tDiff(nSeeTimes)
  real :: y1(Num_WaveLengths_High), y2(Num_WaveLengths_High), x1, x2, x
  real :: m(Num_WaveLengths_High), k(Num_WaveLengths_High)
  character(len=2) :: dday, dhour, dminute
  character(len=7) :: dtime
  integer :: N, NN, iMin(1) = 0, iError
  real, dimension(Num_WaveLengths_High) :: Timed_Flux
  real :: r

  call read_fism_file(iError)

  if (iError /= 0) then
    call stop_gitm("Stopping in get_FISM_fluxes. "// &
                  "Error in EUV data. Check times!")
  endif

  call start_timing("EUV_Data")
  SeeTime(:) = 0

  do N = 1, nSeeTimes
    SeeTime(N) = TimeSee(N)
  enddo

  tDiff = CurrentTime - SeeTime

  where (tDiff .lt. 0) tDiff = 1.e20
  iMin = minloc(tDiff)

  if (iMin(1) == nSeeTimes .or. &
       tDiff(iMin(1)) == 0) then
     Timed_Flux = SeeFlux(:, iMin(1))
  else
     ! Do a simple linear interpolation:
     r = (CurrentTime - SeeTime(iMin(1))) / &
          (SeeTime(iMin(1)+1) - SeeTime(iMin(1)))
     Timed_Flux = &
          (1-r) * SeeFlux(:, iMin(1)) + &
          (r) * SeeFlux(:, iMin(1) + 1)
  endif

  if (CurrentTime .ge. FlareTimes(iFlare) .and. &
      CurrentTime - dt .le. FlareTimes(iFlare)) then
    FlareStartIndex = iMin(1) + 1
    FlareEndIndex = iMin(1) + FlareLength
    Timed_Flux = SeeFlux(:, FlareStartIndex)
    FlareEndTime = SeeTime(FlareEndIndex)
    DuringFlare = .true.
    iFlare = iFlare + 1
  else
     
    if (DuringFlare) then
      if (Seetime(iMin(1) + 1) .lt. FlareTimes(iFlare) .or. &
          FlareTimes(iFlare) .eq. 0) then
        if (CurrentTime .lt. SeeTime(FlareStartIndex)) then
          ! We may not be to the point where the flare has begun
          ! in the SEE data yet...
          Timed_Flux = SeeFlux(:, FlareStartIndex)
        else
          if (CurrentTime .le. FlareEndTime) then
            !Exponentially interpolate between last seetime and
            !next seetim using y = kexp(-mx)

            y1 = SeeFlux(:, iMin(1))
            y2 = SeeFlux(:, iMin(1) + 1)
            x1 = 0
            x2 = SeeTime(iMin(1) + 1) - SeeTime(iMin(1))
            x = CurrentTime - SeeTime(iMin(1))

            m = ALOG(y2/y1)/(x1 - x2)
            k = y1*exp(m*x1)
            Timed_Flux = k*exp(-1*m*x)
          else
            DuringFlare = .False.
          endif
        endif
      endif
    endif
  endif

  !!need to convert from W/m^2 to photons/m^2/s
  do N = 1, Num_WaveLengths_High
    FISM_Flux(N) = Timed_Flux(N)*wavelengthmid(N)*1.0e-10/ &
                   (Planck_Constant*Speed_Light)
  enddo
  call end_timing("EUV_Data")

end subroutine get_FISM_fluxes

! -------------------------------------------------------------------------------
! Read FISM file.  This works with flares and daily FISM files
! -------------------------------------------------------------------------------

subroutine read_fism_file(iError)

  use ModKind
  use ModEUV
  use ModInputs
  use ModTime

  implicit none

  integer, intent(out)     :: iError

  character(len=4) :: autotest, cYear
  character(len=20) :: line, cline
  character(len=iCharLen_) :: cFileToOpen
  integer, dimension(7) :: TimeOfFlare, TimeArray
  real, dimension(6 + Num_Wavelengths_High) :: temp

  real(Real8_) :: TimeDelay, BufferTime = 86400.0*3

  logical :: NotDone, IsThere
  integer :: i, iline, ioerror, nline, nILine = 1, iErrortemp = 0

  iError = 0

  ! If we have been here before and we read the entire file, leave
  if (.not. ReReadEUVFile .and. nSeeTimes > 0) return

  ! This if statement makes sure that we have been here before and
  ! want to be here again
  if (ReReadEUVFile) then
    ! If we still have a lot of data in memory, then don't bother
    ! reading more.
    if (CurrentTime < TimeSee(nSeeTimes)) return
  endif

  ! Assume that we can read the entire file
  ReReadEUVFile = .false.

  cline = ' '

  ! Alter behaviour here a bit....
  ! Can set cEUVFile to 'auto' and it will find it "automatically"
  autotest = cEUVFile
  call lower_case(autotest)
  if (autotest == 'auto') then
     ! Find the current year:
     call i2s(iTimeArray(1), cYear, 4)
     ! Set the file name
     cFileToOpen = 'UA/DataIn/FISM/fismflux_daily_' // cYear // '.dat'
  else
     cFileToOpen = cEUVFile
  endif
  
  inquire(file=cFileToOpen, EXIST=IsThere)
  if (.not. IsThere) &
    call stop_gitm(trim(cFileToOpen)//" cannot be found by read_inputs")

  if (iDebugLevel >= 0) then
     write(*,*) '> Reading FISM file : ', trim(cFileToOpen)
     if (EUV_Ratio_Empirical > 0.0) &
          write(*,*) '-> EUV_Ratio_Empirical is set to ', EUV_Ratio_Empirical
     write(*,*) '-> FISM Scaling Factors set to : '
     write(*,*) '   EuvScaleBase : ', EuvScaleBase
     write(*,*) '   EuvScaleSlope : ', EuvScaleSlope
     write(*,*) '   EuvScaleF107aRef : ', EuvScaleF107aRef
  endif
  open(unit=iInputUnit_, file=cFileToOpen, status="old", IOSTAT=iError)

  if (iError /= 0) then
    write(*, *) "Error in opening EUV file  Is this set?"
    write(*, *) "Code : ", iError, cFileToOpen
    call stop_gitm("Stopping in calc_euv")
  endif

  iline = 1

  NotDone = .true.
  do while (NotDone)
    read(iInputUnit_, '(a)', iostat=iError) cLine
    if (cline(1:1) .eq. '#') then

      ! Remove anything after a space or TAB

      i = index(cLine, ' '); if (i > 0) cLine(i:len(cLine)) = ' '
      i = index(cLine, char(9)); if (i > 0) cLine(i:len(cLine)) = ' '
      select case (cLine)

      case ("#FLARES")
        read(iInputUnit_, *) nFlares
        do i = 1, nFlares
          read(iInputUnit_, *) TimeOfFlare(1:6)
          TimeOfFlare(7) = 0
          call time_int_to_real(TimeOfFlare, FlareTimes(i))
        enddo

      case ('#START')
        NotDone = .false.

      end select
    endif

    iline = iline + 1

  enddo

  TimeSee(:) = 0
  iline = 1
  read(iInputUnit_, *, iostat=iError) temp

  iErrortemp = 0
  do while (iErrortemp .eq. 0)

    TimeArray(1:6) = temp(1:6)
    TimeArray(7) = 0
    call time_int_to_real(TimeArray, TimeSee(iLine))
    SeeFlux(:, iline) = temp(7:6 + Num_WaveLengths_High)

    if (TimeSee(iLine) > StartTime - BufferTime .and. &
        TimeSee(iLine) < EndTime + BufferTime) &
      iline = iline + 1

    read(iInputUnit_, *, iostat=iErrortemp) temp

    if (iline > nSeeLinesMax - 1) then
      iErrortemp = 1
      ReReadEUVFile = .true.
    endif

  enddo

  close(iInputUnit_)
  nSeeTimes = iline - 1

  ! Let's check to see if we may need to move on to the next file:
  if (autotest == 'auto') then
     if (TimeSee(nSeeTimes) < EndTime + BufferTime) &
          rereadeuvfile = .true.
  endif
    
  if (nSeeTimes > 2) iError = 0
  if (nSeeTimes < 1) iError = 1

end subroutine read_fism_file

! -------------------------------------------------------------------------------
! Old EUV model
! -------------------------------------------------------------------------------

subroutine read_euv_waves(iError)

  use ModKind
  use ModEUV
  use ModInputs

  implicit none

  integer, intent(out)     :: iError
  character(len=80)        :: cline

  integer :: nFiles, nWaves, iWave, i
  real :: tmp(6)

  call report("read_euv_waves", 2)

  RidleyF107Subtract = 0.0

  open(unit=iInputUnit_, file=cRidleyEUVFile, status="old", IOSTAT=iError)

  if (iError /= 0) then
    write(*, *) "Error in opening RidleyEUVFile  Is this set?"
    write(*, *) "Code : ", iError, cRidleyEUVFile
    call stop_gitm("Stopping in read_euv_waves")
  endif

  read(iInputUnit_, *, iostat=iError) nFiles

  do i = 1, nFiles
    read(iInputUnit_, *, iostat=iError) cline
  enddo

  read(iInputUnit_, *, iostat=iError) nWaves

  if (nWaves /= Num_WaveLengths_High) then
    write(*, *) "Error in RidleyEUVFile - nWaves /= Num_wavelengths_high"
    write(*, *) "Values : ", nWaves, Num_WaveLengths_High
    call stop_gitm("Stopping in read_euv_waves")
  endif

  read(iInputUnit_, *, iostat=iError) cline

  do i = 1, nWaves

    read(iInputUnit_, *, iostat=iError) iWave, tmp
    RidleySlopes(1, i) = tmp(1)
    RidleySlopes(2, i) = tmp(2)
    RidleySlopes(3, i) = tmp(3)
    RidleyIntercepts(i) = tmp(4)
    RidleyPowers(1, i) = tmp(5)
    RidleyPowers(2, i) = tmp(6)

  enddo

  close(iInputUnit_)

end subroutine read_euv_waves

! -------------------------------------------------------------------------------
! This is an updated model of the EUV that was made in Nov 2026
! It uses the form of the "Ridley EUV Model", but is slightly
! different.  There are two forms:
! model1 = base + Ca * (F107A-50) ^ Pa + Cd * (F107-50) ^ Pd
! model2 = base + Ca * (F107A-50) ^ Pa + Cd * (F107-F107A)
! It seems like the old model assumes:
! model = Ri + Rs1 * F107 ^ Rp1 + Rs2 * F107A ^ Rp1 + Rs3 * (F107-F107A)
! So, we need to fit to this form.

subroutine read_euv2026_model(iError)

  use ModKind
  use ModEUV
  use ModInputs

  implicit none

  integer, intent(out)     :: iError
  character(len=80)        :: cline
  character(len=iCharLen_) :: cFilename

  integer :: nFiles, nWaves, iWave, i
  real :: tmp1(8), tmp2(7)
  logical :: isDone

  call report("read_euv2026_model", 0)

  RidleyF107Subtract = 50.0

  if (RidleyEUVModelNumber == 1) then
    cFilename = cEUVFileModel1
  endif
  if (RidleyEUVModelNumber == 2) then
    cFilename = cEUVFileModel2
  endif

  open(unit=iInputUnit_, file=cFilename, status="old", IOSTAT=iError)

  if (iError /= 0) then
    write(*, *) "Error in opening EUV model file  Is this set?"
    write(*, *) "Code : ", iError, cFilename
    call stop_gitm("Stopping in read_euv2026_model")
  endif

  cline = ''
  isDone = .false.
  do while (.not.isDone)
    read(iInputUnit_, *, iostat=iError) cLine
    if (cLine(1:7) == "#PARAMS") isDone = .true.
    if (iError /= 0) isDone = .true.
  enddo

  if (iError /= 0) then
    write(*, *) "Error in reading EUV model file"
    write(*, *) "Code : ", iError, cFilename
    call stop_gitm("Stopping in read_euv2026_model")
  endif

  do i = 1, Num_WaveLengths_High

    if (RidleyEUVModelNumber == 1) then
      read(iInputUnit_, *, iostat=iError) iWave, tmp1
      RidleySlopes(1, i) = tmp1(6)
      RidleySlopes(2, i) = tmp1(4)
      RidleySlopes(3, i) = 0.0
      RidleyIntercepts(i) = tmp1(3)
      RidleyPowers(1, i) = tmp1(5)
      RidleyPowers(2, i) = tmp1(7)
    else
      read(iInputUnit_, *, iostat=iError) iWave, tmp2
      RidleySlopes(1, i) = 0.0
      RidleySlopes(2, i) = tmp2(4)
      RidleySlopes(3, i) = tmp2(6)
      RidleyIntercepts(i) = tmp2(3)
      RidleyPowers(1, i) = 1.0
      RidleyPowers(2, i) = tmp2(5)
    endif

  enddo

  close(iInputUnit_)

end subroutine read_euv2026_model

