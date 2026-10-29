! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

subroutine get_temperature(lon, lat, alt, t, h)

  use ModInputs
  use ModPlanet

  implicit none

  real, intent(in) :: lon, lat, alt
  real, intent(out) :: t, h

  real    :: tAve, tDiff, n, r, g, m
  integer :: iSpecies
  !---------------------------------------------------------------------------
  if (UseMsis) then

    call initialize_msis_routines
    call get_msis_temperature(lon, lat, alt, t, h)

  else

    tAve = (TempMax + TempMin)/2
    tDiff = (TempMax - TempMin)/2

    t = tAve + tDiff*tanh((alt/1000.0 - TempHeight)/TempWidth)

    r = RBody + alt
    g = Gravitational_Constant*(RBody/r)**2

    h = 0.0
    n = 0.0
    m = 0.0
    do iSpecies = 1, nSpecies
      m = m + exp(LogNS0(iSpecies))*mass(iSpecies)
      n = n + exp(LogNS0(iSpecies))
    enddo

    m = m/n
    h = Boltzmanns_Constant*t/(m*g)

  endif

end subroutine get_temperature

!=============================================================================

subroutine fill_scale_heights(TrialdHFactor, ScaleHeights, AltTop)

  ! Build the scale heights a given dHFactor yields, and report the altitude
  ! of the top physical level (iAlt = nAlts).

  use ModGITM
  use ModInputs
  use ModTime

  implicit none

  real, intent(in)  :: TrialdHFactor
  real, intent(out) :: ScaleHeights(nAlts)
  real, intent(out) :: AltTop

  integer :: iAlt
  real    :: geo_lat, geo_lst, geo_lon, geo_alt, h, t
  !----------------------------------------------------------------------------

  ScaleHeights = 0.0

  do iAlt = 1, nAlts

    geo_lat = 0.0
    geo_lst = 12.0
    geo_lon = mod(geo_lst*15.0 - utime/3600.0*15.0 + 360.0, 360.0)
    geo_alt = AltMin
    if (iAlt > 1) geo_alt = AltMin + sum(ScaleHeights(1:iAlt - 1))*TrialdHFactor

    geo_lon = geo_lon*pi/180.0

    call get_temperature(geo_lon, geo_lat, geo_alt, t, h)
    ScaleHeights(iAlt) = h

  enddo

  AltTop = geo_alt

end subroutine fill_scale_heights

!=============================================================================

subroutine init_altitude

  !---------------------------------------------------------------------------
  !  Build the stretched altitude grid.  Levels start at AltMin and are spaced
  !  dHFactor scale heights apart, so nAlts, AltMin and MSIS set where the top
  !  lands; AltMax is a ceiling, and dHFactor is reduced until the top fits
  !  under it.  Setting both AltMax and #DHFACTOR skips that and uses them as
  !  given, which is how to exceed AltMaxLimit or dHFactorLimit.
  !---------------------------------------------------------------------------

  use ModGITM
  use ModInputs
  use ModTime

  implicit none

  integer, parameter :: nIterMax = 60
  real, parameter :: AltTolerance = 100.0    ! m
  real, parameter :: dHFactorFloor = 0.01

  integer :: iAlt, iLoop, nAltsFits

  real :: ScaleHeights(nAlts)
  real :: dHFactorCoarsest, dHFactorUsed, dHFactorLow, dHFactorHigh
  real :: AltCeiling, AltTop
  logical :: IsTrustedGrid, DoSolve
  !----------------------------------------------------------------------------

  if (IsAltMaxSet .and. AltMax <= AltMin) then
    write(*, *) 'AltMax must be above AltMin in #ALTITUDE.'
    write(*, *) 'AltMin, AltMax (km) : ', AltMin/1000.0, AltMax/1000.0
    call stop_gitm('Incorrect altitude range in init_altitude')
  endif

  ! With both AltMax & dHFactor, trust the user
  IsTrustedGrid = IsAltMaxSet .and. IsDHFactorSet

  ! The spacing to aim for, and the altitude it may not pass
  dHFactorCoarsest = dHFactorLimit
  if (IsDHFactorSet) dHFactorCoarsest = dHFactor

  AltCeiling = AltMaxLimit
  if (IsAltMaxSet) AltCeiling = min(AltMax, AltMaxLimit)

  if (dHFactorCoarsest > dHFactorLimit .and. iProc == 0) then
    write(*, '(a,f7.4,a,f4.2,a)') &
      'WARNING!!  init_altitude :  dHFactor=', dHFactorCoarsest, &
      ' is coarser than the recommended ', &
      dHFactorLimit, '. Results may be unreliable'
  endif

  if (IsAltMaxSet .and. AltMax > AltMaxLimit .and. iProc == 0) then
    write(*, '(a,f0.1,a)') &
      'WARNING!!  init_altitude : AltMax=', AltMax/1000.0, &
      " km, is beyond GITM's validated range. Results may be unreliable."
  endif

  ! The requested spacing stands unless it overshoots a ceiling we enforce
  call fill_scale_heights(dHFactorCoarsest, ScaleHeights, AltTop)

  DoSolve = (AltTop > AltCeiling) .and. (.not. IsTrustedGrid)

  if (DoSolve) then

    ! How many levels would have fit at the requested spacing.  This is the
    ! nAlts to recompile with to keep the resolution instead of the range.
    nAltsFits = nAlts
    do iAlt = 2, nAlts
      if (AltMin + sum(ScaleHeights(1:iAlt - 1))*dHFactorCoarsest > AltCeiling) then
        nAltsFits = iAlt - 1
        exit
      endif
    enddo

    ! Bisect between a spacing that certainly fits and the coarsest allowed.
    dHFactorLow = dHFactorFloor
    dHFactorHigh = dHFactorCoarsest

    call fill_scale_heights(dHFactorLow, ScaleHeights, AltTop)
    if (AltTop > AltCeiling) then
      write(*, '(a,i0,a,f0.1,a,f0.1,a)') &
        ' init_altitude : ', nAlts, ' levels cannot fit between AltMin ', &
        AltMin/1000.0, ' km and ', AltCeiling/1000.0, ' km at any spacing.'
      write(*, '(a)') &
        '   Raise AltMax, lower AltMin, or recompile with fewer levels.'
      call stop_gitm('Cannot fit the vertical grid under the ceiling')
    endif

    do iLoop = 1, nIterMax
      dHFactorUsed = 0.5*(dHFactorLow + dHFactorHigh)
      call fill_scale_heights(dHFactorUsed, ScaleHeights, AltTop)
      if (AltTop > AltCeiling) then
        dHFactorHigh = dHFactorUsed
      else
        dHFactorLow = dHFactorUsed
        if (AltCeiling - AltTop < AltTolerance) exit
      endif
    enddo

    ! End on a spacing known to fit, whichever side the loop stopped on
    dHFactorUsed = dHFactorLow
    call fill_scale_heights(dHFactorUsed, ScaleHeights, AltTop)

  else

    dHFactorUsed = dHFactorCoarsest

  endif

  if (iProc == 0) then
    write(*, '(a,i6)') '    nAlts (compile-time) : ', nAlts
    write(*, '(a,f9.2)') '    AltMin          (km) : ', AltMin/1000.0
    write(*, '(a,f9.2)') '    AltMax          (km) : ', AltTop/1000.0
    write(*, '(a,f7.4)') '    dHFactor             : ', dHFactorUsed

    ! Only reachable on a trusted grid; every other path caps at AltMaxLimit
    if (AltTop > AltMaxLimit) &
      write(*, '(a,f0.1,a)') '   WARNING!!  this top is beyond the ', &
      AltMaxLimit/1000.0, ' km GITM is validated to. Results may be unreliable.'

    if (DoSolve) then
      write(*, '(a,f7.4,a,f0.1,a)') '   reduced from ', dHFactorCoarsest, &
        ' to fit under ', AltCeiling/1000.0, ' km'
      write(*, '(a,f7.4,a,i0)') '   to keep ', dHFactorCoarsest, &
        ', recompile with nAlts = ', nAltsFits
    endif
  endif

  ! Record the grid that was built, not the one that was asked for
  dHFactor = dHFactorUsed
  AltMax = AltTop

  ! Fill min & lower 2 ghost cells' altitudes
  Altitude_GB(:, :, 1, 1:nBlocks) = AltMin
  Altitude_GB(:, :, 0, 1:nBlocks) = AltMin - dHFactorUsed*ScaleHeights(1)
  Altitude_GB(:, :, -1, 1:nBlocks) = AltMin - 2*dHFactorUsed*ScaleHeights(1)

  do iAlt = 2, nAlts + 1
    Altitude_GB(:, :, iAlt, 1:nBlocks) = Altitude_GB(:, :, iAlt - 1, 1:nBlocks) &
                                         + dHFactorUsed*ScaleHeights(iAlt - 1)
    if (iDebugLevel > 3) write(*, *) "Altitude, dHFactor, ScaleHeight : ", &
      Altitude_GB(1, 1, iAlt, 1), dHFactorUsed, ScaleHeights(iAlt - 1)
  enddo

  Altitude_GB(:, :, nAlts + 2, 1:nBlocks) = Altitude_GB(:, :, nAlts + 1, 1:nBlocks) &
                                            + dHFactorUsed*ScaleHeights(nAlts)

end subroutine init_altitude
