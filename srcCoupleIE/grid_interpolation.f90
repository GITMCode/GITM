! --------------------------------------------------------------------
!\
! This routine finds a point on in the spherical file system, given
! a Theta, Phi:
! LocIn(1) = Phi
! LocIn(2) = Theta
! It returns a 5 element array:
! LocOut(1) = Index of Longitude
! LocOut(2) = Index of Latitude
! LocOut(4) = Multiplication factor for Longitude
! LocOut(5) = Multiplication factor for Latitude
!/
! --------------------------------------------------------------------

subroutine find_ua_point(this, LocIn, LocOut)

  use ModErrors

  implicit none

  class(ieModel) :: this
  real, dimension(2), intent(in)  :: LocIn
  real, dimension(4), intent(out) :: LocOut
  real :: MLTIn, LatIn, MLTUp, MLTDown
  integer :: LatIndex, MltIndex, iPt = 0
  real :: MLTMin = 0, MLTMax = 24.0, MLThalf = 12.0, LatMin = -90.0, LatMax = 90.0

  logical :: IsFound

  character(len=*), parameter :: NameSub = "find_point"

  ! this was taken from AMIE, not sure will be correct for IE
  ! Might need to convert MLTs to Lons or vice versa
  LocOut = -1.0

  !\
  ! Check to see if the point is even on the grid.
  !/

  MLTIn = mod(LocIn(1) + MLTMax, MLTMax)

  LatIn = LocIn(2)
  if (LatIn > LatMax) then
    LatIn = 2*LatMax - LatIn
    MLTIn = mod(MLTIn + MLThalf, MLTMax)
  endif
  if (LatIn < LatMin) then
    LatIn = 2*LatMin - LatIn
    MLTIn = mod(MLTIn + MLThalf, MLTMax)
  endif

  if (MLTIn > MLTMax .or. MLTIn < MLTMin .or. &
      LatIn > LatMax .or. LatIn < LatMin) then
    ! should update to be var based
    call set_error("Input lat / mlt is outside of -90-90 and 0-24 range! " &
                   //NameSub)
    return
  endif

  ! only works for regular grid, should eventually check in case IE uses
  ! something else

  MLTs: do iPt = 1, this%havenMlts
    MLTUp = this%haveMLTs(iPt + 1, 2)
    MLTDown = this%haveMLTs(iPt, 2)
    if (MLTUp == 0.0 .and. MLTDown >= 23.0) MLTUp = 24.0
    if ((MltUp > MLTIn) .and. (MltDown <= MLTIn)) then
      MltIndex = iPt
      exit MLTs
    endif
  enddo MLTs

  do iPt = 1, this%havenLats - 1
    if ((this%haveLats(2, iPt + 1) < LatIn) .and. (this%haveLats(2, iPt) >= LatIn)) &
      LatIndex = iPt
  enddo

  ! check my work to be sure!
  !MLTUp = this%haveMLTs(MltIndex+1, LatIndex)
  !MLTDown = this%haveMLTs(MltIndex, LatIndex)

  ! store needed values
  LocOut(1) = MltIndex
  LocOut(2) = LatIndex

  LocOut(3) = (MLTIn - MLTDown)/(MLTUp - MLTDown)
  ! To keep the MLT and LAT ratios uses consistent, we need to calculate
  ! the ratios backward, since lat shrinks with increasing index:
  LocOut(4) = (this%haveLats(MltIndex, LatIndex) - LatIn)/ &
              (this%haveLats(MltIndex, LatIndex) - &
               this%haveLats(MltIndex, LatIndex + 1))

  if (this%iDebugLevel > 4) &
    write(*, *) 'file point!', LatIn, MltIn, LocOut
  if (this%iDebugLevel > 7) then
    if (LatIn > this%haveLats(MltIndex, LatIndex)) &
      write(*, *) "LatIn, LatUp", LatIn, this%haveLats(MltIndex, LatIndex)
    if (LatIn <= this%haveLats(MltIndex, LatIndex + 1)) &
      write(*, *) "LatIn, LatDown", LatIn, this%haveLats(MltIndex, LatIndex + 1)
    if (MLTIn >= MLTUp) write(*, *) "MLTIn, MLTUp", MltIn, MltUP
    if (MltIn < MltDown) write(*, *) "MltIn, MltDown", MltIn, MltDown
  endif

end subroutine find_ua_point
!==============================================================================
subroutine set_ie_ua_interpolation_indices(this, mltsIn, latsIn)

  class(ieModel) :: this
  real, intent(in) :: mltsIn(this%neednMlts, this%neednLats), &
                      latsIn(this%neednMlts, this%neednLats)
  real, dimension(2) :: mlt_and_lat
  real, dimension(4) :: interpolation_info
  integer :: iError, iLat, iMlt, im, il

  if (this%iDebugLevel > 2) &
    write(*, *) "=> Getting IE->UA interpolation indices", &
    this%neednMlts, this%neednLats
  if (allocated(this%IeUaInterpolationIndices)) then
    deallocate(this%IeUaInterpolationIndices)
    deallocate(this%IeUaInterpolationRatios)
  endif
  allocate(this%IeUaInterpolationIndices(this%neednMlts, this%neednLats, 2), &
           stat=iError)
  if (iError /= 0) then
    call set_error("Error allocating IeUaInterpolationIndices!")
    return
  endif
  allocate(this%IeUaInterpolationRatios(this%neednMlts, this%neednLats, 2), &
           stat=iError)
  if (iError /= 0) then
    call set_error("Error allocating IeUaInterpolationRatios!")
    return
  endif

  this%IeUaInterpolationIndices = -1

  do iMlt = 1, this%neednMlts
    do iLat = 1, this%neednLats
      mlt_and_lat(1) = mltsIn(iMlt, iLat)
      mlt_and_lat(2) = latsIn(iMlt, iLat)
      call this%find_ua_point(mlt_and_lat, interpolation_info)
      if (iError == 0) then
        this%IeUaInterpolationIndices(iMlt, iLat, 1:2) = &
          interpolation_info(1:2)
        this%IeUaInterpolationRatios(iMlt, iLat, 1:2) = &
          interpolation_info(3:4)

        iM = this%IeUaInterpolationIndices(iMLT, iLat, 1)
        iL = this%IeUaInterpolationIndices(iMLT, iLat, 2)
      else
        this%IeUaInterpolationIndices(iMlt, iLat, 1:2) = -1
      endif
    enddo
  enddo

end subroutine set_ie_ua_interpolation_indices
!==============================================================================

subroutine get_ie_values_for_ua(this, iVarToGetIn, valueOut)

  class(ieModel) :: this
  character(len=*), intent(in) :: iVarToGetIn
  real, intent(inout) ::  valueOut(this%neednMlts, this%neednLats)
  real :: current_var(this%havenMlts, this%havenLats)
  character(len=*), parameter :: NameSub = "get_ie_values_for_ua"

  valueOut = 0.0

  select case (iVarToGetIn)
  case ("pot")
    current_var = this%havePotential
  case ("def")
    current_var = this%haveDiffuseEeFlux
  case ("dae")
    current_var = this%haveDiffuseEAveE
  case ("mef")
    current_var = this%haveMonoEeFlux
  case ("mae")
    current_var = this%haveMonoEAveE
  case ("wef")
    current_var = this%haveWaveEeFlux
  case ("wae")
    current_var = this%haveWaveEAveE
  case ("ief")
    current_var = this%haveDiffuseIeFlux
  case ("iae")
    current_var = this%haveDiffuseIAveE
  case default
    call CON_stop(NameSub//": "//iVarToGetIn//" is not a valid variable to get")
  end select

  call this%remap_conservative(this%haveMlts, this%haveLats, current_var, &
                                this%needMlts, this%needLats, valueOut)

end subroutine get_ie_values_for_ua

!==============================================================================
subroutine get_ie_spec_for_ua(this, iVarToGetIn, valueOut)

  class(ieModel) :: this
  character(len=*), intent(in) :: iVarToGetIn
  real, intent(inout) ::  valueOut(this%neednMlts, this%neednLats, this%needNEnergyBins)
  integer :: iBin
  real :: current_var(this%havenMlts, this%havenLats, this%havenEnergyBins)
  real :: slice_fine(this%havenMlts, this%havenLats)
  real :: slice_coarse(this%neednMlts, this%neednLats)
  character(len=*), parameter :: NameSub = "get_ie_spec_for_ua"

  valueOut = 0.0

  select case (iVarToGetIn)
  case ("ele")
    current_var = this%haveElecNflux
  case ("hyd")
    current_var = this%haveHydrNflux
  case default
    call CON_stop(NameSub//": "//iVarToGetIn//" is not a valid variable to get")
  end select

  do iBin = 1, this%needNEnergyBins
    slice_fine = current_var(:,:,iBin)
    call this%remap_conservative(this%haveMlts, this%haveLats, slice_fine, &
                                  this%needMlts, this%needLats, slice_coarse)
    valueOut(:,:,iBin) = slice_coarse
  enddo

end subroutine get_ie_spec_for_ua

! --------------------------------------------------------------------
!\
! Remap a field from the (high resolution) IE grid onto the (lower
! resolution) grid the caller needs, by area-weighted averaging.
! Cell edges are the midpoints between neighboring centers; in MLT the
! grid is periodic over 24 hours.  Each output cell gets the average of
! the input cells it overlaps, weighted by the overlap area
! dLon*dSin(Lat).  An output cell the input grid never reaches is left
! at zero.
!
! Both grids are separable, MLT varying only with the first index and
! latitude only with the second, which is what lets the overlap sum
! factor into a longitude weight times a latitude weight.  Nothing here
! works without that, so it is checked on entry.
!
! Watch the closing duplicate column.  haveMlts(havenMlts) is a copy of
! haveMlts(1) (grid_routines.f90, set_ie_mlts), and MagLonMC runs
! 0..360 inclusive over nMagLons+1 columns (calc_electrodynamics.f90).
! Count that column as a cell of its own and the last cell comes out
! ~24 h wide, dragging every output value halfway toward it.
!/
! --------------------------------------------------------------------

subroutine remap_high_to_low_conservative(this, fine_mlts, fine_lats, fine_vals, &
                                          coarse_mlts, coarse_lats, coarse_vals)

  implicit none

  class(ieModel) :: this
  real, intent(in)  :: fine_mlts(this%havenMlts, this%havenLats)
  real, intent(in)  :: fine_lats(this%havenMlts, this%havenLats)
  real, intent(in)  :: fine_vals(this%havenMlts, this%havenLats)
  real, intent(in)  :: coarse_mlts(this%neednMlts, this%neednLats)
  real, intent(in)  :: coarse_lats(this%neednMlts, this%neednLats)
  real, intent(out) :: coarse_vals(this%neednMlts, this%neednLats)

  character(len=*), parameter :: NameSub = "remap_high_to_low_conservative"

  real, parameter :: cPi = 3.14159265358979323846
  real, parameter :: cDegToRad = cPi/180.0
  real, parameter :: cMltTol = 1.0e-4    ! hours
  real, parameter :: cLatTol = 1.0e-4    ! degrees, separability check only

  integer :: nfM, nfL, ncM, ncL      ! full sizes, duplicate column included
  integer :: nfMu, ncMu              ! distinct MLT columns
  integer :: fm, fl, cm, cl, iShift

  real, allocatable :: fLonLo(:), fLonHi(:), cLonLo(:), cLonHi(:)
  real, allocatable :: fSinLo(:), fSinHi(:), cSinLo(:), cSinHi(:)
  real, allocatable :: wLon(:, :), wLat(:, :)
  real, allocatable :: lonTot(:), latTot(:)
  real, allocatable :: partial(:, :)
  real, allocatable :: fLatC(:), cLatC(:)   ! contiguous copies of the lat rows

  real :: lo, hi, shift, over

  !--------------------------------------------------------------------------
  coarse_vals = 0.0

  nfM = this%havenMlts
  nfL = this%havenLats
  ncM = this%neednMlts
  ncL = this%neednLats

  if (nfM < 3 .or. nfL < 2 .or. ncM < 3 .or. ncL < 2) &
    call CON_stop(NameSub//": grid too small to build cell edges")

  ! Both grids must be full lat/MLT rectangles: every MLT column covers the
  ! same latitudes, and no cells are masked or missing.  The denominator
  ! below is factored as lonTot*latTot, which is exact only under that
  ! assumption.  A ragged or masked input would still be counted in the
  ! area while contributing nothing to the sum, so the result would come
  ! out diluted toward zero.  A gap that spans every column, of the kind a
  ! hemispheric IE grid leaves, is fine: lat_cell_bounds leaves it
  ! uncovered and latTot goes to zero there.
  call check_separable(nfM, nfL, fine_mlts, fine_lats, "IE (have)")
  call check_separable(ncM, ncL, coarse_mlts, coarse_lats, "UA (need)")

  ! ---- distinct MLT columns (drop a closing duplicate) --------------------
  nfMu = nfM
  if (abs(wrap_signed(fine_mlts(nfM, 1) - fine_mlts(1, 1))) < cMltTol) &
    nfMu = nfM - 1
  ncMu = ncM
  if (abs(wrap_signed(coarse_mlts(ncM, 1) - coarse_mlts(1, 1))) < cMltTol) &
    ncMu = ncM - 1

  allocate (fLonLo(nfMu), fLonHi(nfMu), cLonLo(ncMu), cLonHi(ncMu))
  allocate (fSinLo(nfL), fSinHi(nfL), cSinLo(ncL), cSinHi(ncL))
  allocate (wLon(ncMu, nfMu), wLat(ncL, nfL))
  allocate (lonTot(ncMu), latTot(ncL), partial(ncMu, nfL))
  allocate (fLatC(nfL), cLatC(ncL))

  fLatC = fine_lats(1, 1:nfL)
  cLatC = coarse_lats(1, 1:ncL)

  call mlt_cell_bounds(nfMu, fine_mlts(1:nfMu, 1), fLonLo, fLonHi)
  call mlt_cell_bounds(ncMu, coarse_mlts(1:ncMu, 1), cLonLo, cLonHi)
  call lat_cell_bounds(nfL, fLatC, fSinLo, fSinHi)
  call lat_cell_bounds(ncL, cLatC, cSinLo, cSinHi)

  ! ---- longitude overlap, in hours ----------------------------------------
  ! The two bound sets each span exactly 24 h but start at different
  ! offsets, so a cell can overlap its partner in the window below, at or
  ! above.  Hence three shifts.  They cannot double count, because
  ! mlt_cell_bounds guarantees no cell is wider than 24 h.
  wLon = 0.0
  do cm = 1, ncMu
    do fm = 1, nfMu
      over = 0.0
      do iShift = -1, 1
        shift = real(iShift)*24.0
        lo = max(fLonLo(fm) + shift, cLonLo(cm))
        hi = min(fLonHi(fm) + shift, cLonHi(cm))
        if (hi > lo) over = over + (hi - lo)
      enddo
      wLon(cm, fm) = over
    enddo
    lonTot(cm) = sum(wLon(cm, 1:nfMu))
  enddo

  ! ---- latitude overlap, in sin(lat) --------------------------------------
  wLat = 0.0
  do cl = 1, ncL
    do fl = 1, nfL
      lo = max(fSinLo(fl), cSinLo(cl))
      hi = min(fSinHi(fl), cSinHi(cl))
      if (hi > lo) wLat(cl, fl) = hi - lo
    enddo
    latTot(cl) = sum(wLat(cl, 1:nfL))
  enddo

  ! ---- apply: area weight factors, so sum over fm then over fl ------------
  do cm = 1, ncMu
    do fl = 1, nfL
      partial(cm, fl) = sum(wLon(cm, 1:nfMu)*fine_vals(1:nfMu, fl))
    enddo
  enddo

  ! Geometric overlap is the only domain test: an output cell the input
  ! grid does not reach has latTot == 0 and stays at the zero set above.
  ! One straddling the edge of the input grid keeps the average of the
  ! part that is covered, rather than being discarded.
  do cl = 1, ncL
    if (latTot(cl) <= 0.0) cycle
    do cm = 1, ncMu
      if (lonTot(cm) <= 0.0) cycle
      coarse_vals(cm, cl) = sum(wLat(cl, 1:nfL)*partial(cm, 1:nfL)) &
                            /(lonTot(cm)*latTot(cl))
    enddo
  enddo

  ! ---- refill the output's own closing duplicate column -------------------
  if (ncMu < ncM) coarse_vals(ncM, 1:ncL) = coarse_vals(1, 1:ncL)

  deallocate (fLonLo, fLonHi, cLonLo, cLonHi)
  deallocate (fSinLo, fSinHi, cSinLo, cSinHi)
  deallocate (wLon, wLat, lonTot, latTot, partial, fLatC, cLatC)

contains

  ! Signed difference of two MLTs, wrapped into (-12, 12]
  real function wrap_signed(dIn)
    real, intent(in) :: dIn
    wrap_signed = dIn - 24.0*nint(dIn/24.0)
  end function wrap_signed

  ! Forward difference of two MLTs, wrapped into (0, 24]
  real function wrap_forward(dIn)
    real, intent(in) :: dIn
    wrap_forward = dIn - 24.0*floor(dIn/24.0)
    if (wrap_forward <= 0.0) wrap_forward = wrap_forward + 24.0
  end function wrap_forward

  ! Periodic MLT cells from centers: edges at the midpoints, bounds
  ! strictly increasing over exactly 24 h (no mod(), so no cell ever
  ! comes back inverted across the 0/24 seam).
  subroutine mlt_cell_bounds(n, center, boundLo, boundHi)
    integer, intent(in) :: n
    real, intent(in) :: center(n)
    real, intent(out) :: boundLo(n), boundHi(n)
    real :: dFwd(n)
    integer :: i, ip
    do i = 1, n
      ip = i + 1
      if (ip > n) ip = 1
      dFwd(i) = wrap_forward(center(ip) - center(i))
    enddo
    boundLo(1) = center(1) - 0.5*dFwd(n)
    boundHi(1) = boundLo(1) + 0.5*(dFwd(n) + dFwd(1))
    do i = 2, n
      boundLo(i) = boundHi(i - 1)
      boundHi(i) = boundLo(i) + 0.5*(dFwd(i - 1) + dFwd(i))
    enddo
  end subroutine mlt_cell_bounds

  ! Latitude cells from centers, returned as sin(lat) bounds so the
  ! caller can multiply by a longitude width to get an area.  Works
  ! for centers ordered either way.
  !
  ! An interval far wider than its neighbors is a gap in the grid, not a
  ! cell boundary.  set_ie_lats leaves one if the IE grid covers a single
  ! hemisphere.  Splitting it at the midpoint would grow the two rows
  ! beside it across the gap and hand the caller polar-cap values at the
  ! equator, so those rows keep their interior width and the gap is left
  ! uncovered.
  subroutine lat_cell_bounds(n, center, sinLo, sinHi)
    integer, intent(in) :: n
    real, intent(in) :: center(n)
    real, intent(out) :: sinLo(n), sinHi(n)
    real, parameter :: cGapFac = 3.0
    real :: d(n - 1), dEff(n - 1), dm, dp, sgn, nbr, lo, hi
    logical :: IsGap(n - 1)
    integer :: j, jj

    do j = 1, n - 1
      d(j) = center(j + 1) - center(j)
    enddo
    sgn = sign(1.0, d(1))

    IsGap = .false.
    do j = 1, n - 1
      nbr = -1.0
      if (j > 1) nbr = abs(d(j - 1))
      if (j < n - 1) then
        if (nbr < 0.0) then
          nbr = abs(d(j + 1))
        else
          nbr = min(nbr, abs(d(j + 1)))
        endif
      endif
      if (nbr > 0.0 .and. abs(d(j)) > cGapFac*nbr) IsGap(j) = .true.
    enddo

    ! width to use beside a gap: the nearest ordinary spacing
    do j = 1, n - 1
      dEff(j) = abs(d(j))
      if (.not. IsGap(j)) CYCLE
      do jj = 1, n - 1
        if (j - jj >= 1) then
          if (.not. IsGap(j - jj)) then
            dEff(j) = abs(d(j - jj)); EXIT
          endif
        endif
        if (j + jj <= n - 1) then
          if (.not. IsGap(j + jj)) then
            dEff(j) = abs(d(j + jj)); EXIT
          endif
        endif
      enddo
    enddo

    do j = 1, n
      if (j > 1) then
        dm = dEff(j - 1)
      else
        dm = dEff(1)
      endif
      if (j < n) then
        dp = dEff(j)
      else
        dp = dEff(n - 1)
      endif
      lo = center(j) - sgn*0.5*dm
      hi = center(j) + sgn*0.5*dp
      lo = max(-90.0, min(90.0, lo))
      hi = max(-90.0, min(90.0, hi))
      sinLo(j) = sin(min(lo, hi)*cDegToRad)
      sinHi(j) = sin(max(lo, hi)*cDegToRad)
    enddo
  end subroutine lat_cell_bounds

  ! MLT must not vary down a column, latitude must not vary across a
  ! row.  Both hold for the IE grid (set_ie_mlts / set_ie_lats) and for
  ! the GITM magnetic grid (MagLonMC depends on i only, MagLatMC on j
  ! only).  If that ever stops being true, every value this routine
  ! returns is wrong, so stop rather than log it.
  subroutine check_separable(nM, nL, mlts, lats, cLabel)
    integer, intent(in) :: nM, nL
    real, intent(in) :: mlts(nM, nL), lats(nM, nL)
    character(len=*), intent(in) :: cLabel
    integer :: i, j
    do j = 2, nL
      do i = 1, nM
        if (abs(wrap_signed(mlts(i, j) - mlts(i, 1))) > cMltTol) &
          call CON_stop(NameSub//": "//cLabel// &
                        " MLT varies with latitude index; remap assumes it does not")
      enddo
    enddo
    do i = 2, nM
      do j = 1, nL
        if (abs(lats(i, j) - lats(1, j)) > cLatTol) &
          call CON_stop(NameSub//": "//cLabel// &
                        " latitude varies with MLT index; remap assumes it does not")
      enddo
    enddo
  end subroutine check_separable

end subroutine remap_high_to_low_conservative
