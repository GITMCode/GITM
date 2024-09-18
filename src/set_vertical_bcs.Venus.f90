! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

!\
! ------------------------------------------------------------
! set_boundary
! ------------------------------------------------------------
!/

subroutine set_vertical_bcs(LogRho, LogNS, Vel_GD, Temp, LogINS, iVel, VertVel)

  ! Fill in ghost cells at the top and bottom

  use ModSizeGitm, only: nAlts
  use ModPlanet, only: nSpecies, Mass, nIons, &
                       IsEarth, iN2_, iNO_, iN4S_, iO_
  use ModGITM, only: TempUnit, iEast_, iNorth_, iUp_
  use ModInputs
  use ModConstants
  use ModTime, only: UTime, iJulianDay, currenttime
  use ModVertical, only: Lat, Lon, Gravity_G, Altitude_G, dAlt_F

  use EUA_ModMsis90, ONLY: meter6

  implicit none

  real, intent(inout) :: &
    LogRho(-1:nAlts + 2), &
    LogNS(-1:nAlts + 2, nSpecies), &
    LogINS(-1:nAlts + 2, nIons), &
    Vel_GD(-1:nAlts + 2, 3), &
    IVel(-1:nAlts + 2, 3), &
    Temp(-1:nAlts + 2), &
    VertVel(-1:nAlts + 2, nSpecies)

  integer :: iSpecies, iAlt
  real    :: InvScaleHeightS, InvScaleHgt, Alt, Lst, Ap = 4.0, dn, dt
  logical :: IsFirstTime = .true., UseMsisBCs = .false.
  real    :: HP
  integer :: ierror

  integer, dimension(25) :: sw

  !-----------------------------------------------------------
  ! Bottom
  !-----------------------------------------------------------

  ! Zero winds at the boundary:

  Vel_GD(-1:0, iEast_) = 0.0
  Vel_GD(-1:0, iNorth_) = 0.0

  Vel_GD(-1:0, iUp_) = 0.0
  VertVel(-1:0, :) = 0.0
  IVel(-1:0, iUp_) = 0.0

  ! Constant gradient at bottom of atmosphere for ions:
  do iSpecies = 1, nIons - 1
    ! dn = (LogINS(2,iSpecies) - LogINS(1,iSpecies))
    ! Set a zero gradient at the bottom
    dn = 0.0
    LogINS(0, iSpecies) = LogINS(1, iSpecies) - dn
    LogINS(-1, iSpecies) = LogINS(0, iSpecies) - dn
  end do

  LogINS(0, nIons) = alog(sum(exp(LogINS(0, 1:nIons - 1))))
  LogINS(1, nIons) = alog(sum(exp(LogINS(1, 1:nIons - 1))))

  if (nSpecies >= iN4S_) then

    dn = (LogNS(2, iN4S_) - LogNS(1, iN4S_))
    if (dn >= 0) then
      LogNS(0, iN4S_) = LogNS(1, iN4S_) - dn
      LogNS(-1, iN4S_) = LogNS(0, iN4S_) - dn
    else
      LogNS(0, iN4S_) = LogNS(1, iN4S_) + dn
      LogNS(-1, iN4S_) = LogNS(0, iN4S_) + dn
    end if

    if (VertVel(1, iN4S_) .gt. 0.0) then
      ! Don't allow upwelling of N4S
      VertVel(0, iN4S_) = -1.0*VertVel(1, iN4S_)
      VertVel(-1, iN4S_) = -1.0*VertVel(2, iN4S_)
    else
      VertVel(0, iN4S_) = 1.0*VertVel(1, iN4S_)
      VertVel(-1, iN4S_) = 1.0*VertVel(1, iN4S_)
    end if

  end if

  ! Lower boundary for O on Mars
  if (nSpecies >= iO_) then

    dn = (LogNS(2, iO_) - LogNS(1, iO_))
    if (dn >= 0) then
      LogNS(0, iO_) = LogNS(1, iO_) - dn
      LogNS(-1, iO_) = LogNS(0, iO_) - dn
    else
      LogNS(0, iO_) = LogNS(1, iO_) + dn
      LogNS(-1, iO_) = LogNS(0, iO_) + dn
    end if

    if (VertVel(1, iO_) .gt. 0.0) then
      ! Don't allow upwelling of O
      VertVel(0, iO_) = -1.0*VertVel(1, iO_)
      VertVel(-1, iO_) = -1.0*VertVel(2, iO_)
    else
      VertVel(0, iO_) = 1.0*VertVel(1, iO_)
      VertVel(-1, iO_) = 1.0*VertVel(1, iO_)
    end if

  end if

  !-----------------------------------------------------------
  ! Top
  !-----------------------------------------------------------

  ! Slip flow at the top

  Vel_GD(nAlts + 1:nAlts + 2, iEast_) = Vel_GD(nAlts, iEast_)
  Vel_GD(nAlts + 1:nAlts + 2, iNorth_) = Vel_GD(nAlts, iNorth_)

  IVel(nAlts + 1:nAlts + 2, iEast_) = IVel(nAlts, iEast_)
  IVel(nAlts + 1:nAlts + 2, iNorth_) = IVel(nAlts, iNorth_)

  ! Things can go up or down in the ions

  IVel(nAlts + 1, iUp_) = IVel(nAlts, iUp_)
  IVel(nAlts + 2, iUp_) = IVel(nAlts - 1, iUp_)

  ! We only let stuff flow out in the neutrals

  if (Vel_GD(nAlts, iUp_) > 0.) then
    Vel_GD(nAlts + 1:nAlts + 2, iUp_) = Vel_GD(nAlts, iUp_)
    VertVel(nAlts + 1, :) = VertVel(nAlts, :)
    VertVel(nAlts + 2, :) = VertVel(nAlts, :)
  else
    ! Vel_GD(nAlts+1:nAlts+2,iUp_) = 0.0 ! -Vel(nAlts)
    Vel_GD(nAlts + 1, iUp_) = -Vel_GD(nAlts, iUp_)
    Vel_GD(nAlts + 2, iUp_) = -Vel_GD(nAlts - 1, iUp_)
    VertVel(nAlts + 1, :) = -VertVel(nAlts, :)
    VertVel(nAlts + 2, :) = -VertVel(nAlts - 1, :)
  end if

  ! Constant temperature (zero gradient)

  Temp(nAlts + 1) = Temp(nAlts)
  Temp(nAlts + 2) = Temp(nAlts)

  dn = (LogRho(nAlts) - LogRho(nAlts - 1))
  LogRho(nAlts + 1) = LogRho(nAlts) + dn
  LogRho(nAlts + 2) = LogRho(nAlts + 1) + dn

  ! Limit the slope of the ion density

  do iSpecies = 1, nIons - 1
    dn = (exp(LogINS(nAlts, iSpecies)) - exp(LogINS(nAlts - 1, iSpecies)))
    if (dn > 0) dn = -0.01*exp(LogINS(nAlts, iSpecies))
    LogINS(nAlts + 1, iSpecies) = alog(exp(LogINS(nAlts, iSpecies)) + dn)
    LogINS(nAlts + 2, iSpecies) = alog(exp(LogINS(nAlts + 1, iSpecies)) + dn)
  end do

  LogINS(nAlts + 1, nIons) = alog(sum(exp(LogINS(nAlts + 1, 1:nIons - 1))))
  LogINS(nAlts + 2, nIons) = alog(sum(exp(LogINS(nAlts + 2, 1:nIons - 1))))

  ! Hydrostatic pressure for the neutrals

  do iSpecies = 1, nSpecies
    do iAlt = nAlts + 1, nAlts + 2
      InvScaleHeightS = -Gravity_G(iAlt)* &
                        Mass(iSpecies)/(Temp(iAlt)*Boltzmanns_Constant)
      LogNS(iAlt, iSpecies) = &
        LogNS(iAlt - 1, iSpecies) - dAlt_F(iAlt)*InvScaleHeightS
      if (LogNS(nAlts + 1, iSpecies) > 75.0 .or. &
          LogNS(nAlts + 2, iSpecies) > 75.0) then
        write(*, *) "======> bcs : ", iSpecies, 1.0e-3/InvScaleHeightS, &
          Gravity_G(nAlts), Mass(iSpecies), Temp(nAlts), &
          LogNS(nAlts, iSpecies), LogNS(nAlts + 1, iSpecies), &
          dAlt_F(nAlts), LogNS(nAlts + 2, iSpecies)
      end if
    end do
  end do

end subroutine set_vertical_bcs

