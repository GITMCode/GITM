! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

subroutine calc_electron_temperature(iBlock)

  !  Take Tion and Telec values from empirical datasets from Viking (settei.F)
  !  -- Pre-MAVEN
  !  Take Telec values from Ergun et al (2015) plus Linear Fit to Tn at 130 km
  !  -- Post-MAVEN  (dayside low SZA experiment for dayside orbit application)

  use ModGITM
  use ModPlanet, only: ialtminiono

  implicit none

  integer, intent(in) :: iBlock
  integer :: iLon, iLat, iAlt, iminiono, k130
  real :: Alt, TN130
! real,parameter :: TL = 510.
  real, parameter :: TL2 = 510.
  real, parameter :: TH = 3140.
  real, parameter :: Z0 = 241.
  real, parameter :: H0 = 60.
  real, parameter :: TE180 = 810.

  call report("Electron Density", 2)
  !Electron Temperature
  !FOX[93] FORMULATION FOR DY TE FROM ROHRBAUGH ET. AL. [79]

! do iLon = 1, nLons
!    do iLat = 1, nLats
!       iMinIono = iAltMinIono(iLon,iLat,iBlock)
!       do ialt = iminiono, nAlts
!
!          Alt = Altitude_GB(iLon,iLat,iAlt,iBlock) /1000.0
!          if (Alt < 130.0) eTemperature(iLon,iLat,iAlt,iBlock) = &
!               Temperature(iLon,iLat,iAlt,iBlock) * TempUnit(iLon,iLat,iAlt)
!
!          if (Alt > 180.0) eTemperature(iLon,iLat,iAlt,iBlock) = &
!               4200.0 - 3750.0*exp((180-Alt)/89.6)
!
!          if (Alt >= 130.0 .and. Alt <= 180.0) eTemperature(iLon,iLat,iAlt,iBlock) = &
!               700.0-536.0*exp((130.0 - Alt)/65.4)
!
!       enddo
!    enddo
! enddo

  !Electron Temperature
  !Bougher[2017] FORMULATION FOR low SZA Dayside TE:
  !    Combo of Ergun ea (2015) tanh above 180 km and Linear Fit to Tn(130 km) below.

  do iLon = 1, nLons
    do iLat = 1, nLats
      iMinIono = iAltMinIono(iLon, iLat, iBlock)
      do ialt = iminiono, nAlts

        Alt = Altitude_GB(iLon, iLat, iAlt, iBlock)/1000.0
        if (Alt < 130.0) Then
          eTemperature(iLon, iLat, iAlt, iBlock) = &
            Temperature(iLon, iLat, iAlt, iBlock)*TempUnit(iLon, iLat, iAlt)
          k130 = ialt
        end if

        if (Alt >= 130.0 .and. Alt <= 180.0) Then
          TN130 = Temperature(iLon, iLat, k130, iBlock)*TempUnit(iLon, iLat, k130)
          eTemperature(iLon, iLat, iAlt, iBlock) = TN130 + (TE180 - TN130)*(Alt - 130.)/50.
        end if

        !        if (Alt > 180.0) eTemperature(iLon,iLat,iAlt,iBlock) = &
        !           0.5*(TH+TL) + 0.5*(TH-TL)*tanh((Alt-Z0)/H0)
        if (Alt > 180.0) eTemperature(iLon, iLat, iAlt, iBlock) = &
          0.5*(TH + TL2) + 0.5*(TH - TL2)*tanh((Alt - Z0)/H0)

      end do
    end do
  end do

  !Ion Temperature
  !FOX[93] FORMULATION FOR DY TI FROM ROHRBAUGH ET. AL. [79]

  do iLon = 1, nLons
    do iLat = 1, nLats
      iMinIono = iAltMinIono(iLon, iLat, iBlock)
      do ialt = iminiono, nAlts

        Alt = Altitude_GB(iLon, iLat, iAlt, iBlock)/1000.0
        if (Alt < 180.0) iTemperature(iLon, iLat, iAlt, iBlock) = &
          Temperature(iLon, iLat, iAlt, iBlock)*TempUnit(iLon, iLat, iAlt)

        if (Alt > 300.0) iTemperature(iLon, iLat, iAlt, iBlock) = &
          eTemperature(iLon, iLat, iAlt, iBlock)

        if (Alt >= 180.0 .and. Alt <= 300.0) iTemperature(iLon, iLat, iAlt, iBlock) = &
          10.0**(2.243 + (Alt - 180.0)/95.0)

      end do
    end do
  end do

end subroutine calc_electron_temperature

