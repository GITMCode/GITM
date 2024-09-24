! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

subroutine check_for_nans_ions(cMarker)

  use ModSizeGitm
  use ModGITM
  use ModPlanet
  use ieee_arithmetic
  
  implicit none

  character(LEN=*), intent(in) :: cMarker
  integer :: iLon, iLat, iAlt, iIon
  logical :: IsFound

  IsFound = .false.
  do iLon = -1, nLons + 2
    do iLat = -1, nLats + 2
      do iAlt = -1, nAlts + 2
        do iIon = 1, nIons
          if (ieee_is_nan(iDensityS(iLon, iLat, iAlt, iIon, 1))) then
            write(*, *) 'Nan found in iDensityS : '
            write(*, *) cMarker
            write(*, *) iLon, iLat, iAlt, iProc, iIon
            IsFound = .true.
          end if
          if (iDensityS(iLon, iLat, iAlt, iIon, 1) < 0.0) then
            write(*, *) 'Negative density found in iDensityS : '
            write(*, *) cMarker
            write(*, *) iLon, iLat, iAlt, iProc, iIon
            IsFound = .true.
          end if
        end do
      end do
    end do
  end do

  if (IsFound) call stop_gitm("Stopping...")

end subroutine check_for_nans_ions

subroutine check_for_nans_neutrals(cMarker)

  use ModSizeGitm
  use ModGITM
  use ModPlanet
  use ieee_arithmetic

  implicit none

  character(LEN=*), intent(in) :: cMarker
  integer :: iLon, iLat, iAlt, iNeu
  logical :: IsFound

  IsFound = .false.
  do iLon = -1, nLons + 2
    do iLat = -1, nLats + 2
      do iAlt = -1, nAlts + 2
        do iNeu = 1, nSpecies
          if (ieee_is_nan(nDensityS(iLon, iLat, iAlt, iNeu, 1))) then
            write(*, *) 'Nan found in nDensityS : '
            write(*, *) cMarker
            write(*, *) iLon, iLat, iAlt, iProc, iNeu
            IsFound = .true.
          end if
          if (nDensityS(iLon, iLat, iAlt, iNeu, 1) < 0.0) then
            write(*, *) 'Negative density found in nDensityS : '
            write(*, *) cMarker
            write(*, *) iLon, iLat, iAlt, iProc, iNeu
            IsFound = .true.
          end if
        end do
      end do
    end do
  end do

  if (IsFound) call stop_gitm("Stopping...")

end subroutine check_for_nans_neutrals

subroutine check_for_nans_temps(cMarker)

  use ModSizeGitm
  use ModGITM
  use ModPlanet
  use ieee_arithmetic

  implicit none

  character(LEN=*), intent(in) :: cMarker
  integer :: iLon, iLat, iAlt, iNeu
  logical :: IsFound

  IsFound = .false.
  do iLon = -1, nLons + 2
    do iLat = -1, nLats + 2
      do iAlt = -1, nAlts + 2
        if (ieee_is_nan(Temperature(iLon, iLat, iAlt, 1))) then
          write(*, *) 'Nan found in Temperature : '
          write(*, *) cMarker
          write(*, *) iLon, iLat, iAlt, iProc
          IsFound = .true.
        end if
        if (ieee_is_nan(iTemperature(iLon, iLat, iAlt, 1))) then
          write(*, *) 'Nan found in iTemperature : '
          write(*, *) cMarker
          write(*, *) iLon, iLat, iAlt, iProc
          IsFound = .true.
        end if
        if (ieee_is_nan(eTemperature(iLon, iLat, iAlt, 1))) then
          write(*, *) 'Nan found in eTemperature : '
          write(*, *) cMarker
          write(*, *) iLon, iLat, iAlt, iProc
          IsFound = .true.
        end if
        ! Check for negative Temperatures:
        if (Temperature(iLon, iLat, iAlt, 1) < 0.0) then
          write(*, *) 'Negative found in Temperature : '
          write(*, *) cMarker
          write(*, *) iLon, iLat, iAlt, iProc
          IsFound = .true.
        end if
        if (iTemperature(iLon, iLat, iAlt, 1) < 0.0) then
          write(*, *) 'Negative found in iTemperature : '
          write(*, *) cMarker
          write(*, *) iLon, iLat, iAlt, iProc
          IsFound = .true.
        end if
        if (eTemperature(iLon, iLat, iAlt, 1) < 0.0) then
          write(*, *) 'Negative found in eTemperature : '
          write(*, *) cMarker
          write(*, *) iLon, iLat, iAlt, iProc
          IsFound = .true.
        end if
      end do
    end do
  end do

  if (IsFound) call stop_gitm("Stopping...")

end subroutine check_for_nans_temps
