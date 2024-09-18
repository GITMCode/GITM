! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

subroutine set_bcs

  use ModGITM
  use ModInputs

  implicit none

  integer :: iBlock, iLon, iLat, iSpecies

  do iBlock = 1, nBlocks

    ! Bottom
    Velocity(:, :, -1, iUp_, iBlock) = -Velocity(:, :, 2, iUp_, iBlock)
    Velocity(:, :, 0, iUp_, iBlock) = -Velocity(:, :, 1, iUp_, iBlock)

    ! Fixed LogRho and Temp

    ! if we don't touch LogRho(:,:,-1:0,iBlock) it will never change

    ! Top

    do iLon = 1, nLons
      do iLat = 1, nLats

        if (Velocity(iLon, iLat, nAlts, iUp_, iBlock) > 0.) then

          Velocity(iLon, iLat, nAlts + 1, iUp_, iBlock) = &
            Velocity(iLon, iLat, nAlts, iUp_, iBlock)*0.0
          Velocity(iLon, iLat, nAlts + 2, iUp_, iBlock) = &
            Velocity(iLon, iLat, nAlts, iUp_, iBlock)*0.0

        else

          Velocity(iLon, iLat, nAlts + 1, iUp_, iBlock) = 0.0 ! -Vel(nAlts)
          Velocity(iLon, iLat, nAlts + 2, iUp_, iBlock) = 0.0 ! -Vel(nAlts-1)

        end if
      end do
    end do

!     Temperature(:,:,nAlts+1,iBlock) = TempMax/TempUnit
!     Temperature(:,:,nAlts+2,iBlock) = TempMax/TempUnit

    Temperature(:, :, nAlts + 1, iBlock) = Temperature(:, :, nAlts, iBlock)
    Temperature(:, :, nAlts + 2, iBlock) = Temperature(:, :, nAlts, iBlock)

    do iSpecies = 1, nSpecies
      LogNS(:, :, nAlts + 1, iSpecies, iBlock) = &
        LogNS(:, :, nAlts, iSpecies, iBlock) + &
        dAlt(nAlts)*Gravity(nAlts)/Temperature(:, :, nAlts, iBlock)
      LogNS(:, :, nAlts + 1, iSpecies, iBlock) = &
        LogNS(:, :, nAlts, iSpecies, iBlock) + &
        2*dAlt(nAlts)*Gravity(nAlts)/Temperature(:, :, nAlts, iBlock)
    end do

    LogRho(:, :, nAlts + 1, iBlock) = &
      LogRho(:, :, nAlts, iBlock) + &
      dAlt(nAlts)*Gravity(nAlts)/Temperature(:, :, nAlts, iBlock)
    LogRho(:, :, nAlts + 1, iBlock) = &
      LogRho(:, :, nAlts, iBlock) + &
      2*dAlt(nAlts)*Gravity(nAlts)/Temperature(:, :, nAlts, iBlock)

  end do

end subroutine set_bcs
