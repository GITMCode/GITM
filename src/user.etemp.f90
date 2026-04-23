! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

subroutine user_create_perturbation

  use ModGITM
  use ModInputs

  NDensityS(nLons/2, nLats/2, 2, 1:3, 1) = NDensityS(nLons/2, nLats/2, 2, 1:3, 1)*50.0
  temperature(nLons/2, nLats/2, 2, 1) = temperature(nLons/2, nLats/2, 2, 1)*100.0

end subroutine user_create_perturbation

subroutine user_perturbation

  use ModGITM
  use ModInputs
  use ModNumConst
  use ModKind, ONLY: Real8_
  use ModTime
  use ModSources, only: UserHeatingRate

  implicit none

  real         :: latcenter, loncenter, amp
  real         :: latwidth, lonwidth
  real         :: f(-1:nLons + 2, -1:nLats + 2)
  integer      :: iBlock, iSpecies, iLat, iLon, iAlt
  real(Real8_) :: PerturbTimeStart, PerturbTimeEnd, MidTime
  real         :: tsave = 0.0
  real         :: dla, dlo, lac, loc

  UserHeatingRate = 0.0

  ! Start 1 hour1 after start of the simulation
  PerturbTimeStart = StartTime + 1.0*3600.0 + 60.0
  ! End 1 minute after that
  PerturbTimeEnd = PerturbTimeStart + 60.0

  dla = latitude(2, 1) - latitude(1, 1)
  lac = (LatEnd + LatStart)/2.0
  dlo = longitude(2, 1) - longitude(1, 1)
  loc = (LonEnd + LonStart)/2.0

  latcenter = lac + dla/2.0
  loncenter = loc + dlo/2.0
  latwidth = dla/2.0
  lonwidth = dlo/2.0

  if (CurrentTime < PerturbTimeStart) then
    tsave = sum(temperature(1:nLons, 1:nLats, 1, 1:nBlocks))/ &
            (nLons*nLats*nBlocks)
  endif

  DuringPerturb = .false.

  if (CurrentTime >= PerturbTimeStart .and. &
      CurrentTime < PerturbTimeEnd) then

    MidTime = (PerturbTimeStart + PerturbTimeEnd)/2.0

    DuringPerturb = .true.

    amp = exp(-((CurrentTime - MidTime)/(PerturbTimeEnd - PerturbTimeStart)*5)**2)
!     write(*,*) "Perturbing!", tsave, amp,latcenter/3.1415*180.0,&
!          loncenter/3.1415*180.0, latwidth, lonwidth

    do iBlock = 1, nBlocks
      do iLon = 1, nLons
        do iLat = 1, nLats
          if ((abs(latitude(iLat, iBlock) - latcenter) < 4*latwidth) .and. &
              (abs(longitude(iLon, iBlock) - loncenter) < 4*lonwidth)) then
            f(iLon, iLat) = & !amp*&
              exp(-((latitude(iLat, iBlock) - latcenter)/latwidth)**2)* &
              exp(-((longitude(iLon, iBlock) - loncenter)/lonwidth)**2)
          else
            f(iLon, iLat) = 0.0
          endif

          iAlt = 1
          UserHeatingRate(iLon, iLat, iAlt, iBlock) = &
            1000.0*4.184e6*f(iLon, iLat)/ &
            cellvolume(iLon, iLat, iAlt, iBlock)/ &
            TempUnit(iLon, iLat, iAlt)/ &
            cp(iLon, iLat, iAlt, iBlock)/ &
            rho(iLon, iLat, iAlt, iBlock)
        enddo
      enddo

!        temperature(:,:,-1,iBlock) = tsave + 5.0 * f * tsave
!        temperature(:,:, 0,iBlock) = tsave + 5.0 * f * tsave
!        temperature(:,:, 1,iBlock) = tsave + 5.0 * f * tsave
!        do iSpecies = 1, nSpecies
!           VerticalVelocity(:, :, -1, iSpecies, iBlock) = 2*f
!           VerticalVelocity(:, :,  0, iSpecies, iBlock) = 2*f
!           VerticalVelocity(:, :,  1, iSpecies, iBlock) = 2*f
!        enddo
!        Velocity(:,:,-1, iUp_, iBlock) = 2*f
!        Velocity(:,:, 0, iUp_, iBlock) = 2*f
!        Velocity(:,:, 1, iUp_, iBlock) = 2*f
    enddo

  endif

!  NDensityS(nLons/2,nLats/2,2,1:3,1) = NDensityS(nLons/2,nLats/2,2,1:3,1)*50.0
!  temperature(nLons/2,nLats/2,2,1) = temperature(nLons/2,nLats/2,2,1)*100.0

end subroutine user_perturbation



