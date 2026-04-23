! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE
!
! ModOutputGather.f90
!
! Data gathering subroutines for GITM output. Each output type has a
! gather_XXXX routine that fills a contiguous buffer from the model
! state arrays. The buffer is then handed to the I/O backend for writing.
!
! Buffer layout: buffer(nVars, nX, nY, nZ) with nVars as leading dimension.

module ModOutputGather

  implicit none

contains

  ! ==================================================================
  ! Dispatcher: call the correct gather routine based on cType
  ! ==================================================================
  subroutine gather_output(cType, iBlock, buffer, nV, nX, nY, nZ, &
                           iiLon, iiLat, iiAlt, rLon, rLat, rAlt)
    character(len=5), intent(in) :: cType
    integer, intent(in) :: iBlock, nV, nX, nY, nZ
    real, intent(out) :: buffer(nV, nX, nY, nZ)
    integer, intent(in) :: iiLon, iiLat, iiAlt
    real, intent(in) :: rLon, rLat, rAlt

    select case (cType)
    case ('3DALL'); call gather_3dall(iBlock, buffer, nV, nX, nY, nZ)
    case ('3DLST'); call gather_3dlst(iBlock, buffer, nV, nX, nY, nZ)
    case ('3DNEU'); call gather_3dneu(iBlock, buffer, nV, nX, nY, nZ)
    case ('3DION'); call gather_3dion(iBlock, buffer, nV, nX, nY, nZ)
    case ('3DTHM'); call gather_3dthm(iBlock, buffer, nV, nX, nY, nZ)
    case ('3DCHM'); call gather_3dchm(iBlock, buffer, nV, nX, nY, nZ)
    case ('3DGLO'); call gather_3dglo(iBlock, buffer, nV, nX, nY, nZ)
    case ('3DMAG'); call gather_3dmag(iBlock, buffer, nV, nX, nY, nZ)
    case ('3DHME'); call gather_3dhme(iBlock, buffer, nV, nX, nY, nZ)
    case ('3DMOH'); call gather_3dmoh(iBlock, buffer, nV, nX, nY, nZ)
    case ('3DMOV'); call gather_3dmov(iBlock, buffer, nV, nX, nY, nZ)
    case ('2DGEL'); call gather_2dgel(iBlock, buffer, nV, nX, nY, nZ)
    case ('2DMEL'); call gather_2dmel(iBlock, buffer, nV, nX, nY, nZ)
    case ('2DTEC'); call gather_2dtec(iBlock, buffer, nV, nX, nY, nZ)
    case ('2DANC'); call gather_2danc(iBlock, buffer, nV, nX, nY, nZ)
    case ('2DHME'); call gather_2dhme(iBlock, buffer, nV, nX, nY, nZ)
    case ('1DALL')
      call gather_1dall(iBlock, buffer, nV, nZ, iiLon, iiLat, rLon, rLat)
    case ('1DGLO'); call gather_1dglo(iBlock, buffer, nV, nZ)
    case ('1DTHM'); call gather_1dthm(iBlock, buffer, nV, nZ)
    case ('1DCHM'); call gather_1dchm(iBlock, buffer, nV, nZ)
    case ('1DNEW')
      call gather_1dnew(iBlock, buffer, nV, nZ, iiLon, iiLat, rLon, rLat)
    case ('0DALL')
      call gather_0dall(iBlock, buffer, nV, iiLon, iiLat, iiAlt, rLon, rLat, rAlt)
    case ('3DUSR')
      call gather_3dusr(iBlock, buffer, nV, nX, nY, nZ)
    case ('2DUSR')
      call gather_2dusr(iBlock, buffer, nV, nX, nY, nZ)
    case ('1DUSR')
      call gather_1dusr(iBlock, buffer, nV, nZ, iiLon, iiLat, rLon, rLat)
    case ('0DUSR')
      call gather_0dusr(iBlock, buffer, nV, iiLon, iiLat, iiAlt, rLon, rLat, rAlt)
    case default
      ! Unhandled types: buffer stays unset
      buffer = 0.0
    end select
  end subroutine gather_output

  ! ==================================================================
  ! Diagnostic: report gather/registry variable count mismatch
  ! ==================================================================

  subroutine gather_error(routine_name, got, expected)
    character(len=*), intent(in) :: routine_name
    integer, intent(in) :: got, expected
    write(*, *) "FATAL: "//trim(routine_name)//": gathered ", &
      got, " vars but registry expects ", expected
    stop
  end subroutine gather_error

  ! ==================================================================
  ! Shared fill helpers: neutral and ion variable sets
  !
  ! These fill the standard variable set shared by multiple output
  ! types (3DALL, 3DNEU, 3DHME, 1DALL, 1DNEW), ensuring the variable
  ! order is defined in exactly one place per access pattern.
  ! ==================================================================

  ! Fill neutral variables for one 3D grid point (direct array access):
  !   Rho, NDensityS(1:nSpeciesTotal), Temperature*TempUnit,
  !   Velocity(1:3), VerticalVelocity(1:nSpecies)
  subroutine fill_neutral_vars_3d(iv, buffer, nV, nX, nY, nZ, &
                                  jx, jy, jz, iLon, iLat, iAlt, iBlock)
    use ModGITM
    use ModInputs
    integer, intent(inout) :: iv
    integer, intent(in) :: nV, nX, nY, nZ, jx, jy, jz
    integer, intent(in) :: iLon, iLat, iAlt, iBlock
    real, intent(inout) :: buffer(nV, nX, nY, nZ)
    integer :: i

    iv = iv + 1; buffer(iv, jx, jy, jz) = Rho(iLon, iLat, iAlt, iBlock)
    do i = 1, nSpeciesTotal
      iv = iv + 1; buffer(iv, jx, jy, jz) = NDensityS(iLon, iLat, iAlt, i, iBlock)
    enddo
    iv = iv + 1; buffer(iv, jx, jy, jz) = &
      Temperature(iLon, iLat, iAlt, iBlock) * TempUnit(iLon, iLat, iAlt)
    do i = 1, 3
      iv = iv + 1; buffer(iv, jx, jy, jz) = Velocity(iLon, iLat, iAlt, i, iBlock)
    enddo
    do i = 1, nSpecies
      iv = iv + 1; buffer(iv, jx, jy, jz) = VerticalVelocity(iLon, iLat, iAlt, i, iBlock)
    enddo
  end subroutine fill_neutral_vars_3d

  ! Fill ion variables for one 3D grid point (direct array access):
  !   IDensityS(1:nIons), eTemperature, ITemperature, IVelocity(1:3)
  subroutine fill_ion_vars_3d(iv, buffer, nV, nX, nY, nZ, &
                              jx, jy, jz, iLon, iLat, iAlt, iBlock)
    use ModGITM
    use ModInputs
    integer, intent(inout) :: iv
    integer, intent(in) :: nV, nX, nY, nZ, jx, jy, jz
    integer, intent(in) :: iLon, iLat, iAlt, iBlock
    real, intent(inout) :: buffer(nV, nX, nY, nZ)
    integer :: i

    do i = 1, nIons
      iv = iv + 1; buffer(iv, jx, jy, jz) = IDensityS(iLon, iLat, iAlt, i, iBlock)
    enddo
    iv = iv + 1; buffer(iv, jx, jy, jz) = eTemperature(iLon, iLat, iAlt, iBlock)
    iv = iv + 1; buffer(iv, jx, jy, jz) = ITemperature(iLon, iLat, iAlt, iBlock)
    do i = 1, 3
      iv = iv + 1; buffer(iv, jx, jy, jz) = IVelocity(iLon, iLat, iAlt, i, iBlock)
    enddo
  end subroutine fill_ion_vars_3d

  ! Fill neutral variables for one altitude using bilinear interpolation.
  ! Same variable set as fill_neutral_vars_3d but interpolated at (rLon, rLat).
  subroutine fill_neutral_vars_interp2d(iv, buffer, nV, nZ, jz, &
                                        iAlt, iBlock, iiLon, iiLat, rLon, rLat)
    use ModGITM
    use ModInputs
    integer, intent(inout) :: iv
    integer, intent(in) :: nV, nZ, jz, iAlt, iBlock, iiLon, iiLat
    real, intent(in) :: rLon, rLat
    real, intent(inout) :: buffer(nV, 1, 1, nZ)
    real :: Tmp(0:nLons + 1, 0:nLats + 1)
    integer :: i

    Tmp = Rho(0:nLons + 1, 0:nLats + 1, iAlt, iBlock)
    iv = iv + 1; buffer(iv, 1, 1, jz) = inter2d(Tmp, iiLon, iiLat, rLon, rLat)

    do i = 1, nSpeciesTotal
      Tmp = NDensityS(0:nLons + 1, 0:nLats + 1, iAlt, i, iBlock)
      iv = iv + 1; buffer(iv, 1, 1, jz) = inter2d(Tmp, iiLon, iiLat, rLon, rLat)
    enddo

    Tmp = Temperature(0:nLons + 1, 0:nLats + 1, iAlt, iBlock) * &
          TempUnit(0:nLons + 1, 0:nLats + 1, iAlt)
    iv = iv + 1; buffer(iv, 1, 1, jz) = inter2d(Tmp, iiLon, iiLat, rLon, rLat)

    do i = 1, 3
      Tmp = Velocity(0:nLons + 1, 0:nLats + 1, iAlt, i, iBlock)
      iv = iv + 1; buffer(iv, 1, 1, jz) = inter2d(Tmp, iiLon, iiLat, rLon, rLat)
    enddo

    do i = 1, nSpecies
      Tmp = VerticalVelocity(0:nLons + 1, 0:nLats + 1, iAlt, i, iBlock)
      iv = iv + 1; buffer(iv, 1, 1, jz) = inter2d(Tmp, iiLon, iiLat, rLon, rLat)
    enddo
  end subroutine fill_neutral_vars_interp2d

  ! Fill ion variables for one altitude using bilinear interpolation.
  ! Same variable set as fill_ion_vars_3d but interpolated at (rLon, rLat).
  subroutine fill_ion_vars_interp2d(iv, buffer, nV, nZ, jz, &
                                    iAlt, iBlock, iiLon, iiLat, rLon, rLat)
    use ModGITM
    use ModInputs
    integer, intent(inout) :: iv
    integer, intent(in) :: nV, nZ, jz, iAlt, iBlock, iiLon, iiLat
    real, intent(in) :: rLon, rLat
    real, intent(inout) :: buffer(nV, 1, 1, nZ)
    real :: Tmp(0:nLons + 1, 0:nLats + 1)
    integer :: i

    do i = 1, nIons
      Tmp = IDensityS(0:nLons + 1, 0:nLats + 1, iAlt, i, iBlock)
      iv = iv + 1; buffer(iv, 1, 1, jz) = inter2d(Tmp, iiLon, iiLat, rLon, rLat)
    enddo

    Tmp = eTemperature(0:nLons + 1, 0:nLats + 1, iAlt, iBlock)
    iv = iv + 1; buffer(iv, 1, 1, jz) = inter2d(Tmp, iiLon, iiLat, rLon, rLat)

    Tmp = iTemperature(0:nLons + 1, 0:nLats + 1, iAlt, iBlock)
    iv = iv + 1; buffer(iv, 1, 1, jz) = inter2d(Tmp, iiLon, iiLat, rLon, rLat)

    Tmp = IVelocity(0:nLons + 1, 0:nLats + 1, iAlt, iEast_, iBlock)
    iv = iv + 1; buffer(iv, 1, 1, jz) = inter2d(Tmp, iiLon, iiLat, rLon, rLat)
    Tmp = IVelocity(0:nLons + 1, 0:nLats + 1, iAlt, iNorth_, iBlock)
    iv = iv + 1; buffer(iv, 1, 1, jz) = inter2d(Tmp, iiLon, iiLat, rLon, rLat)
    Tmp = IVelocity(0:nLons + 1, 0:nLats + 1, iAlt, iUp_, iBlock)
    iv = iv + 1; buffer(iv, 1, 1, jz) = inter2d(Tmp, iiLon, iiLat, rLon, rLat)
  end subroutine fill_ion_vars_interp2d

  ! ==================================================================
  ! 3D output gather routines
  ! ==================================================================

  subroutine gather_3dall(iBlock, buffer, nV, nX, nY, nZ)
    use ModGITM
    use ModInputs
    integer, intent(in) :: iBlock, nV, nX, nY, nZ
    real, intent(out) :: buffer(nV, nX, nY, nZ)
    integer :: iAlt, iLat, iLon, jx, jy, jz, iv

    do iAlt = -1, nAlts + 2
      jz = iAlt + 2
      do iLat = -1, nLats + 2
        jy = iLat + 2
        do iLon = -1, nLons + 2
          jx = iLon + 2
          iv = 0
          iv = iv + 1; buffer(iv, jx, jy, jz) = Longitude(iLon, iBlock)
          iv = iv + 1; buffer(iv, jx, jy, jz) = Latitude(iLat, iBlock)
          iv = iv + 1; buffer(iv, jx, jy, jz) = Altitude_GB(iLon, iLat, iAlt, iBlock)
          call fill_neutral_vars_3d(iv, buffer, nV, nX, nY, nZ, &
                                    jx, jy, jz, iLon, iLat, iAlt, iBlock)
          call fill_ion_vars_3d(iv, buffer, nV, nX, nY, nZ, &
                                jx, jy, jz, iLon, iLat, iAlt, iBlock)
        enddo
      enddo
    enddo
    if (iv /= nV) call gather_error('gather_3dall', iv, nV)
  end subroutine gather_3dall

  subroutine gather_3dlst(iBlock, buffer, nV, nX, nY, nZ)
    use ModGITM
    use ModInputs
    integer, intent(in) :: iBlock, nV, nX, nY, nZ
    real, intent(out) :: buffer(nV, nX, nY, nZ)
    integer :: iAlt, iLat, iLon, jx, jy, jz, iv, i

    do iAlt = -1, nAlts + 2
      jz = iAlt + 2
      do iLat = -1, nLats + 2
        jy = iLat + 2
        do iLon = -1, nLons + 2
          jx = iLon + 2
          iv = 0
          iv = iv + 1; buffer(iv, jx, jy, jz) = Longitude(iLon, iBlock)
          iv = iv + 1; buffer(iv, jx, jy, jz) = Latitude(iLat, iBlock)
          iv = iv + 1; buffer(iv, jx, jy, jz) = Altitude_GB(iLon, iLat, iAlt, iBlock)
          if (iRhoOutputList) then
            iv = iv + 1; buffer(iv, jx, jy, jz) = Rho(iLon, iLat, iAlt, iBlock)
          endif
          do i = 1, nSpeciesTotal
            if (iNeutralDensityOutputList(i)) then
              iv = iv + 1; buffer(iv, jx, jy, jz) = NDensityS(iLon, iLat, iAlt, i, iBlock)
            endif
          enddo
          do i = 1, 3
            if (iNeutralWindOutputList(i)) then
              iv = iv + 1; buffer(iv, jx, jy, jz) = Velocity(iLon, iLat, iAlt, i, iBlock)
            endif
          enddo
          do i = 1, nIons
            if (iIonDensityOutputList(i)) then
              iv = iv + 1; buffer(iv, jx, jy, jz) = IDensityS(iLon, iLat, iAlt, i, iBlock)
            endif
          enddo
          do i = 1, 3
            if (iIonWindOutputList(i)) then
              iv = iv + 1; buffer(iv, jx, jy, jz) = IVelocity(iLon, iLat, iAlt, i, iBlock)
            endif
          enddo
          if (iTemperatureOutputList(1)) then
            iv = iv + 1; buffer(iv, jx, jy, jz) = &
              Temperature(iLon, iLat, iAlt, iBlock) * TempUnit(iLon, iLat, iAlt)
          endif
          if (iTemperatureOutputList(2)) then
            iv = iv + 1; buffer(iv, jx, jy, jz) = ITemperature(iLon, iLat, iAlt, iBlock)
          endif
          if (iTemperatureOutputList(3)) then
            iv = iv + 1; buffer(iv, jx, jy, jz) = eTemperature(iLon, iLat, iAlt, iBlock)
          endif
        enddo
      enddo
    enddo
    if (iv /= nV) call gather_error('gather_3dlst', iv, nV)
  end subroutine gather_3dlst

  subroutine gather_3dneu(iBlock, buffer, nV, nX, nY, nZ)
    use ModGITM
    use ModInputs
    integer, intent(in) :: iBlock, nV, nX, nY, nZ
    real, intent(out) :: buffer(nV, nX, nY, nZ)
    integer :: iAlt, iLat, iLon, jx, jy, jz, iv

    do iAlt = -1, nAlts + 2
      jz = iAlt + 2
      do iLat = -1, nLats + 2
        jy = iLat + 2
        do iLon = -1, nLons + 2
          jx = iLon + 2
          iv = 0
          iv = iv + 1; buffer(iv, jx, jy, jz) = Longitude(iLon, iBlock)
          iv = iv + 1; buffer(iv, jx, jy, jz) = Latitude(iLat, iBlock)
          iv = iv + 1; buffer(iv, jx, jy, jz) = Altitude_GB(iLon, iLat, iAlt, iBlock)
          call fill_neutral_vars_3d(iv, buffer, nV, nX, nY, nZ, &
                                    jx, jy, jz, iLon, iLat, iAlt, iBlock)
        enddo
      enddo
    enddo
    if (iv /= nV) call gather_error('gather_3dneu', iv, nV)
  end subroutine gather_3dneu

  subroutine gather_3dion(iBlock, buffer, nV, nX, nY, nZ)
    use ModGITM
    use ModInputs
    use ModElectrodynamics
    integer, intent(in) :: iBlock, nV, nX, nY, nZ
    real, intent(out) :: buffer(nV, nX, nY, nZ)
    integer :: iAlt, iLat, iLon, jx, jy, jz, iv, i

    do iAlt = -1, nAlts + 2
      jz = iAlt + 2
      do iLat = -1, nLats + 2
        jy = iLat + 2
        do iLon = -1, nLons + 2
          jx = iLon + 2
          iv = 0
          iv = iv + 1; buffer(iv, jx, jy, jz) = Longitude(iLon, iBlock)
          iv = iv + 1; buffer(iv, jx, jy, jz) = Latitude(iLat, iBlock)
          iv = iv + 1; buffer(iv, jx, jy, jz) = Altitude_GB(iLon, iLat, iAlt, iBlock)
          do i = 1, nIons
            iv = iv + 1; buffer(iv, jx, jy, jz) = IDensityS(iLon, iLat, iAlt, i, iBlock)
          enddo
          iv = iv + 1; buffer(iv, jx, jy, jz) = eTemperature(iLon, iLat, iAlt, iBlock)
          iv = iv + 1; buffer(iv, jx, jy, jz) = ITemperature(iLon, iLat, iAlt, iBlock)
          do i = 1, 3
            iv = iv + 1; buffer(iv, jx, jy, jz) = IVelocity(iLon, iLat, iAlt, i, iBlock)
          enddo
          iv = iv + 1; buffer(iv, jx, jy, jz) = ed1(iLon, iLat, iAlt)
          iv = iv + 1; buffer(iv, jx, jy, jz) = ed2(iLon, iLat, iAlt)
          iv = iv + 1; buffer(iv, jx, jy, jz) = je1(iLon, iLat, iAlt)
          iv = iv + 1; buffer(iv, jx, jy, jz) = je2(iLon, iLat, iAlt)
          iv = iv + 1; buffer(iv, jx, jy, jz) = mLatitude(iLon, iLat, iAlt, iBlock)
          iv = iv + 1; buffer(iv, jx, jy, jz) = mLongitude(iLon, iLat, iAlt, iBlock)
          do i = 1, 4
            iv = iv + 1; buffer(iv, jx, jy, jz) = B0(iLon, iLat, iAlt, i, iBlock)
          enddo
          iv = iv + 1; buffer(iv, jx, jy, jz) = Potential(iLon, iLat, iAlt, iBlock)
          do i = 1, 3
            iv = iv + 1; buffer(iv, jx, jy, jz) = EField(iLon, iLat, iAlt, i)
          enddo
          iv = iv + 1; buffer(iv, jx, jy, jz) = sqrt(sum(EField(iLon, iLat, iAlt, :)**2))
          iv = iv + 1; buffer(iv, jx, jy, jz) = Collisions(iLon, iLat, iAlt, iVIN_)
          do i = 1, 3
            iv = iv + 1; buffer(iv, jx, jy, jz) = IonPressureGradient(iLon, iLat, iAlt, i, iBlock)
          enddo
        enddo
      enddo
    enddo
    if (iv /= nV) call gather_error('gather_3dion', iv, nV)
  end subroutine gather_3dion

  subroutine gather_3dthm(iBlock, buffer, nV, nX, nY, nZ)
    use ModGITM
    use ModInputs
    use ModSources
    use ModElectrodynamics
    use ModEUV, only: EuvTotal
    integer, intent(in) :: iBlock, nV, nX, nY, nZ
    real, intent(out) :: buffer(nV, nX, nY, nZ)
    integer :: iAlt, iLat, iLon, jx, jy, jz, iv
    integer :: iiAlt, iiLat, iiLon

    do iAlt = -1, nAlts + 2
      jz = iAlt + 2
      iiAlt = max(min(iAlt, nAlts), 1)
      do iLat = -1, nLats + 2
        jy = iLat + 2
        iiLat = min(max(iLat, 1), nLats)
        do iLon = -1, nLons + 2
          jx = iLon + 2
          iiLon = min(max(iLon, 1), nLons)
          iv = 0
          iv = iv + 1; buffer(iv, jx, jy, jz) = Longitude(iLon, iBlock)
          iv = iv + 1; buffer(iv, jx, jy, jz) = Latitude(iLat, iBlock)
          iv = iv + 1; buffer(iv, jx, jy, jz) = Altitude_GB(iLon, iLat, iAlt, iBlock)
          iv = iv + 1; buffer(iv, jx, jy, jz) = &
            EuvHeating(iiLon, iiLat, iiAlt, iBlock) * TempUnit(iiLon, iiLat, iiAlt)
          iv = iv + 1; buffer(iv, jx, jy, jz) = &
            Conduction(iiLon, iiLat, iiAlt) * TempUnit(iiLon, iiLat, iiAlt) / dt
          iv = iv + 1; buffer(iv, jx, jy, jz) = MoleConduction(iiLon, iiLat, iiAlt)
          iv = iv + 1; buffer(iv, jx, jy, jz) = EddyCond(iiLon, iiLat, iiAlt)
          iv = iv + 1; buffer(iv, jx, jy, jz) = EddyCondAdia(iiLon, iiLat, iiAlt)
          iv = iv + 1; buffer(iv, jx, jy, jz) = &
            ChemicalHeatingRate(iiLon, iiLat, iiAlt) * TempUnit(iiLon, iiLat, iiAlt) / dt
          iv = iv + 1; buffer(iv, jx, jy, jz) = &
            JouleHeating(iiLon, iiLat, iiAlt) * TempUnit(iiLon, iiLat, iiAlt)
          iv = iv + 1; buffer(iv, jx, jy, jz) = &
            -NOCooling(iiLon, iiLat, iiAlt) * TempUnit(iiLon, iiLat, iiAlt)
          iv = iv + 1; buffer(iv, jx, jy, jz) = &
            -OCooling(iiLon, iiLat, iiAlt) * TempUnit(iiLon, iiLat, iiAlt)
          iv = iv + 1; buffer(iv, jx, jy, jz) = EuvTotal(iiLon, iiLat, iiAlt, iBlock)
          iv = iv + 1; buffer(iv, jx, jy, jz) = cp(iiLon, iiLat, iiAlt, iBlock)
          iv = iv + 1; buffer(iv, jx, jy, jz) = Rho(iiLon, iiLat, iiAlt, iBlock)
          iv = iv + 1; buffer(iv, jx, jy, jz) = sqrt(sum(EField(iLon, iLat, iAlt, :)**2))
          iv = iv + 1; buffer(iv, jx, jy, jz) = Sigma_Pedersen(iLon, iLat, iAlt)
          iv = iv + 1; buffer(iv, jx, jy, jz) = &
            AuroralIonRateS(iiLon, iiLat, iiAlt, iO_3P_, iBlock)
          iv = iv + 1; buffer(iv, jx, jy, jz) = &
            AuroralIonRateS(iiLon, iiLat, iiAlt, iO2_, iBlock)
          iv = iv + 1; buffer(iv, jx, jy, jz) = &
            AuroralIonRateS(iiLon, iiLat, iiAlt, iN2_, iBlock)
        enddo
      enddo
    enddo
    if (iv /= nV) call gather_error('gather_3dthm', iv, nV)
  end subroutine gather_3dthm

  subroutine gather_3dchm(iBlock, buffer, nV, nX, nY, nZ)
    use ModGITM
    use ModInputs
    use ModSources
    use ModConstants
    integer, intent(in) :: iBlock, nV, nX, nY, nZ
    real, intent(out) :: buffer(nV, nX, nY, nZ)
    integer :: iAlt, iLat, iLon, jx, jy, jz, iv, iReact
    integer :: iiAlt, iiLat, iiLon

    do iAlt = -1, nAlts + 2
      jz = iAlt + 2
      iiAlt = max(min(iAlt, nAlts), 1)
      do iLat = -1, nLats + 2
        jy = iLat + 2
        iiLat = min(max(iLat, 1), nLats)
        do iLon = -1, nLons + 2
          jx = iLon + 2
          iiLon = min(max(iLon, 1), nLons)
          iv = 0
          iv = iv + 1; buffer(iv, jx, jy, jz) = Longitude(iLon, iBlock)
          iv = iv + 1; buffer(iv, jx, jy, jz) = Latitude(iLat, iBlock)
          iv = iv + 1; buffer(iv, jx, jy, jz) = Altitude_GB(iLon, iLat, iAlt, iBlock)
          do iReact = 1, nReactions
            iv = iv + 1; buffer(iv, jx, jy, jz) = &
              ChemicalHeatingSpecies(iiLon, iiLat, iiAlt, iReact) / Element_Charge
          enddo
          iv = iv + 1; buffer(iv, jx, jy, jz) = &
            ChemicalHeatingRate(iiLon, iiLat, iiAlt) * &
            cp(iiLon, iiLat, iiAlt, iBlock) * &
            Rho(iiLon, iiLat, iiAlt, iBlock) * TempUnit(iiLon, iiLat, iiAlt) / &
            Element_Charge
        enddo
      enddo
    enddo
    if (iv /= nV) call gather_error('gather_3dchm', iv, nV)
  end subroutine gather_3dchm

  subroutine gather_3dglo(iBlock, buffer, nV, nX, nY, nZ)
    use ModGITM
    use ModInputs
    integer, intent(in) :: iBlock, nV, nX, nY, nZ
    real, intent(out) :: buffer(nV, nX, nY, nZ)

    ! STUB: 3DGLO variables are registered in ModOutputRegistry but data
    ! gathering is not yet implemented.  Buffer is zeroed; implement this
    ! subroutine when real glow emission data is available.
    buffer = 0.0
  end subroutine gather_3dglo

  subroutine gather_3dmag(iBlock, buffer, nV, nX, nY, nZ)
    use ModGITM
    use ModInputs
    integer, intent(in) :: iBlock, nV, nX, nY, nZ
    real, intent(out) :: buffer(nV, nX, nY, nZ)
    integer :: iAlt, iLat, iLon, jx, jy, jz, iv, i

    do iAlt = -1, nAlts + 2
      jz = iAlt + 2
      do iLat = -1, nLats + 2
        jy = iLat + 2
        do iLon = -1, nLons + 2
          jx = iLon + 2
          iv = 0
          iv = iv + 1; buffer(iv, jx, jy, jz) = Longitude(iLon, iBlock)
          iv = iv + 1; buffer(iv, jx, jy, jz) = Latitude(iLat, iBlock)
          iv = iv + 1; buffer(iv, jx, jy, jz) = Altitude_GB(iLon, iLat, iAlt, iBlock)
          iv = iv + 1; buffer(iv, jx, jy, jz) = mLatitude(iLon, iLat, iAlt, iBlock)
          iv = iv + 1; buffer(iv, jx, jy, jz) = mLongitude(iLon, iLat, iAlt, iBlock)
          do i = 1, 4
            iv = iv + 1; buffer(iv, jx, jy, jz) = B0(iLon, iLat, iAlt, i, iBlock)
          enddo
        enddo
      enddo
    enddo
    if (iv /= nV) call gather_error('gather_3dmag', iv, nV)
  end subroutine gather_3dmag

  subroutine gather_3dhme(iBlock, buffer, nV, nX, nY, nZ)
    use ModGITM
    use ModInputs
    use ModSources
    integer, intent(in) :: iBlock, nV, nX, nY, nZ
    real, intent(out) :: buffer(nV, nX, nY, nZ)
    integer :: iAlt, iLat, iLon, jx, jy, jz, iv, i, iiAlt

    do iAlt = -1, nAlts + 2
      jz = iAlt + 2
      iiAlt = max(min(iAlt, nAlts), 1)
      do iLat = 1, nLats
        jy = iLat
        do iLon = 1, nLons
          jx = iLon
          iv = 0
          iv = iv + 1; buffer(iv, jx, jy, jz) = Longitude(iLon, iBlock)
          iv = iv + 1; buffer(iv, jx, jy, jz) = Latitude(iLat, iBlock)
          iv = iv + 1; buffer(iv, jx, jy, jz) = Altitude_GB(iLon, iLat, iAlt, iBlock)
          call fill_neutral_vars_3d(iv, buffer, nV, nX, nY, nZ, &
                                    jx, jy, jz, iLon, iLat, iAlt, iBlock)
          call fill_ion_vars_3d(iv, buffer, nV, nX, nY, nZ, &
                                jx, jy, jz, iLon, iLat, iAlt, iBlock)
          ! HME-specific variables (source terms use clamped iiAlt)
          iv = iv + 1; buffer(iv, jx, jy, jz) = &
            PhotoElectronHeating(iLon, iLat, iiAlt, iBlock) * dt * TempUnit(iLon, iLat, iiAlt)
          iv = iv + 1; buffer(iv, jx, jy, jz) = &
            JouleHeating(iLon, iLat, iiAlt) * dt * TempUnit(iLon, iLat, iiAlt)
          iv = iv + 1; buffer(iv, jx, jy, jz) = cp(iLon, iLat, iiAlt, iBlock)
          iv = iv + 1; buffer(iv, jx, jy, jz) = mLatitude(iLon, iLat, iAlt, iBlock)
          iv = iv + 1; buffer(iv, jx, jy, jz) = mLongitude(iLon, iLat, iAlt, iBlock)
          do i = 1, 4
            iv = iv + 1; buffer(iv, jx, jy, jz) = B0(iLon, iLat, iAlt, i, iBlock)
          enddo
          iv = iv + 1; buffer(iv, jx, jy, jz) = Potential(iLon, iLat, iAlt, iBlock)
          iv = iv + 1; buffer(iv, jx, jy, jz) = PotentialY(iLon, iLat, iAlt, iBlock)
          do i = 1, 3
            iv = iv + 1; buffer(iv, jx, jy, jz) = EField(iLon, iLat, iAlt, i)
          enddo
        enddo
      enddo
    enddo
    if (iv /= nV) call gather_error('gather_3dhme', iv, nV)
  end subroutine gather_3dhme

  subroutine gather_3dmoh(iBlock, buffer, nV, nX, nY, nZ)
    use ModGITM
    use ModInputs
    use ModSources
    integer, intent(in) :: iBlock, nV, nX, nY, nZ
    real, intent(out) :: buffer(nV, nX, nY, nZ)
    integer :: iAlt, iLat, iLon, jx, jy, jz, iv
    integer :: iiAlt, iiLat, iiLon

    do iAlt = -1, nAlts + 2
      jz = iAlt + 2
      iiAlt = max(min(iAlt, nAlts), 1)
      do iLat = -1, nLats + 2
        jy = iLat + 2
        iiLat = max(min(iLat, nLats), 1)
        do iLon = -1, nLons + 2
          jx = iLon + 2
          iiLon = max(min(iLon, nLons), 1)
          iv = 0
          iv = iv + 1; buffer(iv, jx, jy, jz) = Longitude(iLon, iBlock)
          iv = iv + 1; buffer(iv, jx, jy, jz) = Latitude(iLat, iBlock)
          iv = iv + 1; buffer(iv, jx, jy, jz) = Altitude_GB(iLon, iLat, iAlt, iBlock)
          iv = iv + 1; buffer(iv, jx, jy, jz) = Rho(iLon, iLat, iAlt, iBlock)
          iv = iv + 1; buffer(iv, jx, jy, jz) = Velocity(iLon, iLat, iAlt, iEast_, iBlock)
          iv = iv + 1; buffer(iv, jx, jy, jz) = Velocity(iLon, iLat, iAlt, iNorth_, iBlock)
          iv = iv + 1; buffer(iv, jx, jy, jz) = Viscosity(iiLon, iiLat, iiAlt, iEast_) / dt
          iv = iv + 1; buffer(iv, jx, jy, jz) = Viscosity(iiLon, iiLat, iiAlt, iNorth_) / dt
          iv = iv + 1; buffer(iv, jx, jy, jz) = IonDrag(iiLon, iiLat, iiAlt, iEast_)
          iv = iv + 1; buffer(iv, jx, jy, jz) = IonDrag(iiLon, iiLat, iiAlt, iNorth_)
          iv = iv + 1; buffer(iv, jx, jy, jz) = HorizAdvection(iiLon, iiLat, iiAlt, 1)
          iv = iv + 1; buffer(iv, jx, jy, jz) = HorizAdvection(iiLon, iiLat, iiAlt, 2)
          iv = iv + 1; buffer(iv, jx, jy, jz) = HorizPressureGrad(iiLon, iiLat, iiAlt, 1)
          iv = iv + 1; buffer(iv, jx, jy, jz) = HorizPressureGrad(iiLon, iiLat, iiAlt, 2)
          iv = iv + 1; buffer(iv, jx, jy, jz) = HorizCoriolis(iiLon, iiLat, iiAlt, 1)
          iv = iv + 1; buffer(iv, jx, jy, jz) = HorizCoriolis(iiLon, iiLat, iiAlt, 2)
          iv = iv + 1; buffer(iv, jx, jy, jz) = Centrifugal(iiLon, iiLat, iiAlt, 2)
        enddo
      enddo
    enddo
    if (iv /= nV) call gather_error('gather_3dmoh', iv, nV)
  end subroutine gather_3dmoh

  subroutine gather_3dmov(iBlock, buffer, nV, nX, nY, nZ)
    use ModGITM
    use ModInputs
    use ModSources
    integer, intent(in) :: iBlock, nV, nX, nY, nZ
    real, intent(out) :: buffer(nV, nX, nY, nZ)
    integer :: iAlt, iLat, iLon, jx, jy, jz, iv, iSpecies
    integer :: iiAlt, iiLat, iiLon

    do iAlt = -1, nAlts + 2
      jz = iAlt + 2
      iiAlt = max(min(iAlt, nAlts), 1)
      do iLat = -1, nLats + 2
        jy = iLat + 2
        iiLat = max(min(iLat, nLats), 1)
        do iLon = -1, nLons + 2
          jx = iLon + 2
          iiLon = max(min(iLon, nLons), 1)
          iv = 0
          iv = iv + 1; buffer(iv, jx, jy, jz) = Longitude(iLon, iBlock)
          iv = iv + 1; buffer(iv, jx, jy, jz) = Latitude(iLat, iBlock)
          iv = iv + 1; buffer(iv, jx, jy, jz) = Altitude_GB(iLon, iLat, iAlt, iBlock)
          do iSpecies = 1, nSpecies
            iv = iv + 1; buffer(iv, jx, jy, jz) = &
              VerticalVelocity(iLon, iLat, iAlt, iSpecies, iBlock)
            iv = iv + 1; buffer(iv, jx, jy, jz) = &
              VerticalIonDrag(iiLon, iiLat, iiAlt, iSpecies)
            iv = iv + 1; buffer(iv, jx, jy, jz) = &
              VertAdvection(iiLon, iiLat, iiAlt, iSpecies)
          enddo
          iv = iv + 1; buffer(iv, jx, jy, jz) = VertCoriolis(iiLon, iiLat, iiAlt)
          iv = iv + 1; buffer(iv, jx, jy, jz) = VertCentrifugal(iiLon, iiLat, iiAlt)
          iv = iv + 1; buffer(iv, jx, jy, jz) = EffectiveGravity(iiLon, iiLat, iiAlt)
        enddo
      enddo
    enddo
    if (iv /= nV) call gather_error('gather_3dmov', iv, nV)
  end subroutine gather_3dmov

  ! ==================================================================
  ! 2D output gather routines
  ! ==================================================================

  subroutine gather_2dgel(iBlock, buffer, nV, nX, nY, nZ)
    use ModGITM
    use ModInputs
    use ModElectrodynamics
    integer, intent(in) :: iBlock, nV, nX, nY, nZ
    real, intent(out) :: buffer(nV, nX, nY, nZ)
    integer :: iLat, iLon, iv

    do iLat = 1, nLats
      do iLon = 1, nLons
        iv = 0
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = Longitude(iLon, iBlock)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = Latitude(iLat, iBlock)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = Altitude_GB(iLon, iLat, 1, iBlock)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = Potential(iLon, iLat, 1, iBlock)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = PedersenConductance(iLon, iLat, iBlock)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = HallConductance(iLon, iLat, iBlock)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = ElectronAverageEnergyDiffuse(iLon, iLat)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = ElectronEnergyFluxDiffuse(iLon, iLat)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = ElectronAverageEnergyWave(iLon, iLat)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = ElectronEnergyFluxWave(iLon, iLat)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = ElectronAverageEnergyMono(iLon, iLat)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = ElectronEnergyFluxMono(iLon, iLat)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = IonAverageEnergy(iLon, iLat)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = IonEnergyFlux(iLon, iLat)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = DivJuAlt(iLon, iLat)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = PedersenFieldLine(iLon, iLat)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = HallFieldLine(iLon, iLat)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = DivJuFieldLine(iLon, iLat)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = LengthFieldLine(iLon, iLat)
      enddo
    enddo
    if (iv /= nV) call gather_error('gather_2dgel', iv, nV)
  end subroutine gather_2dgel

  subroutine gather_2dmel(iBlock, buffer, nV, nX, nY, nZ)
    use ModGITM
    use ModInputs
    use ModElectrodynamics
    use ModConstants, only: Pi
    integer, intent(in) :: iBlock, nV, nX, nY, nZ
    real, intent(out) :: buffer(nV, nX, nY, nZ)
    integer :: iLat, iLon, iv

    do iLat = 1, nMagLats
      do iLon = 1, nMagLons + 1
        iv = 0
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = MagLonMC(iLon, iLat) * Pi / 180.0
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = MagLatMC(iLon, iLat) * Pi / 180.0
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = Altitude_GB(iLon, iLat, 1, iBlock)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = MagLocTimeMC(iLon, iLat) * Pi / 180.0
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = GeoLatMC(iLon, iLat)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = GeoLonMC(iLon, iLat)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = SigmaPedersenMC(iLon, iLat)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = SigmaHallMC(iLon, iLat)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = DivJuAltMC(iLon, iLat)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = LengthMC(iLon, iLat)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = SigmaPPMC(iLon, iLat)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = SigmaLLMC(iLon, iLat)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = SigmaHHMC(iLon, iLat)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = SigmaCCMC(iLon, iLat)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = SigmaPLMC(iLon, iLat)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = SigmaLPMC(iLon, iLat)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = KDpmMC(iLon, iLat)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = KdlmMC(iLon, iLat)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = solver_a_mc(iLon, iLat)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = solver_b_mc(iLon, iLat)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = solver_c_mc(iLon, iLat)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = solver_d_mc(iLon, iLat)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = solver_e_mc(iLon, iLat)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = solver_s_mc(iLon, iLat)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = DynamoPotentialMC(iLon, iLat)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = Ed1new(iLon, iLat)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = Ed2new(iLon, iLat)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = kpmMC(iLon, iLat)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = klmMC(iLon, iLat)
      enddo
    enddo
    if (iv /= nV) call gather_error('gather_2dmel', iv, nV)
  end subroutine gather_2dmel

  subroutine gather_2dtec(iBlock, buffer, nV, nX, nY, nZ)
    use ModGITM
    use ModInputs
    use ModEUV, only: Sza
    integer, intent(in) :: iBlock, nV, nX, nY, nZ
    real, intent(out) :: buffer(nV, nX, nY, nZ)
    integer :: iLat, iLon, iv

    call calc_vtec(iBlock)

    do iLat = 1, nLats
      do iLon = 1, nLons
        iv = 0
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = Longitude(iLon, iBlock)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = Latitude(iLat, iBlock)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = Altitude_GB(iLon, iLat, 1, iBlock)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = Sza(iLon, iLat, iBlock)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = VTEC(iLon, iLat, iBlock)
      enddo
    enddo
    if (iv /= nV) call gather_error('gather_2dtec', iv, nV)
  end subroutine gather_2dtec

  subroutine gather_2danc(iBlock, buffer, nV, nX, nY, nZ)
    use ModGITM
    use ModInputs
    use ModElectrodynamics
    use ModEUV, only: Sza
    use ModSources
    integer, intent(in) :: iBlock, nV, nX, nY, nZ
    real, intent(out) :: buffer(nV, nX, nY, nZ)
    integer :: iLat, iLon, iv

    call calc_vtec(iBlock)

    do iLat = 1, nLats
      do iLon = 1, nLons
        iv = 0
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = Longitude(iLon, iBlock)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = Latitude(iLat, iBlock)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = Altitude_GB(iLon, iLat, nAlts, iBlock)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = LocalTime(iLon)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = Sza(iLon, iLat, iBlock)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = VTEC(iLon, iLat, iBlock)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = JouleHeating2d(iLon, iLat)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = HeatTransfer2d(iLon, iLat)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = EuvHeating2d(iLon, iLat)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = PhotoElectronHeating2d(iLon, iLat)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = ChemicalHeating2d(iLon, iLat)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = RadiativeCooling2d(iLon, iLat)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = CO2Cooling2d(iLon, iLat)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = NOCooling2d(iLon, iLat)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = OCooling2d(iLon, iLat)
      enddo
    enddo
    if (iv /= nV) call gather_error('gather_2danc', iv, nV)
  end subroutine gather_2danc

  subroutine gather_2dhme(iBlock, buffer, nV, nX, nY, nZ)
    use ModGITM
    use ModInputs
    integer, intent(in) :: iBlock, nV, nX, nY, nZ
    real, intent(out) :: buffer(nV, nX, nY, nZ)
    integer :: iLat, iLon, iv

    call calc_vtec(iBlock)

    do iLat = 1, nLats
      do iLon = 1, nLons
        iv = 0
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = Longitude(iLon, iBlock)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = Latitude(iLat, iBlock)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = Altitude_GB(iLon, iLat, 1, iBlock)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = LocalTime(iLon)
        iv = iv + 1; buffer(iv, iLon, iLat, 1) = VTEC(iLon, iLat, iBlock)
      enddo
    enddo
    if (iv /= nV) call gather_error('gather_2dhme', iv, nV)
  end subroutine gather_2dhme

  ! ==================================================================
  ! 1D output gather routines
  ! ==================================================================

  subroutine gather_1dall(iBlock, buffer, nV, nZ, iiLon, iiLat, rLon, rLat)
    use ModGITM
    use ModInputs
    integer, intent(in) :: iBlock, nV, nZ, iiLon, iiLat
    real, intent(in) :: rLon, rLat
    real, intent(out) :: buffer(nV, 1, 1, nZ)
    integer :: iAlt, jz, iv

    do iAlt = -1, nAlts + 2
      jz = iAlt + 2
      iv = 0

      iv = iv + 1; buffer(iv, 1, 1, jz) = &
        rLon * Longitude(iiLon, iBlock) + (1 - rLon) * Longitude(iiLon + 1, iBlock)
      iv = iv + 1; buffer(iv, 1, 1, jz) = &
        rLat * Latitude(iiLat, iBlock) + (1 - rLat) * Latitude(iiLat + 1, iBlock)
      iv = iv + 1; buffer(iv, 1, 1, jz) = Altitude_GB(iiLon, iiLat, iAlt, iBlock)

      call fill_neutral_vars_interp2d(iv, buffer, nV, nZ, jz, &
                                      iAlt, iBlock, iiLon, iiLat, rLon, rLat)
      call fill_ion_vars_interp2d(iv, buffer, nV, nZ, jz, &
                                  iAlt, iBlock, iiLon, iiLat, rLon, rLat)
    enddo
    if (iv /= nV) call gather_error('gather_1dall', iv, nV)
  end subroutine gather_1dall

  subroutine gather_1dglo(iBlock, buffer, nV, nZ)
    use ModGITM
    use ModInputs
    integer, intent(in) :: iBlock, nV, nZ
    real, intent(out) :: buffer(nV, 1, 1, nZ)

    ! STUB: 1DGLO variables are registered in ModOutputRegistry but data
    ! gathering is not yet implemented.  Buffer is zeroed; implement this
    ! subroutine when real glow emission data is available.
    buffer = 0.0
  end subroutine gather_1dglo

  subroutine gather_1dthm(iBlock, buffer, nV, nZ)
    use ModGITM
    use ModInputs
    use ModSources
    use ModEUV, only: EuvTotal
    integer, intent(in) :: iBlock, nV, nZ
    real, intent(out) :: buffer(nV, 1, 1, nZ)
    integer :: iAlt, iiAlt, jz, iv, i

    do iAlt = -1, nAlts + 2
      jz = iAlt + 2
      iiAlt = max(min(iAlt, nAlts), 1)
      iv = 0
      iv = iv + 1; buffer(iv, 1, 1, jz) = Longitude(1, 1)
      iv = iv + 1; buffer(iv, 1, 1, jz) = Latitude(1, 1)
      iv = iv + 1; buffer(iv, 1, 1, jz) = Altitude_GB(1, 1, iAlt, 1)
      iv = iv + 1; buffer(iv, 1, 1, jz) = &
        EuvHeating(1, 1, iiAlt, 1) * dt * TempUnit(1, 1, iiAlt)
      iv = iv + 1; buffer(iv, 1, 1, jz) = &
        Conduction(1, 1, iiAlt) * TempUnit(1, 1, iiAlt)
      iv = iv + 1; buffer(iv, 1, 1, jz) = MoleConduction(1, 1, iiAlt)
      iv = iv + 1; buffer(iv, 1, 1, jz) = EddyCond(1, 1, iiAlt)
      iv = iv + 1; buffer(iv, 1, 1, jz) = EddyCondAdia(1, 1, iiAlt)
      iv = iv + 1; buffer(iv, 1, 1, jz) = &
        ChemicalHeatingRate(1, 1, iiAlt) * TempUnit(1, 1, iiAlt)
      iv = iv + 1; buffer(iv, 1, 1, jz) = &
        JouleHeating(1, 1, iiAlt) * dt * TempUnit(1, 1, iiAlt)
      iv = iv + 1; buffer(iv, 1, 1, jz) = &
        -RadCooling(1, 1, iiAlt, 1) * dt * TempUnit(1, 1, iiAlt)
      iv = iv + 1; buffer(iv, 1, 1, jz) = &
        -OCooling(1, 1, iiAlt) * dt * TempUnit(1, 1, iiAlt)
      iv = iv + 1; buffer(iv, 1, 1, jz) = EuvTotal(1, 1, iiAlt, 1) * dt
      do i = 1, nSpeciesTotal
        iv = iv + 1; buffer(iv, 1, 1, jz) = NeutralSourcesTotal(iiAlt, i)
      enddo
      do i = 1, nSpeciesTotal
        iv = iv + 1; buffer(iv, 1, 1, jz) = NeutralLossesTotal(iiAlt, i)
      enddo
    enddo
    if (iv /= nV) call gather_error('gather_1dthm', iv, nV)
  end subroutine gather_1dthm

  subroutine gather_1dchm(iBlock, buffer, nV, nZ)
    use ModGITM
    use ModInputs
    use ModSources
    use ModConstants
    integer, intent(in) :: iBlock, nV, nZ
    real, intent(out) :: buffer(nV, 1, 1, nZ)
    integer :: iAlt, iiAlt, jz, iv, iReact

    do iAlt = -1, nAlts + 2
      jz = iAlt + 2
      iiAlt = max(min(iAlt, nAlts), 1)
      iv = 0
      iv = iv + 1; buffer(iv, 1, 1, jz) = Longitude(1, iBlock)
      iv = iv + 1; buffer(iv, 1, 1, jz) = Latitude(1, iBlock)
      iv = iv + 1; buffer(iv, 1, 1, jz) = Altitude_GB(1, 1, iAlt, iBlock)
      do iReact = 1, nReactions
        iv = iv + 1; buffer(iv, 1, 1, jz) = &
          ChemicalHeatingSpecies(1, 1, iiAlt, iReact) / Element_Charge
      enddo
      iv = iv + 1; buffer(iv, 1, 1, jz) = &
        ChemicalHeatingRate(1, 1, iiAlt) * &
        cp(1, 1, iiAlt, iBlock) * &
        Rho(1, 1, iiAlt, iBlock) * TempUnit(1, 1, iiAlt) / Element_Charge
    enddo
    if (iv /= nV) call gather_error('gather_1dchm', iv, nV)
  end subroutine gather_1dchm

  subroutine gather_1dnew(iBlock, buffer, nV, nZ, iiLon, iiLat, rLon, rLat)
    use ModGITM
    use ModInputs
    use ModEUV, only: Sza
    integer, intent(in) :: iBlock, nV, nZ, iiLon, iiLat
    real, intent(in) :: rLon, rLat
    real, intent(out) :: buffer(nV, 1, 1, nZ)
    real :: Tmp(0:nLons + 1, 0:nLats + 1)
    integer :: iAlt, jz, iv, i

    do iAlt = -1, nAlts + 2
      jz = iAlt + 2
      iv = 0

      ! 1DNEW has different coords: Lon, LocalTime, Lat, SZA, Alt
      iv = iv + 1; buffer(iv, 1, 1, jz) = &
        rLon * Longitude(iiLon, iBlock) + (1 - rLon) * Longitude(iiLon + 1, iBlock)
      iv = iv + 1; buffer(iv, 1, 1, jz) = &
        rLon * LocalTime(iiLon) + (1 - rLon) * LocalTime(iiLon + 1)
      iv = iv + 1; buffer(iv, 1, 1, jz) = &
        rLat * Latitude(iiLat, iBlock) + (1 - rLat) * Latitude(iiLat + 1, iBlock)
      iv = iv + 1; buffer(iv, 1, 1, jz) = Sza(iiLon, iiLat, iBlock)
      iv = iv + 1; buffer(iv, 1, 1, jz) = Altitude_GB(iiLon, iiLat, iAlt, iBlock)

      call fill_neutral_vars_interp2d(iv, buffer, nV, nZ, jz, &
                                      iAlt, iBlock, iiLon, iiLat, rLon, rLat)
      call fill_ion_vars_interp2d(iv, buffer, nV, nZ, jz, &
                                  iAlt, iBlock, iiLon, iiLat, rLon, rLat)

      ! 1DNEW-specific: mixing ratios
      do i = 1, nSpecies
        Tmp = NDensityS(0:nLons + 1, 0:nLats + 1, iAlt, i, iBlock) / &
              NDensity(0:nLons + 1, 0:nLats + 1, iAlt, iBlock)
        iv = iv + 1; buffer(iv, 1, 1, jz) = inter2d(Tmp, iiLon, iiLat, rLon, rLat)
      enddo
    enddo
    if (iv /= nV) call gather_error('gather_1dnew', iv, nV)
  end subroutine gather_1dnew

  ! ==================================================================
  ! 0D output gather routines
  ! ==================================================================

  subroutine gather_0dall(iBlock, buffer, nV, iiLon, iiLat, iiAlt, rLon, rLat, rAlt)
    use ModGITM
    use ModInputs
    use ModEUV, only: HeatingEfficiency_CB
    use ModSources, only: JouleHeating, RadCooling, EuvHeating, Conduction
    integer, intent(in) :: iBlock, nV, iiLon, iiLat, iiAlt
    real, intent(in) :: rLon, rLat, rAlt
    real, intent(out) :: buffer(nV, 1, 1, 1)
    real :: Tmp(0:nLons + 1, 0:nLats + 1, 0:nAlts + 1)
    integer :: jAlt, iv, i

    jAlt = max(min(iiAlt, nAlts), 1)
    iv = 0

    iv = iv + 1; buffer(iv, 1, 1, 1) = &
      rLon * Longitude(iiLon, iBlock) + (1 - rLon) * Longitude(iiLon + 1, iBlock)
    iv = iv + 1; buffer(iv, 1, 1, 1) = &
      rLat * Latitude(iiLat, iBlock) + (1 - rLat) * Latitude(iiLat + 1, iBlock)
    iv = iv + 1; buffer(iv, 1, 1, 1) = &
      rAlt * Altitude_GB(iiLon, iiLat, iiAlt, iBlock) + &
      (1 - rAlt) * Altitude_GB(iiLon + 1, iiLat + 1, iiAlt + 1, iBlock)

    Tmp = Rho(0:nLons + 1, 0:nLats + 1, 0:nAlts + 1, iBlock)
    iv = iv + 1; buffer(iv, 1, 1, 1) = inter3d(Tmp, iiLon, iiLat, iiAlt, rLon, rLat, rAlt)

    do i = 1, nSpeciesTotal
      Tmp = NDensityS(0:nLons + 1, 0:nLats + 1, 0:nAlts + 1, i, iBlock)
      iv = iv + 1; buffer(iv, 1, 1, 1) = inter3d(Tmp, iiLon, iiLat, iiAlt, rLon, rLat, rAlt)
    enddo

    Tmp = Temperature(0:nLons + 1, 0:nLats + 1, 0:nAlts + 1, iBlock) * &
          TempUnit(0:nLons + 1, 0:nLats + 1, 0:nAlts + 1)
    iv = iv + 1; buffer(iv, 1, 1, 1) = inter3d(Tmp, iiLon, iiLat, iiAlt, rLon, rLat, rAlt)

    do i = 1, 3
      Tmp = Velocity(0:nLons + 1, 0:nLats + 1, 0:nAlts + 1, i, iBlock)
      iv = iv + 1; buffer(iv, 1, 1, 1) = inter3d(Tmp, iiLon, iiLat, iiAlt, rLon, rLat, rAlt)
    enddo

    do i = 1, nSpecies
      Tmp = VerticalVelocity(0:nLons + 1, 0:nLats + 1, 0:nAlts + 1, i, iBlock)
      iv = iv + 1; buffer(iv, 1, 1, 1) = inter3d(Tmp, iiLon, iiLat, iiAlt, rLon, rLat, rAlt)
    enddo

    do i = 1, nIons
      Tmp = IDensityS(0:nLons + 1, 0:nLats + 1, 0:nAlts + 1, i, iBlock)
      iv = iv + 1; buffer(iv, 1, 1, 1) = inter3d(Tmp, iiLon, iiLat, iiAlt, rLon, rLat, rAlt)
    enddo

    Tmp = eTemperature(0:nLons + 1, 0:nLats + 1, 0:nAlts + 1, iBlock)
    iv = iv + 1; buffer(iv, 1, 1, 1) = inter3d(Tmp, iiLon, iiLat, iiAlt, rLon, rLat, rAlt)
    Tmp = iTemperature(0:nLons + 1, 0:nLats + 1, 0:nAlts + 1, iBlock)
    iv = iv + 1; buffer(iv, 1, 1, 1) = inter3d(Tmp, iiLon, iiLat, iiAlt, rLon, rLat, rAlt)

    Tmp = IVelocity(0:nLons + 1, 0:nLats + 1, 0:nAlts + 1, iEast_, iBlock)
    iv = iv + 1; buffer(iv, 1, 1, 1) = inter3d(Tmp, iiLon, iiLat, iiAlt, rLon, rLat, rAlt)
    Tmp = IVelocity(0:nLons + 1, 0:nLats + 1, 0:nAlts + 1, iNorth_, iBlock)
    iv = iv + 1; buffer(iv, 1, 1, 1) = inter3d(Tmp, iiLon, iiLat, iiAlt, rLon, rLat, rAlt)
    Tmp = IVelocity(0:nLons + 1, 0:nLats + 1, 0:nAlts + 1, iUp_, iBlock)
    iv = iv + 1; buffer(iv, 1, 1, 1) = inter3d(Tmp, iiLon, iiLat, iiAlt, rLon, rLat, rAlt)

    ! Mixing ratios
    do i = 1, nSpecies
      Tmp = NDensityS(0:nLons + 1, 0:nLats + 1, 0:nAlts + 1, i, iBlock) / &
            NDensity(0:nLons + 1, 0:nLats + 1, 0:nAlts + 1, iBlock)
      iv = iv + 1; buffer(iv, 1, 1, 1) = inter3d(Tmp, iiLon, iiLat, iiAlt, rLon, rLat, rAlt)
    enddo

    ! Heating terms
    iv = iv + 1; buffer(iv, 1, 1, 1) = dt * RadCooling(1, 1, jAlt, iBlock) * TempUnit(1, 1, jAlt)
    iv = iv + 1; buffer(iv, 1, 1, 1) = dt * EuvHeating(1, 1, jAlt, iBlock) * TempUnit(1, 1, jAlt)
    iv = iv + 1; buffer(iv, 1, 1, 1) = Conduction(1, 1, jAlt) * TempUnit(1, 1, jAlt)
    iv = iv + 1; buffer(iv, 1, 1, 1) = &
      dt * EuvHeating(1, 1, jAlt, iBlock) * TempUnit(1, 1, jAlt) - &
      dt * RadCooling(1, 1, jAlt, iBlock) * TempUnit(1, 1, jAlt) + &
      Conduction(1, 1, jAlt) * TempUnit(1, 1, jAlt)
    iv = iv + 1; buffer(iv, 1, 1, 1) = HeatingEfficiency_CB(1, 1, jAlt, iBlock)
    if (iv /= nV) call gather_error('gather_0dall', iv, nV)
  end subroutine gather_0dall

  ! ==================================================================
  ! USR output gather routines
  ! ==================================================================

  subroutine gather_3dusr(iBlock, buffer, nV, nX, nY, nZ)
    use ModGITM
    use ModUserGITM
    integer, intent(in) :: iBlock, nV, nX, nY, nZ
    real, intent(out) :: buffer(nV, nX, nY, nZ)
    integer :: iAlt, iLat, iLon, jx, jy, jz

    do iAlt = -1, nAlts + 2
      jz = iAlt + 2
      do iLat = -1, nLats + 2
        jy = iLat + 2
        do iLon = -1, nLons + 2
          jx = iLon + 2
          buffer(1, jx, jy, jz) = Longitude(iLon, iBlock)
          buffer(2, jx, jy, jz) = Latitude(iLat, iBlock)
          buffer(3, jx, jy, jz) = Altitude_GB(iLon, iLat, iAlt, iBlock)
          if (nV > 3) &
            buffer(4:nV, jx, jy, jz) = UserData3D(iLon, iLat, iAlt, 1:nV - 3, iBlock)
        enddo
      enddo
    enddo
  end subroutine gather_3dusr

  subroutine gather_2dusr(iBlock, buffer, nV, nX, nY, nZ)
    use ModGITM
    use ModUserGITM
    integer, intent(in) :: iBlock, nV, nX, nY, nZ
    real, intent(out) :: buffer(nV, nX, nY, nZ)
    integer :: iLat, iLon

    do iLat = 1, nLats
      do iLon = 1, nLons
        buffer(1, iLon, iLat, 1) = Longitude(iLon, iBlock)
        buffer(2, iLon, iLat, 1) = Latitude(iLat, iBlock)
        buffer(3, iLon, iLat, 1) = Altitude_GB(iLon, iLat, 1, iBlock)
        if (nV > 3) &
          buffer(4:nV, iLon, iLat, 1) = UserData2D(iLon, iLat, 1, 1:nV - 3, iBlock)
      enddo
    enddo
  end subroutine gather_2dusr

  ! 1DUSR: mirrors output_1dUser — interpolated Lon/Lat, Altitude, then
  ! whatever the user put in UserData1D rows 1:nV-3 at each altitude.
  subroutine gather_1dusr(iBlock, buffer, nV, nZ, iiLon, iiLat, rLon, rLat)
    use ModGITM
    use ModUserGITM
    integer, intent(in) :: iBlock, nV, nZ, iiLon, iiLat
    real, intent(in) :: rLon, rLat
    real, intent(out) :: buffer(nV, 1, 1, nZ)
    integer :: iAlt, jz

    do iAlt = -1, nAlts + 2
      jz = iAlt + 2
      buffer(1, 1, 1, jz) = &
        rLon * Longitude(iiLon, iBlock) + (1 - rLon) * Longitude(iiLon + 1, iBlock)
      buffer(2, 1, 1, jz) = &
        rLat * Latitude(iiLat, iBlock) + (1 - rLat) * Latitude(iiLat + 1, iBlock)
      buffer(3, 1, 1, jz) = Altitude_GB(iiLon, iiLat, iAlt, iBlock)
      if (nV > 3) &
        buffer(4:nV, 1, 1, jz) = UserData1D(1, 1, iAlt, 1:nV - 3)
    enddo
  end subroutine gather_1dusr

  ! 0DUSR: single point at (iiLon, iiLat, iiAlt), interpolated.
  subroutine gather_0dusr(iBlock, buffer, nV, iiLon, iiLat, iiAlt, rLon, rLat, rAlt)
    use ModGITM
    use ModUserGITM
    integer, intent(in) :: iBlock, nV, iiLon, iiLat, iiAlt
    real, intent(in) :: rLon, rLat, rAlt
    real, intent(out) :: buffer(nV, 1, 1, 1)

    buffer(1, 1, 1, 1) = &
      rLon * Longitude(iiLon, iBlock) + (1 - rLon) * Longitude(iiLon + 1, iBlock)
    buffer(2, 1, 1, 1) = &
      rLat * Latitude(iiLat, iBlock) + (1 - rLat) * Latitude(iiLat + 1, iBlock)
    buffer(3, 1, 1, 1) = &
      rAlt * Altitude_GB(iiLon, iiLat, iiAlt, iBlock) + &
      (1 - rAlt) * Altitude_GB(iiLon, iiLat, iiAlt + 1, iBlock)
    ! 0D user data is a scalar per variable; no UserData0D array exists —
    ! the user populates UserData1D(1,1,1,1:nVarsUser0d-3) at a fixed point.
    if (nV > 3) &
      buffer(4:nV, 1, 1, 1) = UserData1D(1, 1, 1, 1:nV - 3)
  end subroutine gather_0dusr

  ! ==================================================================
  ! Interpolation helper functions
  ! ==================================================================

  ! Assumed-shape arrays: no dependency on nLons/nLats/nAlts at module scope
  real function inter2d(variable, iiLon, iiLat, rLon, rLat)
    real, intent(in) :: variable(:, :), rLon, rLat
    integer, intent(in) :: iiLon, iiLat

    inter2d = &
      rLon * rLat * variable(iiLon, iiLat) + &
      (1 - rLon) * rLat * variable(iiLon + 1, iiLat) + &
      rLon * (1 - rLat) * variable(iiLon, iiLat + 1) + &
      (1 - rLon) * (1 - rLat) * variable(iiLon + 1, iiLat + 1)
  end function inter2d

  real function inter3d(variable, iiLon, iiLat, iiAlt, rLon, rLat, rAlt)
    real, intent(in) :: variable(:, :, :)
    real, intent(in) :: rLon, rLat, rAlt
    integer, intent(in) :: iiLon, iiLat, iiAlt

    inter3d = &
      rLon * rLat * rAlt * variable(iiLon, iiLat, iiAlt) + &
      (1 - rLon) * rLat * rAlt * variable(iiLon + 1, iiLat, iiAlt) + &
      rLon * (1 - rLat) * rAlt * variable(iiLon, iiLat + 1, iiAlt) + &
      (1 - rLon) * (1 - rLat) * rAlt * variable(iiLon + 1, iiLat + 1, iiAlt) + &
      rLon * rLat * (1 - rAlt) * variable(iiLon, iiLat, iiAlt + 1) + &
      (1 - rLon) * rLat * (1 - rAlt) * variable(iiLon + 1, iiLat, iiAlt + 1) + &
      rLon * (1 - rLat) * (1 - rAlt) * variable(iiLon, iiLat + 1, iiAlt + 1) + &
      (1 - rLon) * (1 - rLat) * (1 - rAlt) * variable(iiLon + 1, iiLat + 1, iiAlt + 1)
  end function inter3d

end module ModOutputGather
