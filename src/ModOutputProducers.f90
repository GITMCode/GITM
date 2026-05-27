! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE
!
! ModOutputProducers.f90
!
! Container-based output producers.  Each output type has two procedures:
!   define_schema_<type>(c)   — register variable names/units/shapes once
!   fill_<type>(c, iBlock)    — put data into the container each timestep
!
! Module-level containers(:) array lives for the whole run; allocated once in
! init_output_containers, reused across timesteps.
!
! Currently implemented: 2DMEL, 2DTEC, 2DGEL.

module ModOutputProducers

  use ModOutputContainer

  implicit none

  type(OutputContainer), allocatable, target :: containers(:)

  ! Index constants — one per output type handled here.
  integer, parameter :: iCont_2DMEL = 1
  integer, parameter :: iCont_2DTEC = 2
  integer, parameter :: iCont_2DGEL = 3
  integer, parameter :: iCont_3DNEU = 4
  integer, parameter :: iCont_3DALL = 5
  integer, parameter :: iCont_3DION = 6
  integer, parameter :: iCont_3DTHM = 7
  integer, parameter :: iCont_3DCHM = 8
  integer, parameter :: nContainers = 8

  logical :: containers_initialized = .false.

contains

  ! ---------------------------------------------------------------------------
  ! init_output_containers — allocate and define schemas.  Called once.
  ! ---------------------------------------------------------------------------
  subroutine init_output_containers()
    if (containers_initialized) return
    allocate(containers(nContainers))
    call define_schema_2dmel(containers(iCont_2DMEL))
    call define_schema_2dtec(containers(iCont_2DTEC))
    call define_schema_2dgel(containers(iCont_2DGEL))
    call define_schema_3dneu(containers(iCont_3DNEU))
    call define_schema_3dall(containers(iCont_3DALL))
    call define_schema_3dion(containers(iCont_3DION))
    call define_schema_3dthm(containers(iCont_3DTHM))
    call define_schema_3dchm(containers(iCont_3DCHM))
    containers_initialized = .true.
  end subroutine init_output_containers

  ! ---------------------------------------------------------------------------
  ! define_schema_2dmel — register all 2DMEL variables.
  !
  ! Grid: magnetic (nMagLons+1) x nMagLats, no altitude.
  ! Axis variables are 1D; all others are 2D (degenerate 3rd dim = 1).
  ! Units: degrees for lon/lat (deg, not rad), hours for MLT.
  ! GeoLat/GeoLon are stored in radians in the model; conversion at put time.
  ! ---------------------------------------------------------------------------
  subroutine define_schema_2dmel(c)
    use ModElectrodynamics, only: nMagLats, nMagLons
    type(OutputContainer), intent(inout) :: c

    c%cType    = '2DMEL'
    ! gridKind drives BOTH participation (prepare: only iProc==0 writes) AND
    ! file layout (mpiio_write_container: MAG_2D writes at offset 0, not
    ! striped by global block index). Changing the participation rule for
    ! MAG_2D requires updating the offset logic in lockstep.
    c%gridKind = GRID_MAG_2D

    ! Axis variables (1D coordinate arrays).
    call c%define_var('mlon', units='degrees_east', &
                      shape3=[nMagLons + 1, 1, 1], is_axis=.true., &
                      longName='Magnetic longitude')
    call c%define_var('mlat', units='degrees_north', &
                      shape3=[nMagLats, 1, 1], is_axis=.true., &
                      longName='Magnetic latitude')

    ! 2D data variables (nMagLons+1 x nMagLats x 1).
    call c%define_var('MLT', units='hr', &
                      shape3=[nMagLons + 1, nMagLats, 1], &
                      longName='Magnetic local time')
    call c%define_var('GeoLat', units='degrees_north', &
                      shape3=[nMagLons + 1, nMagLats, 1], &
                      longName='Geographic latitude')
    call c%define_var('GeoLon', units='degrees_east', &
                      shape3=[nMagLons + 1, nMagLats, 1], &
                      longName='Geographic longitude')
    call c%define_var('PedCond', units='S', &
                      shape3=[nMagLons + 1, nMagLats, 1], &
                      longName='Height-integrated Pedersen conductance')
    call c%define_var('HalCond', units='S', &
                      shape3=[nMagLons + 1, nMagLats, 1], &
                      longName='Height-integrated Hall conductance')
    call c%define_var('DivJuAlt', units='A/m2', &
                      shape3=[nMagLons + 1, nMagLats, 1], &
                      longName='Height-integrated divergence of horizontal current')
    call c%define_var('FL_length', units='m', &
                      shape3=[nMagLons + 1, nMagLats, 1], &
                      longName='Field line length')
    call c%define_var('Sigma_PP', units='S', &
                      shape3=[nMagLons + 1, nMagLats, 1], &
                      longName='Sigma PP conductance component')
    call c%define_var('Sigma_LL', units='S', &
                      shape3=[nMagLons + 1, nMagLats, 1], &
                      longName='Sigma LL conductance component')
    call c%define_var('Sigma_H', units='S', &
                      shape3=[nMagLons + 1, nMagLats, 1], &
                      longName='Hall conductance component')
    call c%define_var('Sigma_C', units='S', &
                      shape3=[nMagLons + 1, nMagLats, 1], &
                      longName='Cowling conductance component')
    call c%define_var('Sigma_PL', units='S', &
                      shape3=[nMagLons + 1, nMagLats, 1], &
                      longName='Sigma PL conductance component')
    call c%define_var('Sigma_LP', units='S', &
                      shape3=[nMagLons + 1, nMagLats, 1], &
                      longName='Sigma LP conductance component')
    call c%define_var('K_mphi', units='', &
                      shape3=[nMagLons + 1, nMagLats, 1], &
                      longName='Dynamo K_mphi coefficient')
    call c%define_var('K_mlambda', units='', &
                      shape3=[nMagLons + 1, nMagLats, 1], &
                      longName='Dynamo K_mlambda coefficient')
    call c%define_var('Solver_A', units='', &
                      shape3=[nMagLons + 1, nMagLats, 1], &
                      longName='Dynamo solver coefficient A')
    call c%define_var('Solver_B', units='', &
                      shape3=[nMagLons + 1, nMagLats, 1], &
                      longName='Dynamo solver coefficient B')
    call c%define_var('Solver_C', units='', &
                      shape3=[nMagLons + 1, nMagLats, 1], &
                      longName='Dynamo solver coefficient C')
    call c%define_var('Solver_D', units='', &
                      shape3=[nMagLons + 1, nMagLats, 1], &
                      longName='Dynamo solver coefficient D')
    call c%define_var('Solver_E', units='', &
                      shape3=[nMagLons + 1, nMagLats, 1], &
                      longName='Dynamo solver coefficient E')
    call c%define_var('Solver_S', units='', &
                      shape3=[nMagLons + 1, nMagLats, 1], &
                      longName='Dynamo solver coefficient S')
    call c%define_var('Phi_dynamo', units='V', &
                      shape3=[nMagLons + 1, nMagLats, 1], &
                      longName='Dynamo electric potential')
    call c%define_var('Ed1new', units='V/m', &
                      shape3=[nMagLons + 1, nMagLats, 1], &
                      longName='Dynamo electric field component 1')
    call c%define_var('Ed2new', units='V/m', &
                      shape3=[nMagLons + 1, nMagLats, 1], &
                      longName='Dynamo electric field component 2')
    call c%define_var('Kphi', units='', &
                      shape3=[nMagLons + 1, nMagLats, 1], &
                      longName='Dynamo Kphi coefficient')
    call c%define_var('Klamda', units='', &
                      shape3=[nMagLons + 1, nMagLats, 1], &
                      longName='Dynamo Klambda coefficient')
  end subroutine define_schema_2dmel

  ! ---------------------------------------------------------------------------
  ! fill_2dmel — put 2DMEL data into the container for iBlock.
  !
  ! Only iBlock == 1 writes (GRID_MAG_2D); all others return after prepare().
  ! GeoLat/GeoLon are in radians in the model; multiplied by cRadToDeg here.
  ! MagLonMC/MagLatMC are already in degrees.
  ! ---------------------------------------------------------------------------
  subroutine fill_2dmel(c, iBlock)
    use ModElectrodynamics
    use ModConst, only: cRadToDeg
    type(OutputContainer), intent(inout) :: c
    integer, intent(in) :: iBlock

    call c%prepare(iBlock)
    if (.not. c%this_rank_writes) return

    ! Guard: arrays may not be allocated when dynamo is disabled.
    if (.not. (allocated(MagLonMC) .and. allocated(MagLatMC))) return

    ! Axis arrays: 1D slices from the 2D grid (all rows are identical for axes).
    call c%put('mlon', real(MagLonMC(1:nMagLons + 1, 1), output_kind))
    call c%put('mlat', real(MagLatMC(1, 1:nMagLats), output_kind))

    ! 2D data variables.
    call c%put('MLT',      real(MagLocTimeMC, output_kind))
    call c%put('GeoLat',   real(GeoLatMC * cRadToDeg, output_kind))
    call c%put('GeoLon',   real(GeoLonMC * cRadToDeg, output_kind))
    call c%put('PedCond',  real(SigmaPedersenMC, output_kind))
    call c%put('HalCond',  real(SigmaHallMC, output_kind))
    call c%put('DivJuAlt', real(DivJuAltMC, output_kind))
    call c%put('FL_length',real(LengthMC, output_kind))
    call c%put('Sigma_PP', real(SigmaPPMC, output_kind))
    call c%put('Sigma_LL', real(SigmaLLMC, output_kind))
    call c%put('Sigma_H',  real(SigmaHHMC, output_kind))
    call c%put('Sigma_C',  real(SigmaCCMC, output_kind))
    call c%put('Sigma_PL', real(SigmaPLMC, output_kind))
    call c%put('Sigma_LP', real(SigmaLPMC, output_kind))
    call c%put('K_mphi',    real(KDpmMC, output_kind))
    call c%put('K_mlambda', real(KDlmMC, output_kind))
    call c%put('Solver_A',  real(solver_a_mc, output_kind))
    call c%put('Solver_B',  real(solver_b_mc, output_kind))
    call c%put('Solver_C',  real(solver_c_mc, output_kind))
    call c%put('Solver_D',  real(solver_d_mc, output_kind))
    call c%put('Solver_E',  real(solver_e_mc, output_kind))
    call c%put('Solver_S',  real(solver_s_mc, output_kind))
    call c%put('Phi_dynamo',real(DynamoPotentialMC, output_kind))
    call c%put('Ed1new',    real(Ed1new, output_kind))
    call c%put('Ed2new',    real(Ed2new, output_kind))
    call c%put('Kphi',      real(kpmMC, output_kind))
    call c%put('Klamda',    real(klmMC, output_kind))
  end subroutine fill_2dmel

  ! ---------------------------------------------------------------------------
  ! define_schema_2dtec — register all 2DTEC variables.
  !
  ! Grid: geographic 2D, per-block, nLons x nLats.
  ! Longitude/Latitude stored in radians in model; converted to degrees at put time.
  ! ---------------------------------------------------------------------------
  subroutine define_schema_2dtec(c)
    use ModGITM, only: nLons, nLats
    type(OutputContainer), intent(inout) :: c

    c%cType    = '2DTEC'
    c%gridKind = GRID_GEO_2D

    ! Axis variables (1D coordinate arrays).
    call c%define_var('Longitude', units='degrees_east', &
                      shape3=[nLons, 1, 1], is_axis=.true., &
                      longName='Geographic longitude')
    call c%define_var('Latitude', units='degrees_north', &
                      shape3=[nLats, 1, 1], is_axis=.true., &
                      longName='Geographic latitude')

    ! 2D data variables (nLons x nLats x 1).
    call c%define_var('Altitude', units='m', &
                      shape3=[nLons, nLats, 1], &
                      longName='Altitude above surface')
    call c%define_var('SZA', units='rad', &
                      shape3=[nLons, nLats, 1], &
                      longName='Solar zenith angle')
    call c%define_var('TEC', units='TECU', &
                      shape3=[nLons, nLats, 1], &
                      longName='Vertical total electron content')
  end subroutine define_schema_2dtec

  ! ---------------------------------------------------------------------------
  ! fill_2dtec — put 2DTEC data into the container for iBlock.
  !
  ! All blocks participate (GRID_GEO_2D).
  ! Longitude/Latitude are in radians in the model; converted to degrees here.
  ! ---------------------------------------------------------------------------
  subroutine fill_2dtec(c, iBlock)
    use ModGITM
    use ModEUV, only: Sza
    use ModConst, only: cRadToDeg
    type(OutputContainer), intent(inout) :: c
    integer, intent(in) :: iBlock

    call c%prepare(iBlock)
    if (.not. c%this_rank_writes) return

    call calc_vtec(iBlock)

    call c%put('Longitude', real(Longitude(1:nLons, iBlock) * cRadToDeg, output_kind))
    call c%put('Latitude',  real(Latitude(1:nLats, iBlock) * cRadToDeg, output_kind))
    call c%put('Altitude',  real(Altitude_GB(1:nLons, 1:nLats, 1, iBlock), output_kind))
    call c%put('SZA',       real(Sza(1:nLons, 1:nLats, iBlock), output_kind))
    call c%put('TEC',       real(VTEC(1:nLons, 1:nLats, iBlock), output_kind))
  end subroutine fill_2dtec

  ! ---------------------------------------------------------------------------
  ! define_schema_2dgel — register all 2DGEL variables.
  !
  ! Grid: geographic 2D, per-block, nLons x nLats.
  ! Longitude/Latitude stored in radians in model; converted to degrees at put time.
  ! ---------------------------------------------------------------------------
  subroutine define_schema_2dgel(c)
    use ModGITM, only: nLons, nLats
    type(OutputContainer), intent(inout) :: c

    c%cType    = '2DGEL'
    c%gridKind = GRID_GEO_2D

    ! Axis variables (1D coordinate arrays).
    call c%define_var('Longitude', units='degrees_east', &
                      shape3=[nLons, 1, 1], is_axis=.true., &
                      longName='Geographic longitude')
    call c%define_var('Latitude', units='degrees_north', &
                      shape3=[nLats, 1, 1], is_axis=.true., &
                      longName='Geographic latitude')

    ! 2D data variables (nLons x nLats x 1).
    call c%define_var('Altitude', units='m', &
                      shape3=[nLons, nLats, 1], &
                      longName='Altitude above surface')
    call c%define_var('pot', units='V', &
                      shape3=[nLons, nLats, 1], &
                      longName='Electric potential')
    call c%define_var('PedCond', units='S', &
                      shape3=[nLons, nLats, 1], &
                      longName='Height-integrated Pedersen conductance')
    call c%define_var('HalCond', units='S', &
                      shape3=[nLons, nLats, 1], &
                      longName='Height-integrated Hall conductance')
    call c%define_var('AveE', units='eV', &
                      shape3=[nLons, nLats, 1], &
                      longName='Diffuse electron precipitation average energy')
    call c%define_var('eFlux', units='ergs/cm2/s', &
                      shape3=[nLons, nLats, 1], &
                      longName='Diffuse electron precipitation energy flux')
    call c%define_var('AveE_W', units='eV', &
                      shape3=[nLons, nLats, 1], &
                      longName='Wave-driven electron precipitation average energy')
    call c%define_var('eFlux_W', units='ergs/cm2/s', &
                      shape3=[nLons, nLats, 1], &
                      longName='Wave-driven electron precipitation energy flux')
    call c%define_var('AveE_M', units='eV', &
                      shape3=[nLons, nLats, 1], &
                      longName='Monoenergetic electron precipitation average energy')
    call c%define_var('eFlux_M', units='ergs/cm2/s', &
                      shape3=[nLons, nLats, 1], &
                      longName='Monoenergetic electron precipitation energy flux')
    call c%define_var('AveE_I', units='eV', &
                      shape3=[nLons, nLats, 1], &
                      longName='Ion precipitation average energy')
    call c%define_var('eFlux_I', units='ergs/cm2/s', &
                      shape3=[nLons, nLats, 1], &
                      longName='Ion precipitation energy flux')
    call c%define_var('DivJuAlt', units='A/m2', &
                      shape3=[nLons, nLats, 1], &
                      longName='Height-integrated divergence of horizontal current')
    call c%define_var('PedFLCond', units='S', &
                      shape3=[nLons, nLats, 1], &
                      longName='Field-line integrated Pedersen conductance')
    call c%define_var('HalFLCond', units='S', &
                      shape3=[nLons, nLats, 1], &
                      longName='Field-line integrated Hall conductance')
    call c%define_var('DivJuFL', units='A/m2', &
                      shape3=[nLons, nLats, 1], &
                      longName='Field-line divergence of current')
    call c%define_var('FL_length', units='m', &
                      shape3=[nLons, nLats, 1], &
                      longName='Field line length')
  end subroutine define_schema_2dgel

  ! ---------------------------------------------------------------------------
  ! fill_2dgel — put 2DGEL data into the container for iBlock.
  !
  ! All blocks participate (GRID_GEO_2D).
  ! Longitude/Latitude are in radians in the model; converted to degrees here.
  ! ModElectrodynamics 2D arrays (ElectronAverageEnergyDiffuse etc.) are
  ! block-local: no iBlock index, hold the current block's data at call time.
  ! ---------------------------------------------------------------------------
  subroutine fill_2dgel(c, iBlock)
    use ModGITM
    use ModElectrodynamics
    use ModConst, only: cRadToDeg
    type(OutputContainer), intent(inout) :: c
    integer, intent(in) :: iBlock

    call c%prepare(iBlock)
    if (.not. c%this_rank_writes) return

    call c%put('Longitude', real(Longitude(1:nLons, iBlock) * cRadToDeg, output_kind))
    call c%put('Latitude',  real(Latitude(1:nLats, iBlock) * cRadToDeg, output_kind))
    call c%put('Altitude',  real(Altitude_GB(1:nLons, 1:nLats, 1, iBlock), output_kind))
    call c%put('pot',       real(Potential(1:nLons, 1:nLats, 1, iBlock), output_kind))
    call c%put('PedCond',   real(PedersenConductance(1:nLons, 1:nLats, iBlock), output_kind))
    call c%put('HalCond',   real(HallConductance(1:nLons, 1:nLats, iBlock), output_kind))
    call c%put('AveE',      real(ElectronAverageEnergyDiffuse(1:nLons, 1:nLats), output_kind))
    call c%put('eFlux',     real(ElectronEnergyFluxDiffuse(1:nLons, 1:nLats), output_kind))
    call c%put('AveE_W',    real(ElectronAverageEnergyWave(1:nLons, 1:nLats), output_kind))
    call c%put('eFlux_W',   real(ElectronEnergyFluxWave(1:nLons, 1:nLats), output_kind))
    call c%put('AveE_M',    real(ElectronAverageEnergyMono(1:nLons, 1:nLats), output_kind))
    call c%put('eFlux_M',   real(ElectronEnergyFluxMono(1:nLons, 1:nLats), output_kind))
    call c%put('AveE_I',    real(IonAverageEnergy(1:nLons, 1:nLats), output_kind))
    call c%put('eFlux_I',   real(IonEnergyFlux(1:nLons, 1:nLats), output_kind))
    call c%put('DivJuAlt',  real(DivJuAlt(1:nLons, 1:nLats), output_kind))
    call c%put('PedFLCond', real(PedersenFieldLine(1:nLons, 1:nLats), output_kind))
    call c%put('HalFLCond', real(HallFieldLine(1:nLons, 1:nLats), output_kind))
    call c%put('DivJuFL',   real(DivJuFieldLine(1:nLons, 1:nLats), output_kind))
    call c%put('FL_length', real(LengthFieldLine(1:nLons, 1:nLats), output_kind))
  end subroutine fill_2dgel

  ! ---------------------------------------------------------------------------
  ! define_schema_3dneu — register all 3DNEU variables.
  ! ---------------------------------------------------------------------------
  subroutine define_schema_3dneu(c)
    use ModGITM, only: nLons, nLats, nAlts
    use ModPlanet, only: nSpeciesTotal, nSpecies, cSpecies
    type(OutputContainer), intent(inout) :: c
    integer :: i

    c%cType    = '3DNEU'
    c%gridKind = GRID_GEO_3D

    call c%define_var('Longitude', units='degrees_east', &
                      shape3=[nLons, 1, 1], is_axis=.true., &
                      longName='Geographic longitude')
    call c%define_var('Latitude', units='degrees_north', &
                      shape3=[nLats, 1, 1], is_axis=.true., &
                      longName='Geographic latitude')
    call c%define_var('Altitude', units='m', &
                      shape3=[nLons, nLats, nAlts], &
                      longName='Altitude above surface')
    call c%define_var('Rho', units='kg/m3', &
                      shape3=[nLons, nLats, nAlts], &
                      longName='Total neutral mass density')
    do i = 1, nSpeciesTotal
      call c%define_var('['//trim(cSpecies(i))//']', units='/m3', &
                        shape3=[nLons, nLats, nAlts], &
                        longName='Number density of '//trim(cSpecies(i)))
    end do
    call c%define_var('Temperature', units='K', &
                      shape3=[nLons, nLats, nAlts], &
                      longName='Neutral temperature')
    call c%define_var('Vn (east)', units='m/s', &
                      shape3=[nLons, nLats, nAlts], &
                      longName='Eastward neutral wind')
    call c%define_var('Vn (north)', units='m/s', &
                      shape3=[nLons, nLats, nAlts], &
                      longName='Northward neutral wind')
    call c%define_var('Vn (up)', units='m/s', &
                      shape3=[nLons, nLats, nAlts], &
                      longName='Vertical neutral wind')
    do i = 1, nSpecies
      call c%define_var('Vn (up,'//trim(cSpecies(i))//')', units='m/s', &
                        shape3=[nLons, nLats, nAlts], &
                        longName='Vertical neutral wind for '//trim(cSpecies(i)))
    end do
  end subroutine define_schema_3dneu

  ! ---------------------------------------------------------------------------
  ! fill_3dneu — put 3DNEU data into the container for iBlock.
  ! ---------------------------------------------------------------------------
  subroutine fill_3dneu(c, iBlock)
    use ModGITM
    use ModConst, only: cRadToDeg
    use ModPlanet, only: nSpeciesTotal, nSpecies, cSpecies
    type(OutputContainer), intent(inout) :: c
    integer, intent(in) :: iBlock
    integer :: i

    call c%prepare(iBlock)
    if (.not. c%this_rank_writes) return

    call c%put('Longitude', real(Longitude(1:nLons, iBlock) * cRadToDeg, output_kind))
    call c%put('Latitude',  real(Latitude(1:nLats, iBlock) * cRadToDeg, output_kind))
    call c%put('Altitude',  real(Altitude_GB(1:nLons, 1:nLats, 1:nAlts, iBlock), output_kind))
    call c%put('Rho',       real(Rho(1:nLons, 1:nLats, 1:nAlts, iBlock), output_kind))
    do i = 1, nSpeciesTotal
      call c%put('['//trim(cSpecies(i))//']', &
                 real(NDensityS(1:nLons, 1:nLats, 1:nAlts, i, iBlock), output_kind))
    end do
    call c%put('Temperature', &
               real(Temperature(1:nLons, 1:nLats, 1:nAlts, iBlock) * &
                    TempUnit(1:nLons, 1:nLats, 1:nAlts), output_kind))
    call c%put('Vn (east)',  real(Velocity(1:nLons, 1:nLats, 1:nAlts, 1, iBlock), output_kind))
    call c%put('Vn (north)', real(Velocity(1:nLons, 1:nLats, 1:nAlts, 2, iBlock), output_kind))
    call c%put('Vn (up)',    real(Velocity(1:nLons, 1:nLats, 1:nAlts, 3, iBlock), output_kind))
    do i = 1, nSpecies
      call c%put('Vn (up,'//trim(cSpecies(i))//')', &
                 real(VerticalVelocity(1:nLons, 1:nLats, 1:nAlts, i, iBlock), output_kind))
    end do
  end subroutine fill_3dneu

  ! ---------------------------------------------------------------------------
  ! define_schema_3dall — register all 3DALL variables.
  ! ---------------------------------------------------------------------------
  subroutine define_schema_3dall(c)
    use ModGITM, only: nLons, nLats, nAlts
    use ModPlanet, only: nSpeciesTotal, nSpecies, nIons, cSpecies, cIons
    type(OutputContainer), intent(inout) :: c
    integer :: i

    c%cType    = '3DALL'
    c%gridKind = GRID_GEO_3D

    call c%define_var('Longitude', units='degrees_east', &
                      shape3=[nLons, 1, 1], is_axis=.true., &
                      longName='Geographic longitude')
    call c%define_var('Latitude', units='degrees_north', &
                      shape3=[nLats, 1, 1], is_axis=.true., &
                      longName='Geographic latitude')
    call c%define_var('Altitude', units='m', &
                      shape3=[nLons, nLats, nAlts], &
                      longName='Altitude above surface')
    call c%define_var('Rho', units='kg/m3', &
                      shape3=[nLons, nLats, nAlts], &
                      longName='Total neutral mass density')
    do i = 1, nSpeciesTotal
      call c%define_var('['//trim(cSpecies(i))//']', units='/m3', &
                        shape3=[nLons, nLats, nAlts], &
                        longName='Number density of '//trim(cSpecies(i)))
    end do
    call c%define_var('Temperature', units='K', &
                      shape3=[nLons, nLats, nAlts], &
                      longName='Neutral temperature')
    call c%define_var('Vn (east)', units='m/s', &
                      shape3=[nLons, nLats, nAlts], &
                      longName='Eastward neutral wind')
    call c%define_var('Vn (north)', units='m/s', &
                      shape3=[nLons, nLats, nAlts], &
                      longName='Northward neutral wind')
    call c%define_var('Vn (up)', units='m/s', &
                      shape3=[nLons, nLats, nAlts], &
                      longName='Vertical neutral wind')
    do i = 1, nSpecies
      call c%define_var('Vn (up,'//trim(cSpecies(i))//')', units='m/s', &
                        shape3=[nLons, nLats, nAlts], &
                        longName='Vertical neutral wind for '//trim(cSpecies(i)))
    end do
    do i = 1, nIons
      call c%define_var('['//trim(cIons(i))//']', units='/m3', &
                        shape3=[nLons, nLats, nAlts], &
                        longName='Number density of '//trim(cIons(i)))
    end do
    call c%define_var('eTemperature', units='K', &
                      shape3=[nLons, nLats, nAlts], &
                      longName='Electron temperature')
    call c%define_var('iTemperature', units='K', &
                      shape3=[nLons, nLats, nAlts], &
                      longName='Ion temperature')
    call c%define_var('Vi (east)', units='m/s', &
                      shape3=[nLons, nLats, nAlts], &
                      longName='Eastward ion drift')
    call c%define_var('Vi (north)', units='m/s', &
                      shape3=[nLons, nLats, nAlts], &
                      longName='Northward ion drift')
    call c%define_var('Vi (up)', units='m/s', &
                      shape3=[nLons, nLats, nAlts], &
                      longName='Vertical ion drift')
  end subroutine define_schema_3dall

  ! ---------------------------------------------------------------------------
  ! fill_3dall — put 3DALL data into the container for iBlock.
  ! ---------------------------------------------------------------------------
  subroutine fill_3dall(c, iBlock)
    use ModGITM
    use ModConst, only: cRadToDeg
    use ModPlanet, only: nSpeciesTotal, nSpecies, nIons, cSpecies, cIons
    type(OutputContainer), intent(inout) :: c
    integer, intent(in) :: iBlock
    integer :: i

    call c%prepare(iBlock)
    if (.not. c%this_rank_writes) return

    call c%put('Longitude', real(Longitude(1:nLons, iBlock) * cRadToDeg, output_kind))
    call c%put('Latitude',  real(Latitude(1:nLats, iBlock) * cRadToDeg, output_kind))
    call c%put('Altitude',  real(Altitude_GB(1:nLons, 1:nLats, 1:nAlts, iBlock), output_kind))
    call c%put('Rho',       real(Rho(1:nLons, 1:nLats, 1:nAlts, iBlock), output_kind))
    do i = 1, nSpeciesTotal
      call c%put('['//trim(cSpecies(i))//']', &
                 real(NDensityS(1:nLons, 1:nLats, 1:nAlts, i, iBlock), output_kind))
    end do
    call c%put('Temperature', &
               real(Temperature(1:nLons, 1:nLats, 1:nAlts, iBlock) * &
                    TempUnit(1:nLons, 1:nLats, 1:nAlts), output_kind))
    call c%put('Vn (east)',  real(Velocity(1:nLons, 1:nLats, 1:nAlts, 1, iBlock), output_kind))
    call c%put('Vn (north)', real(Velocity(1:nLons, 1:nLats, 1:nAlts, 2, iBlock), output_kind))
    call c%put('Vn (up)',    real(Velocity(1:nLons, 1:nLats, 1:nAlts, 3, iBlock), output_kind))
    do i = 1, nSpecies
      call c%put('Vn (up,'//trim(cSpecies(i))//')', &
                 real(VerticalVelocity(1:nLons, 1:nLats, 1:nAlts, i, iBlock), output_kind))
    end do
    do i = 1, nIons
      call c%put('['//trim(cIons(i))//']', &
                 real(IDensityS(1:nLons, 1:nLats, 1:nAlts, i, iBlock), output_kind))
    end do
    call c%put('eTemperature', real(eTemperature(1:nLons, 1:nLats, 1:nAlts, iBlock), output_kind))
    call c%put('iTemperature', real(ITemperature(1:nLons, 1:nLats, 1:nAlts, iBlock), output_kind))
    call c%put('Vi (east)',    real(IVelocity(1:nLons, 1:nLats, 1:nAlts, 1, iBlock), output_kind))
    call c%put('Vi (north)',   real(IVelocity(1:nLons, 1:nLats, 1:nAlts, 2, iBlock), output_kind))
    call c%put('Vi (up)',      real(IVelocity(1:nLons, 1:nLats, 1:nAlts, 3, iBlock), output_kind))
  end subroutine fill_3dall

  ! ---------------------------------------------------------------------------
  ! define_schema_3dion — register all 3DION variables.
  ! ---------------------------------------------------------------------------
  subroutine define_schema_3dion(c)
    use ModGITM, only: nLons, nLats, nAlts
    use ModPlanet, only: nIons, cIons
    type(OutputContainer), intent(inout) :: c
    integer :: i

    c%cType    = '3DION'
    c%gridKind = GRID_GEO_3D

    call c%define_var('Longitude', units='degrees_east', &
                      shape3=[nLons, 1, 1], is_axis=.true., &
                      longName='Geographic longitude')
    call c%define_var('Latitude', units='degrees_north', &
                      shape3=[nLats, 1, 1], is_axis=.true., &
                      longName='Geographic latitude')
    call c%define_var('Altitude', units='m', &
                      shape3=[nLons, nLats, nAlts], &
                      longName='Altitude above surface')
    do i = 1, nIons
      call c%define_var('['//trim(cIons(i))//']', units='/m3', &
                        shape3=[nLons, nLats, nAlts], &
                        longName='Number density of '//trim(cIons(i)))
    end do
    call c%define_var('eTemperature', units='K', &
                      shape3=[nLons, nLats, nAlts], &
                      longName='Electron temperature')
    call c%define_var('iTemperature', units='K', &
                      shape3=[nLons, nLats, nAlts], &
                      longName='Ion temperature')
    call c%define_var('Vi (east)',  units='m/s', shape3=[nLons, nLats, nAlts], &
                      longName='Eastward ion drift')
    call c%define_var('Vi (north)', units='m/s', shape3=[nLons, nLats, nAlts], &
                      longName='Northward ion drift')
    call c%define_var('Vi (up)',    units='m/s', shape3=[nLons, nLats, nAlts], &
                      longName='Vertical ion drift')
    call c%define_var('Ed1', units='V/m', shape3=[nLons, nLats, nAlts], &
                      longName='Electric field component 1 (magnetic east)')
    call c%define_var('Ed2', units='V/m', shape3=[nLons, nLats, nAlts], &
                      longName='Electric field component 2 (magnetic north)')
    call c%define_var('Je1', units='A/m2', shape3=[nLons, nLats, nAlts], &
                      longName='Current density component 1 (magnetic east)')
    call c%define_var('Je2', units='A/m2', shape3=[nLons, nLats, nAlts], &
                      longName='Current density component 2 (magnetic north)')
    call c%define_var('Magnetic Latitude',  units='deg', shape3=[nLons, nLats, nAlts], &
                      longName='Magnetic latitude')
    call c%define_var('Magnetic Longitude', units='deg', shape3=[nLons, nLats, nAlts], &
                      longName='Magnetic longitude')
    call c%define_var('B.F. East',     units='T', shape3=[nLons, nLats, nAlts], &
                      longName='Eastward magnetic field component')
    call c%define_var('B.F. North',    units='T', shape3=[nLons, nLats, nAlts], &
                      longName='Northward magnetic field component')
    call c%define_var('B.F. Vertical', units='T', shape3=[nLons, nLats, nAlts], &
                      longName='Vertical magnetic field component')
    call c%define_var('B.F. Magnitude', units='T', shape3=[nLons, nLats, nAlts], &
                      longName='Magnetic field magnitude')
    call c%define_var('Potential', units='V', shape3=[nLons, nLats, nAlts], &
                      longName='Electric potential')
    call c%define_var('E.F. East',      units='V/m', shape3=[nLons, nLats, nAlts], &
                      longName='Eastward electric field component')
    call c%define_var('E.F. North',     units='V/m', shape3=[nLons, nLats, nAlts], &
                      longName='Northward electric field component')
    call c%define_var('E.F. Vertical',  units='V/m', shape3=[nLons, nLats, nAlts], &
                      longName='Vertical electric field component')
    call c%define_var('E.F. Magnitude', units='V/m', shape3=[nLons, nLats, nAlts], &
                      longName='Electric field magnitude')
    call c%define_var('IN Collision Freq', units='Hz', shape3=[nLons, nLats, nAlts], &
                      longName='Ion-neutral collision frequency')
    call c%define_var('PressGrad (east)',  units='m/s2', shape3=[nLons, nLats, nAlts], &
                      longName='Pressure gradient acceleration eastward')
    call c%define_var('PressGrad (north)', units='m/s2', shape3=[nLons, nLats, nAlts], &
                      longName='Pressure gradient acceleration northward')
    call c%define_var('PressGrad (up)',   units='m/s2', shape3=[nLons, nLats, nAlts], &
                      longName='Pressure gradient acceleration upward')
  end subroutine define_schema_3dion

  ! ---------------------------------------------------------------------------
  ! fill_3dion — put 3DION data into the container for iBlock.
  ! ---------------------------------------------------------------------------
  subroutine fill_3dion(c, iBlock)
    use ModGITM
    use ModConst, only: cRadToDeg
    use ModPlanet, only: nIons, cIons
    use ModElectrodynamics
    type(OutputContainer), intent(inout) :: c
    integer, intent(in) :: iBlock
    integer :: i

    call c%prepare(iBlock)
    if (.not. c%this_rank_writes) return

    call c%put('Longitude', real(Longitude(1:nLons, iBlock) * cRadToDeg, output_kind))
    call c%put('Latitude',  real(Latitude(1:nLats, iBlock) * cRadToDeg, output_kind))
    call c%put('Altitude',  real(Altitude_GB(1:nLons, 1:nLats, 1:nAlts, iBlock), output_kind))
    do i = 1, nIons
      call c%put('['//trim(cIons(i))//']', &
                 real(IDensityS(1:nLons, 1:nLats, 1:nAlts, i, iBlock), output_kind))
    end do
    call c%put('eTemperature', real(eTemperature(1:nLons, 1:nLats, 1:nAlts, iBlock), output_kind))
    call c%put('iTemperature', real(ITemperature(1:nLons, 1:nLats, 1:nAlts, iBlock), output_kind))
    call c%put('Vi (east)',    real(IVelocity(1:nLons, 1:nLats, 1:nAlts, 1, iBlock), output_kind))
    call c%put('Vi (north)',   real(IVelocity(1:nLons, 1:nLats, 1:nAlts, 2, iBlock), output_kind))
    call c%put('Vi (up)',      real(IVelocity(1:nLons, 1:nLats, 1:nAlts, 3, iBlock), output_kind))
    call c%put('Ed1', real(ed1(1:nLons, 1:nLats, 1:nAlts), output_kind))
    call c%put('Ed2', real(ed2(1:nLons, 1:nLats, 1:nAlts), output_kind))
    call c%put('Je1', real(je1(1:nLons, 1:nLats, 1:nAlts), output_kind))
    call c%put('Je2', real(je2(1:nLons, 1:nLats, 1:nAlts), output_kind))
    call c%put('Magnetic Latitude',  real(MLatitude(1:nLons, 1:nLats, 1:nAlts, iBlock), output_kind))
    call c%put('Magnetic Longitude', real(MLongitude(1:nLons, 1:nLats, 1:nAlts, iBlock), output_kind))
    call c%put('B.F. East',     real(B0(1:nLons, 1:nLats, 1:nAlts, 1, iBlock), output_kind))
    call c%put('B.F. North',    real(B0(1:nLons, 1:nLats, 1:nAlts, 2, iBlock), output_kind))
    call c%put('B.F. Vertical', real(B0(1:nLons, 1:nLats, 1:nAlts, 3, iBlock), output_kind))
    call c%put('B.F. Magnitude',real(B0(1:nLons, 1:nLats, 1:nAlts, 4, iBlock), output_kind))
    call c%put('Potential', real(Potential(1:nLons, 1:nLats, 1:nAlts, iBlock), output_kind))
    call c%put('E.F. East',      real(EField(1:nLons, 1:nLats, 1:nAlts, 1), output_kind))
    call c%put('E.F. North',     real(EField(1:nLons, 1:nLats, 1:nAlts, 2), output_kind))
    call c%put('E.F. Vertical',  real(EField(1:nLons, 1:nLats, 1:nAlts, 3), output_kind))
    call c%put('E.F. Magnitude', &
               real(sqrt(sum(EField(1:nLons, 1:nLats, 1:nAlts, :)**2, dim=4)), output_kind))
    call c%put('IN Collision Freq', real(Collisions(1:nLons, 1:nLats, 1:nAlts, iVIN_), output_kind))
    call c%put('PressGrad (east)',  real(IonPressureGradient(1:nLons, 1:nLats, 1:nAlts, 1, iBlock), output_kind))
    call c%put('PressGrad (north)', real(IonPressureGradient(1:nLons, 1:nLats, 1:nAlts, 2, iBlock), output_kind))
    call c%put('PressGrad (up)',    real(IonPressureGradient(1:nLons, 1:nLats, 1:nAlts, 3, iBlock), output_kind))
  end subroutine fill_3dion

  ! ---------------------------------------------------------------------------
  ! define_schema_3dthm — register all 3DTHM variables.
  ! ---------------------------------------------------------------------------
  subroutine define_schema_3dthm(c)
    use ModGITM, only: nLons, nLats, nAlts
    type(OutputContainer), intent(inout) :: c

    c%cType    = '3DTHM'
    c%gridKind = GRID_GEO_3D

    call c%define_var('Longitude', units='degrees_east', &
                      shape3=[nLons, 1, 1], is_axis=.true., &
                      longName='Geographic longitude')
    call c%define_var('Latitude', units='degrees_north', &
                      shape3=[nLats, 1, 1], is_axis=.true., &
                      longName='Geographic latitude')
    call c%define_var('Altitude', units='m', &
                      shape3=[nLons, nLats, nAlts], &
                      longName='Altitude above surface')
    call c%define_var('EUV Heating (K/s)',              units='K/s', shape3=[nLons, nLats, nAlts], longName='EUV heating rate')
    call c%define_var('Conduction (K/s)',               units='K/s', shape3=[nLons, nLats, nAlts], longName='Thermal conduction heating rate')
    call c%define_var('Molecular Conduction (K/s)',     units='K/s', shape3=[nLons, nLats, nAlts], longName='Molecular conduction heating rate')
    call c%define_var('Eddy Conduction (K/s)',          units='K/s', shape3=[nLons, nLats, nAlts], longName='Eddy conduction heating rate')
    call c%define_var('Eddy Adiabatic Conduction (K/s)', units='K/s', shape3=[nLons, nLats, nAlts], longName='Eddy adiabatic conduction heating rate')
    call c%define_var('Chemical Heating (K/s)',         units='K/s', shape3=[nLons, nLats, nAlts], longName='Chemical heating rate')
    call c%define_var('Joule Heating (K/s)',            units='K/s', shape3=[nLons, nLats, nAlts], longName='Joule heating rate')
    call c%define_var('NO Cooling (K/s)',               units='K/s', shape3=[nLons, nLats, nAlts], longName='NO infrared cooling rate')
    call c%define_var('O Cooling (K/s)',                units='K/s', shape3=[nLons, nLats, nAlts], longName='O infrared cooling rate')
    call c%define_var('Total Abs EUV',                 units='W/m2', shape3=[nLons, nLats, nAlts], longName='Total absorbed EUV flux')
    call c%define_var('Cp',                            units='J/kg/K', shape3=[nLons, nLats, nAlts], longName='Specific heat at constant pressure')
    call c%define_var('Rho',                           units='kg/m3', shape3=[nLons, nLats, nAlts], longName='Total neutral mass density')
    call c%define_var('E-Field Mag',                   units='V/m', shape3=[nLons, nLats, nAlts], longName='Electric field magnitude')
    call c%define_var('Sigma Ped',                     units='S/m', shape3=[nLons, nLats, nAlts], longName='Pedersen conductivity')
    call c%define_var('Ionization Rate O_3P',          units='/s',  shape3=[nLons, nLats, nAlts], longName='Ionization rate for O(3P)')
    call c%define_var('Ionization Rate O2',            units='/s',  shape3=[nLons, nLats, nAlts], longName='Ionization rate for O2')
    call c%define_var('Ionization Rate N2',            units='/s',  shape3=[nLons, nLats, nAlts], longName='Ionization rate for N2')
  end subroutine define_schema_3dthm

  ! ---------------------------------------------------------------------------
  ! fill_3dthm — put 3DTHM data into the container for iBlock.
  ! ---------------------------------------------------------------------------
  subroutine fill_3dthm(c, iBlock)
    use ModGITM
    use ModConst, only: cRadToDeg
    use ModSources
    use ModElectrodynamics, only: Sigma_Pedersen
    use ModEUV, only: EuvTotal
    type(OutputContainer), intent(inout) :: c
    integer, intent(in) :: iBlock

    call c%prepare(iBlock)
    if (.not. c%this_rank_writes) return

    call c%put('Longitude', real(Longitude(1:nLons, iBlock) * cRadToDeg, output_kind))
    call c%put('Latitude',  real(Latitude(1:nLats, iBlock) * cRadToDeg, output_kind))
    call c%put('Altitude',  real(Altitude_GB(1:nLons, 1:nLats, 1:nAlts, iBlock), output_kind))
    call c%put('EUV Heating (K/s)', &
               real(EuvHeating(1:nLons, 1:nLats, 1:nAlts, iBlock) * &
                    TempUnit(1:nLons, 1:nLats, 1:nAlts), output_kind))
    call c%put('Conduction (K/s)', &
               real(Conduction(1:nLons, 1:nLats, 1:nAlts) * &
                    TempUnit(1:nLons, 1:nLats, 1:nAlts) / dt, output_kind))
    call c%put('Molecular Conduction (K/s)',      real(MoleConduction(1:nLons, 1:nLats, 1:nAlts), output_kind))
    call c%put('Eddy Conduction (K/s)',           real(EddyCond(1:nLons, 1:nLats, 1:nAlts), output_kind))
    call c%put('Eddy Adiabatic Conduction (K/s)', real(EddyCondAdia(1:nLons, 1:nLats, 1:nAlts), output_kind))
    call c%put('Chemical Heating (K/s)', &
               real(ChemicalHeatingRate(1:nLons, 1:nLats, 1:nAlts) * &
                    TempUnit(1:nLons, 1:nLats, 1:nAlts) / dt, output_kind))
    call c%put('Joule Heating (K/s)', &
               real(JouleHeating(1:nLons, 1:nLats, 1:nAlts) * &
                    TempUnit(1:nLons, 1:nLats, 1:nAlts), output_kind))
    call c%put('NO Cooling (K/s)', &
               real(-NOCooling(1:nLons, 1:nLats, 1:nAlts) * &
                    TempUnit(1:nLons, 1:nLats, 1:nAlts), output_kind))
    call c%put('O Cooling (K/s)', &
               real(-OCooling(1:nLons, 1:nLats, 1:nAlts) * &
                    TempUnit(1:nLons, 1:nLats, 1:nAlts), output_kind))
    call c%put('Total Abs EUV', real(EuvTotal(1:nLons, 1:nLats, 1:nAlts, iBlock), output_kind))
    call c%put('Cp',            real(cp(1:nLons, 1:nLats, 1:nAlts, iBlock), output_kind))
    call c%put('Rho',           real(Rho(1:nLons, 1:nLats, 1:nAlts, iBlock), output_kind))
    call c%put('E-Field Mag',   &
               real(sqrt(sum(EField(1:nLons, 1:nLats, 1:nAlts, :)**2, dim=4)), output_kind))
    call c%put('Sigma Ped',     real(Sigma_Pedersen(1:nLons, 1:nLats, 1:nAlts), output_kind))
    call c%put('Ionization Rate O_3P', real(AuroralIonRateS(1:nLons, 1:nLats, 1:nAlts, iO_3P_, iBlock), output_kind))
    call c%put('Ionization Rate O2',   real(AuroralIonRateS(1:nLons, 1:nLats, 1:nAlts, iO2_,   iBlock), output_kind))
    call c%put('Ionization Rate N2',   real(AuroralIonRateS(1:nLons, 1:nLats, 1:nAlts, iN2_,   iBlock), output_kind))
  end subroutine fill_3dthm

  ! ---------------------------------------------------------------------------
  ! define_schema_3dchm — register all 3DCHM variables.
  ! ---------------------------------------------------------------------------
  subroutine define_schema_3dchm(c)
    use ModGITM, only: nLons, nLats, nAlts
    type(OutputContainer), intent(inout) :: c

    c%cType    = '3DCHM'
    c%gridKind = GRID_GEO_3D

    call c%define_var('Longitude', units='degrees_east', &
                      shape3=[nLons, 1, 1], is_axis=.true., &
                      longName='Geographic longitude')
    call c%define_var('Latitude', units='degrees_north', &
                      shape3=[nLats, 1, 1], is_axis=.true., &
                      longName='Geographic latitude')
    call c%define_var('Altitude', units='m', &
                      shape3=[nLons, nLats, nAlts], &
                      longName='Altitude above surface')
    call c%define_var('N2+ + e',    units='/m3/s', shape3=[nLons,nLats,nAlts], longName='N2+ + e recombination rate')
    call c%define_var('O2+ + e',    units='/m3/s', shape3=[nLons,nLats,nAlts], longName='O2+ + e recombination rate')
    call c%define_var('N2+ + O',    units='/m3/s', shape3=[nLons,nLats,nAlts], longName='N2+ + O reaction rate')
    call c%define_var('NO+ + e',    units='/m3/s', shape3=[nLons,nLats,nAlts], longName='NO+ + e recombination rate')
    call c%define_var('N+ + O2',    units='/m3/s', shape3=[nLons,nLats,nAlts], longName='N+ + O2 reaction rate')
    call c%define_var('NO + N',     units='/m3/s', shape3=[nLons,nLats,nAlts], longName='NO + N reaction rate')
    call c%define_var('O+ + O2',    units='/m3/s', shape3=[nLons,nLats,nAlts], longName='O+ + O2 reaction rate')
    call c%define_var('N + O2',     units='/m3/s', shape3=[nLons,nLats,nAlts], longName='N + O2 reaction rate')
    call c%define_var('O2+ + N',    units='/m3/s', shape3=[nLons,nLats,nAlts], longName='O2+ + N reaction rate')
    call c%define_var('O2+ + NO',   units='/m3/s', shape3=[nLons,nLats,nAlts], longName='O2+ + NO reaction rate')
    call c%define_var('O2+ + N2',   units='/m3/s', shape3=[nLons,nLats,nAlts], longName='O2+ + N2 reaction rate')
    call c%define_var('N2+ + O2',   units='/m3/s', shape3=[nLons,nLats,nAlts], longName='N2+ + O2 reaction rate')
    call c%define_var('N+ + O',     units='/m3/s', shape3=[nLons,nLats,nAlts], longName='N+ + O reaction rate')
    call c%define_var('O+ + N2',    units='/m3/s', shape3=[nLons,nLats,nAlts], longName='O+ + N2 reaction rate')
    call c%define_var('O(1D) + N2', units='/m3/s', shape3=[nLons,nLats,nAlts], longName='O(1D) + N2 quench rate')
    call c%define_var('O(1D) + O2', units='/m3/s', shape3=[nLons,nLats,nAlts], longName='O(1D) + O2 quench rate')
    call c%define_var('O(1D) + O',  units='/m3/s', shape3=[nLons,nLats,nAlts], longName='O(1D) + O quench rate')
    call c%define_var('O(1D) + e',  units='/m3/s', shape3=[nLons,nLats,nAlts], longName='O(1D) + e quench rate')
    call c%define_var('N(2D) + O2', units='/m3/s', shape3=[nLons,nLats,nAlts], longName='N(2D) + O2 reaction rate')
    call c%define_var('O+(2D)+e',   units='/m3/s', shape3=[nLons,nLats,nAlts], longName='O+(2D) + e recombination rate')
    call c%define_var('N(2D) + O',  units='/m3/s', shape3=[nLons,nLats,nAlts], longName='N(2D) + O quench rate')
    call c%define_var('N(2D) + e',  units='/m3/s', shape3=[nLons,nLats,nAlts], longName='N(2D) + e quench rate')
    call c%define_var('O+(2D + N2', units='/m3/s', shape3=[nLons,nLats,nAlts], longName='O+(2D) + N2 reaction rate')
    call c%define_var('O+(2P) + e', units='/m3/s', shape3=[nLons,nLats,nAlts], longName='O+(2P) + e recombination rate')
    call c%define_var('O+(2P) + O', units='/m3/s', shape3=[nLons,nLats,nAlts], longName='O+(2P) + O reaction rate')
    call c%define_var('O+(2P) + N2', units='/m3/s', shape3=[nLons,nLats,nAlts], longName='O+(2P) + N2 reaction rate')
    call c%define_var('Chemical Heating Rate', units='K/s', shape3=[nLons,nLats,nAlts], longName='Chemical heating rate')
  end subroutine define_schema_3dchm

  ! ---------------------------------------------------------------------------
  ! fill_3dchm — put 3DCHM data into the container for iBlock.
  ! ---------------------------------------------------------------------------
  subroutine fill_3dchm(c, iBlock)
    use ModGITM
    use ModConst, only: cRadToDeg
    use ModSources
    use ModConstants, only: Element_Charge
    type(OutputContainer), intent(inout) :: c
    integer, intent(in) :: iBlock
    integer :: i
    character(len=20), parameter :: chem_names(nReactions) = [ &
      'N2+ + e             ', 'O2+ + e             ', 'N2+ + O             ', &
      'NO+ + e             ', 'N+ + O2             ', 'NO + N              ', &
      'O+ + O2             ', 'N + O2              ', 'O2+ + N             ', &
      'O2+ + NO            ', 'O2+ + N2            ', 'N2+ + O2            ', &
      'N+ + O              ', 'O+ + N2             ', 'O(1D) + N2          ', &
      'O(1D) + O2          ', 'O(1D) + O           ', 'O(1D) + e           ', &
      'N(2D) + O2          ', 'O+(2D)+e            ', 'N(2D) + O           ', &
      'N(2D) + e           ', 'O+(2D + N2          ', 'O+(2P) + e          ', &
      'O+(2P) + O          ', 'O+(2P) + N2         ']

    call c%prepare(iBlock)
    if (.not. c%this_rank_writes) return

    call c%put('Longitude', real(Longitude(1:nLons, iBlock) * cRadToDeg, output_kind))
    call c%put('Latitude',  real(Latitude(1:nLats, iBlock) * cRadToDeg, output_kind))
    call c%put('Altitude',  real(Altitude_GB(1:nLons, 1:nLats, 1:nAlts, iBlock), output_kind))
    do i = 1, nReactions
      call c%put(trim(chem_names(i)), &
                 real(ChemicalHeatingSpecies(1:nLons, 1:nLats, 1:nAlts, i) / Element_Charge, output_kind))
    end do
    call c%put('Chemical Heating Rate', &
               real(ChemicalHeatingRate(1:nLons, 1:nLats, 1:nAlts) * &
                    cp(1:nLons, 1:nLats, 1:nAlts, iBlock) * &
                    Rho(1:nLons, 1:nLats, 1:nAlts, iBlock) * &
                    TempUnit(1:nLons, 1:nLats, 1:nAlts) / Element_Charge, output_kind))
  end subroutine fill_3dchm

end module ModOutputProducers
