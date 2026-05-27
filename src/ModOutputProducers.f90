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
  integer, parameter :: nContainers = 3

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

end module ModOutputProducers
