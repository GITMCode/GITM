! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE
!
! ModOutputNetCDF.F90  (uppercase .F90: compiled with C preprocessor)
!
! PnetCDF parallel output backend for GITM.
!
! When GITM is built in an environment where pnetcdf-config is available,
! src/Makefile automatically detects it and adds -DHavePNetCDF to the
! compiler flags.  If pnetcdf-config is not found, this module compiles as
! a stub: all routines are no-ops except netcdf_open_file, which aborts.
!
! Each output type is written as a single NetCDF file collectively by all ranks.
! Ghost cells are stripped; each rank writes its interior slice to the correct
! position in the global (lon, lat[, alt]) array.
!
! Global grid layout (block ordering is lon-inner, lat-outer):
!   iBlockLon_0 = mod(iBLK-1, nBlocksLon)     ! 0-based lon block index
!   iBlockLat_0 =    (iBLK-1) / nBlocksLon    ! 0-based lat block index
!   lon_start   = iBlockLon_0 * nLons_i + 1   ! 1-based start in global array
!   lat_start   = iBlockLat_0 * nLats_i + 1
!
! For 3D types (nGhostCells=2): nX=nLons+4, nY=nLats+4, nZ=nAlts+4
!   Interior stripped: buffer(:, 3:nX-2, 3:nY-2, 3:nZ-2)
! For 2D types (nGhostCells=0): nX=nLons, nY=nLats, nZ=1
!   No stripping: full buffer used
!
! All registered variables (including the full Longitude/Latitude/Altitude
! fields) are written as data variables via varids(:).  Separate 1-D
! coordinate variables (time, and lon/lat/z in degrees/km) are also defined
! in define_file_structure; data variables additionally reference the full
! coordinate fields through their CF coordinates= attribute.
!
! Build requirements (auto-detected by src/Makefile via pnetcdf-config):
!   Compile: -I$(PNETCDF_INC)  and  -DHavePNetCDF
!   Link:    -L$(PNETCDF_LIB) -lpnetcdf

module ModOutputNetCDF

#ifdef HavePNetCDF
  use ModMpi
  use pnetcdf
#endif

  implicit none

#ifdef HavePNetCDF
  integer :: ncid = -1       ! PnetCDF file handle; -1 = no file open
  integer, allocatable :: varids(:), & ! variable IDs for all data variables
                          coordids(:) ! ID's for coordinates (1-D data)
  ! Cached from open_file for use in write_container
  integer :: nc_nBlocksLon = 0
  integer :: nc_nDims = 3    ! spatial dims of open file: 2 or 3
  logical :: nc_needsIndepWrite = .false.  ! whether open file uses magnetic grid

  ! Multi-time state (UseNetcdfMultiTime mode)
  integer, parameter :: nc_max_types = 50
  integer :: nc_time_record = 1        ! current time record being written (1-based)
  integer :: nc_varid_time = -1        ! varid of the time coordinate variable
  integer :: nc_type_records(nc_max_types)  ! records written per output type so far
  data nc_type_records/nc_max_types*0/
  character(len=600) :: nc_type_filenames(nc_max_types)  ! last filename opened per type
  data nc_type_filenames/nc_max_types*' '/

  ! Leg naming (#NETCDFOPTIONS), computed once per process by build_leg_names:
  ! nc_name_part   = '_<runname>' when netCDFFileName is set (all modes)
  ! nc_restart_part= '_restart_yymmdd_hhmmss' at the restart-file simulation
  !                  time (leg start, NOT wall clock), only for named
  !                  single/daily restart legs, and kept for the whole leg so
  !                  daily rollover files carry it too
  ! nc_guard_active: single/daily runs that are restarts or carry a run name
  !                  refuse (stop_gitm) to overwrite an existing target file;
  !                  unnamed fresh runs keep the status-quo clobber
  logical :: nc_names_built = .false.
  character(len=200) :: nc_name_part = ''
  character(len=40) :: nc_restart_part = ''
  logical :: nc_guard_active = .false.

#endif

contains

  ! ==================================================================
  ! Open a PnetCDF file collectively and define all dims/vars.
  ! Called from write_output.f90 before the block loop.
  ! Supports 2D (lon×lat) and 3D (lon×lat×alt) output types.
  !
  ! UseNetcdfMultiTime=F (default): one file per output type per timestep.
  ! UseNetcdfMultiTime=T: one file per output type for all timesteps;
  !   data variables gain an NF_UNLIMITED time dimension as their last axis.
  ! ==================================================================
  subroutine netcdf_open_file(dir, cType, cTime, cL)
    use ModGITM, only: iCommGITM, iProc
#ifdef HavePNetCDF
    use ModInputs, only: nBlocksLon, nBlocksLat, UseNetcdfMultiTime, NetCdfAppendOption
    use ModOutputProducers, only: get_container_idx, containers, init_output_containers, &
                                  nContainers
    use ModOutputContainer, only: GRID_MAG_2D, GRID_HIME, GRID_GEO_3D
    use ModTime, only: iTimeArray, CurrentTime
    ! Not compiled in standalone. Change 1965 when it is!
    ! use ModTimeConvert, only: iYearBase
#endif
    character(len=*), intent(in) :: dir, cType, cTime
    integer, intent(in) :: cL

#ifdef HavePNetCDF
    character(len=600) :: filename
    integer :: iTypeIdx, nVars, nDims, iV, ierr, ndimid
    integer :: nX, nY, nZ, nGC
    integer(kind=MPI_OFFSET_KIND) :: nLons_g, nLats_g, nAlts_g
    integer(kind=MPI_OFFSET_KIND) :: tstart(1), tcnt(1)
    character(len=80) :: varname
    character(len=20) :: cDate
    ! Double, not the container's real*4: float32 seconds-since-1965 can only
    ! resolve ~128 s at modern dates, which visibly corrupts the time axis.
    real(kind=8) :: time_val(1)

    ! The open path runs before the block loop, so containers may not be built
    ! yet (init_output_containers normally runs during the first fill). Ensure
    ! the schemas exist before reading metadata; the call is idempotent.
    call init_output_containers()

    ! The per-type filename/record caches are fixed-size; overflow would make
    ! types share a cache slot, re-creating (and with the guard below, refusing)
    ! files mid-run. Fail loudly at the first open instead.
    if (nContainers > nc_max_types) &
      call stop_gitm("ModOutputNetCDF: nContainers exceeds nc_max_types"// &
                     "; raise nc_max_types in ModOutputNetCDF")

    call build_leg_names()

    ! Look up variable metadata from the output container
    iTypeIdx = get_container_idx(cType)
    ! An unregistered type reaching open_file is a programming error (e.g. the
    ! Is1D 3D->1D rename producing a name with no container); silently skipping
    ! it would drop the user's requested output.
    if (iTypeIdx < 1) &
      call stop_gitm("ModOutputNetCDF: unknown output type '"//trim(cType)//"'")
    nVars = containers(iTypeIdx)%nVars
    ! Probe types (1DTHM, 0DALL, ...) have no netcdf file layout yet: forcing
    ! them through the 2D path would write a degenerate lon x lat file whose
    ! x-dimension is really altitude. Until a probe layout is designed, halt
    ! loudly (the pre-fix behavior was worse: no file, no error).
    if (cType(1:1) == '0' .or. cType(1:1) == '1') &
      call stop_gitm("ModOutputNetCDF: probe type "//trim(cType)// &
                     " is not supported by the netcdf backend yet; use legacy or mpiio")
    ! On-disk rank comes from gridKind, not the type-name digit: 2D types on
    ! a 3D grid would otherwise mis-parse. Only the HIME pair needs the digit
    ! (2DHME/3DHME share GRID_HIME).
    if (containers(iTypeIdx)%gridKind == GRID_GEO_3D) then
      nDims = 3
    elseif (containers(iTypeIdx)%gridKind == GRID_HIME) then
      read(cType(1:1), '(i1)') nDims
    else
      nDims = 2
    endif
    nc_nBlocksLon = nBlocksLon
    nc_nDims = nDims
    ! Use independent I/O mode for types where not all blocks participate:
    ! magnetic grids (only block 1) or regional outputs (only in region)
    nc_needsIndepWrite = containers(iTypeIdx)%gridKind == GRID_MAG_2D .or. &
                         containers(iTypeIdx)%gridKind == GRID_HIME

    ! Global on-disk dims from the container's own variable shapes (max over
    ! non-axis vars, ghost padding removed) — the same rule the write path
    ! uses for start/cnt — instead of hardcoded per-gridKind grid sizes.
    ! MAG_2D holds one full grid written by a single rank (for 2DMEL this
    ! reproduces the old nMagLons+1 x nMagLats sizing exactly); geographic
    ! grids tile the per-block interior across the block decomposition.
    nX = 0; nY = 0; nZ = 0
    do iV = 1, nVars
      if (.not. containers(iTypeIdx)%vars(iV)%is_axis) then
        nX = max(nX, containers(iTypeIdx)%vars(iV)%shape3(1))
        nY = max(nY, containers(iTypeIdx)%vars(iV)%shape3(2))
        nZ = max(nZ, containers(iTypeIdx)%vars(iV)%shape3(3))
      endif
    enddo
    if (nX == 0 .or. nY == 0 .or. nZ == 0) &
      call stop_gitm("ModOutputNetCDF: type "//trim(cType)// &
                     " has no non-axis variables to size its file from")
    nGC = 0
    if (containers(iTypeIdx)%gridKind /= GRID_MAG_2D) &
      nGC = containers(iTypeIdx)%nGhostCells
    if (containers(iTypeIdx)%gridKind == GRID_MAG_2D) then
      nLons_g = int(nX, MPI_OFFSET_KIND)
      nLats_g = int(nY, MPI_OFFSET_KIND)
    else
      nLons_g = int(nBlocksLon*(nX - 2*nGC), MPI_OFFSET_KIND)
      nLats_g = int(nBlocksLat*(nY - 2*nGC), MPI_OFFSET_KIND)
    endif
    nAlts_g = int(max(nZ - 2*nGC, 1), MPI_OFFSET_KIND)

    if (UseNetcdfMultiTime) then
      ! ----------------------------------------------------------------
      ! Multi-time mode! Either:
      ! - single: one file per output type, accumulates all times.
      ! - daily: one output file, per type, per day
      ! ----------------------------------------------------------------
      ! Build filename first; create-vs-reopen is driven by whether this filename
      ! has been opened in this session.
      ! This correctly handles daily rollover (new date -> new file) and fresh
      ! runs in the same directory (not-yet-used filename -> always create).
      if (NetCdfAppendOption == 'single') then
        filename = trim(dir)//"/"//trim(cType)//trim(nc_name_part)// &
                   trim(nc_restart_part)//".nc"
      else ! Add date to 'daily' files
        write(cDate, '(i4.4,"-",i2.2,"-",i2.2)') iTimeArray(1), iTimeArray(2), iTimeArray(3)
        filename = trim(dir)//"/"//trim(cType)//trim(nc_name_part)// &
                   trim(nc_restart_part)//"_"//trim(cDate)//".nc"
      endif

      if (trim(filename) /= trim(nc_type_filenames(iTypeIdx))) then
        ! ---- New file: create and define everything ----
        nc_type_records(iTypeIdx) = 1
        nc_type_filenames(iTypeIdx) = trim(filename)
        if (nc_guard_active) then
          ! Restart legs and named runs must never overwrite prior data: an
          ! existing target means a same-leg re-run or a scenario branch in
          ! place. NF_NOCLOBBER makes the existence check atomic.
          ierr = nfmpi_create(iCommGITM, trim(filename), &
                              IOR(NF_NOCLOBBER, NF_64BIT_DATA), MPI_INFO_NULL, ncid)
          if (ierr == NF_EEXIST) then
            if (iProc == 0) then
              write(*, *) ""
              write(*, *) "ModOutputNetCDF: output file already exists: ", trim(filename)
              write(*, *) "Refusing to overwrite it. To proceed:"
              write(*, *) " - re-running the same leg: move/delete the file, or"
              write(*, *) "   change netCDFFileName in #NETCDFOPTIONS"
              write(*, *) " - restarting an unnamed run: set netCDFFileName in"
              write(*, *) "   #NETCDFOPTIONS (each leg then gets its own file)"
              write(*, *) " - branching a scenario: use a fresh run directory"
              write(*, *) ""
            endif
            call stop_gitm("netcdf output target exists; see message above")
          endif
        else
          ierr = nfmpi_create(iCommGITM, trim(filename), &
                              IOR(NF_CLOBBER, NF_64BIT_DATA), MPI_INFO_NULL, ncid)
        endif
        if (ierr /= NF_NOERR) then
          ! Collective create over iCommGITM: ierr is consistent across ranks,
          ! and stop_gitm calls MPI_Abort, so halting here cannot deadlock.
          write(*, *) "ModOutputNetCDF: nfmpi_create failed for ", trim(filename), &
            ": ", nfmpi_strerror(ierr)
          call stop_gitm("ModOutputNetCDF: cannot create output file; see message above")
        endif
        call define_file_structure(nDims, nVars, iTypeIdx, nLons_g, nLats_g, nAlts_g, &
                                   iTimeArray, isMultiTime=.true.)
      else
        ! ---- Same file: reopen and append next record ----
        nc_type_records(iTypeIdx) = nc_type_records(iTypeIdx) + 1
        ierr = nfmpi_open(iCommGITM, trim(filename), NF_WRITE, MPI_INFO_NULL, ncid)
        if (ierr /= NF_NOERR) then
          ! Collective reopen over iCommGITM: consistent ierr, MPI_Abort backstop.
          write(*, *) "ModOutputNetCDF: nfmpi_open failed for ", trim(filename), &
            ": ", nfmpi_strerror(ierr)
          call stop_gitm("ModOutputNetCDF: cannot reopen output file; see message above")
        endif
        if (allocated(varids)) deallocate(varids, coordids)
        ndimid = merge(3, 2, nDims == 3)
        allocate(varids(nVars), coordids(ndimid))
        do iV = 1, nVars
          varname = trim(containers(iTypeIdx)%vars(iV)%name)
          call sanitize_nc_name(varname)
          ierr = nfmpi_inq_varid(ncid, trim(varname), varids(iV))
          if (ierr /= NF_NOERR) &
            write(*, *) "netcdf_open: inq_varid '", trim(varname), "' failed: ", &
            nfmpi_strerror(ierr)
        enddo
        ierr = nfmpi_inq_varid(ncid, "time", nc_varid_time)
        ierr = nfmpi_inq_varid(ncid, "lon", coordids(1))
        ierr = nfmpi_inq_varid(ncid, "lat", coordids(2))
        if (nDims == 3) ierr = nfmpi_inq_varid(ncid, "z", coordids(3))
      endif
      nc_time_record = nc_type_records(iTypeIdx)

    else
      ! ----------------------------------------------------------------
      ! Single-time mode (default): one file per output type per timestep.
      ! The run name applies here too (same-directory ensembles); no restart
      ! stamp and no guard; colliding epochs are recomputed and overwritten
      ! in place like the legacy/mpiio backends (not guaranteed bit-identical
      ! on restart; known GITM behavior).
      ! ----------------------------------------------------------------
      nc_time_record = 1
      filename = trim(dir)//"/"//trim(cType)//trim(nc_name_part)// &
                 "_"//cTime(1:cL)//".nc"
      ierr = nfmpi_create(iCommGITM, trim(filename), &
                          IOR(NF_CLOBBER, NF_64BIT_DATA), MPI_INFO_NULL, ncid)
      if (ierr /= NF_NOERR) then
        ! Collective create over iCommGITM: consistent ierr, MPI_Abort backstop.
        write(*, *) "ModOutputNetCDF: nfmpi_create failed for ", trim(filename), &
          ": ", nfmpi_strerror(ierr)
        call stop_gitm("ModOutputNetCDF: cannot create output file; see message above")
      endif
      call define_file_structure(nDims, nVars, iTypeIdx, nLons_g, nLats_g, nAlts_g, &
                                 iTimeArray, isMultiTime=.false.)
    endif

    ! Write this timestep's time value (both modes; nc_time_record=1 for single-time)
    tstart(1) = int(nc_time_record, MPI_OFFSET_KIND)
    tcnt(1) = 1_MPI_OFFSET_KIND
    time_val(1) = CurrentTime
    ierr = nfmpi_put_vara_double_all(ncid, nc_varid_time, tstart, tcnt, time_val)
    if (ierr /= NF_NOERR) &
      write(*, *) "ModOutputNetCDF: failed to write time variable: ", nfmpi_strerror(ierr)

    ! For types with conditional block participation, switch to independent data mode.
    ! This allows only some blocks/processors to write without deadlocking collective I/O.
    ! Examples: magnetic grid types (only block 1) or regional types (only blocks in region).
    ! All processes call this together (it's a collective operation).
    if (nc_needsIndepWrite) then
      ierr = nfmpi_begin_indep_data(ncid)
      if (ierr /= NF_NOERR) &
        write(*, *) "ModOutputNetCDF: nfmpi_begin_indep_data failed: ", nfmpi_strerror(ierr)
    endif

#else
    ! Every rank stops (a rank-0-only stop would leave the others running
    ! with no output); rank 0 prints the message.
    if (iProc == 0) then
      write(*, *) ""
      write(*, *) "ERROR: GITM was built without PnetCDF support."
      write(*, *) "       Recompile in an environment where pnetcdf-config is available."
    endif
    call stop_gitm("NetCDF cannot be used")
#endif
  end subroutine netcdf_open_file

#ifdef HavePNetCDF
  ! ==================================================================
  ! Compute the leg's filename parts once per process (#NETCDFOPTIONS).
  !
  ! The restart stamp uses the restart-file SIMULATION time (= leg start),
  ! not wall clock: filenames stay reproducible, and re-running the same
  ! leg collides with its own earlier file so the NF_NOCLOBBER guard in
  ! netcdf_open_file fires instead of silently minting a second file.
  ! Computed once and never rebuilt, so every file this leg opens
  ! (including daily-rollover files on later days) carries the same stamp.
  ! ==================================================================
  subroutine build_leg_names()
    use ModInputs, only: NetCdfFileName, DoRestart, UseNetcdfMultiTime
    use ModTime, only: RestartTime

    integer :: iLegTime(7)

    if (nc_names_built) return
    nc_names_built = .true.

    ! A name too long for nc_name_part would silently truncate the filename
    ! (and could alias two runs onto one file); refuse instead.
    if (len_trim(NetCdfFileName) > len(nc_name_part) - 1) &
      call stop_gitm("ModOutputNetCDF: netCDFFileName is too long")
    if (len_trim(NetCdfFileName) > 0) &
      nc_name_part = '_'//trim(NetCdfFileName)

    ! Restart stamp: only named single/daily legs get one. Unnamed restarted
    ! legs target the prior leg's file and are stopped by the guard instead
    ! (multi-leg append runs must opt in with a run name).
    if (UseNetcdfMultiTime .and. DoRestart .and. &
        len_trim(NetCdfFileName) > 0) then
      call time_real_to_int(RestartTime, iLegTime)
      write(nc_restart_part, &
            '("_restart_",3i2.2,"_",3i2.2)') &
        mod(iLegTime(1), 100), iLegTime(2), iLegTime(3), &
        iLegTime(4), iLegTime(5), iLegTime(6)
    endif

    nc_guard_active = UseNetcdfMultiTime .and. &
                      (DoRestart .or. len_trim(NetCdfFileName) > 0)
  end subroutine build_leg_names
#endif

  ! ==================================================================
  ! Close the PnetCDF file (collective).
  ! If in independent data mode (magnetic grid types), exit that mode first.
  ! ==================================================================
  subroutine netcdf_close_file()
#ifdef HavePNetCDF
    integer :: ierr

    if (ncid == -1) return
    ! If we entered independent data mode for magnetic grids, exit before close.
    if (nc_needsIndepWrite) then
      ierr = nfmpi_end_indep_data(ncid)
      if (ierr /= NF_NOERR) &
        write(*, *) "ModOutputNetCDF: nfmpi_end_indep_data failed: ", nfmpi_strerror(ierr)
    endif
    ierr = nfmpi_close(ncid)
    if (ierr /= NF_NOERR) &
      write(*, *) "ModOutputNetCDF: nfmpi_close failed: ", nfmpi_strerror(ierr)
    ncid = -1
    if (allocated(varids)) deallocate(varids)
    if (allocated(coordids)) deallocate(coordids)
#endif
  end subroutine netcdf_close_file

  ! netcdf_write_header is intentionally absent.
  ! ActiveBackend%writes_header = .false. for the netcdf backend, so
  ! output_common.f90 never calls write_header and no .header file is created.
  ! All metadata (time, variable names, units) is embedded in the .nc file.

#ifdef HavePNetCDF
  ! ==================================================================
  ! Define all dims, variables, and coordinate metadata for a newly
  ! created NetCDF file. Called from netcdf_open_file after nfmpi_create.
  !
  ! isMultiTime=T: time dim is unlimited, data vars include time axis
  ! isMultiTime=F: time dim is size 1, data vars have spatial dims only,
  !                ISO timestamp written as global "time" attribute
  ! ==================================================================
  subroutine define_file_structure(nDims, nVars, iTypeIdx, nLons_g, nLats_g, nAlts_g, &
                                   iTimeArray, isMultiTime)
    use ModOutputProducers, only: containers
    use ModOutputContainer, only: GRID_MAG_2D, GRID_HIME
    integer, intent(in) :: nDims, nVars, iTypeIdx
    integer(kind=MPI_OFFSET_KIND), intent(in) :: nLons_g, nLats_g, nAlts_g
    integer, intent(in) :: iTimeArray(7)
    logical, intent(in) :: isMultiTime

    integer :: dimid_lon, dimid_lat, dimid_alt, dimid_time
    integer :: dimid_vars(4), ndimid, ndimid_var, iV, ierr
    integer(kind=MPI_OFFSET_KIND) :: attlen
    character(len=80) :: varname
    character(len=80) :: coord_units
    character(len=30) :: coord_str
    character(len=40) :: coordTag
    logical :: isCoordField
    character(len=8)  :: dstr
    character(len=10) :: tstr
    ! Types where only some blocks write (magnetic-grid, regional) leave the
    ! rest of the grid unwritten; fill those elements with NaN instead of 0 so
    ! consumers don't read them as real data.
    logical :: doFillNaN
    real(8) :: fillNaN

    doFillNaN = containers(iTypeIdx)%gridKind == GRID_MAG_2D .or. &
                containers(iTypeIdx)%gridKind == GRID_HIME
    fillNaN = transfer(int(z'7FF8000000000000', 8), 1.0_8)

    ! CF global attributes
    attlen = 6_MPI_OFFSET_KIND
    coord_units = trim(containers(iTypeIdx)%cType)//" output from GITM"
    attlen = int(len_trim(coord_units), MPI_OFFSET_KIND)
    ierr = nfmpi_put_att_text(ncid, NF_GLOBAL, "title", attlen, trim(coord_units))
    coord_units = "GITM (Global Ionosphere-Thermosphere Model)"
    attlen = int(len_trim(coord_units), MPI_OFFSET_KIND)
    ierr = nfmpi_put_att_text(ncid, NF_GLOBAL, "source", attlen, trim(coord_units))
    call date_and_time(date=dstr, time=tstr)
    coord_units = dstr(1:4)//"-"//dstr(5:6)//"-"//dstr(7:8)//" "// &
                  tstr(1:2)//":"//tstr(3:4)//":"//tstr(5:6)//" Created by GITM"
    attlen = int(len_trim(coord_units), MPI_OFFSET_KIND)
    ierr = nfmpi_put_att_text(ncid, NF_GLOBAL, "history", attlen, trim(coord_units))
    coord_units = "CF-1.8"
    attlen = int(len_trim(coord_units), MPI_OFFSET_KIND)
    ierr = nfmpi_put_att_text(ncid, NF_GLOBAL, "Conventions", attlen, trim(coord_units))

    ! Spatial dims
    ierr = nfmpi_def_dim(ncid, "lon", nLons_g, dimid_lon)
    ierr = nfmpi_def_dim(ncid, "lat", nLats_g, dimid_lat)
    if (nDims == 3) ierr = nfmpi_def_dim(ncid, "z", nAlts_g, dimid_alt)

    ! Time dim: unlimited for multi-time, fixed-size-1 for single-time
    if (isMultiTime) then
      ierr = nfmpi_def_dim(ncid, "time", NFMPI_UNLIMITED, dimid_time)
    else
      ierr = nfmpi_def_dim(ncid, "time", 1_MPI_OFFSET_KIND, dimid_time)
    endif

    ! Build dim arrays; data vars include time axis only in multi-time mode
    ndimid = merge(3, 2, nDims == 3)
    ndimid_var = merge(ndimid + 1, ndimid, isMultiTime)
    dimid_vars(1) = dimid_lon
    dimid_vars(2) = dimid_lat
    if (nDims == 3) dimid_vars(3) = dimid_alt
    dimid_vars(ndimid + 1) = dimid_time

    if (allocated(varids)) deallocate(varids, coordids)
    allocate(varids(nVars), coordids(ndimid))

    ! CF auxiliary coordinates: the full Longitude/Latitude/Altitude fields
    ! are stored as data variables; every other data variable references them
    ! via a coordinates= attribute (CF ch. 5) so consumers can associate them
    ! without guessing. Types without those fields (e.g. magnetic grids) get
    ! no tag.
    coordTag = ''
    do iV = 1, nVars
      varname = trim(containers(iTypeIdx)%vars(iV)%name)
      if (varname == 'Longitude' .or. varname == 'Latitude' .or. &
          varname == 'Altitude') coordTag = trim(coordTag)//' '//trim(varname)
    enddo
    coordTag = adjustl(coordTag)

    ! Data variables
    do iV = 1, nVars
      varname = trim(containers(iTypeIdx)%vars(iV)%name)
      isCoordField = varname == 'Longitude' .or. varname == 'Latitude' .or. &
                     varname == 'Altitude'
      call sanitize_nc_name(varname)
      ierr = nfmpi_def_var(ncid, trim(varname), NF_FLOAT, ndimid_var, &
                           dimid_vars(1:ndimid_var), varids(iV))
      if (ierr /= NF_NOERR) &
        write(*, *) "netcdf_open: def_var '", trim(varname), "' failed: ", nfmpi_strerror(ierr)
      if (doFillNaN) then
        ierr = nfmpi_def_var_fill(ncid, varids(iV), 0, fillNaN)
        ierr = nfmpi_put_att_real(ncid, varids(iV), "_FillValue", NF_FLOAT, &
                                  1_MPI_OFFSET_KIND, [real(fillNaN, kind=4)])
      endif
      if (len_trim(containers(iTypeIdx)%vars(iV)%units) > 0) then
        attlen = int(len_trim(containers(iTypeIdx)%vars(iV)%units), MPI_OFFSET_KIND)
        ierr = nfmpi_put_att_text(ncid, varids(iV), "units", attlen, &
                                  trim(containers(iTypeIdx)%vars(iV)%units))
      endif
      if (len_trim(containers(iTypeIdx)%vars(iV)%longName) > 0) then
        attlen = int(len_trim(containers(iTypeIdx)%vars(iV)%longName), MPI_OFFSET_KIND)
        ierr = nfmpi_put_att_text(ncid, varids(iV), "long_name", attlen, &
                                  trim(containers(iTypeIdx)%vars(iV)%longName))
      endif
      ! 'gitm_name' holds the raw (unsanitized) name as written in .bin files,
      ! for matching variables across the two formats
      attlen = int(len_trim(containers(iTypeIdx)%vars(iV)%name), MPI_OFFSET_KIND)
      ierr = nfmpi_put_att_text(ncid, varids(iV), "gitm_name", attlen, &
                                trim(containers(iTypeIdx)%vars(iV)%name))
      if (len_trim(coordTag) > 0 .and. .not. isCoordField) then
        attlen = int(len_trim(coordTag), MPI_OFFSET_KIND)
        ierr = nfmpi_put_att_text(ncid, varids(iV), "coordinates", attlen, &
                                  trim(coordTag))
      endif
    enddo

    ! Time coordinate variable
    ! 1965 is left out in case we need to change the base-year
    write(coord_units, '("seconds since ",i4.4,"-01-01 00:00:00")') 1965
    ierr = nfmpi_def_var(ncid, "time", NF_DOUBLE, 1, [dimid_time], nc_varid_time)
    attlen = int(len_trim(coord_units), MPI_OFFSET_KIND)
    ierr = nfmpi_put_att_text(ncid, nc_varid_time, "units", attlen, trim(coord_units))
    attlen = 8_MPI_OFFSET_KIND
    ierr = nfmpi_put_att_text(ncid, nc_varid_time, "calendar", attlen, "standard")
    attlen = 1_MPI_OFFSET_KIND
    ierr = nfmpi_put_att_text(ncid, nc_varid_time, "axis", attlen, "T")
    attlen = 4_MPI_OFFSET_KIND
    ierr = nfmpi_put_att_text(ncid, nc_varid_time, "standard_name", attlen, "time")

    ! Spatial coordinate variables. Magnetic grids (2DMEL) hold magnetic
    ! lon/lat in these variables: the geographic CF identifiers (axis X/Y,
    ! standard_name longitude/latitude) are reserved for geographic
    ! coordinates and would make CF-aware tools georeference the magnetic
    ! frame, so they are skipped there (units stay as-is).
    ierr = nfmpi_def_var(ncid, "lon", NF_FLOAT, 1, [dimid_lon], coordids(1))
    attlen = 12_MPI_OFFSET_KIND
    ierr = nfmpi_put_att_text(ncid, coordids(1), "units", attlen, "degrees_east")
    if (containers(iTypeIdx)%gridKind /= GRID_MAG_2D) then
      attlen = 1_MPI_OFFSET_KIND
      ierr = nfmpi_put_att_text(ncid, coordids(1), "axis", attlen, "X")
      attlen = 9_MPI_OFFSET_KIND
      ierr = nfmpi_put_att_text(ncid, coordids(1), "standard_name", attlen, "longitude")
    endif

    ierr = nfmpi_def_var(ncid, "lat", NF_FLOAT, 1, [dimid_lat], coordids(2))
    attlen = 13_MPI_OFFSET_KIND
    ierr = nfmpi_put_att_text(ncid, coordids(2), "units", attlen, "degrees_north")
    if (containers(iTypeIdx)%gridKind /= GRID_MAG_2D) then
      attlen = 1_MPI_OFFSET_KIND
      ierr = nfmpi_put_att_text(ncid, coordids(2), "axis", attlen, "Y")
      attlen = 8_MPI_OFFSET_KIND
      ierr = nfmpi_put_att_text(ncid, coordids(2), "standard_name", attlen, "latitude")
    endif

    if (nDims == 3) then
      ierr = nfmpi_def_var(ncid, "z", NF_FLOAT, 1, [dimid_alt], coordids(3))
      attlen = 2_MPI_OFFSET_KIND
      ierr = nfmpi_put_att_text(ncid, coordids(3), "units", attlen, "km")
      attlen = 1_MPI_OFFSET_KIND
      ierr = nfmpi_put_att_text(ncid, coordids(3), "axis", attlen, "Z")
      attlen = 8_MPI_OFFSET_KIND
      ierr = nfmpi_put_att_text(ncid, coordids(3), "standard_name", attlen, "altitude")
      attlen = 2_MPI_OFFSET_KIND
      ierr = nfmpi_put_att_text(ncid, coordids(3), "positive", attlen, "up")
    endif

    ! Single-time: write current time as ISO global attribute
    if (.not. isMultiTime) then
      write(coord_str, '(i4.4,"-",i2.2,"-",i2.2,"T",i2.2,":",i2.2,":",i2.2,"Z")') &
        iTimeArray(1), iTimeArray(2), iTimeArray(3), &
        iTimeArray(4), iTimeArray(5), iTimeArray(6)
      attlen = int(len_trim(coord_str), MPI_OFFSET_KIND)
      ierr = nfmpi_put_att_text(ncid, NF_GLOBAL, "time", attlen, trim(coord_str))
    endif

    ierr = nfmpi_enddef(ncid)
    if (ierr /= NF_NOERR) &
      write(*, *) "ModOutputNetCDF: nfmpi_enddef failed: ", nfmpi_strerror(ierr)
  end subroutine define_file_structure

  ! ==================================================================
  ! Sanitize a string for use as a NetCDF variable name:
  !   ()[]        -> removed
  !   +           -> 'p'
  !   -           -> 'm'
  !   spaces/misc -> '_'
  ! After substitution: collapse '__' runs, strip leading/trailing '_',
  ! then prefix with 'var_' if the result doesn't start with a letter.
  ! ==================================================================
  subroutine sanitize_nc_name(name)
    character(len=*), intent(inout) :: name
    character(len=len(name)) :: tmp
    integer :: i, j, n
    character :: prev

    ! Pass 1: build tmp with substitutions; remove '(' ')' '[' ']'
    j = 0
    do i = 1, len_trim(name)
      select case (name(i:i))
      case ('(', ')', '[', ']')
        ! remove: do not copy
      case ('+')
        j = j + 1
        if (j <= len(tmp)) tmp(j:j) = 'p'
      case ('-')
        j = j + 1
        if (j <= len(tmp)) tmp(j:j) = 'm'
      case (' ', '/', '#', '%', '!', '*', '^', '.', ',', '{', '}', '\', &
            '<', '>', '=', '&', '|', '~', '`', '"', "'", '?', '@', '$', ';', ':')
        j = j + 1
        if (j <= len(tmp)) tmp(j:j) = '_'
      case default
        j = j + 1
        if (j <= len(tmp)) tmp(j:j) = name(i:i)
      end select
    enddo
    if (j < len(tmp)) tmp(j + 1:) = ' '
    n = j

    ! Pass 2: copy tmp -> name, collapsing consecutive '_' into one
    j = 0
    name = ' '
    prev = ' '
    do i = 1, n
      if (tmp(i:i) == '_' .and. prev == '_') cycle
      j = j + 1
      name(j:j) = tmp(i:i)
      prev = tmp(i:i)
    enddo

    ! Strip leading '_'
    do while (len_trim(name) > 0 .and. name(1:1) == '_')
      name = name(2:)
    enddo

    ! Strip trailing '_'
    n = len_trim(name)
    do while (n > 0 .and. name(n:n) == '_')
      name(n:n) = ' '
      n = n - 1
    enddo

    ! NetCDF classic names must start with a letter; prefix with 'var_' if not
    if (len_trim(name) > 0) then
      select case (name(1:1))
      case ('a':'z', 'A':'Z')
        ! ok
      case default
        name = 'var_'//name
      end select
    endif

  end subroutine sanitize_nc_name
#endif

  ! ---------------------------------------------------------------------------
  ! netcdf_write_container: write a filled OutputContainer to the open PnetCDF file.
  !
  ! Called from output_common.f90 after fill_<type>().  The file is already open
  ! (netcdf_open_file ran before the block loop).
  !
  ! Two I/O modes, driven by nc_needsIndepWrite (set at open_file time):
  !
  !   nc_needsIndepWrite=T (MAG_2D / usesMagGrid):
  !     File is in independent data mode.  Only the participating rank writes;
  !     others return early.  start=(1,1): single block covers the full grid.
  !
  !   nc_needsIndepWrite=F (GEO_2D block-partitioned):
  !     File is in collective data mode.  All ranks participate; each writes
  !     its block at the correct global offset.
  !     iBlockLon_0 = mod(iBlock-1, nc_nBlocksLon)
  !     iBlockLat_0 = (iBlock-1) / nc_nBlocksLon
  !     start(1) = iBlockLon_0 * nX + 1
  !     start(2) = iBlockLat_0 * nY + 1
  !
  ! Axis vars (is_axis=T) are 1D in the container but 2D in the data variables:
  !   lon axis (shape3(1)==nX): replicate along Y -> buf2d(:, iy)
  !   lat axis (shape3(1)==nY): replicate along X -> buf2d(ix, :)
  !
  ! 1D coord variables (lon, lat) are written once per file (nc_time_record==1)
  ! from the first axis var encountered for each axis, at the block's offset.
  ! ---------------------------------------------------------------------------
  subroutine netcdf_write_container(c, iBlock, iTimeArray)
    use ModOutputContainer, only: OutputContainer, output_kind, &
                                  GRID_GEO_3D, GRID_GEO_2D, GRID_MAG_2D, GRID_HIME
#ifdef HavePNetCDF
    use ModInputs, only: UseNetcdfMultiTime
#endif
    type(OutputContainer), intent(inout) :: c
    integer, intent(in) :: iBlock
    integer, intent(in) :: iTimeArray(7)

#ifdef HavePNetCDF
    ! start/cnt are sized for 3D+time; 2D paths use start(1:3)/cnt(1:3) slices.
    integer(kind=MPI_OFFSET_KIND) :: start(4), cnt(4), cstart(1), ccnt(1)
    integer :: nX, nY, nZ, i, ix, iy, iz, iAxis, ierr
    integer :: iBlockLon_0, iBlockLat_0
    integer :: nGC, ax_lo, ay_lo, az_lo
    logical :: is3d, doStrip
    real(kind=4), allocatable :: buf2d(:, :), buf3d(:, :, :), coord1d(:)

    if (.not. c%this_rank_writes) return
    if (ncid == -1) return
    if (c%nVars == 0) return

    ! Ghost cell stripping: netcdf writes interior only for geographic grids.
    ! Magnetic grids (GRID_MAG_2D) have no ghost concept and are never stripped.
    doStrip = (c%nGhostCells > 0) .and. &
              (c%gridKind == GRID_GEO_2D .or. c%gridKind == GRID_GEO_3D &
               .or. c%gridKind == GRID_HIME)
    nGC = merge(c%nGhostCells, 0, doStrip)

    ! Derive grid dims from non-axis vars.  Container shapes include any ghost
    ! padding; subtract it here so nX/nY/nZ are interior sizes (the on-disk shape).
    nX = 0; nY = 0; nZ = 0
    do i = 1, c%nVars
      if (.not. c%vars(i)%is_axis) then
        nX = max(nX, c%vars(i)%shape3(1))
        nY = max(nY, c%vars(i)%shape3(2))
        nZ = max(nZ, c%vars(i)%shape3(3))
      endif
    enddo
    if (nX == 0 .or. nY == 0 .or. nZ == 0) return
    nX = nX - 2*nGC
    nY = nY - 2*nGC
    nZ = nZ - 2*nGC

    ! Lower-bound offsets into container storage when stripping (1-based slice start).
    ax_lo = nGC + 1
    ay_lo = nGC + 1
    az_lo = nGC + 1

    ! On-disk rank, not gridKind: GRID_HIME carries both 2DHME (rank 2) and
    ! 3DHME (rank 3), so keying on GEO_3D alone would push 3DHME through the
    ! 2-D path — 2-element start/cnt against 3-dim variables, alt level 1 only.
    is3d = (nc_nDims == 3)

    if (nc_needsIndepWrite .and. c%gridKind == GRID_MAG_2D) then
      ! MAG_2D: single rank writes the full magnetic grid starting at (1,1).
      start(1) = 1_MPI_OFFSET_KIND
      start(2) = 1_MPI_OFFSET_KIND
    else
      ! GEO_2D / GEO_3D / GRID_HIME: each block writes at its global lon/lat offset.
      ! For GRID_HIME (regional), only participating blocks call this, but they still
      ! write at their correct position in the full global grid.
      iBlockLon_0 = mod(iBlock - 1, nc_nBlocksLon)
      iBlockLat_0 = (iBlock - 1)/nc_nBlocksLon
      start(1) = int(iBlockLon_0*nX, MPI_OFFSET_KIND) + 1_MPI_OFFSET_KIND
      start(2) = int(iBlockLat_0*nY, MPI_OFFSET_KIND) + 1_MPI_OFFSET_KIND
    endif
    cnt(1) = int(nX, MPI_OFFSET_KIND)
    cnt(2) = int(nY, MPI_OFFSET_KIND)

    if (is3d) then
      start(3) = 1_MPI_OFFSET_KIND            ! alt always starts at 1
      start(4) = int(nc_time_record, MPI_OFFSET_KIND)
      cnt(3) = int(nZ, MPI_OFFSET_KIND)
      cnt(4) = 1_MPI_OFFSET_KIND
    else
      start(3) = int(nc_time_record, MPI_OFFSET_KIND)
      cnt(3) = 1_MPI_OFFSET_KIND
    endif

    ! Axis vars: first is_axis var = lon, second = lat (position-based, safe
    ! when nLons==nLats).  Non-axis vars are written as full 2D or 3D arrays.
    if (is3d) then
      allocate(buf3d(nX, nY, nZ))
    else
      allocate(buf2d(nX, nY))
    endif

    iAxis = 0
    do i = 1, c%nVars
      if (c%vars(i)%is_axis) then
        iAxis = iAxis + 1
        if (is3d) then
          if (iAxis == 1) then
            ! Longitude axis: replicate along Y and Z.
            do iz = 1, nZ
              do iy = 1, nY
                buf3d(:, iy, iz) = real(c%vars(i)%data(ax_lo:ax_lo + nX - 1, 1, 1), kind=4)
              enddo
            enddo
          else
            ! Latitude axis: replicate along X and Z.
            do iz = 1, nZ
              do ix = 1, nX
                buf3d(ix, :, iz) = real(c%vars(i)%data(ay_lo:ay_lo + nY - 1, 1, 1), kind=4)
              enddo
            enddo
          endif
        else
          if (iAxis == 1) then
            ! Longitude axis: replicate along Y.
            do iy = 1, nY
              buf2d(:, iy) = real(c%vars(i)%data(ax_lo:ax_lo + nX - 1, 1, 1), kind=4)
            enddo
          else
            ! Latitude axis: replicate along X.
            do ix = 1, nX
              buf2d(ix, :) = real(c%vars(i)%data(ay_lo:ay_lo + nY - 1, 1, 1), kind=4)
            enddo
          endif
        endif
      else
        if (is3d) then
          buf3d = real(c%vars(i)%data(ax_lo:ax_lo + nX - 1, &
                                      ay_lo:ay_lo + nY - 1, &
                                      az_lo:az_lo + nZ - 1), kind=4)
        else
          buf2d = real(c%vars(i)%data(ax_lo:ax_lo + nX - 1, &
                                      ay_lo:ay_lo + nY - 1, 1), kind=4)
        endif
      endif

      if (is3d) then
        if (nc_needsIndepWrite) then
          if (UseNetcdfMultiTime) then
            ierr = nfmpi_put_vara_real(ncid, varids(i), start(1:4), cnt(1:4), buf3d)
          else
            ierr = nfmpi_put_vara_real(ncid, varids(i), start(1:3), cnt(1:3), buf3d)
          endif
        else
          if (UseNetcdfMultiTime) then
            ierr = nfmpi_put_vara_real_all(ncid, varids(i), start(1:4), cnt(1:4), buf3d)
          else
            ierr = nfmpi_put_vara_real_all(ncid, varids(i), start(1:3), cnt(1:3), buf3d)
          endif
        endif
      else
        if (nc_needsIndepWrite) then
          if (UseNetcdfMultiTime) then
            ierr = nfmpi_put_vara_real(ncid, varids(i), start(1:3), cnt(1:3), buf2d)
          else
            ierr = nfmpi_put_vara_real(ncid, varids(i), start(1:2), cnt(1:2), buf2d)
          endif
        else
          if (UseNetcdfMultiTime) then
            ierr = nfmpi_put_vara_real_all(ncid, varids(i), start(1:3), cnt(1:3), buf2d)
          else
            ierr = nfmpi_put_vara_real_all(ncid, varids(i), start(1:2), cnt(1:2), buf2d)
          endif
        endif
      endif
      if (ierr /= NF_NOERR) &
        write(*, *) "netcdf_write_container: put_vara failed for ", &
        trim(c%vars(i)%name), ": ", nfmpi_strerror(ierr)
    enddo

    if (is3d) then
      deallocate(buf3d)
    else
      deallocate(buf2d)
    endif

    ! 1D coordinate variables written once per file (nc_time_record==1).
    ! lon/lat use block offset; z uses absolute [1,nZ] (all blocks span full alt).
    if (nc_time_record == 1) then
      ! lon coord (var 1 by convention)
      ccnt(1) = int(nX, MPI_OFFSET_KIND)
      cstart(1) = start(1)
      allocate(coord1d(nX))
      coord1d = real(c%vars(1)%data(ax_lo:ax_lo + nX - 1, 1, 1), kind=4)
      if (nc_needsIndepWrite) then
        ierr = nfmpi_put_vara_real(ncid, coordids(1), cstart, ccnt, coord1d)
      else
        ierr = nfmpi_put_vara_real_all(ncid, coordids(1), cstart, ccnt, coord1d)
      endif
      deallocate(coord1d)

      ! lat coord (var 2 by convention)
      ccnt(1) = int(nY, MPI_OFFSET_KIND)
      cstart(1) = start(2)
      allocate(coord1d(nY))
      coord1d = real(c%vars(2)%data(ay_lo:ay_lo + nY - 1, 1, 1), kind=4)
      if (nc_needsIndepWrite) then
        ierr = nfmpi_put_vara_real(ncid, coordids(2), cstart, ccnt, coord1d)
      else
        ierr = nfmpi_put_vara_real_all(ncid, coordids(2), cstart, ccnt, coord1d)
      endif
      deallocate(coord1d)

      ! z coord (var 3 by convention for 3D; altitude in km)
      if (is3d) then
        ccnt(1) = int(nZ, MPI_OFFSET_KIND)
        cstart(1) = 1_MPI_OFFSET_KIND
        allocate(coord1d(nZ))
        coord1d = real(c%vars(3)%data(ax_lo, ay_lo, az_lo:az_lo + nZ - 1)/1000.0d0, kind=4)
        if (nc_needsIndepWrite) then
          ierr = nfmpi_put_vara_real(ncid, coordids(3), cstart, ccnt, coord1d)
        else
          ierr = nfmpi_put_vara_real_all(ncid, coordids(3), cstart, ccnt, coord1d)
        endif
        deallocate(coord1d)
      endif
    endif
#endif
  end subroutine netcdf_write_container

end module ModOutputNetCDF
