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
! All registered variables (including Longitude/Latitude/Altitude from
! add_coord_vars) are written as data variables via varids(:).  Separate
! 1-D NetCDF coordinate variables are not defined to avoid name conflicts.
!
! Build requirements (auto-detected by src/Makefile via pnetcdf-config):
!   Compile: -I$(PNETCDF_INC)  and  -DHavePNetCDF
!   Link:    -L$(PNETCDF_LIB) -lpnetcdf

module ModOutputNetCDF

  use ModOutputRegistry, only: OutputTypeInfo
#ifdef HavePNetCDF
  use ModMpi
  use pnetcdf
#endif

  implicit none

#ifdef HavePNetCDF
  integer :: ncid = -1       ! PnetCDF file handle; -1 = no file open
  integer, allocatable :: varids(:), & ! variable IDs for all data variables
                          coordids(:) ! ID's for coordinates (1-D data)
  ! Cached from open_file for use in write_block
  integer :: nc_nBlocksLon = 0
  integer :: nc_nDims = 3    ! spatial dims of open file: 2 or 3
  logical :: nc_needsIndepWrite = .false.  ! whether open file uses magnetic grid

  ! Multi-time state (UseNetcdfMultiTime mode)
  integer, parameter :: nc_max_types = 50
  integer :: nc_time_record = 1        ! current time record being written (1-based)
  integer :: nc_varid_time = -1        ! varid of the time coordinate variable
  integer :: nc_type_records(nc_max_types)  ! records written per output type so far
  data nc_type_records/nc_max_types*0/
  character(len=300) :: nc_multitime_filename = ''  ! filename for current type
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
    use ModInputs, only: nBlocksLon, nBlocksLat, UseNetcdfMultiTime
    use ModSizeGitm, only: nLons, nLats, nAlts
    use ModElectrodynamics, only: nMagLons, nMagLats
    use ModOutputRegistry, only: find_output_type, RegisteredTypes
    use ModTime, only: iTimeArray, CurrentTime
    ! Not compiled in standalone. Change 1965 when it is!
    ! use ModTimeConvert, only: iYearBase
#endif
    character(len=*), intent(in) :: dir, cType, cTime
    integer, intent(in) :: cL

#ifdef HavePNetCDF
    character(len=300) :: filename
    integer :: iTypeIdx, nVars, nDims, iV, ierr, ndimid
    integer(kind=MPI_OFFSET_KIND) :: nLons_g, nLats_g, nAlts_g
    integer(kind=MPI_OFFSET_KIND) :: tstart(1), tcnt(1)
    character(len=80) :: varname
    real(kind=8) :: time_val(1)

    ! Look up variable metadata from registry
    iTypeIdx = find_output_type(cType)
    if (iTypeIdx < 1) then
      write(*, *) "ModOutputNetCDF: unknown output type '", trim(cType), "'"
      ncid = -1
      return
    endif
    nVars = RegisteredTypes(iTypeIdx)%nVars
    nDims = RegisteredTypes(iTypeIdx)%nDims
    nc_nBlocksLon = nBlocksLon
    nc_nDims = nDims
    ! Use independent I/O mode for types where not all blocks participate:
    ! magnetic grids (only block 1) or regional outputs (only in region)
    nc_needsIndepWrite = RegisteredTypes(iTypeIdx)%usesMagGrid .or. &
                         RegisteredTypes(iTypeIdx)%isRegional

    ! Only 2D and 3D output supported in PnetCDF backend
    if (nDims < 2) then
      ncid = -1
      return
    endif

    if (RegisteredTypes(iTypeIdx)%usesMagGrid) then
      nLons_g = int(nMagLons + 1, MPI_OFFSET_KIND)
      nLats_g = int(nMagLats, MPI_OFFSET_KIND)
    else
      nLons_g = int(nBlocksLon*nLons, MPI_OFFSET_KIND)
      nLats_g = int(nBlocksLat*nLats, MPI_OFFSET_KIND)
    endif
    nAlts_g = int(nAlts, MPI_OFFSET_KIND)

    if (UseNetcdfMultiTime) then
      ! ----------------------------------------------------------------
      ! Multi-time mode: one file per output type accumulates all times.
      ! ----------------------------------------------------------------
      if (iTypeIdx > nc_max_types) iTypeIdx = nc_max_types
      nc_type_records(iTypeIdx) = nc_type_records(iTypeIdx) + 1
      nc_time_record = nc_type_records(iTypeIdx)
      filename = trim(dir)//"/"//trim(cType)//".nc"
      nc_multitime_filename = filename

      if (nc_time_record == 1) then
        ! ---- First record: create file and define everything ----
        ierr = nfmpi_create(iCommGITM, trim(filename), &
                            IOR(NF_CLOBBER, NF_64BIT_DATA), MPI_INFO_NULL, ncid)
        if (ierr /= NF_NOERR) then
          write(*, *) "ModOutputNetCDF: nfmpi_create failed: ", nfmpi_strerror(ierr)
          ncid = -1; return
        endif
        call define_file_structure(nDims, nVars, iTypeIdx, nLons_g, nLats_g, nAlts_g, &
                                   iTimeArray, isMultiTime=.true.)
      else
        ! ---- Subsequent records: reopen existing file and re-fetch varids ----
        ierr = nfmpi_open(iCommGITM, trim(filename), NF_WRITE, MPI_INFO_NULL, ncid)
        if (ierr /= NF_NOERR) then
          write(*, *) "ModOutputNetCDF: nfmpi_open failed: ", nfmpi_strerror(ierr)
          ncid = -1; return
        endif
        if (allocated(varids)) deallocate(varids, coordids)
        ndimid = merge(3, 2, nDims == 3)
        allocate(varids(nVars), coordids(ndimid))
        do iV = 1, nVars
          varname = trim(RegisteredTypes(iTypeIdx)%vars(iV)%shortName)
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

    else
      ! ----------------------------------------------------------------
      ! Single-time mode (default): one file per output type per timestep
      ! ----------------------------------------------------------------
      nc_time_record = 1
      filename = trim(dir)//"/"//trim(cType)//"_"//cTime(1:cL)//".nc"
      ierr = nfmpi_create(iCommGITM, trim(filename), &
                          IOR(NF_CLOBBER, NF_64BIT_DATA), MPI_INFO_NULL, ncid)
      if (ierr /= NF_NOERR) then
        write(*, *) "ModOutputNetCDF: nfmpi_create failed: ", nfmpi_strerror(ierr)
        ncid = -1; return
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
    if (RegisteredTypes(iTypeIdx)%usesMagGrid .or. RegisteredTypes(iTypeIdx)%isRegional) then
      ierr = nfmpi_begin_indep_data(ncid)
      if (ierr /= NF_NOERR) &
        write(*, *) "ModOutputNetCDF: nfmpi_begin_indep_data failed: ", nfmpi_strerror(ierr)
    endif

#else
    if (iProc == 1) then
      write(*, *) ""
      write(*, *) "ERROR: GITM was built without PnetCDF support."
      write(*, *) "       Recompile in an environment where pnetcdf-config is available."
      call stop_gitm("NetCDF cannot be used")
    endif
#endif
  end subroutine netcdf_open_file

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

  ! ==================================================================
  ! Write one block's interior data to the open PnetCDF file.
  !
  ! buffer(nV, nX, nY, nZ) — full buffer for this block
  ! iBLK                   — 1-based global block index
  !
  ! For types with conditionalIO flag (magnetic grids, regional outputs):
  !   - File opened in independent data mode (via nfmpi_begin_indep_data)
  !   - Uses non-collective nfmpi_put_vara_double; only block 1 (or blocks in
  !     region) write; others skip without deadlock
  ! For regular block-partitioned types:
  !   - File in collective mode; all blocks call collective nfmpi_put_vara_double_all
  !
  ! Dispatches on nc_nDims cached from open_file:
  !   3D: strip 2 ghost cells each side; write to (lon×lat×alt[×time]) var
  !   2D: no ghost cells (nZ=1); write to (lon×lat[×time]) var
  ! ==================================================================
  subroutine netcdf_write_block(iUnit, buffer, nV, nX, nY, nZ, iBLK)
    use ModGITM, only: Longitude, Latitude, Altitude_GB
    use ModConst, only: cRadToDeg
    use ModInputs, only: UseNetcdfMultiTime
    integer, intent(in) :: iUnit, nV, nX, nY, nZ, iBLK
    real, intent(in) :: buffer(nV, nX, nY, nZ)

#ifdef HavePNetCDF
    integer :: iBlockLon_0, iBlockLat_0
    integer(kind=MPI_OFFSET_KIND) :: start(4), cnt(4)
    integer(kind=MPI_OFFSET_KIND) :: cstart(1), ccnt(1)
    integer :: nLi, nLai, nAi, lo_x, hi_x, lo_y, hi_y, lo_z, hi_z
    integer :: ndimid, iV, ierr
    real, allocatable :: slice(:, :, :)
    real(kind=8), allocatable :: coord1d(:)

    if (ncid == -1) return

    ! 0-based block spatial indices
    iBlockLon_0 = mod(iBLK - 1, nc_nBlocksLon)
    iBlockLat_0 = (iBLK - 1)/nc_nBlocksLon

    if (nc_nDims == 3) then
      ! 3D type: buffer has 2 ghost cells each side
      nLi = nX - 4; lo_x = 3; hi_x = nX - 2
      nLai = nY - 4; lo_y = 3; hi_y = nY - 2
      nAi = nZ - 4; lo_z = 3; hi_z = nZ - 2
      ndimid = 3
    else
      ! 2D type: buffer has no ghost cells (nX=nLons, nY=nLats, nZ=1)
      nLi = nX; lo_x = 1; hi_x = nX
      nLai = nY; lo_y = 1; hi_y = nY
      nAi = nZ; lo_z = 1; hi_z = nZ
      ndimid = 2
    endif

    start(1) = int(iBlockLon_0*nLi, MPI_OFFSET_KIND) + 1_MPI_OFFSET_KIND
    start(2) = int(iBlockLat_0*nLai, MPI_OFFSET_KIND) + 1_MPI_OFFSET_KIND
    start(3) = 1_MPI_OFFSET_KIND
    start(4) = int(nc_time_record, MPI_OFFSET_KIND)

    cnt(1) = int(nLi, MPI_OFFSET_KIND)
    cnt(2) = int(nLai, MPI_OFFSET_KIND)
    cnt(3) = int(nAi, MPI_OFFSET_KIND)
    cnt(4) = 1_MPI_OFFSET_KIND

    allocate(slice(nLi, nLai, nAi))

    ! ndimid_var: spatial dims only (single-time) or spatial+time (multi-time)
    ! For 2D: ndimid=2, so start/cnt(3) = time; for 3D: ndimid=3, start/cnt(4) = time.
    ! We always pass start(1:ndimid+1) / cnt(1:ndimid+1) in multi-time mode
    ! so that start(ndimid+1)=nc_time_record, cnt(ndimid+1)=1.

    if (nc_needsIndepWrite) then
      ! For magnetic grid types, use non-collective I/O to avoid MPI deadlock
      ! (only block 1 has data to write, but all processes must reach close_file() together).
      do iV = 1, nV
        slice = buffer(iV, lo_x:hi_x, lo_y:hi_y, lo_z:hi_z)
        if (UseNetcdfMultiTime) then
          ierr = nfmpi_put_vara_double(ncid, varids(iV), &
                                       start(1:ndimid + 1), cnt(1:ndimid + 1), slice)
        else
          ierr = nfmpi_put_vara_double(ncid, varids(iV), &
                                       start(1:ndimid), cnt(1:ndimid), slice)
        endif
        if (ierr /= NF_NOERR) &
          write(*, *) "netcdf_write_block: put_vara_double failed for var ", iV, &
          ": ", nfmpi_strerror(ierr)
      enddo
    else
      ! Collective put for regular block-partitioned types
      do iV = 1, nV
        slice = buffer(iV, lo_x:hi_x, lo_y:hi_y, lo_z:hi_z)
        if (UseNetcdfMultiTime) then
          ierr = nfmpi_put_vara_double_all(ncid, varids(iV), &
                                           start(1:ndimid + 1), cnt(1:ndimid + 1), slice)
        else
          ierr = nfmpi_put_vara_double_all(ncid, varids(iV), &
                                           start(1:ndimid), cnt(1:ndimid), slice)
        endif
        if (ierr /= NF_NOERR) &
          write(*, *) "netcdf_write_block: put_vara_double_all failed for var ", iV, &
          ": ", nfmpi_strerror(ierr)
      enddo

      ! --- 1-D coordinate variables --- (Not present in magnetic/needsIndepWrite outputs)
      ! Written only on the first time record; coordinates don't change between timesteps.
      if (nc_time_record == 1) then
        ! lon: each block writes its own interior lon range (radians east)
        cstart(1) = start(1)
        ccnt(1) = cnt(1)
        allocate(coord1d(nLi))
        coord1d = real(Longitude(1:nLi, 1), kind=8)!*cRadToDeg
        ierr = nfmpi_put_vara_double_all(ncid, coordids(1), cstart, ccnt, coord1d)
        deallocate(coord1d)

        ! lat: each block writes its own interior lat range (rad north)
        cstart(1) = start(2)
        ccnt(1) = cnt(2)
        allocate(coord1d(nLai))
        coord1d = real(Latitude(1:nLai, 1), kind=8)!*cRadToDeg
        ierr = nfmpi_put_vara_double_all(ncid, coordids(2), cstart, ccnt, coord1d)
        deallocate(coord1d)

        ! z (3D only): all blocks cover the same altitude levels (m)
        if (nc_nDims == 3) then
          cstart(1) = 1_MPI_OFFSET_KIND
          ccnt(1) = cnt(3)
          allocate(coord1d(nAi))
          coord1d = real(Altitude_GB(1, 1, 1:nAi, 1), kind=8)!/1000.0d0
          ierr = nfmpi_put_vara_double_all(ncid, coordids(3), cstart, ccnt, coord1d)
          deallocate(coord1d)
        endif
      endif

    endif

    deallocate(slice)
#endif
  end subroutine netcdf_write_block

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
    use ModOutputRegistry, only: RegisteredTypes
    integer, intent(in) :: nDims, nVars, iTypeIdx
    integer(kind=MPI_OFFSET_KIND), intent(in) :: nLons_g, nLats_g, nAlts_g
    integer, intent(in) :: iTimeArray(7)
    logical, intent(in) :: isMultiTime

    integer :: dimid_lon, dimid_lat, dimid_alt, dimid_time
    integer :: dimid_vars(4), ndimid, ndimid_var, iV, ierr
    integer(kind=MPI_OFFSET_KIND) :: attlen
    character(len=80) :: varname
    character(len=40) :: coord_units
    character(len=30) :: coord_str
    character(len=8) :: dstr, tstr

    ! CF global attributes
    attlen = 6_MPI_OFFSET_KIND
    coord_units = trim(RegisteredTypes(iTypeIdx)%code)//" output from GITM"
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

    ! Data variables
    do iV = 1, nVars
      varname = trim(RegisteredTypes(iTypeIdx)%vars(iV)%shortName)
      call sanitize_nc_name(varname)
      ierr = nfmpi_def_var(ncid, trim(varname), NF_DOUBLE, ndimid_var, &
                           dimid_vars(1:ndimid_var), varids(iV))
      if (ierr /= NF_NOERR) &
        write(*, *) "netcdf_open: def_var '", trim(varname), "' failed: ", nfmpi_strerror(ierr)
      if (len_trim(RegisteredTypes(iTypeIdx)%vars(iV)%units) > 0) then
        attlen = int(len_trim(RegisteredTypes(iTypeIdx)%vars(iV)%units), MPI_OFFSET_KIND)
        ierr = nfmpi_put_att_text(ncid, varids(iV), "units", attlen, &
                                  trim(RegisteredTypes(iTypeIdx)%vars(iV)%units))
      endif
      if (len_trim(RegisteredTypes(iTypeIdx)%vars(iV)%longName) > 0) then
        attlen = int(len_trim(RegisteredTypes(iTypeIdx)%vars(iV)%longName), MPI_OFFSET_KIND)
        ierr = nfmpi_put_att_text(ncid, varids(iV), "long_name", attlen, &
                                  trim(RegisteredTypes(iTypeIdx)%vars(iV)%longName))
      endif
    enddo

    ! Time coordinate variable
    write(coord_units, '("seconds since ",i4.4,"-01-01 00:00:00")') 1965
    ierr = nfmpi_def_var(ncid, "time", NF_DOUBLE, 1, [dimid_time], nc_varid_time)
    attlen = int(len_trim(coord_units), MPI_OFFSET_KIND)
    ierr = nfmpi_put_att_text(ncid, nc_varid_time, "units", attlen, trim(coord_units))
    attlen = 8_MPI_OFFSET_KIND
    ierr = nfmpi_put_att_text(ncid, nc_varid_time, "calendar", attlen, "standard")

    ! Spatial coordinate variables
    ierr = nfmpi_def_var(ncid, "lon", NF_DOUBLE, 1, [dimid_lon], coordids(1))
    attlen = 12_MPI_OFFSET_KIND
    ierr = nfmpi_put_att_text(ncid, coordids(1), "units", attlen, "degrees_east")

    ierr = nfmpi_def_var(ncid, "lat", NF_DOUBLE, 1, [dimid_lat], coordids(2))
    attlen = 13_MPI_OFFSET_KIND
    ierr = nfmpi_put_att_text(ncid, coordids(2), "units", attlen, "degrees_north")

    if (nDims == 3) then
      ierr = nfmpi_def_var(ncid, "z", NF_DOUBLE, 1, [dimid_alt], coordids(3))
      attlen = 2_MPI_OFFSET_KIND
      ierr = nfmpi_put_att_text(ncid, coordids(3), "units", attlen, "km")
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
  ! Replace characters invalid in NetCDF variable names with '_'
  ! ==================================================================
  subroutine sanitize_nc_name(name)
    character(len=*), intent(inout) :: name
    integer :: i
    do i = 1, len_trim(name)
      select case (name(i:i))
      case (' ', '/', '(', ')', '[', ']', '#', '%', '!', '-', '*', '^', '.', &
            ',', '{', '}', '\', '<', '>', '=', '&', '|', '~', '`', '"', "'", '?', '@', '$', ';', ':')
        name(i:i) = '_'
      case ('+')
        name(i:i) = '_plus'
      end select
    enddo
    ! NetCDF classic names must start with a letter; prefix with 'v_' if not
    if (len_trim(name) > 0) then
      select case (name(1:1))
      case ('a':'z', 'A':'Z')
        ! ok
      case default
        name = 'v_'//name
      end select
    endif
  end subroutine sanitize_nc_name
#endif

end module ModOutputNetCDF
