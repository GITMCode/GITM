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
  integer, allocatable :: varids(:)  ! variable IDs for all data variables

  ! Cached from open_file for use in write_block
  integer :: nc_nBlocksLon = 0
  integer :: nc_nDims = 3    ! spatial dims of open file: 2 or 3
#endif

contains

  ! ==================================================================
  ! Open a PnetCDF file collectively and define all dims/vars.
  ! Called from write_output.f90 before the block loop.
  ! Supports 2D (lon×lat) and 3D (lon×lat×alt) output types.
  ! ==================================================================
  subroutine netcdf_open_file(dir, cType, cTime, cL)
    use ModGITM, only: iCommGITM, iProc
#ifdef HavePNetCDF
    use ModInputs, only: nBlocksLon, nBlocksLat
    use ModSizeGitm, only: nLons, nLats, nAlts
    use ModOutputRegistry, only: find_output_type, RegisteredTypes
    use ModTime, only: iTimeArray
#endif
    character(len=*), intent(in) :: dir, cType, cTime
    integer, intent(in) :: cL

#ifdef HavePNetCDF
    character(len=300) :: filename
    integer :: iTypeIdx, nVars, nDims, iV, ierr, ndimid
    integer :: dimid_lon, dimid_lat, dimid_alt, dimid_vars(3)
    integer(kind=MPI_OFFSET_KIND) :: nLons_g, nLats_g, nAlts_g, attlen
    character(len=80) :: varname
    character(len=30) :: time_str

    ! Build filename
    filename = trim(dir)//"/"//trim(cType)//"_"//cTime(1:cL)//".nc"

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

    ! Only 2D and 3D output supported in PnetCDF backend
    if (nDims < 2) then
      ncid = -1
      return
    endif

    nLons_g = int(nBlocksLon * nLons, MPI_OFFSET_KIND)
    nLats_g = int(nBlocksLat * nLats, MPI_OFFSET_KIND)
    nAlts_g = int(nAlts, MPI_OFFSET_KIND)

    ! --- Create file (collective across all ranks) ---
    ierr = nfmpi_create(iCommGITM, trim(filename), &
                        IOR(NF_CLOBBER, NF_64BIT_DATA), &
                        MPI_INFO_NULL, ncid)
    if (ierr /= NF_NOERR) then
      write(*, *) "ModOutputNetCDF: nfmpi_create failed: ", nfmpi_strerror(ierr)
      ncid = -1
      return
    endif

    ! --- Define dimensions ---
    ierr = nfmpi_def_dim(ncid, "lon", nLons_g, dimid_lon)
    ierr = nfmpi_def_dim(ncid, "lat", nLats_g, dimid_lat)
    if (nDims == 3) then
      ierr = nfmpi_def_dim(ncid, "alt", nAlts_g, dimid_alt)
    endif

    ! --- Global time attribute ---
    write(time_str, '(i4.4,"-",i2.2,"-",i2.2,"T",i2.2,":",i2.2,":",i2.2,"Z")') &
      iTimeArray(1), iTimeArray(2), iTimeArray(3), &
      iTimeArray(4), iTimeArray(5), iTimeArray(6)
    attlen = int(len_trim(time_str), MPI_OFFSET_KIND)
    ierr = nfmpi_put_att_text(ncid, NF_GLOBAL, "time", attlen, trim(time_str))

    ! --- Dimension IDs for data variables ---
    ! All registered vars (including Longitude/Latitude/Altitude from add_coord_vars)
    ! are defined as nDims-dimensional data variables to avoid name conflicts with
    ! any separately defined 1-D coordinate variables.
    if (nDims == 3) then
      ndimid = 3
      dimid_vars(1) = dimid_lon
      dimid_vars(2) = dimid_lat
      dimid_vars(3) = dimid_alt
    else
      ndimid = 2
      dimid_vars(1) = dimid_lon
      dimid_vars(2) = dimid_lat
    endif

    if (allocated(varids)) deallocate(varids)
    allocate(varids(nVars))

    do iV = 1, nVars
      varname = trim(RegisteredTypes(iTypeIdx)%vars(iV)%shortName)
      call sanitize_nc_name(varname)
      ierr = nfmpi_def_var(ncid, trim(varname), NF_DOUBLE, ndimid, &
                           dimid_vars(1:ndimid), varids(iV))
      if (ierr /= NF_NOERR) &
        write(*, *) "netcdf_open: def_var '", trim(varname), "' failed: ", &
                    nfmpi_strerror(ierr)
      if (len_trim(RegisteredTypes(iTypeIdx)%vars(iV)%units) > 0) then
        attlen = int(len_trim(RegisteredTypes(iTypeIdx)%vars(iV)%units), MPI_OFFSET_KIND)
        ierr = nfmpi_put_att_text(ncid, varids(iV), "units", attlen, &
                                  trim(RegisteredTypes(iTypeIdx)%vars(iV)%units))
      endif
    enddo

    ! --- End define mode (collective) ---
    ierr = nfmpi_enddef(ncid)
    if (ierr /= NF_NOERR) &
      write(*, *) "ModOutputNetCDF: nfmpi_enddef failed: ", nfmpi_strerror(ierr)

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
  ! ==================================================================
  subroutine netcdf_close_file()
#ifdef HavePNetCDF
    integer :: ierr

    if (ncid == -1) return
    ierr = nfmpi_close(ncid)
    if (ierr /= NF_NOERR) &
      write(*, *) "ModOutputNetCDF: nfmpi_close failed: ", nfmpi_strerror(ierr)
    ncid = -1
    if (allocated(varids)) deallocate(varids)
#endif
  end subroutine netcdf_close_file

  ! ==================================================================
  ! Write one block's interior data to the open PnetCDF file.
  !
  ! buffer(nV, nX, nY, nZ) — full buffer for this block
  ! iBLK                   — 1-based global block index
  !
  ! Dispatches on nc_nDims cached from open_file:
  !   3D: strip 2 ghost cells each side; write to (lon×lat×alt) var
  !   2D: no ghost cells (nZ=1); write to (lon×lat) var
  ! ==================================================================
  subroutine netcdf_write_block(iUnit, buffer, nV, nX, nY, nZ, iBLK)
    integer, intent(in) :: iUnit, nV, nX, nY, nZ, iBLK
    real, intent(in) :: buffer(nV, nX, nY, nZ)

#ifdef HavePNetCDF
    integer :: iBlockLon_0, iBlockLat_0
    integer(kind=MPI_OFFSET_KIND) :: start(3), cnt(3)
    integer :: nLi, nLai, nAi, lo_x, hi_x, lo_y, hi_y, lo_z, hi_z
    integer :: ndimid, iV, ierr
    real, allocatable :: slice(:, :, :)

    if (ncid == -1) return

    ! 0-based block spatial indices
    iBlockLon_0 = mod(iBLK - 1, nc_nBlocksLon)
    iBlockLat_0 = (iBLK - 1) / nc_nBlocksLon

    if (nc_nDims == 3) then
      ! 3D type: buffer has 2 ghost cells each side
      nLi  = nX - 4;  lo_x = 3;  hi_x = nX - 2
      nLai = nY - 4;  lo_y = 3;  hi_y = nY - 2
      nAi  = nZ - 4;  lo_z = 3;  hi_z = nZ - 2
      ndimid = 3
    else
      ! 2D type: buffer has no ghost cells (nX=nLons, nY=nLats, nZ=1)
      nLi  = nX;  lo_x = 1;  hi_x = nX
      nLai = nY;  lo_y = 1;  hi_y = nY
      nAi  = nZ;  lo_z = 1;  hi_z = nZ
      ndimid = 2
    endif

    start(1) = int(iBlockLon_0 * nLi,  MPI_OFFSET_KIND) + 1_MPI_OFFSET_KIND
    start(2) = int(iBlockLat_0 * nLai, MPI_OFFSET_KIND) + 1_MPI_OFFSET_KIND
    start(3) = 1_MPI_OFFSET_KIND

    cnt(1) = int(nLi,  MPI_OFFSET_KIND)
    cnt(2) = int(nLai, MPI_OFFSET_KIND)
    cnt(3) = int(nAi,  MPI_OFFSET_KIND)

    allocate(slice(nLi, nLai, nAi))

    do iV = 1, nV
      slice = buffer(iV, lo_x:hi_x, lo_y:hi_y, lo_z:hi_z)
      ierr = nfmpi_put_vara_double_all(ncid, varids(iV), &
                                       start(1:ndimid), cnt(1:ndimid), slice)
      if (ierr /= NF_NOERR) &
        write(*, *) "netcdf_write_block: put_vara_double_all failed for var ", iV, &
                    ": ", nfmpi_strerror(ierr)
    enddo

    deallocate(slice)
#endif
  end subroutine netcdf_write_block

  ! netcdf_write_header is intentionally absent.
  ! ActiveBackend%writes_header = .false. for the netcdf backend, so
  ! output_common.f90 never calls write_header and no .header file is created.
  ! All metadata (time, variable names, units) is embedded in the .nc file.

#ifdef HavePNetCDF
  ! ==================================================================
  ! Replace characters invalid in NetCDF variable names with '_'
  ! ==================================================================
  subroutine sanitize_nc_name(name)
    character(len=*), intent(inout) :: name
    integer :: i
    do i = 1, len_trim(name)
      select case (name(i:i))
      case (' ', '/', '(', ')', '[', ']', '#', '%', '!', '+', '-', '*', '^', '.', &
            ',', '{', '}', '\', '<', '>', '=', '&', '|', '~', '`', '"', "'", '?', '@', '$', ';', ':')
        name(i:i) = '_'
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
