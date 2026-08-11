! Copyright 2026, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

module ModCO2Fomichev

  use ModGITM
  use ModSources, only: CO2Cooling
  use ModInputs
  use ModMpi
  use ModSizeGitm
  
  implicit None

  integer :: nFomichevLines
  character(len=iCharLen_) :: cFomichevText(nInputMaxLines) = ''

  logical :: isFirst
  
  ! 43 = [2 - 12.50] in 0.25 bins
  integer, parameter :: nX0s = 43
  ! 73 = [2 - 20.00] in 0.25 bins
  integer, parameter :: nX0sWhole = 83
  ! 44 = next cell after region 1
  integer, parameter :: iR2Start = 43
  ! 60 = 44 + (16.5 - 12.5) * 4
  integer, parameter :: iR2End = 60
  
  integer, parameter :: nColsL = 10
  integer, parameter :: nLowerTables = 4
  real :: x0(nX0s), dx0
  real :: aTables(nLowerTables, nColsL, nX0s)
  real :: bTables(nLowerTables, nColsL, nX0s)
  real :: aTable(nColsL, nX0s)
  real :: bTable(nColsL, nX0s)
  integer, parameter :: nRows10 = 52
  real :: table10(2, nRows10)
  integer, parameter :: nRows11 = 24
  real :: table11(3, nRows11)
  real :: co2values(nLowerTables)
  data co2values/150.0, 360.0, 540.0, 720.0/

  real :: x0Whole(nX0sWhole)
  real :: pressureWhole(nX0sWhole)
  real :: tempWhole(nX0sWhole), tempReference(nX0sWhole)
  real :: ndenWhole(nX0sWhole), ndenReference(nX0sWhole)
  real :: co2Whole(nX0sWhole), co2denWhole(nX0sWhole)
  real :: n2Whole(nX0sWhole)
  real :: o2Whole(nX0sWhole)
  real :: oWhole(nX0sWhole)
  real :: psiWhole(nX0sWhole)
  real :: capPsi165
  ! mean mass
  real :: muWhole(nX0sWhole)
  ! altitude integral of CO2 in cm2
  real :: loguWhole(nX0sWhole)
  real :: uWhole(nX0sWhole)
  real :: epsilonWhole(nX0sWhole)
  real :: dzWhole(nX0sWhole)
  real :: smallDWhole(nX0sWhole)
  real :: capDWhole(nX0sWhole)
  real :: lambdaWhole(nX0sWhole)
  
  ! These two are real:
  real :: n2Base = 0.7808
  real :: o2Base = 0.2095
  ! This is just a minimal value:
  real :: oBase = 0.0005
  
  ! GITM values:
  real :: gitm_x0(nLons, nLats, -1:nAlts+2)
  real :: gitm_temp(-1:nLons+2, -1:nLats+2, -1:nAlts+2)
  real :: gitm_co2_int(-1:nLons+2, -1:nLats+2, -1:nAlts+2)
  integer :: gitm_iX(nLons, nLats, -1:nAlts+2)
  real :: gitm_t(-1:nAlts+2)

contains

  ! ----------------------------------------------------------------
  ! Read in the whole file and store it in a string
  ! ----------------------------------------------------------------
  
  subroutine read_and_store_file(cFile)

    character(len=iCharLen_) :: line  
    character(len=*), intent(in) :: cFile
    integer :: iError
    logical :: IsThere

    if (iProc == 0) then

       call report("Reading Fomichev File", 0)

       nFomichevLines = 1

       inquire(file=cFile, EXIST=IsThere)
       if (.not. IsThere) &
            call stop_gitm(trim(cFile)//" cannot be found by ModCO2Fomichev")

       open(iInputUnit_, file=cFile, status="old")

       iError = 0
       do while (iError == 0)

          read(iInputUnit_, '(a)', iostat=iError) line

          if (nFomichevLines > nInputMaxLines) &
               call stop_gitm("Too many lines of input in read_inputs")

          cFomichevText(nFomichevLines) = line
          nFomichevLines = nFomichevLines + 1

       enddo

       close(iInputUnit_)

       if (nFomichevLines == 0) &
            call stop_gitm("No lines of input read by ModCO2Fomichev")

    endif

    ! Broadcast the number of lines and the text itself to all processors
    call MPI_Bcast(nFomichevLines, 1, MPI_Integer, 0, iCommGITM, ierror)

    if (iError > 0) &
         call stop_gitm("nFomichevLines could not be broadcast ModCO2Fomichev")

    call MPI_Bcast(cFomichevText, len(cFomichevText(1))*nFomichevLines, &
         MPI_Character, 0, iCommGITM, ierror)
    
    if (iError > 0) &
         call stop_gitm("cFomichevText could not be broadcast by read_mpi")

  end subroutine read_and_store_file

  ! ----------------------------------------------------------------
  ! Read in a generic table
  ! ----------------------------------------------------------------
  
  subroutine parse_general_table(nRows, nCols, iStartLine, table)

    integer, intent(in) :: nRows, nCols, iStartLine
    real, intent(out) :: table(nCols, nRows)
    real :: oneRow(nCols)
    integer :: iRow, iCol, iLine, iError

    do iRow = 1, nRows
       iLine = iStartLine + iRow - 1
       read(cFomichevText(iLine), *, iostat=iError) (oneRow(iCol), iCol=1,nCols)
       table(:, iRow) = oneRow
    enddo
    
  end subroutine parse_general_table

  ! ----------------------------------------------------------------
  ! find table
  ! ----------------------------------------------------------------

  subroutine find_table(cTableName, iRow)

    character(len=*), intent(in) :: cTableName
    integer, intent(out) :: iRow
    logical :: isDone
    integer :: iLine
    
    iRow = -1
    isDone = .false.
    iLine = 1
    do while (.not. isDone)
       if (index(cFomichevText(iLine), cTableName) > 0) then
          isDone = .true.
          iRow = iLine
       endif
       iLine = iLine + 1
       if (iLine >= nFomichevLines) isDone = .true.
    enddo
    
  end subroutine find_table
  
  ! ----------------------------------------------------------------
  ! Read tables
  ! ----------------------------------------------------------------

  subroutine read_all_tables()

    integer :: iRow, iX
    real :: lowerTable(nColsL, nX0s)

    call find_table('table 2', iRow)
    call parse_general_table(nX0s, nColsL, iRow + 2, lowerTable)
    aTables(1, :, :) = lowerTable
    call find_table('table 3', iRow)
    call parse_general_table(nX0s, nColsL, iRow + 2, lowerTable)
    bTables(1, :, :) = lowerTable
    call find_table('table 4', iRow)
    call parse_general_table(nX0s, nColsL, iRow + 2, lowerTable)
    aTables(2, :, :) = lowerTable
    call find_table('table 5', iRow)
    call parse_general_table(nX0s, nColsL, iRow + 2, lowerTable)
    bTables(2, :, :) = lowerTable
    call find_table('table 6', iRow)
    call parse_general_table(nX0s, nColsL, iRow + 2, lowerTable)
    aTables(3, :, :) = lowerTable
    call find_table('table 7', iRow)
    call parse_general_table(nX0s, nColsL, iRow + 2, lowerTable)
    bTables(3, :, :) = lowerTable
    call find_table('table 8', iRow)
    call parse_general_table(nX0s, nColsL, iRow + 2, lowerTable)
    aTables(4, :, :) = lowerTable
    call find_table('table 9', iRow)
    call parse_general_table(nX0s, nColsL, iRow + 2, lowerTable)
    bTables(4, :, :) = lowerTable

    call find_table('table 10', iRow)
    call parse_general_table(nRows10, 2, iRow + 2, table10)
    call find_table('table 11', iRow)
    call parse_general_table(nRows11, 3, iRow + 2, table11)

  end subroutine read_all_tables
  
  ! ----------------------------------------------------------------
  ! Calculate a and b table
  ! ----------------------------------------------------------------
  subroutine calculate_a_and_b_tables(co2value)
    real :: co2value
    integer :: iX
    real :: rX
    
    if (co2value < co2values(1)) then
       iX = 2
       rX = 1.0
    else
       if (co2value > co2values(nLowerTables)) then
          iX = nLowerTables
          rX = 0.0
       else
          if (co2value <= co2values(2)) iX = 2
          if (co2value <= co2values(3)) iX = 3
          if (co2value <= co2values(4)) iX = 4
          rX = (co2values(iX) - co2Value) / (co2values(iX) - co2Values(iX-1))
       endif
    endif

    aTable = rX * aTables(iX-1, :, :) + (1.0 - rX) * aTables(iX, :, :)
    bTable = rX * bTables(iX-1, :, :) + (1.0 - rX) * bTables(iX, :, :)
    
  end subroutine calculate_a_and_b_tables

  ! ----------------------------------------------------------------
  ! Initialize X0s (both small and large)
  ! ----------------------------------------------------------------
  subroutine initialize_x0s
    integer :: iRow, iX
    dx0 = 0.25
    do iX = 1, nX0s
       x0(iX) = aTables(1, 1, iX) 
    enddo
    x0Whole(1) = x0(1)
    do iX = 2, nX0sWhole
       x0Whole(iX) = x0Whole(iX-1) + dx0
    enddo
  end subroutine initialize_x0s

  ! ----------------------------------------------------------------
  ! Set Temperature
  !   - This is based on Figure 4
  ! ----------------------------------------------------------------
  subroutine set_reference_temperature

    integer, parameter :: nControls = 9
    real :: controlTemps(nControls)
    real :: controlX0s(nControls)
    integer :: iUpper, iX
    real :: xR
    
    ! SAW:
    ! 0.0 - 270.0
    ! 1.0 - 210
    ! 3.5 - 200
    ! 4.0 - 185
    ! 8.0 - 260
    ! 11.5 - 220
    ! 12.75 - 220
    ! 14.5 = 200
    ! 17.0 - 260
    data controlTemps/270,210,200,185,260,220,220,200,260/
    data controlX0s/0.0,1.0,3.5,4.0,8.0,11.5,12.75,14.5,17.0/
    iUpper = 1
    do iX = 1, nX0sWhole
       if (x0Whole(iX) > controlX0s(nControls)) then
          tempReference(iX) = controlTemps(nControls)
       else
          do while (x0Whole(iX) > controlX0s(iUpper))
             iUpper = iUpper + 1
          enddo
          xR = (controlX0s(iUpper) - x0Whole(iX)) / &
               (controlX0s(iUpper) - controlX0s(iUpper-1))
          tempReference(iX) = (1.0 - xR) * controlTemps(iUpper) + &
               xR * controlTemps(iUpper - 1)
       endif
    enddo

    ! 1e5 vs 1e3 for pascal vs mbar
    pressureWhole = 1.0e5 / exp(x0Whole)
    ndenReference = pressureWhole / (tempReference * Boltzmanns_Constant)
    
  end subroutine set_reference_temperature
  
  ! ----------------------------------------------------------------
  !
  ! ----------------------------------------------------------------

  subroutine calc_gitm_x0

    ! Equation 1: (pressure is supposed to be in mbar,
    !              but gitm is in pascals -> 1e3 goes to 1e5)    
    gitm_x0 = log(1e5 / pressure(1:nLons, 1:nLats, :, 1))
    gitm_iX = ceiling((gitm_x0 - x0(1)) / dx0) + 1
   
  end subroutine calc_gitm_x0
  
  ! ----------------------------------------------------------------
  !
  ! ----------------------------------------------------------------

  subroutine calc_gitm_co2integral

    integer :: iAlt

    gitm_co2_int(:, :, nAlts + 2) = &
         NDensityS(:, :, nAlts + 2, iCO2_, 1) * dAlt_GB(:, :, nAlts + 2, 1)
    do iAlt = nAlts+1, 0, -1
       gitm_co2_int(:, :, iAlt) = &
            gitm_co2_int(:, :, iAlt + 1) + &
            (NDensityS(:, :, iAlt-1, iCO2_, 1) + &
            NDensityS(:, :, iAlt, iCO2_, 1))/2 * dAlt_GB(:, :, iAlt, 1)
    enddo
    iAlt = -1
    gitm_co2_int(:, :, iAlt) = &
         gitm_co2_int(:, :, iAlt + 1) + &
         NDensityS(:, :, iAlt, iCO2_, 1) * dAlt_GB(:, :, iAlt, 1)

    ! Convert to /cm2:
    gitm_co2_int = gitm_co2_int / 1e4

  end subroutine calc_gitm_co2integral
  
  ! ----------------------------------------------------------------
  !
  ! ----------------------------------------------------------------

  subroutine complete_co2integral
    
    integer :: j

    do j = nX0sWhole-1, 1, -1
       if (uWhole(j) == 0) then
          ! integrate, /1e4 is to convert to /cm2
          uWhole(j) = uWhole(j + 1) + co2denWhole(j) * dzWhole(j) / 1e4
       endif
    enddo
    
  end subroutine complete_co2integral
  
  ! ----------------------------------------------------------------
  !
  ! ----------------------------------------------------------------

  subroutine push_gitm_with_r(iLon, iLat, iAlt, iAltX)
    
    integer, intent(in) :: iLon, iLat, iAlt, iAltX
    real :: r

    integer :: iAltG

    iAltG = iAlt
    r = -1
    do while (r < 0)
       r = (gitm_x0(iLon, iLat, iAltG + 1) - x0Whole(iAltX)) / &
            (gitm_x0(iLon, iLat, iAltG + 1) - gitm_x0(iLon, iLat, iAltG))
       if (r < 0) iAltG = iAltG + 1
    enddo

    tempWhole(iAltX) = &
         r * gitm_temp(iLon, iLat, iAlt) + &
         (1 - r) * gitm_temp(iLon, iLat, iAlt+1)
    ndenWhole(iAltX) = &
         r * nDensity(iLon, iLat, iAlt, 1) + &
         (1 - r) * nDensity(iLon, iLat, iAlt+1, 1)
    n2Whole(iAltX) = ( &
         r * nDensityS(iLon, iLat, iAlt, iN2_, 1) + &
         (1 - r) * nDensityS(iLon, iLat, iAlt+1, iN2_, 1)) / ndenWhole(iAltX)
    o2Whole(iAltX) = ( &
         r * nDensityS(iLon, iLat, iAlt, iO2_, 1) + &
         (1 - r) * nDensityS(iLon, iLat, iAlt+1, iO2_, 1)) / ndenWhole(iAltX)
    co2Whole(iAltX) = ( &
         r * nDensityS(iLon, iLat, iAlt, iCO2_, 1) + &
         (1 - r) * nDensityS(iLon, iLat, iAlt+1, iCO2_, 1)) / ndenWhole(iAltX)
    oWhole(iAltX) = ( &
         r * nDensityS(iLon, iLat, iAlt, iO_3P_, 1) + &
         (1 - r) * nDensityS(iLon, iLat, iAlt+1, iO_3P_, 1) + &
         r * nDensityS(iLon, iLat, iAlt, iO_1D_, 1) + &
         (1 - r) * nDensityS(iLon, iLat, iAlt+1, iO_1D_, 1)) / ndenWhole(iAltX)
    uWhole(iAltX) = ( &
         r * gitm_co2_int(iLon, iLat, iAlt) + &
         (1 - r) * gitm_co2_int(iLon, iLat, iAlt+1))
    
  end subroutine push_gitm_with_r
  
  ! ----------------------------------------------------------------
  !
  ! ----------------------------------------------------------------

  subroutine move_to_x0_grid(iLon, iLat)

    integer, intent(in) :: iLon, iLat
    integer :: iAltG, iAltX
    real :: r
    
    tempWhole = tempReference
    ndenWhole = ndenReference
    co2Whole = CO2ppm / 1e6
    n2Whole = n2Base
    o2Whole = o2Base
    oWhole = oBase
    uWhole = 0.0
    
    do iAltG = 1, nAlts-1
       iAltX = gitm_iX(iLon, iLat, iAltG)
       if (iAltX <= nX0sWhole) then
          call push_gitm_with_r(iLon, iLat, iAltG, iAltX)
          do while (iAltX < gitm_iX(iLon, iLat, iAltG + 1))
             iAltX = iAltX + 1
             if (iAltX <= nX0sWhole) then
                call push_gitm_with_r(iLon, iLat, iAltG, iAltX)
             endif
          enddo
       endif

    enddo

    ! co2 number density in /m3 (for calculation of u)
    co2denWhole = co2Whole * ndenWhole
    
    ! This is the mean weight in g/mol:
    muWhole = &
         28.014 * n2Whole + &
         32.00 * o2Whole + &
         16.00 * oWhole + &
         44.01 * co2Whole

    ! Convert to /cm3 instead of /m3:
    ndenWhole = ndenWhole / 1e6
    
  end subroutine move_to_x0_grid

  ! ----------------------------------------------------------------
  !
  ! ----------------------------------------------------------------

  subroutine move_to_gitm_grid(iLon, iLat)

    integer, intent(in) :: iLon, iLat
    integer :: iAltG, iAltX
    real :: r

    CO2Cooling(iLon, iLat, :) = 0.0
    do iAltG = 1, nAlts-1
       iAltX = gitm_iX(iLon, iLat, iAltG)
       if (iAltX <= nX0sWhole) then
          r = (gitm_x0(iLon, iLat, iAltG) - x0Whole(iAltX-1)) / (x0Whole(iAltX) - x0Whole(iAltX-1))
          CO2Cooling(iLon, iLat, iAltG) = &
               (1.0 - r) * epsilonWhole(iAltX - 1) + r * epsilonWhole(iAltX) 
       endif
    enddo

    ! Epsilon was in ergs / g / s and we need it in W/m
    ! GITM expects cooling to be positive, since it subtracts things
    CO2Cooling(iLon, iLat, :) = - CO2Cooling(iLon, iLat, :) * 1e-7 * (rho(iLon, iLat, 1:nAlts, 1)*1000.0)
    
  end subroutine move_to_gitm_grid

  ! ----------------------------------------------------------------
  !
  ! ----------------------------------------------------------------

  subroutine calc_dz

    integer :: j
    real :: radius, alt, g, h
    real, parameter :: re = 6371.0

    ! Guess initial altitude (x0 = 2, from figure 1):
    alt = 15.0
    
    do j = 1, nX0sWhole - 1
       g = 9.8 * (re / (re+alt))**2
       h = Boltzmanns_Constant * tempWhole(j) / (muWhole(j) * amu * g)
       dzWhole(j) = -log(pressureWhole(j+1)/pressureWhole(j)) * h
       alt = alt + dzWhole(j)/1000.0
    enddo
    
  end subroutine calc_dz
  
  ! ----------------------------------------------------------------
  ! This is below and including x0 = 12.5, using tables 2-9
  ! ----------------------------------------------------------------

  subroutine calc_region1_cooling

    integer :: iX, iXr, j, jr
    real :: dx(-5:3)
    ! table 1:
    data dx/-6.25, -3, -1.75, -0.75, -0.25, 0, 0.25, 0.75, 1.50/    
    
    do iX = 1, nX0s
       epsilonWhole(iX) = 0.0
       do j = -5, 3
          jr = j + 6
          iXr = iX + dx(j)/dx0
          if (iXr < 1) iXr = 1
          epsilonWhole(iX) = epsilonWhole(iX) + &
               (aTable(jr, iX) + bTable(jr, iX) * psiWhole(iX)) * psiWhole(iXr)
       enddo
    enddo
    
  end subroutine calc_region1_cooling

  ! ----------------------------------------------------------------
  ! Use Table 11 to get alpha from log(u) and xj
  ! ----------------------------------------------------------------

  subroutine get_alpha_from_table11(logu, xj, alpha)
    real, intent(in) :: logu, xj
    real, intent(out) :: alpha
    integer :: iRow, iOff
    real :: r

    iRow = (xj - 12.5) / 0.25 * 4 + 1
    iOff = 0
    if ((iRow < 1) .or. (iRow > nRows11)) then
       alpha = 1
    else
       do while (logu > table11(2, iRow + iOff) .and. iOff < 3)
          iOff = iOff + 1
       enddo
       ! this means it is lower than the bottom, so it needs to be extrapolated:
       if (iOff == 0) iOff = 1
       iRow = iRow + iOff
       r = (table11(2, iRow) - logu) / (table11(2, iRow) - table11(2, iRow - 1))
       alpha = (1 - r) * table11(3, iRow) + r * table11(3, iRow - 1)
    endif
    if (alpha < 1) alpha = 1
  end subroutine get_alpha_from_table11
  

  ! ----------------------------------------------------------------
  ! Use Table 10 to convert log(u) to Lu
  ! ----------------------------------------------------------------

  subroutine get_Lu_from_table10(logu, logLu)
    real, intent(in) :: logu
    real, intent(out) :: logLu
    integer :: iRow
    real :: r
    if (logu < table10(1,1)) then
       logLu = table10(2,1)
    else
       if (logu > table10(1,nRows10)) then
          logLu = table10(2,nRows10)
       else
          iRow = 2
          do while (logu > table10(1,iRow))
             iRow = iRow + 1
          enddo
          r = (table10(1, iRow) - logu) / (table10(1, iRow) - table10(1, iRow - 1))
          logLu = (1 - r) * table10(2, iRow) + r * table10(2, iRow - 1)
       endif
    endif
  end subroutine get_Lu_from_table10
  
  ! ----------------------------------------------------------------
  !
  ! ----------------------------------------------------------------
  
  subroutine get_Lu_and_alpha(logu, Lu, alpha)

    real, intent(in) :: logu(nX0sWhole)
    real, intent(out) :: Lu(iR2Start : iR2End)
    real, intent(out) :: alpha(iR2Start : iR2End)

    integer :: j
    real :: logLu, alpha1
    
    do j = iR2Start, iR2End
       call get_Lu_from_table10(logu(j), logLu)
       Lu(j) = exp(-logLu)
       call get_alpha_from_table11(logu(j), x0Whole(j), alpha1)
       alpha(j) = alpha1
    enddo
        
  end subroutine get_Lu_and_alpha

  ! ----------------------------------------------------------------
  ! Region 2 = 12.75 - 16.5, using tables 10 and 11
  ! ----------------------------------------------------------------

  subroutine calc_region2_cooling

    integer :: j

    real :: epsilontilde(iR2Start : iR2End)
    real :: Lu(iR2Start : iR2End)
    real :: alpha(iR2Start : iR2End)
    real :: n2part, o2part, opart
    real :: kn2, ko2, ko
    real :: Dj, Dj_1, left, right
    
    epsilontilde = 0.0

    ! Calculate Lambda from equation 8:
    ko = 3e-12
    do j = iR2Start, nX0sWhole
       ! break apart equation 8 into parts:
       kn2 = 5.5e-17 * sqrt(tempWhole(j)) + &
            6.7e-10 * exp(-83.8 * (tempWhole(j) ** (-0.33333)))
       ko2 = 1e-15 * exp( &
            23.37 - &
            230.9 * (tempWhole(j) ** (-0.33333)) + &
            564.0 * (tempWhole(j) ** (-0.66667)))
       n2part = n2Whole(j) * kn2
       o2part = o2Whole(j) * ko2
       opart = oWhole(j) * ko
       lambdaWhole(j) = 1.5988 / (1.5988 + ndenWhole(j) * (n2part + o2part + opart))
    enddo

    ! This is needed for equation 12:
    loguWhole = log(uWhole)
    call get_Lu_and_alpha(loguWhole, Lu, alpha)

    ! Equation 12 (alpha is 1.0 when Xj is 14.0 and above:
    do j = iR2Start, iR2End
       smallDWhole(j) = alpha(j) * Lu(j)
    enddo

    ! Equation 10:
    epsilontilde(iR2Start) = 1.10036e-10 * epsilonWhole(iR2Start) / &
         (co2Whole(iR2Start) * (1.0 - lambdaWhole(iR2Start)))

    ! This starts at 12.75 (instead of 12.5):
    do j = iR2Start + 1, iR2End
       ! Equation 11:
       Dj = 0.25 * (smallDWhole(j - 1) + 3.0 * smallDWhole(j))
       Dj_1 = 0.25 * (3.0 * smallDWhole(j - 1) + smallDWhole(j))
       ! Equation 9:
       left = (1.0 - lambdaWhole(j) * (1.0 - Dj))
       right = (1.0 - lambdaWhole(j - 1) * (1.0 - Dj_1))
       epsilontilde(j) = &
            (right * epsilontilde(j - 1) + &
            Dj_1 * psiWhole(j - 1) - &
            Dj * psiWhole(j)) / left
       ! Equation 7:
       epsilonWhole(j) =  2.63187e11 * co2Whole(j) * (1.0 - lambdaWhole(j)) / muWhole(j) * epsilontilde(j)
    enddo
    capPsi165 = epsilontilde(iR2End) + psiWhole(iR2End) 
    
  end subroutine calc_region2_cooling

  ! ----------------------------------------------------------------
  ! Region 3 = 16.75 +
  ! ----------------------------------------------------------------

  subroutine calc_region3_cooling

    integer :: j

    ! Equation 13
    do j = iR2End + 1, nX0sWhole
       epsilonWhole(j) = &
            2.63187e11 * co2Whole(j) * (1.0 - lambdaWhole(j)) / muWhole(j) * (capPsi165 - psiWhole(j))
    enddo

  end subroutine calc_region3_cooling
    
  ! ----------------------------------------------------------------
  ! This are the main subroutines
  ! ----------------------------------------------------------------
  
  ! ----------------------------------------------------------------
  ! Initialize - this needs to be called once and only once
  ! ----------------------------------------------------------------

  subroutine initialize_fomichev_cooling

    call read_and_store_file(cFomichevFile)
    call read_all_tables
    call initialize_x0s
    call calculate_a_and_b_tables(CO2ppm)
    call set_reference_temperature
    isFirst = .true.
    
  end subroutine initialize_fomichev_cooling

  ! ----------------------------------------------------------------
  ! Calculate the CO2 cooling
  !   This is done in three regions:
  !    - Region 1: below (and including x0 = 12.5)
  !    - Region 2: 12.75 - 16.5
  !    - Region 3: above 16.5
  !   Use 1D columns from GITM
  !    - transform pressure to x0 (3D)
  !    - map GITM values onto x0 1D grid (temp, N2, O2, O, CO2, num den)
  !    - calculate cooling on x0 grid
  !    - map back to GITM grid
  ! ----------------------------------------------------------------

  subroutine calc_co2fomichev_cooling

    integer :: iLon, iLat, j
    
    gitm_temp = temperature(:,:,:,1) * tempunit
    
    call calc_gitm_x0
    call calc_gitm_co2integral

    do iLat = 1, nLats
       do iLon = 1, nLons
    
          call move_to_x0_grid(iLon, iLat)

          call calc_dz
          call complete_co2integral
    
          psiWhole = exp(-960.217 / tempWhole)
          call calc_region1_cooling
          call calc_region2_cooling
          call calc_region3_cooling

          call move_to_gitm_grid(iLon, iLat)

       enddo
    enddo

    ! I need to figure out why it is sometime positive.
    where (co2cooling < 0)
       co2cooling = 0.0
    endwhere

  end subroutine calc_co2fomichev_cooling
  
  
end module ModCO2Fomichev
  

