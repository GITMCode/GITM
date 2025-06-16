# Common Inputs {#indices.sec}

This only touches on the most frequently changed input options. For a full
reference of all available inputs, please see [All Inputs](set_inputs.md)

!!!note
    All of the auxiliary input (data) files can be kept in the same directory
    as the GITM executable and the `UAM.in` file. Otherwise, the path should be
    specified relative to the GITM executable.

## IMF and Solar Wind {#imf.sec}

This file controls the high-latitude electric field and aurora when
using models that depend on the solar wind and interplanetary magnetic
field (IMF). It allows GITM to dynamically control these quantities. You
can create either realistic IMF files or hypothetical ones.

A script is provided to automatically download these files, and can be called
using the `srcPython/omni_download_write_swmf.py` python file.

Here is an example file:

    This file was created by Aaron Ridley to do some wicked cool science thing.

    The format is:
     Year MM DD HH Mi SS mS   Bx  By   Bz     Vx   Vy   Vz    N        T

    Year=year
    MM = Month
    DD = Day
    HH = Hour
    Mi = Minute
    SS = Second
    mS = Millisecond
    Bx = IMF Bx GSM Component (nT)
    By = IMF By GSM Component (nT)
    Bz = IMF Bz GSM Component (nT)
    Vx = Solar Wind Vx (km/s)
    Vy = Solar Wind Vy (km/s)
    Vz = Solar Wind Vz (km/s)
    N  = Solar Wind Density (/cm3)
    T  = Solar Wind Temperature (K)

    #DELAY
    900.0

    #START
     2000  3 20  2 53  0  0  0.0 0.0  2.0 -400.0  0.0  0.0  5.0  50000.0
     2000  3 20  2 54  0  0  0.0 0.0  2.0 -400.0  0.0  0.0  5.0  50000.0
     2000  3 20  2 55  0  0  0.0 0.0  2.0 -400.0  0.0  0.0  5.0  50000.0
     2000  3 20  2 56  0  0  0.0 0.0  2.0 -400.0  0.0  0.0  5.0  50000.0
     2000  3 20  2 57  0  0  0.0 0.0  2.0 -400.0  0.0  0.0  5.0  50000.0
     2000  3 20  2 58  0  0  0.0 0.0  2.0 -400.0  0.0  0.0  5.0  50000.0
     2000  3 20  2 59  0  0  0.0 0.0  2.0 -400.0  0.0  0.0  5.0  50000.0
     2000  3 20  3  0  0  0  0.0 0.0 -2.0 -400.0  0.0  0.0  5.0  50000.0
     2000  3 20  3  1  0  0  0.0 0.0 -2.0 -400.0  0.0  0.0  5.0  50000.0
     2000  3 20  3  2  0  0  0.0 0.0 -2.0 -400.0  0.0  0.0  5.0  50000.0
     2000  3 20  3  3  0  0  0.0 0.0 -2.0 -400.0  0.0  0.0  5.0  50000.0
     2000  3 20  3  4  0  0  0.0 0.0 -2.0 -400.0  0.0  0.0  5.0  50000.0

This file is provided to GITM with:

    #MHD_INDICES
    imf_file_name.dat

## Hemispheric Power {#hp.sec}

The hemispheric power files describe the dynamic variation of the
auroral power going into each hemisphere. Models such as FRE[^FRE] use
the Hemispheric Power to determine which level of the model it should
use. The Hemispheric Power is converted to a Hemispheric Power Index
using the formula:

```math
\begin{equation}
HPI = 2.09 \log(HP)^{1.0475}
\end{equation}
```

[^FRE]: Fuller-Rowell, T. J., and D. S. Evans (1987), Height-integrated Pedersen
    and Hall conductivity patterns inferred from the TIROS-NOAA satellite data,
    J. Geophys. Res., 92(A7), 7606–7618, 
    doi:[10.1029/JA092iA07p07606](https://doi.org/10.1029/JA092iA07p07606).

The National Oceanic and Atmospheric Administration (NOAA) provides
these hemispheric power files for public use online at
<http://www.swpc.noaa.gov/ftpmenu/lists/hpi.html>. There are two types
of formats used for hemispheric power files (due to a change in the NOAA
output format in 2007). Both file formats can be used by GITM, and are
shown in the examples below.

!!! attention 
    GITM can read a NOAA HPI file, however the recommended way to provide
    hemispheric power is to have GITM derive HPI from the AE-index.
    
    See [SME Indices](#sme-indices) for more information. 

Example file 1 for data prior to 2007:

    # Prepared by the U.S. Dept. of Commerce, NOAA, Space Environment Center.
    # Please send comments and suggestions to sec@sec.noaa.gov 
    # 
    # Source: NOAA POES (Whatever is aloft)
    # Units: gigawatts

    # Format:

    # The first line of data contains the four-digit year of the data.
    # Each following line is formatted as in this example:

    # NOAA-12(S)  10031     9.0  4    .914

    # Please note that if the first line of data in the file has a
    # day-of-year of 365 (or 366) and a HHMM of greater than 2300, 
    # that polar pass started at the end of the previous year and
    # ended on day-of-year 001 of the current year.

    # A7    NOAA POES Satellite number
    # A3    (S) or (N) - hemisphere
    # I3    Day of year
    # I4    UT hour and minute
    # F8.1  Estimated Hemispheric Power in gigawatts
    # I3    Hemispheric Power Index (activity level)
    # F8.3  Normalizing factor

    2000
    NOAA-15(N)  10023    35.5  7    1.085
    NOAA-14(S)  10044    25.3  7     .843
    NOAA-15(S)  10114    29.0  7     .676
    NOAA-14(N)  10135   108.7 10    1.682
    NOAA-15(N)  10204    36.4  7    1.311
    .
    .
    .

Example file 2 for data in and after 2007:

    :Data_list: power_2010.txt
    :Created: Sun Jan  2 10:12:58 UTC 2011


    # Prepared by the U.S. Dept. of Commerce, NOAA, Space Environment Center.
    # Please send comments and suggestions to sec@sec.noaa.gov 
    # 
    # Source: NOAA POES (Whatever is aloft)
    # Units: gigawatts

    # Format:

    # Each line is formatted as in this example:

    # 2006-09-05 00:54:25 NOAA-16 (S)  7  29.67   0.82

    # A19   Date and UT at the center of the polar pass as YYYY-MM-DD hh:mm:ss
    # 1X    (Space)
    # A7    NOAA POES Satellite number
    # 1X    (Space)
    # A3    (S) or (N) - hemisphere
    # I3    Hemispheric Power Index (activity level)
    # F7.2  Estimated Hemispheric Power in gigawatts
    # F7.2  Normalizing factor

    2010-01-01 00:14:37 NOAA-17 (N)  1   1.45   1.16
    2010-01-01 00:44:33 NOAA-19 (N)  1   1.45   1.17
    .
    .
    .

The file type is automatically inferred. To provide an HPI file, use:

    #NOAAHPI_INDICES
    hemi-power-file.txt

## SME Indices

*[SME]: SuperMag auroral Electrojet

To use models such as FTA[^1] to drive the aurora, GUITM must be provided
Auroral Electrojet (AE) indices. Normally this is from SuperMag, but any source
may be used. 

[^1]: Wu, C., Ridley, A. J., DeJong, A. D., & Paxton, L. J. (2021). FTA: A Feature Tracking Empirical Model Of Auroral Precipitation. Space Weather, 19, e2020SW002629. <https://doi.org/10.1029/2020SW002629>.

These files are normally of the format:

    File created by python code using SuperMAGGetIndices

    ============================================================
    <year>  <month>  <day>  <hour>  <min>  <sec>  <SME (nT)>  <SML (nT)>  <SMU (nT)>
    2002  12  21  00  00  00   616.74  -354.47   262.26
    2002  12  21  00  01  00   623.75  -354.72   269.03
    2002  12  21  00  02  00   617.18  -349.28   267.90
    2002  12  21  00  03  00   633.56  -350.01   283.55
    2002  12  21  00  04  00   664.55  -357.88   306.68

A python routine to download these files over a given date range can be found
in `srcPython/supermag_download_ae.py`.


The corresponding section in `UAM.in` is read as:

    #SME_INDICES
    ae_file-name.dat        ae file name
    none                    onset file
    T                       use AE for HP

The lines following the AE file are for the AL-onset file. This can be set to
`none` if you do not have one. The next line tells GITM whether to derive HP
from AE or to use a NOAA HPI file (if one is required). Even if AE is not
required, it is recommended to provide one to derive HP as it is often more
representative of geomagnetic conditions than the NOAA HPI.

The formula to calculate hemispheric power from AE is taken from (Wu et al.,
2021)[^wuetal] and is given as:

```math
\begin{align}
HemisphericPower = 0.102 * AE + 8.953
\end{align}
```

[^wuetal]: Wu, C., Ridley, A. J., DeJong, A. D., & Paxton, L. J. (2021).
FTA: A Feature Tracking Empirical Model Of Auroral Precipitation. 
Space Weather, 19, e2020SW002629. <https://doi.org/10.1029/2020SW002629>


## Solar Irradiance {#solar_irradiance.sec}

To provide GITM with realistic solar irradiance, the solar EUV must be
specified. This can be done through a file containing modeled or
observed solar irradiance data. An example from the FISM model is shown
below.

    #START
        2009       3      20       0       0       0   0.00389548   0.00235693
       0.00127776  0.000907677  0.000652528  0.000372993  0.000250124  0.000194781
      0.000389686  0.000118650   0.00642058   0.00618358  0.000133490  7.67560e-05
      7.80045e-05  0.000145722  5.92577e-05  5.95070e-05  0.000102437  6.48526e-05
      8.94509e-05  0.000101928  5.94333e-05  5.36012e-05  1.51744e-05  1.10265e-05
      1.26937e-05  2.16591e-05  9.57055e-06  1.82608e-05  7.07992e-05  2.55451e-05
      1.12451e-05  6.89255e-05  3.03882e-05  2.33862e-05  2.98026e-05  4.44682e-05
      1.50847e-05  3.00909e-05  8.18379e-05  3.52176e-05  0.000416491  0.000269080
      0.000269080  0.000275734  6.60872e-05  4.46671e-05  0.000220697  0.000512933
      3.85239e-05  9.30928e-05  2.71239e-05  1.23011e-05  1.05722e-05  9.30876e-06
      7.08442e-07  3.54221e-07  1.77110e-07
        2009       3      20       0       1       0   0.00389548   0.00235693
       0.00127776  0.000907677  0.000652528  0.000372993  0.000250124  0.000194781
      0.000389686  0.000118650   0.00642058   0.00618358  0.000133490  7.67560e-05
      7.80045e-05  0.000145722  5.92577e-05  5.95070e-05  0.000102437  6.48526e-05
      8.94509e-05  0.000101928  5.94333e-05  5.36012e-05  1.51744e-05  1.10265e-05
      1.26937e-05  2.16591e-05  9.57055e-06  1.82608e-05  7.07992e-05  2.55451e-05
      1.12451e-05  6.89255e-05  3.03882e-05  2.33862e-05  2.98026e-05  4.44682e-05
      1.50847e-05  3.00909e-05  8.18379e-05  3.52176e-05  0.000416491  0.000269080
      0.000269080  0.000275734  6.60872e-05  4.46671e-05  0.000220697  0.000512933
      3.85239e-05  9.30928e-05  2.71239e-05  1.23011e-05  1.05722e-05  9.30876e-06
      7.08442e-07  3.54221e-07  1.77110e-07
    .
    .
    .

GITM knows to use the provided solar irradiance file through the EUV_DATA input
option specified in the `UAM.in` file:

    #EUV_DATA
    T
    UA/DataIn/FISM/fismflux_daily_[year].dat

These files are included by default, however will not be read unless this option
is added to the `UAM.in` file.

## Satellites {#sat_aux.sec}

GITM can provide output data at a list of times and locations using the
SATELLITE input option. Although this is designed to output data along
a satellite orbit, any list of locations may be used. There is currently
no routine to create a satellite input file, but the format is simple
and may be easily constructed from a satellite ASCII data file using
awk or Python. Here is a sample satellite input file:

    year mm dd hh mm ss msec long lat alt
    #START
    2002 4 16 23 34 25 0 299.16 -2.21 0.00 
    2002 4 16 23 34 25 0 293.63 -1.21 0.00 
    2002 4 16 23 34 25 0 291.28 -0.75 0.00 
    2002 4 16 23 34 25 0 289.83 -0.45 0.00 
    2002 4 16 23 34 25 0 288.79 -0.21 0.00 
    2002 4 16 23 34 25 0 287.98 -0.01 0.00 
    2002 4 16 23 34 25 0 287.32  0.16 0.00 
    2002 4 16 23 34 25 0 286.76  0.31 0.00 
    2002 4 16 23 34 25 0 286.26  0.46 0.00 
    2002 4 16 23 34 25 0 285.81  0.60 0.00 
    2002 4 16 23 34 25 0 285.39  0.74 0.00

Note that the satellite output type is not specified in this sample file. This
is because altitude entry doesn't matter at this time, GITM ignores the altitude
and outputs altitudinal profiles of the atmospheric characteristics at each
geographic location and universal time. Although millisecond accuracy is
provided, GITM should not be output at a resolution smaller than 1 second. The
temporal resolution in the satellite file does not need to match the output
resolution.

The types of outputs are specified in the [`#OUTPUT`](outputs.md) section of
the UAM file.
