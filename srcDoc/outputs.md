# Outputs

#TODO:

- [ ] Fix the links
- [ ] remove idl code
- [ ] probably(?) remove the python stuff too
- [ ] Fix and/or remove images
- [ ] update with other output types:
    - [ ] 3D*, 2D*
    - [ ] satellite outputs
    - [ ] empirical models only!


Now that you have managed to successfully complete a GITM run you've
found yourself with a bunch of output files. All of the GITM output is
in mks units and this data is contained within several files located in
the `UA/data` directory, as was previously discussed in
Chapter [\[quickstart.ch\]](#quickstart.ch){reference-type="ref"
reference="quickstart.ch"}
Section [\[post_process.sec\]](#post_process.sec){reference-type="ref"
reference="post_process.sec"}. You will have found yourself with several
`iriOut_*.dat` files, a `log*.dat` file, and many `.bin` files in
whichever formats you specified in SAVEPLOT (see
Chapter [\[input.ch\]](#input.ch){reference-type="ref"
reference="input.ch"}
Section [\[def_out.sec\]](#def_out.sec){reference-type="ref"
reference="def_out.sec"}). The `iriOut_*.dat` files are required by the
IRI model and not typically used when analyzing the outcome of the GITM
run.

The log file provides useful information about the run, such as whether
a restart was performed, which physical processes were used, and a list
of the universal time, time-step, neutral temperature ranges (T), solar
and geomagnetic indices, and the neutral velocity (VV) ranges for each
iteration. This file can be very useful when sharing runs with other
users, when revisiting an old run, or merely ensuring that GITM
performed as expected. An example log file is provided below:

    ## Inputs from UAM.in
    # Resart= F
    # Eddy coef:   100.000 Eddy P0:     0.020 Eddy P1:     0.003 Eddy Scaling:     1.000
    # Statistical Models Only:  F Apex:  T
    # EUV Data:  TFile: 
    fismflux.dat                                                                                        
    # AMIE: none           
    none                                                                                                
    # Solar Heating:  T Joule Heating:  T Auroral Heating:  T
    # NO Cooling:  T O Cooling:  T
    # Conduction:  T Turbulent Conduction:  T Updated Turbulent Conduction:  T
    # Pressure Grad:  T Ion Drag:  T Neutral Drag:  T
    # Viscosity:  T Coriolis:  T Gravity:  T
    # Ion Chemistry:  T Ion Advection:  T Neutral Chemistry:  T
     
    #START
       iStep yyyy mm dd hh mm ss  ms      dt min(T) max(T)...
       ...mean(T) min(VV) max(VV) mean(VV) F107 F107A By Bz Vx...
       ...HP HPn HPs SubsolarLon SubsolarLat SubsolarVTEC
           2 2011  9 23  0  0  2 297  2.2979  168.75192  1062.87354...
           ...933.09984 -48.19362    524.93645  1.01910 159.3 127.9 -4.6  0.5 406.9...
           ...11.1 14.4  15.5  3.14145  -0.37655  45.73188
           .
           .
           .

The 3DALL output binary files can contain the following atmospheric
quantities:

-   **Altitude:** Altitude from the surface of the planet (m)

-   **Ar:** Argon density (m$`^{-3}`$)

-   **Ar Mixing Ratio:** Argon mixing ratio

-   **CH4 Mixing Ratio:** Methane mixing ratio

-   **Conduction:** Heat conduction

-   **EuvHeating:** EUV Heating rate

-   **H:** Hydrogen density (m$`^{-3}`$)

-   **H!U+!N:** H$`^+`$ density (m$`^{-3}`$)

-   **H2 Mixing Ratio:** Molecular Hydrogen mixing ratio

-   **HCN Mixing Ratio:** Hydrogen Cyanide mixing ratio

-   **He:** Helium density (m$`^{-3}`$)

-   **He!U+!N:** He$`^+`$ density (m$`^{-3}`$)

-   **Heaing Efficiency:** Heating efficiency

-   **Heat Balance Total:** Heat balance total

-   **Latitude:** Geographic latitude (degrees)

-   **Longitude:** Geographic longitude (degrees)

-   **N!D2!N:** N$`_2`$ density (m$`^{-3}`$)

-   **N!D2!U+!N:** N$`_2^+`$ density (m$`^{-3}`$)

-   **N!U+!N:** N$`^+`$ density (m$`^{-3}`$)

-   **N(!U2!ND):** N($`^2`$D) density (m$`^{-3}`$)

-   **N(!U2!NP):** N($`^2`$P) density (m$`^{-3}`$)

-   **N(!U4!NS):** N($`^4`$S) density (m$`^{-3}`$)

-   **N2 Mixing Ratio:** Molecular nitrogen mixing ratio

-   **NO:** Nitrious Oxide density (m$`^{-3}`$)

-   **NO!U+!N:** NO$`^+`$ density (m$`^{-3}`$)

-   **O!D2!N:** O$`_2`$ density (m$`^{-3}`$)

-   **O!D2!U+!N:** O$`_2^+`$ density (m$`^{-3}`$)

-   **O(!U1!ND):** O($`^1`$D) density (m$`^{-3}`$)

-   **O(!U2!ND)!U+!N:** O($`^2`$D) density (m$`^{-3}`$)

-   **O(!U2!NP)!U+!N:** O($`^2`$P) density (m$`^{-3}`$)

-   **O(!U3!NP):** O($`^3`$P) density (m$`^{-3}`$)

-   **O_4SP\_!U+!N:** O($`_4`$SP)$`^+`$ density (m$`^{-3}`$)

-   **RadCooling:** Radiative Cooling rate

-   **Rho:** Neutral density (m$`^{-3}`$)

-   **Temperature:** Neutral temperature (K)

-   **V!Di!N (east):** Ion velocity towards geographic East (m s$`^{-1}`$)

-   **V!Di!N (north):** Ion velocity towards geographic North (m
    s$`^{-1}`$)

-   **V!Di!N (up):** Vertical ion velocity (m s$`^{-1}`$)

-   **V!Dn!N (east):** Neutral velocity towards geographic East (m
    s$`^{-1}`$)

-   **V!Dn!N (north):** Neutral velocity towards geographic North (m
    s$`^{-1}`$)

-   **V!Dn!N (up):** Vertical neutral velocity (m s$`^{-1}`$)

-   **V!Dn!N (up,N!D2!N):** Vertical N$`_2`$ velocity (m s$`^{-1}`$)

-   **V!Dn!N (up,N(!U4!NS)):** Vertical N($`^4`$S) velocity (m s$`^{-1}`$)

-   **V!Dn!N (up,NO):** Vertical NO velocity (m s$`^{-1}`$)

-   **V!Dn!N (up,O!D2!N):** Vertical O$`_2`$ velocity (m s$`^{-1}`$)

-   **V!Dn!N (up,O(!U3!NP)):** Vertical O($`^3`$P) velocity (m s$`^{-1}`$)

-   **e-:** electron density (m$`^{-3}`$)

-   **eTemperature:** electron temperature (K)

-   **iTemperature:** ion temperature (K)

-   **time:** Universal time

There are many routines available to process and analyze the GITM binary
files. The majority of these routines are written in IDL and are
available in the `srcIDL` directory within the GITM model directory.
Currently 50 routines have been saved in this directory and more are
under development. Alternatively, python routines are currently being
developed and these are located in the `srcPython` directory. Please not
that when using the IDL reader the universal time is read in as epoch
seconds from January 1, 1965 00:00 UT, while when using the python
reader, the time is imported as a datetime object.

# IDL {#idl.sec}

Here is an complete list with some description of the IDL processing and
visualization routines currently available. Please feel free to update
this section for other GITM users when you CVS your vetted GITM
processing routines.

### gitm_read_bin

This is a routine to read a GITM bin file into IDL. This is great when
you want to get a handle on the data and experiment with different
visualization methods.

### thermo_plotsat

This is the most commonly used routine to plot the 1D GITM results. It
can also be used to plot satellite files and other 1D simulations. It is
relatively straight forward to use, but experimentation can be help.
This is an actual program, so you have to `.run` it.

### thermo_gui

This is a someone simplistic graphical user interface code for plotting
3D results. The filename has to be entered manually in the upper left.
You then have to press the button for loading the file. Variables appear
on the left side, and you can select which one you want to plot. You
then select which of the available planes you would like to look at
(lat/lon, lat/alt, or lon/alt) or scroll through the options. This
interface allows you to add wind vectors, plot in polar coordinates, and
plot the log of the variable.

### thermo_batch_new

This code will let you look at at 3D files exactly the same way as
thermo_gui, but is all scripted. There are a few features that this has
that thermo_batch doesn't have:

1.  You can use wildcards for the file name, so that a list of files can
    be read. The postscript file names created for each figure will be
    differentiated by appending numbers sequentially so that no figures
    are overwritten.

2.  When plotting a lat/alt plane, you can do a zonal average.

3.  You can do a global average.

### thermo_plotter

All of the above plotting codes will only plot one plot per page. This
code will plot many more than one plot per page. You can plot multiple
variables on the same page, or multiple files with the same variable, or
both.

### Other IDL Routines

Please feel free to provide a description of these routines so that GITM
users do not waste their time rewriting code that already exists.

::: multicols
3

-   **ask**

-   **c_a_to_r**

-   **c_a_to_s**

-   **chopr**

-   **closedevice**

-   **c_r_to_a**

-   **c_s_to_a**

-   **get_position**

-   **makect**

-   **mklower**

-   **mm**

-   **plotct**

-   **plotdumb**

-   **plotmlt**

-   **pos_space**

-   **read_thermosphere_file**

-   **setdevice**

-   **thermo_batch**

-   **thermo_calcforce**

-   **thermo_champ**

-   **thermo_compare**

-   **thermo_compare_time**

-   **thermo_convert_champfiles**

-   **thermo_guvi**

-   **thermo_magequator**

-   **thermo_make_summary**

-   **thermo_mkguvisat**

-   **thermo_mksatsave**

-   **thermo_mksave**

-   **thermo_mktec**

-   **thermo_on2**

-   **thermo_plotdist**

-   **thermo_plotlog**

-   **thermo_plot_new**

-   **thermo_plot**

-   **thermo_plotsat2**

-   **thermo_plotsat_constalt_ON2**

-   **thermo_plotsat_constalt**

-   **thermo_plotvectors**

-   **thermo_readsat**

-   **thermo_sigma**

-   **thermo_superposed**

-   **thermo_tec**

-   **thermo_temp**

-   **tostr**
:::

# Python {#python.sec}

This section provides an almost complete list of the vetted GITM python
routines. These routines require that you use PyBats, a module included
in SpacePy. This is a library developed for space physics applications
by the scientists at Los Alamos and can be downloaded for free at:
`http://spacepy.lanl.gov`

Another library, Basemap, is required for certain plotting routines.
Basemap is a part of the Matplotlib Toolkit and can be installed using
Fink, Macports, or downloaded at: `http://matplotlib.org/basemap/`

Yet another library, Pysolar, is used to calculate the solar position.
You don't *need* to download Pysolar to run any of the GITM plotting
scripts, but it does expand the functionality. Pysolar is documented at
Github and can be downloaded at: `http://www.pysolar.org`

Python scripts that create movies rely on external programs to do so.
The scripts included here use FFmpeg, which can be installed using Fink,
Macports, or downloaded at: `http://http://www.ffmpeg.org/`

If you have questions about these routines or are at the University of
Michigan and want to start using Python, Dr. Welling is the man to see.
The source code behind the PyBats GITM routines are also located in
`GITM2/srcPython`.

The following programs include the vetted python routines. The examples
shown are meant to be run in ipython, an interactive command-line
interface for python. The terminal window that is running ipython should
be located in the `GITM2/srcPython` directory. The test file for the
example code is one of the files output after running the default
`UAM.in` file.

### gitm.py

GITM is a PyBats submodule that handles input and output from GITM. It
can be helpful for those wishing to write their own GITM processing
routines but doesn't contain any analysis or visualization routines.

Once you have downloaded and installed Spacepy, the gitm submodule can
be accessed via:

Though to be sure that you have the latest version it is best to CVS the
latest version of `gitm.py` and load:

while running the `GITM2/srcPhython` directory. This module contains the
following routines:

-   **GitmBin:** A data class and routine to load a GITM output bin
    file, based on the PyBats data container class PbData. There are two
    keywords arguements that may be associated with this input,
    *ionfile*, and *varlist*. The first keyword, *ionfile*, takes a
    3DION file as input, assigns an attribute called 'ionfile' to the
    data class, and runs one of the functions called **calc_magvel** to
    calculate the ion and neutral velocities in magnetic coordinates.
    The second keyword, *varlist*, takes a list of atmospheric
    quantities such as those listed at the beginning of this chapter.
    More variables may always be added from the same file that created a
    data class by calling the function **append_data**. Geographic
    position in degrees and local time are also added to the output,
    regardless of what variables are specified. A 3DION file may also be
    associated with the output GITM class after the fact by assigning a
    filename to the attribute 'ionfile' and the magnetic velocity
    quantities obtained by running **calc_magvel**.

    -   **append_data:** A routine to append variables specified in an
        input list. These are obtained from the file specified in the
        'file' attribute.

    -   **append_units:** A routine to append unit, scale, and name
        attributes to the variable keys in the data class. Runs
        automatically with **GitmBin**.

    -   **calc_deg:** A routine to compute latitude and longitude in
        degrees instead of radians. Runs automatically with **GitmBin**.

    -   **calc_lt:** A routine to compute local time from universal time
        and longitude. Runs automatically with **GitmBin**.

    -   **calc_magdi:** A routine to compute the magnetic inclination
        and declination from the magnetic field expressed in
        East-North-Vertical coordinates. Runs automatically with
        **GitmBin** when the appropriate inputs are present.

    -   **calc_magvel:** A routine to compute ion and neutral velocities
        in magnetic coordinates. Runs automatically with **GitmBin** if
        the *ionfile* keyword is specified.

    -   **calc_tec:** A routine to calculate the VTEC from any GITM data
        structure that has altitude and electron density.

    -   **calc_2dion:** A routine to calculate the VTEC (if it hasn't
        been done already), $`h_mF_2`$, and $`N_mF_2`$ from any GITM data
        structure that has altitude and electron density

You can load a GITM binary file by entering the following commands.
Comments are preceded by a '\#'.

    In [1]: import spacepy
    In [2]: import gitm # Load the local version of gitm.py, which may be more up-to-date
    In [3]: gdata = gitm.GitmBin(`3DALL_t021124_000000.bin') # example binary file
    In [4]: gdata.attrs
    Out[4]:
    {`endian': `big',
     `file': `3DALL_t021124_000000.bin',
     `nAlt': 54,
     `nLat': 22,
     `nLon': 22,
     `nVars': 39,
     `version': 3.13}

The method used by **calc_2dion** to compute the $`h_mF_2`$ (height of the
F$`_2`$ region density peak) does more than just look for the height of
the electron maximum at the specified locations. Instead it performs a
few checks to ensure that the peak density lies within an altitude range
consistent with the F region and is not an edge artifact. This process
is outlined in Figure [1](#hmf2_flowchart.fig){reference-type="ref"
reference="hmf2_flowchart.fig"}.

<figure id="hmf2_flowchart.fig">
<div class="center">
<img src="Figures/determining_hmF2_flowchart.png" />
</div>
<figcaption>Process for finding the <span
class="math inline"><em>h</em><sub><em>m</em></sub><em>F</em><sub>2</sub></span>
from a GITM electron density altitude profile.</figcaption>
</figure>

One of the intermediate steps in
Figure [1](#hmf2_flowchart.fig){reference-type="ref"
reference="hmf2_flowchart.fig"} is to find the inflection points in the
altitude profile of electron density.
Figure [2](#hmf2_profile.fig){reference-type="ref"
reference="hmf2_profile.fig"} shows an example of an electron density
profile with no local maxima. To determine the height of the F$`_2`$ peak,
the height derivative of the electron density profile is computed and
local minima along the profile located. These local minima correspond to
inflection points in the electron density profile, and can be used to
locate a reasonable $`h_mF_2`$.

<figure id="hmf2_profile.fig">
<div class="center">
<img src="Figures/determining_hmF2_saddle_profile.png" />
</div>
<figcaption>Electron density profile (red) and height derivative of the
electron density profile (blue) for an instance where GITM does not show
a clear F<span class="math inline"><sub>2</sub></span> peak. The <span
class="math inline"><em>h</em><sub><em>m</em></sub><em>F</em><sub>2</sub></span>
is identified by choosing the inflection point (local minima from the
height derivative of the electron density profile) with the largest
electron density (dark grey line). The light grey lines show the
remaining locations of inflection points.</figcaption>
</figure>

### gitm_time.py

gitm_time has not yet been incorporated into PyBats, as it is actively
being developed. This module contains the following routines:

-   **GitmTime:** A data class and routine to load multiple GITM output
    binaries into a structure that includes a universal time (UT)
    dimension. All data types and attributes provided by GitmBin are
    provided in this data class.

    -   **appendgitm:** Add another GitmBin object to an existing
        GitmTime object. This routine is also be used to create a
        GitmTime object.

    -   **appendobs:** Add any type of data to an existing GitmTime
        object. The data can be directly appended or matched to the
        existing GitmTime data. Three match options are available:
        nearest neighbor, running average, and running median. The
        latter two options compute a central value and standard
        deviation using a specified location/time window about each
        GitmTime data point.

    -   **sat_dateloc_ticks:** Define axis ticks that include all the
        information necessary to know where measurements lie in
        spacetime. This is most useful when plotting data along a
        satellite orbit.

-   **load_multiple_gitm_bin:** A routine to load a list of GITM output
    binary files into GitmBin data structures. The output is a list of
    the GitmBin structures, where each element in the list contains the
    data from a GITM output binary.

-   **set_sat_dateloc_label:** Create a label for the ticks created by
    GitmTime.sat_dateloc_ticks. Outputs the label onto a plot on the
    right end of the x-axis.

You can load multiple GITM binary files into a data structure with UT
dependence by entering the following commands. Commands outside of
ipython are preceded by a `$`. Again, the output files used are
produced by running the default `UAM.in` file.

    $ ls 3DALL_t021124_000* > test.list
    $ ipython
    In [1]: import spacepy
    In [2]: import gitm # Load the local version of gitm.py, which may be more up-to-date
    In [3]: import gitm_time as gt
    In [4]: gtdata = gt.GitmTime(`test.list')
    In [5]: print gtdata['time'][:]
    [datetime.datetime(2002, 11, 24, 0, 0)
     datetime.datetime(2002, 11, 24, 0, 5)]
    In [6]: print gtdata['dLon'][0,:,10,27]
    [ -30.  -10.   10.   30.   50.   70.   90.  110.  130.  150.  170.  190.
      210.  230.  250.  270.  290.  310.  330.  350.  370.  390.]
    In [7]: print gtdata['dLat'][0,10,:,27]
    [-105.  -95.  -85.  -75.  -65.  -55.  -45.  -35.  -25.  -15.   -5.    5.
       15.   25.   35.   45.   55.   65.   75.   85.   95.  105.]
    In [8]: print gtdata['Altitude'][0,10,10,:]
    [  96660.90047544   98330.45023772  100000.          101669.54976228
      103348.94180562  105050.15497984  106788.53042153  108584.61276626
      110467.61429397  112482.48190265  114698.91401856  117219.51380361
      120200.6752566   123855.40706002  128245.43857411  133349.83464183
      139220.53740507  145894.02746965  153389.42319611  161708.33330952
      170836.59866663  180746.14929878  191399.58798926  202753.75931196
      214763.18484817  227382.89781152  240569.83529396  254283.5502381
      268486.11667059  283141.61144778  298215.46358093  313673.89525175
      329480.19139848  345605.83193525  362018.21624294  378685.73558327
      395578.19775965  412667.22510966  429926.74623085  447332.78448625
      464864.02590453  482501.63589793  500229.1662886   518032.39821451
      535899.14534895  553819.03817939  571783.3052721   589784.56266012
      607816.61836022  625874.29578468  643953.27746286  662049.96890339
      680161.38144647  698272.79398955]

### gitm_plot_rout.py

Common routines used to format and analyze GITM data.

-   **choose_contour_map:** A routine to choose an appropriate color map
    based on whether the plot will be black and white or color, and
    whether data data range will be centered about zero or not.

-   **add_colorbar:** Add a color bar to a contour plot.

-   **find_order_of_magnitude:** Find the order of magnitude of a
    specified number.

-   **center_polar_cap:** Adjust radial coordinates to produce a
    centered polar plot. Necessary for the northern hemisphere, where
    polar plots assume the radial (latitude) coordinates should be
    centered at zero instead of 90$`^\circ`$. This routine does not depend
    on SpacePy.

-   **find_data_limits:** Find the upper and lower limits for a
    specified data key in a list of GITM data structures at a specified
    location (either single indices or the entire range are permitted
    for latitude, longitude, and altitude).

-   **find_data_limits_irange:** Find the upper and lower limits for a
    specified data key in a list of GITM data structures at a specified
    location range (upper and lower limits or the entire range are
    permitted for latitude, longitude, and altitude).

-   **find_data_limits_ivalues:** Establish the appropriate axis limits
    for a list of GitmBin files at a range or specific latitude,
    longitude, or altitude.

-   **glon_to_localtime:** Compute the local time given a longitude and
    universal time.

-   **localtime_to_glon:** Find the longitude at a specified universal
    time and local time.

-   **find_lon_lat_index:** Find the indexes for the location closest to
    a specified latitude and longitude.

-   **retrieve_key_from_web_name:** Find a data key given a
    website-friendly version of data key names.

-   **find_alt_index:** Find the index closest to the specified
    altitude. Altitude may be specified in km or m.

-   **match_cindi_key:** A routine to retrieve a CINDI data key from a
    GITM key or vice versa.

-   **add_geomagnetic_equator:** A routine to add a line showing the
    geomagnetic equator (specified by IGRF-10) to an existing plot. Line
    style and color may be specified.

-   **add_subsolar_point:** A routine to find the location of the
    subsolar point at a specified Universal Time using the Pysolar
    routines. Returns the geographic location of the point and will also
    add a marker to a plot.

-   **add_solar_terminator:** Computes the location of the solar
    terminator using the Pysolar routines. Returns numpy arrays of the
    geographic coordinates of the solar terminator and will also add a
    line denoting the solar terminator to a plot.

-   **find_sunside_twilight_sza:** A routine to find the maximum angular
    distance between the solar terminator and the sunlight side with
    conjugate flux tube feet in darkness using Pysolar routines. The
    solar zenith angle corresponding to the sunlight boundary is
    returned. The entire day is searched to identify the sunlight
    boundary given any degree of magnetic declination.

-   **create_contour_input_array:** Creates contour input at a specified
    location between GITM grid points.

-   **create_linear_input_array:** Creates linear input at a specified
    location between GITM grid points.

-   **get_meq_offset:** Find the offset in degrees between the
    geographic and geomagnetic equators at a specified longitude.

This example shows how the index for a specified altitude can be found.
Note that GITM saves altitude in meters.

    In [1]: import spacepy
    In [2]: import gitm
    In [3]: import gitm_plot_rout as gpr
    In [4]: gdata = gitm.GitmBin(`3DALL_t021124_000000.bin')
    In [5]: ialt = gpr.find_alt_index(gdata, 10, 10, 250.0, "km")
    In [6]: print ialt, gdata['Altitude'][10,10,ialt]
    27 254283.550238

### solar_rout.py

Routines that use Pysolar to find the location of the solar terminator
and subsolar point.

-   **subsolar_point:** Finds the geographic location of the subsolar
    point at a specified Universal Time.

-   **lat_lon2spherical_xyz:** Converts latitude and longitude to
    spherical coordinates. Assumes a spherical earth.

-   **spherical_xyz2lat_lon:** Converts from spherical coordinates to
    geographic latitude and longitude, assuming a spherical earth.

-   **get_solar_terminator_xyz_matrix:** Finds the location of the solar
    terminator in spherical coordinates.

-   **get_terminator_lat_lon_coordinates:** Finds the location of the
    solar terminator in geographic coordinates using input from a
    solar-oriented spherical coordinate system.

-   **get_solar_terminator_lat_lon:** Finds the location of the solar
    terminator at a specified time in geographic coordinates.

### gitm_loc_rout.py

Routines used to find certain locations or values, as well as routines
to align or match points in different data sets.

-   **find_nearest_location:** A routine to find the nearest neighbor in
    a 1, 2, or 3D coordinate system.

-   **find_nearest_value:** A routine to find the nearest neighbor to a
    specified value.

-   **find_nearest_datetime:** A routine to find the nearest neighbor
    between datetime objects.

-   **match_running_average:** Provide running averages at specified
    times and locations.

-   **match_running_median:** Provide running medians at specified times
    and locations.

-   **gitm_inst_loc:** A routine to align GITM and instrument data.

-   **gitm_net_loc:** A routine to align GITM and data from a large
    network of instruments.

-   **gitm_time_obs_loc:** A routine to find, through interpolation, the
    GITM value at a specific observation location when only one spatial
    coordinate needs to be aligned.

### gitm_3D_global_plots.py

Routines to build and output GITM output variable contour (recommended)
or scatter plots over a geographic range. Several different standard
plot formats are available, and routines useful for creating custom
figures are also included. The Earth's continental boundaries may be
included in any output figure. If they are, shading in the night time
region of the globe may also be included.

-   **gitm_single_3D_image:** This is a basic visualization routine that
    creates a filled contour plot of a single output variable from a
    GITM 3D at a specified altitude or 2D bin file. The output variable
    is plotted as a function of latitude and longitude over the entire
    globe, though the latitude range may be limited. The output plot may
    be polar or rectangular and if the rectangular option is chosen, the
    geomagnetic equator may also be included in the output figure.
    Sample output of the electron temperature is shown in
    figure [3](#gitm_3D_global_plots.fig){reference-type="ref"
    reference="gitm_3D_global_plots.fig"} (a) and (b).

-   **gitm_single_nsglobal_3D_image:** A quick way to examine GITM
    output at both poles. This routine creates two polar contour plots
    centered at the geographic northern and southern poles for a single
    output variable from a GITM 3D at a specified altitude or 2D bin
    file. The equatorial and polar latitude boundaries may both be
    specified, though they cannot change between hemispheres. Sample
    output of the electron temperature is shown in
    figure [3](#gitm_3D_global_plots.fig){reference-type="ref"
    reference="gitm_3D_global_plots.fig"} (c)

-   **gitm_global_3D_snapshot:** A snapshot of a single GITM output over
    the entire globe. This routine creates two polar contour plots
    centered at the geographic northern and southern poles and extending
    to a specified latitude and a single rectangular plot containing the
    latitudes equatorward of this point for a single output variable
    from a GITM 3D at a specified altitude or 2D bin file. The
    geomagnetic equator may also be included in the output figure.
    Sample output of the electron temperature is shown in
    figure [3](#gitm_3D_global_plots.fig){reference-type="ref"
    reference="gitm_3D_global_plots.fig"} (d)

-   **gitm_mult_3D_slices:** This routine creates a single plot
    containing multiple global contours of a GITM output variable from a
    3D or 2D bin file at a list of specified altitudes. These plots may
    be either polar or rectangular, with or without the geomagnetic
    equator, and within a specified latitude range. Sample output of the
    electron temperature is shown in
    figure [4](#gitm_3D_mult_plots.fig){reference-type="ref"
    reference="gitm_3D_mult_plots.fig"}.

<figure id="gitm_3D_global_plots.fig">
<div class="center">

</div>
<figcaption>GITM electron temperature at 456.63 <span
class="math inline"><em>k</em><em>m</em></span> altitude for: (a)
northern latitudes, (b) over the entire globe, (c) over the entire
globe, as viewed from the poles, and (d) as a global
snapshot.</figcaption>
</figure>

<figure id="gitm_3D_mult_plots.fig">
<div class="center">

</div>
<figcaption>GITM electron temperature at seven altitude slices for (a)
northern latitudes and (b) the entire globe.</figcaption>
</figure>

This example shows how to reproduce
Figure [3](#gitm_3D_global_plots.fig){reference-type="ref"
reference="gitm_3D_global_plots.fig"} (a).

    In [1]: import spacepy
    In [2]: import gitm
    In [3]: import gitm_3D_global_plots as g3d
    In [4]: import matplotlib.pyplot as plt
    In [5]: plt.ion() # This makes the plotting happen interactively
    In [6]: gdata = gitm.GitmBin(`3DALL_t021124_000000.bin')
    In [7]: title = "%s UT" % (gdata[`time'])
    In [8]: f = g3d.plot_single_3D_image("polar", "eTemperature", gdata, title, 
                                         "example_polar_plot.png", True, 27, 90, 0)

### plot_3D_global.py

Routines to build contour or scatter plots over a geographic range.
Several different standard plot formats are available, and routines
useful for creating custom figures are also included. The Earth's
continental boundaries may be included in any output figure. If they
are, shading in the night time region of the globe may also be included.
Input data must be provided in separate numpy arrays (1D or 2D for
scatter plots, 2D for contour plots).

-   **plot_single_3D_image:** This is a basic visualization routine that
    creates a filled contour or a colored scatter plot of a single
    output variable. The output variable is plotted as a function of
    latitude and longitude over the entire globe, though the latitude
    range may be limited. The output plot may be polar or rectangular
    and if the rectangular option is chosen, the geomagnetic equator may
    also be included in the output figure. If the polar option is used,
    the latitude range cannot extend beyond $`\pm`$ 90$`^\circ`$. If the
    rectangular option is used, the geomagnetic equator may also be
    plotted.

-   **plot_single_nsglobal_3D_image:** A quick way to examine a variable
    at both poles. This routine creates two polar contour or scatter
    plots centered at the geographic northern and southern poles for a
    single output variable. The equatorial and polar latitude boundaries
    may both be specified, though they cannot change between
    hemispheres.

-   **plot_global_3D_snapshot:** A snapshot of an output variable over
    the entire globe. This routine creates two polar contour or scatter
    plots centered at the geographic northern and southern poles and
    extending to a specified latitude as well as a single rectangular
    plot containing the latitudes equatorward of this point. The
    geomagnetic equator may be output over the data.

-   **plot_mult_3D_slices:** This routine creates a single plot
    containing multiple global contour or scatter plots for a single
    variable at specific indices (corresponding to different altitude,
    universal times, *et cetera*). These plots may be either polar or
    rectangular (with or without the geomagnetic equator). Sample output
    of the electron temperature is shown in
    figure [4](#gitm_3D_mult_plots.fig){reference-type="ref"
    reference="gitm_3D_mult_plots.fig"}. The numpy arrays containing the
    data may be 2D or 3D for scatter plots or 3D for contour plots. Any
    dimension may contain the indices to iterate over for the
    subfigures.

-   **plot_nsglobal_subfigure:** This routine creates a subfigure with
    two polar contour or scatter plots centered at the geographic
    northern and southern poles for a single output variable. The
    equatorial and polar latitude boundaries may both be specified,
    though they cannot change between hemispheres. This is used by
    plot_single_nsglobal_3D_image and may also be used to create a
    subplot with this format.

-   **plot_snapshot_subfigure:** This routine creates a subfigure with
    two polar dials and a rectangular region showing the equatorial
    latitudes for a single output variable. This is used by
    plot_global_3D_snapshot and may also be used to create a subplot
    with this format.

-   **plot_rectangular_3D_global:** This routine plots a single
    rectangular filled contour or colored scatter for an output variable
    as a function of latitude and longitude. Options exist to control
    the colorbar, ticks, labels, and more. A handle to the contour plot
    is returned to allow the output to be further manipulated depending
    on what other subplots are included in the output figure.

-   **plot_polar_3D_global:** This routine plots a single polar filled
    contour for a GITM output variable at a specified altitude index as
    a function of latitude and longitude. Title, colorbar, labels, and
    more may be specified using input options. A handle to the contour
    plot is returned to allow the output to be further manipulated
    depending on what other subplots are included in the output figure.
    The longitude at the top of the plot may also be specified, this
    allows one to ensure a specific local time is always located at the
    top of the dial using a routine like **localtime_to_glon**.

### gitm_alt_plots.py

Routines to build and output GITM output variable linear and contour
plots over an altitude range. Several different standard plot formats
are available, and routines useful for creating custom figures are also
included.

-   **gitm_single_alt_image:** Creates a single linear or contour
    altitude plot.

-   **gitm_mult_alt_image:** Creates a figure with multiple linear or
    contour altitude plots.

-   **gitm_alt_slices:** Creates a figure with a contour plot showing
    the altitude dependence of a quantity as a function of latitude or
    longitude with several linear altitude slices at specified
    locations. An example is shown in
    figure [5](#gitm_alt_slices.fig){reference-type="ref"
    reference="gitm_alt_slices.fig"}

<figure id="gitm_alt_slices.fig">
<div class="center">
<img src="Figures/gitm_alt_slice_test_Te.png" />
</div>
<figcaption>GITM electron temperature at a constant longitude with six
latitude slices.</figcaption>
</figure>

This example shows how to reproduce
Figure [5](#gitm_alt_slices.fig){reference-type="ref"
reference="gitm_alt_slices.fig"}.

    In [1]: import spacepy
    In [2]: import gitm
    In [3]: import gitm_alt_plots as gap
    In [4]: import gitm_plot_rout as gpr
    In [5]: import matplotlib.pyplot as plt
    In [6]: plt.ion() # This makes the plotting happen interactively
    In [7]: gdata = gitm.GitmBin(`3DALL_t021124_000000.bin')
    In [8]: title = "%s UT" % (gdata[`time'])
    In [9]: lat_index = list()
    In [10]: lon_index = list()
    In [11]: (ilon, ilat) = gpr.find_lon_lat_index(gdata, 150.0, -65.0, "degrees")
    In [12]: lon_index.append(ilon)
    In [13]: lat_index.append(ilat)
    In [14]: (ilon, ilat) = gpr.find_lon_lat_index(gdata, 150.0, -45.0, "degrees")
    In [15]: lat_index.append(ilat)
    In [16]: (ilon, ilat) = gpr.find_lon_lat_index(gdata, 150.0, -5.0, "degrees")
    In [17]: lat_index.append(ilat)
    In [18]: (ilon, ilat) = gpr.find_lon_lat_index(gdata, 150.0, 5.0, "degrees")
    In [19]: lat_index.append(ilat)
    In [20]: (ilon, ilat) = gpr.find_lon_lat_index(gdata, 150.0, 45.0, "degrees")
    In [21]: lat_index.append(ilat)
    In [22]: (ilon, ilat) = gpr.find_lon_lat_index(gdata, 150.0, 65.0, "degrees")
    In [23]: lat_index.append(ilat)
    In [24]: f = gap.plot_alt_slices("eTemperature", gdata, lat_index, lon_index, 
                                     title, "example_alt_plot.png")

### plot_alt_profiles.py

Routines to build and output linear and contour plots over an altitude
range. Several different standard plot formats are available, and use
any numpy array as input.

-   **plot_single_alt_image:** Creates a single linear or contour
    altitude plot.

-   **plot_mult_alt_image:** Creates a figure with multiple linear or
    contour altitude plots.

-   **plot_alt_slices:** Creates a figure with a contour plot showing
    the altitude dependence of a quantity as a function of latitude or
    longitude with several linear altitude slices at specified
    locations.

-   **plot_linear_alt:** Plots the the linear altitude dependence of a
    quantity, with altitude on the y-axis.

-   **plot_3D_alt:** Plots the altitude dependence of a quantity as the
    function of another spatiotemporal coordinate with the
    spatiotemporal coordinate on the x-axis, altitude on the y-axis, and
    the desired quantity as a color contour.

### gitm_comparison_plots.py

Routines to make plots that compare GITM data with observations. The
observational sources include satellites, ground-based receivers, and
receiver networks.

-   **extract_data_matched_arrays:** Extract points from matched data
    arrays for elements where neither array contains a specified 'bad'
    value.

-   **extract_gitm_time_arrays:** Routine to extract all positions with
    valid data from a GitmTime object and construct numpy arrays. A
    single universal time, longitude, latitude, and/or altitude may be
    specified.

-   **plot_net_gitm_comp:** A routine to create a plot comparing 2D GITM
    data (such as VTEC or $`h_mF_2`$) to observations taken from a network
    of instruments over the globe. The map format may be rectangular,
    polar, or a combination (provided by plot_nsglobal_subfigure or
    plot_snapshot_subfigure). The top subfigure shows the observations
    as a scatter figure, the middle subfigure shows the GITM data as a
    contour, and the bottom subfigure shows the difference between the
    two as a scatter plot. The difference must be computed outside of
    this program, and so may be the difference, absolute difference,
    percent difference, or any other type of comparison. An example is
    shown in Figure [6](#teccomp.fig){reference-type="ref"
    reference="teccomp.fig"}.

-   **plot_sat_gitm_comp:** Routine to plot satellite and GITM data to
    show how a single physical quantity varies over the orbit. Four
    panels are included; the top panel shows the raw satellite data and
    the GITM data along the track. The second panel shows the matched
    GITM/satellite data. The third panel shows the difference between
    the satellite and GITM data. The fourth panel shows the percent
    difference 100\*(sat-GITM)/sat.

<figure id="teccomp.fig">
<div class="center">
<img src="Figures/gitm_tec_comp.png" />
</div>
<figcaption>Madrigal and GITM vertical TEC comparison.</figcaption>
</figure>

### load_files.py

Routines to load certain types of data files into a dictionary of numpy
arrays, where each data type is used to specify the dictionary keys.

-   **loadCINDIorbit_ASCII:** Loads the Coupled Ion Neutral Dynamics
    Investigation ASCII files provided by the UT Dallas website.

-   **loadMadrigalVTEC_ASCII:** Loads the GPS TEC ASCII files provided
    by the Millstone Hill Madrigal site.

-   **loadGITMsat_ASCII:** Loads the satellite file used as input for a
    GITM run.

-   **loadMadrigalVTEC_HDF5:** Loads the GPS TEC HDF5 files provided by
    the Millstone Hill Madrigal site.

### read_files.py

Routines to read certain file formats and load the data into a
dictionary of numpy arrays, where each data type is used to specify the
dictionary keys.

-   **loadASCII_data_header:** Loads an ASCII file with header lines
    denoted by a '\#' .

-   **loadASCII_data_hline:** Loads an ASCII file with a specified
    number of header lines.

-   **loadASCII_index_profile:** Loads an ASCII file with header lines
    denoted by a '\#' that has been broken up into indexed blocks
    (blocks separated by double newlines, or indexes as specified by
    gnuplot). The indexed structure is maintained in the output
    dictionary by providing a list of numpy arrays for each data column.

-   **load_multASCII_data:** Loads multiple ASCII files into a single
    output dictionary.

-   **loadnetCDF_data:** Loads netCDF files into an output dictionary.

-   **combine_data_dictionaries:** Combines multiple data dictionaries
    into a single dictionary, including only the data keys common to all
    of the inputted data dictionaries.

### read_gps_bin.py

A script to read a VMR GPS binary file, providing data at a specified
Universal Time.

-   **GpsFile:** Class containing all of the data from the VMR GPS
    binary file.

    -   **read_header:** Reads the VMR GPS binary file header, allowing
        data to be easily located as desired.

    -   **read_time:** Read in all the GPS data at the specified
        Universal Time. The specified time and the file times do not
        have to be perfectly aligned, as long as the specified time
        falls between the first and last time in the GPS file, the data
        with the closest temporal proximity will be returned.

### write_files.py

Routines to write output files.

-   **writeASCII_file:** A routine to create an ASCII file from a string
    or list of strings. Will overwrite any file of the same name that
    already exists.

-   **writeASCII_data_w_sorttext:** A routine to create an ASCII file
    from a data dictionary of dictionaries. The first layer of keys is
    used to provide the data columns, the second layer is output as
    additional data column(s) where the keys are output as strings.
    Datetime columns are output as two strings, one containing the date
    information and a second one containing the time of day information.

### plot_stats.py

Routines to compute and plot common statistics.

-   **add_stat_to_line:** Computes the moments, first through third
    quartiles, and mode(s) (as desired) for a dataset and outputs the
    statistics as a formatted text line that can be output to a file and
    in a list.

-   **add_stat_box:** Computes the moments, first through third
    quartiles, and mode(s) (as desired) for a dataset and outputs the
    statistics in a text box onto a plot and in a list.

-   **plot_hist_w_stats:** Calculates and plots a histogram of a
    specified dataset, as well as the moments, first through third
    quartiles, and mode(s) (as desired).

-   **plot_lat_lt_stats:** Calculates and plots histograms and
    statistics for a specified dataset broken up into latitude and local
    time regions.

-   **lat_lt_stat_lines:** Calculates statistics for a specified dataset
    broken up into latitude and local time regions, providing a
    formatted string with the statistics for each region taking up a
    line.

### gitm_movie_script.py

This is a python script that can be run either from ipython using the
command `run gitm_movie_script.py` or the command line using the command
`python gitm_movie_script.py`. Input to this program is prompted
interactively, and includes:

-   `Ordered list of GITM binary files:` A list of GITM binary files in
    chronological order (or whatever other order the movie should be
    played in).

-   `GITM plot type (rectangular, polar, nspolar, snapshot):` The
    keyword for the desired plot type. These are the plot types shown in
    figure [3](#gitm_3D_global_plots.fig){reference-type="ref"
    reference="gitm_3D_global_plots.fig"}, where polar corresponds to
    panel (a), rectangular to panel (b), nspolar to panel (c), and
    snapshot to panel (d).

-   At this point, the routine enters a while-loop to allow multiple
    movies to be made for the same list of GITM binary files

    -   `GITM key to plot on z axis (eg Temperature):` The data key
        corresponding to the data type to plot on the z axis. A list of
        data keys can be found by typing '`gdata.keys()`' into ipython
        after loading one of the listed GITM binaries.

    -   `Altitude to plot z value at (eg 250):` Altitude to plot, may be
        specified in km or m. For 2D parameters, a value must be
        entered, but doesn't matter.

    -   `Units of altitude (km or m):` Units of altitude entered above.

    -   `Use map of Earth? (empty for False):` Enter any value to
        include a Basemap plot of the earth, enter a carriage return to
        exclude the map.

    -   The latitude limits needed depend on the plot type

        -   **nspolar** `Polar latitude limit (degrees):` Specify the
            polar latitude limit (positive, same for both hemispheres).

        -   **nspolar** `Equatorial latitude limit (degrees):` Specify
            the equatorial latitude limit (positive, same for both
            hemispheres).

        -   **snapshot** `Polar latitude limit (degrees):` Specify the
            polar latitude limit (positive, same for both hemispheres).

        -   **polar/rectangular** `Northern latitude limit (degrees):`
            Specify the northernmost latitude (may be negative).

        -   **polar/rectangular** `Southern latitude limit (degrees):`
            Specify the southernmost latitude (must be smaller/more
            negative than the northernmost limit)

    -   `Load another z axis key? (empty for False):` Enter any value to
        include another movie or, enter a carriage return finish.

With this information, movies with appropriate z-variable ranges will be
plot as .png files and combined into a movie using FFmpeg. The image and
movie files will be named using the plot type, z parameter, and altitude
to distinguish them.
