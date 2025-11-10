# All Inputs

When running, GITM looks for a file called `UAM.in` in the same directory as the
executable. This is where all of the settings for a run are read from.

This plain-text file can contain comments and other notes which will not be used. Only
sections that begin with a `#` and match one of the GITM's settings are actually read.
GITM will print a message if it finds a line that begins with `#` that was not used.

For example, in the following snippet, only the first block of text will be recognized
as a valid setting:

```
#LOGFILE
10.0

# LOGFILE
10.0       not read, space after # in first line

#LOG_FILE
10.0       not read, key does not match

#LOGFILE
 10.0      will cause error, space before value!

```

This allows us to add comments or descriptors, if we wish.

## General Configuration

### STARTTIME

This sets the starting time of the simulation. Even when you restart,
the starttime should be to the real start time, not the restart time.

    #STARTTIME
    iYear    (integer)
    iMonth   (integer)
    iDay     (integer)
    iHour    (integer)
    iMinute  (integer)
    iSecond  (integer)

### ENDTIME

This sets the ending time of the simulation.

    #ENDTIME
    iYear    (integer)
    iMonth   (integer)
    iDay     (integer)
    iHour    (integer)
    iMinute  (integer)
    iSecond  (integer)


### ALTITUDE

For Earth, the AltMin is the only variable used here. The altitudes are
set to 0.3 times the scale height reported by MSIS, at the equator for
the specified F107 and F107a values.`

    #ALTITUDE
    AltMin                (real, km)
    AltMax                (real, km)
    UseStretchedAltitude  (logical)

### GRID

If LatStart and LatEnd are set to < -90 and > 90, respectively, then
GITM does a whole sphere. If not, it models between the two. If you want
to do 1-D, set nLons=1, nLats=1 in ModSizeGitm.f90, then recompile, then
set LatStart and LonStart to the point on the Globe you want to model.

    #GRID
    nBlocksLon   (integer)
    nBlocksLat   (integer)
    LatStart     (real)
    LatEnd       (real)
    LonStart     (real)
    LonEnd       (real)


### RESTART

    #RESTART
    DoRestart (logical)

When doing a restart, GITM will look for "restart files" in `UA/restartIN/`. These files
are automatically created when running, so you do not have to move anything to restart
a failed run.

A common use case for restarts is where we stop a run after initializing and then
restart it with different settings. 
For example, maybe we scheduled a run to stop at 12:00 after 36 hours of 
initialization, but now want it to continue with some different settings. To do this:

- Create new run directories for each of the different settings. Since GITM outputs are
  named with the output time, the subsequest runs may overwrite existing files. It is
  up to you on how you want to do this.
- Copy the restart files from the initial runs so the new runs can find them. We need
  all of the files in UA/restartOUT/, so one way to do this is to run
  `cp run_initial/UA/restartOUT/* run_setting01/UA/restartIN/`.
- Copy all of the input files (UAM.in, IMF, SME, etc) from the initial run directory to
  the new one: 
  `cp run_initial/*.dat run_initial/UAM.in run_initial/*.txt run_setting01/`
- Modify the new UAM file however you want. Take note of the important notes below
- Run `GITM.exe` from the new directory

!!! note "Some important notes on restarts"
    - Do not change the start time when restarting! The start time in the original and
      "restart" UAM.in file must be the same. This is because of how GITM handles
      temperature and altitude.
    - The end time can be increased in the new UAM.in file. 
    - If you do not specify the `#RESTART` flag GITM will not know to restart, even if the
      files are present.

### SAVEPLOT

The DtRestart variable sets the time in between writing full restart
files to the UA/restartOUT directory. This sets the output files. The most common
type is 3DALL, which outputs all primary state variables. 

Other types include : 3DALL, 3DNEU, 3DION, 3DTHM,
3DCHM, 3DUSR, 3DGLO, 2DGEL, 2DMEL, 2DUSR, 1DALL, 1DGLO, 1DTHM, 1DNEW,
1DCHM, 1DCMS, 1DUSR. DtPlot sets the frequency of output

    #SAVEPLOT
    DtRestart (real, seconds)
    nOutputTypes  (integer)
    Outputtype (string, 3D, 2D, ION, NEUTRAL, ...)
    DtPlot    (real, seconds)

### PLOTTIMECHANGE

This allows you to change the output cadence of the files for a limited
time. If you have an event then you can output much more often during
that event.

Specify the start and end time the new cadence will be used for, followed by the new Dt
for each output type. This example would be for two output types:

    #PLOTTIMECHANGE
    yyyy mm dd hh mm ss ms (start)
    yyyy mm dd hh mm ss ms (end)
    DtOutputType1           (real, seconds)
    DtOutputType2           (real, seconds)


### LOGFILE

You really want a log file. They are very important. It is output in
UA/data/log\*.dat. You can output the log file at whatever frequency you
would like, but if you set dt to some very small value, you will get an
output every iteration, which is probably a good thing.

    #LOGFILE
    DtLogFile   (real, seconds)

### CCMCFILENAME

Typicaly file is named (e.g.) 1DALL_yymmdd_hhmmss.bin With this it will
be named 1DALL_GITM_yyyy-mm-ddThh-mm-ss.bin

    #CCMCFILENAME
    UseCCMCFileName    (logical)


### SAVEHIMEPLOT

    #SAVEHIMEPLOT
    HIMEPlotLonStart (real)
    HIMEPlotLonEnd (real)
    HIMEPlotLatStart (real)
    HIMEPlotLatEnd (real)

### SATELLITES

See [this page](common_inputs.md#satellite-section) for more details.

    #SATELLITES
    nSats     (integer - max = ', nMaxSats, ')
    SatFile1  (string)
    SatOutputtype1  (string, 0DUSR or 1DUSR or other)
    DtPlot1   (real, seconds)
    etc...


### APPENDFILES

For satellite files, you can have one single file per satellite, instead
of one for every output. This makes GITM output significantly less
files. It only works for satellite files now.

    #APPENDFILES
    DoAppendFiles    (logical)

## Drivers

### F107

Sets the F10.7 and 81 day average F10.7. This is used to set the initial
altitude grid, and drive the lower boundary conditions.

    #F107
    f107  (real)
    f107a (real - 81 day average of f107)

### ELECTRODYNAMICS

Sets the sources for, and the time for updating, the high-latitude (and low-latitude)
electrodynamic drivers, such as the potential and the aurora.

See [the electrodynamics page](internals/electrodynamics.md) for more details

    #ELECTRODYNAMICS
    AuroralModel     [fta], fre, pem, ovation, hpi/ihp, amie, etc.
    DtAurora    (real, seconds [60.0])
    PotentialModel   [weimer05], hepmay, amie, etc.
    DtPotential (real, seconds [60.0])'


### HPI

This sets the hemispheric power of the aurora. Typical it ranges from
1-1000, although 20 is a nominal, quiet time value.

    #HPI
    HemisphericPower  (real)

### KP

GITM will not use this unless the Foster electric
field model is used.

    #KP
    kp  (real)


### SOLARWIND

This sets the driving conditions for the high-latitude electric field
models. This is static for the whole run, though. It is better to use
the MHD_INDICES command to have dynamic driving conditions.

    #SOLARWIND
    bx  (real)
    by  (real)
    bz  (real)
    vx  (real)

### MHD_INDICES

Use this for dynamic IMF and solar wind conditions. The exact format of
the file is discussed further in the manual.

    #MHD_INDICES
    filename  (string)

### EUV_DATA

This is for a FISM or some other solar spectrum file.

    #EUV_DATA
    UseEUVData            (logical)
    cEUVFile              (string)

### AURORAMODS

This is for modifying the aurora a bit. The NormalizeAuroraToHP variable
calculates the modeled hemispheric power and then normalizes it the
hemispheric power read in. AveEFactor - changes the aveE of the aurora
by factor IsKappaAurora - use a kappa instead of Maxwellian
AuroraKappa - kappa to use in the distribution

    #AURORAMODS
    NormalizeAuroraToHP     (logical)
    AveEFactor    (real)
    IsKappaAurora     (logical)
    AuroraKappa    (real)

### AURORATYPES

This is for using different types of aurora. Currently this is supported by Ovation & 
AMIE. The first three options (diffuse, mono, wave) are for electrons only.

    #AURORATYPES
    UseDiffuseAurora   (logical)
    UseMonoAurora (logical)
    UseWaveAurora (logical)
    UseIonAurora (logical)

### USECUSP

This is for specifying a cusp.

    #USECUSP
    UseCusp        (logical)
    CuspAveE       (real)
    CuspEFlux      (real)

### AMIEFILES

    #AMIEFILES
    cAMIEFileNorth  (string)
    cAMIEFileSouth  (string)

### AMIENORTH

There have been issues with some filesystems (Lustre) not allowing large files to be
read across large numbers of processors simultaneously. To deal with this, AMIE files
can be split when they reach a certain size when they are created using the Python
routines in Electrodynamics. 

For these to be read by GITM, we first say how many files there are and then specify
the names of these files. This option, and the option below, are replacements for
#AMIEFILES (above) - #AMIENORTH/#AMIESOUTH should not be used at the same time as 
#AMIEFILES.

    #AMIENORTH
    nAMIENorth       (int)
    cAMIEFileNorth1  (string)
    ...              (string)s
    cAMIEFileNorth(n)  (string)


### AMIESOUTH

See above...

    #AMIESOUTH
    nAMIESouth       (int)
    cAMIEFileSouth1  (string)
    ...              (string)s
    cAMIEFileSouth(n)  (string)


### USEREGIONALAMIE

This is to set up a local region with specified potential from AMIE
files. Use Weimer potential elsewhere. AMIEBoundaryWidth is padded
outside of the region with the geographic lon and lat boundaries set
below.

    #USEREGIONALAMIE
    UseRegionalAMIE      (logical)
    UseTwoAMIEPotentials (logical)
    AMIETimeStart        (yyyy mm dd hh mm ss)
    AMIETimeEnd          (yyyy mm dd hh mm ss)
    AMIELonStart         (real)
    AMIELonEnd           (real)
    AMIELatStart         (real)
    AMIELatEnd           (real)
    AMIEBoundaryWidth    (real)

### INPUTTIMEDELAY

Sets the time delay for the high latitude drivers and solar EUV inputs.

    #INPUTTIMEDELAY
    TimeDelayHighLat (real, seconds)
    TimeDelayEUV     (real, seconds)


## Boundary and Initial Conditions

### INITIAL

This specifies the initial conditions and the lower boundary conditions.
For Earth, we typically just use MSIS and IRI for initial conditions.
For other planets, the vertical BCs can be set here.

    #INITIAL
    UseMSIS        (logical)
    UseIRI         (logical)
    If UseMSIS is .false. then :
    TempMin        (real, bottom temperature)
    TempMax        (real, top initial temperature)
    TempHeight     (real, Height of the middle of temp gradient)
    TempWidth      (real, Width of the temperature gradient)

### MSISTIDES

This says how to use msis tides. The first argument is for using diurnal tide, the
second is for using semi-diurnal tide and the third is to use terdiurnal
tide.

    #MSISTIDES
    UseMSISDiurnal        (logical)
    UseMSISSemidiurnal    (logical)
    UseMSISTerdiurnal     (logical)

### MSISOBC

UseOBCExperiment - use MSIS \[O\] BC shifted by 6 months Only applicable
for MSIS00! MsisOblateFactor - alt = alt \* (1.0-f/2 + f\*cos(lat)) -
seems like -0.1 works well

    #MSISOBC
    UseOBCExperiment        (logical)
    MsisOblateFactor           (real)

### MSIS21

This toggles between using MSIS00 (false) and MSIS-2.1 (true)

    #MSIS21
    UseMsis21       (logical)


### TIDES

This says how to use tides. The first one is using MSIS with no tides.
The second is using MSIS with full up tides. The third is using GSWM
tides, while the forth is for experimentation with using WACCM tides.

    #TIDES
    UseMSISOnly        (logical)
    UseMSISTides       (logical)
    UseGSWMTides       (logical)
    UseWACCMTides      (logical)
    UseHmeTides        (logical)

### GSWMCOMP

If you selected to use GSWM tides above, you can specify which
components to use.

    #GSWMCOMP
    GSWMdiurnal(1)        (logical)
    GSWMdiurnal(2)        (logical)
    GSWMsemidiurnal(1)    (logical)
    GSWMsemidiurnal(2)    (logical)

### USEPERTURBATION

    #USEPERTURBATION
    UsePerturbation        (logical)

### USEBCPERTURBATION

Use Boundary Condition Perturbation are a set of arguments for running WP-GITM.

    #USEBCPERTURBATION
    UseBcPerturbation        (logical)
    If UseBcPerturbation = .true. then:
    iTypeBcPerturb         (int) 
    perturbation characteristics ...

### GITMBCS

    #GITMBCS
    UseGitmBCs
    GitmBCsDir

## Numerics


### STATISTICALMODELSONLY

This command will skip pretty much all of the physics of GITM, and
will reinitialize the model with the MSIS, IRI, APEX, and Electrodynamics values at
the interval set in the second variable. If you want to compare a run to MSIS and
IRI, you can run GITM with this command and get output at exactly the
same cadence and locations, thereby allowing easier comparisons. The dt
can be set as low as needed, so you can run satellites through MSIS and
IRI.

    #STATISTICALMODELSONLY
    UseStatisticalModelsOnly    (logical)
    DtStatisticalModels         (real)

### CFL

The CFL condition sets how close to the maximum time step that GITM will
take. 1.0 is the maximum value. A value of about 0.75 is typical. If
instabilities form, a lower value is probably needed.

    #CFL
    cfl  (real)



### LIMITER

The limiter is quite important. It is a value between 1.0 and 2.0, with
1.0 being more diffuse and robust, and 2.0 being less diffuse and less
robust.

    #LIMITER
    TypeLimiter  (string)
    BetaLimiter  (real between 1.0-minmod and 2.0-mc)

### NANCHECK

This will turn on all of the NaN checks in the code!

    #NANCHECK
    DoCheckForNans (logical)

### DEBUG

This will set how much information the code screams at you - set to 0 to
get minimal, set to 10 to get EVERYTHING. You can also change which
processor is shouting the information - PE 0 is the first one. If you
set the iDebugLevel to 0, you can set the dt of the reporting. If you
set it to a big value, you wont get very many outputs. If you set it to
a tiny value, you will get a LOT of outputs. UseBarriers will force the
code to sync up a LOT.

    #DEBUG
    iDebugLevel (integer)
    iDebugProc  (integer)
    DtReport    (real)
    UseBarriers (logical)

### IONLIMITS

    #IONLIMITS
    MaxVParallel     (real, default=100 m/s)
    MaxEField        (real, default=0.1 V/m)
    MinIonDensity    (real, default=100 m^-3)
    MinIonDensityAdvect    (real, default=1e5 m^-3)

### NEUTRALLIMITS

Neutral limits only get applied if the deniity of the species drops below zero.
This is very rare and only occurs a few times when simulating extreme event. This 
likely does not need to be changed.

As with the ion limits, the lower limit for advected vs non-advected species can be set
differently.

    #NEUTRALLIMITS
    MinNeutralDensity    (real, default=100 m^-3)
    MinNeutralDensityAdvect    (real, default=1e5 m^-3)

### PHOTOELECTRON

    #PHOTOELECTRON
    PhotoElectronHeatingEfficiency   (real)

### NEUTRALHEATING

    #NEUTRALHEATING
    NeutralHeatingEfficiency   (real)

### DON4SHACK

    #DON4SHACK
    DoN4SHack       (logical)


### VERTICALSOURCES

    #VERTICALSOURCES
    MaximumVerticalVelocity      (real)

### AUSMSOLVER

This dictates whether to use the AUSM solver (T) or the Rusanov solver (F)

    #AUSMSOLVER
    Use AUSM Solver      (logical)


### USEIMPLICITIONMOMENTUM

    #USEIMPLICITIONMOMENTUM
    UseImplicitFieldAlignedMomentum      (logical)

### USEIMPROVEDIONADVECTION

    #USEIMPROVEDIONADVECTION
    UseImprovedIonAdvection      (logical)
    UseNighttimeIonBCs           (logical)
    MinTEC                       (real)

### USETESTVISCOSITY

    #USETESTVISCOSITY
    TestViscosityFactor      (real)

### FIXEDDT

If you would like to force GITM to take a fixed dt you can use this. It
will try to take that fixed dt, unless the CFL condition is violated.

    #FIXEDDT
    FixedDt  (real)

## Physics


### THERMO

    #THERMO
    UseSolarHeating   (logical)
    UseJouleHeating   (logical)
    UseAuroralHeating (logical)
    UseNOCooling      (logical)
    UseOCooling       (logical)
    UseConduction     (logical)
    UseTurbulentCond  (logical)
    UseIRHeating      (logical)

### THERMALDIFFUSION

    #THERMALDIFFUSION
    KappaTemp0    (thermal conductivity, real)

### THERMALCONDUCTION

    #THERMALCONDUCTION
    ThermalConduction_AO2 (Conduction A(O2): 3.6e-4, real)
    ThermalConduction_AO  (Conduction A(O): 5.6e-4, real)
    ThermalConduction_s   (Conduction s: 0.75, real)


### DIFFUSION

If you use eddy diffusion, you must specify two pressure levels - under
the first, the eddy diffusion is constant. Between the first and the
second, there is a linear drop-off. Therefore The first pressure must be
larger than the second!

    #DIFFUSION
    UseDiffusion (logical)
    EddyDiffusionCoef (real)
    EddyDiffusionPressure0 (real)
    EddyDiffusionPressure1 (real)

### FORCING

    #FORCING
    UsePressureGradient (logical)
    UseIonDrag          (logical)
    UseNeutralFriction  (logical)
    UseViscosity        (logical)
    UseCoriolis         (logical)
    UseGravity          (logical)

### MODIFYPLANET

    #MODIFYPLANET
    RotationPeriodInput        (real)
    DaysPerYearInput           (real)
    DaysPerYearInput           (real)

### DYNAMO

    #DYNAMO
    UseDynamo              (logical)
    DynamoHighLatBoundary  (real)
    nItersMax              (integer)
    MaxResidual            (V,real)
    IncludeCowling         (logical)
    DynamoLonAverage       (real)

### IONFORCING

    #IONFORCING
    UseExB                 (logical)
    UseIonPressureGradient (logical)
    UseIonGravity          (logical)
    UseNeutralDrag         (logical)

### FIXTILT

    #FIXTILT
    IsFixedTilt   (logical)

### DIPOLE

    #DIPOLE
    MagneticPoleRotation   (real)
    MagneticPoleTilt       (real)
    xDipoleCenter          (real)
    yDipoleCenter          (real)
    zDipoleCenter          (real)

### APEX

    #APEX
    UseApex (logical)
            Sets whether to use a realistic magnetic
            field (T) or a dipole (F)


## Misc

### CPUTIMEMAX

This sets the maximum CPU time that the code should run before it starts
to write a restart file and end the simulation. It is very useful on
systems that have a queueing system and has limited time runs.
Typically, set it for a couple of minutes short of the max wall clock,
since it needs some time to write the restart files.

    #CPUTIMEMAX
    CPUTimeMax    (real)


### PAUSETIME

This will set a time for the code to just pause. Really, this should
never be used.

    #PAUSETIME
    iYear iMonth iDay iHour iMinute iSecond

### ISTEP

This is typically only specified in a restart header. If you specify it
in a start UAM.in it will start the counter to whatever you specify.

    #ISTEP
    iStep     (integer)



### TSIMULATION

This is typically only specified in a restart header. It sets the offset
from the starttime to the currenttime. Should really only be used with
caution.

    #TSIMULATION
    tsimulation    (real)

### DUSTDATA

    #DUST
    cDustFile
    cConrathFile

### OVERWRITEIONOSPHERE

    #OVERWRITEIONOSPHERE
    DoOverwriteIonosphere
    DoOverwriteWithIRI
    DoOverwriteWithSami
    SamiInFile



### DUST

    #DUST
    TauTot
    Conrnu

### DAMPING

This is probably for damping vertical wind oscillations that can occur
in the lower atmosphere.

    #DAMPING
    UseDamping        (logical)

### GRAVITYWAVE

    #GRAVITYWAVE
    UseGravityWave        (logical)


### NEWSTRETCH

    #NEWSTRETCH
    Poleward Edge of Stretch Region (real, degrees)
    StretchWidth  (real, 1.0-20.0, deg)
    StretchingPercentage  (real, 0-1)
    Example (auroral zone):
    #NEWSTRETCH
    65.0 ! location of minimum grid spacing
    5.0         ! Width of stretched region
    0.6         ! Amount of stretch 0 (none) to 1 (lots)

### STRETCH

    #STRETCH
    ConcentrationLatitude (real, degrees)
    StretchingPercentage  (real, 0-1)
    StretchingFactor      (real)
    Example (no stretching):
    #STRETCH
    65.0 ! location of minimum grid spacing
    0.0         ! Amount of stretch 0 (none) to 1 (lots)
    1.0  

### TOPOGRAPHY

    #TOPOGRAPHY
    UseTopography (logical)
    AltMinUniform (real)

### RCMR

    #RCMR
    Input data type to assimilate (RHO/VTEC)
    Output variable to drive (F107/PHOTOELECTRON)
    Initial output estimate
    Number satellites to assimilate
    1st satellite index to assimilate
    2nd satellite index to assimilate
    etc...

### DART

    #DART
    useDART (integer, {default 0=no}, 1=master ensemble member, 2=slave ens.)


### LTERadiation

    #LTERadiation
    DtLTERadiation (real)


