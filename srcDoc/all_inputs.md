# All Inputs

This sets the starting time of the simulation. Even when you restart,
the starttime should be to the real start time, not the restart time.

    #STARTTIME
    iYear    (integer)
    iMonth   (integer)
    iDay     (integer)
    iHour    (integer)
    iMinute  (integer)
    iSecond  (integer)

This sets the ending time of the simulation.

    #ENDTIME
    iYear    (integer)
    iMonth   (integer)
    iDay     (integer)
    iHour    (integer)
    iMinute  (integer)
    iSecond  (integer)

This will set a time for the code to just pause. Really, this should
never be used.

    #PAUSETIME
    iYear iMonth iDay iHour iMinute iSecond

This is typically only specified in a restart header. If you specify it
in a start UAM.in it will start the counter to whatever you specify.

    #ISTEP
    iStep     (integer)

This sets the maximum CPU time that the code should run before it starts
to write a restart file and end the simulation. It is very useful on
systems that have a queueing system and has limited time runs.
Typically, set it for a couple of minutes short of the max wall clock,
since it needs some time to write the restart files.

    #CPUTIMEMAX
    CPUTimeMax    (real)

This command will skip all pretty much all of the physics of GITM, and
will reinitialize the model with the MSIS and IRI values at the interval
set in the second variable. If you want to compare a run to MSIS and
IRI, you can run GITM with this command and get output at exactly the
same cadence and locations, thereby allowing easier comparisons. The dt
can be set as low as needed, so you can run satellites through MSIS and
IRI.

    #STATISTICALMODELSONLY
    UseStatisticalModelsOnly    (logical)
    DtStatisticalModels         (real)

This is typically only specified in a restart header. It sets the offset
from the starttime to the currenttime. Should really only be used with
caution.

    #TSIMULATION
    tsimulation    (real)

Sets the F10.7 and 81 day average F10.7. This is used to set the initial
altitude grid, and drive the lower boundary conditions.

    #F107
    f107  (real)
    f107a (real - 81 day average of f107)

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

This says how to use tides. The first one is using MSIS with no tides.
The second is using MSIS with full up tides. The third is using GSWM
tides, while the forth is for experimentation with using WACCM tides.

    #TIDES
    UseMSISOnly        (logical)
    UseMSISTides       (logical)
    UseGSWMTides       (logical)
    UseWACCMTides      (logical)
    UseHmeTides        (logical)

This says how to use msis tides. The first one is using diurnal tide The
first one is using semi-diurnal tide The first one is using terdiurnal
tide

    #MSISTIDES
    UseMSISDiurnal        (logical)
    UseMSISSemidiurnal    (logical)
    UseMSISTerdiurnal     (logical)

UseOBCExperiment - use MSIS \[O\] BC shifted by 6 months Only applicable
for MSIS00! MsisOblateFactor - alt = alt \* (1.0-f/2 + f\*cos(lat)) -
seems like -0.1 works well

    #MSISOBC
    UseOBCExperiment        (logical)
    MsisOblateFactor           (real)

This toggles between using MSIS00 (false) and MSIS-2.1 (true)

    #MSISOBC
    UseMsis21       (logical)

This says whether you want seconds in output file name. F means no
seconds in output file name.

    #DUST
    cDustFile
    cConrathFile

    #OVERWRITEIONOSPHERE
    DoOverwriteIonosphere
    DoOverwriteWithIRI
    DoOverwriteWithSami
    SamiInFile

    #GITMBCS
    UseGitmBCs
    GitmBCsDir

    #DUST
    TauTot
    Conrnu

If you selected to use GSWM tides above, you can specify which
components to use.

    #GSWMCOMP
    GSWMdiurnal(1)        (logical)
    GSWMdiurnal(2)        (logical)
    GSWMsemidiurnal(1)    (logical)
    GSWMsemidiurnal(2)    (logical)

    #USEPERTURBATION
    UsePerturbation        (logical)

    #USEBCPERTURBATION
    UseBcPerturbation        (logical)
    If UseBcPerturbation = .true. then:
    iTypeBcPerturb         (int) 
    perturbation characteristics ...

This is probably for damping vertical wind oscillations that can occur
in the lower atmosphere.

    #DAMPING
    UseDamping        (logical)

I dont know what this is for\...

    #GRAVITYWAVE
    UseGravityWave        (logical)

This sets the hemispheric power of the aurora. Typical it ranges from
1-1000, although 20 is a nominal, quiet time value.

    #HPI
    HemisphericPower  (real)

I dont think that GITM actually uses this unless the Foster electric
field model is used.

    #KP
    kp  (real)

The CFL condition sets how close to the maximum time step that GITM will
take. 1.0 is the maximum value. A value of about 0.75 is typical. If
instabilities form, a lower value is probably needed.

    #CFL
    cfl  (real)

If you would like to force GITM to take a fixed dt you can use this. It
will try to take that fixed dt, unless the CFL condition is violated.

    #FIXEDDT
    FixedDt  (real)

This sets the driving conditions for the high-latitude electric field
models. This is static for the whole run, though. It is better to use
the MHD_INDICES command to have dynamic driving conditions.

    #SOLARWIND
    bx  (real)
    by  (real)
    bz  (real)
    vx  (real)

Use this for dynamic IMF and solar wind conditions. The exact format of
the file is discussed further in the manual.

    #MHD_INDICES
    filename  (string)

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

This is for using Pat Newells aurora (Ovation).

    #NEWELLAURORA
    UseNewellAurora   (logical)
    UseNewellAveraged (logical)
    UseNewellMono (logical)
    UseNewellWave (logical)
    UseNewellRemoveSpikes (logical)
    UseNewellAverage      (logical)

This is for using Betsy Michells aurora (OvationSME).

    #OVATIONSME
    UseOvationSME     (logical)
    UseOvationSMEMono (logical)
    UseOvationSMEWave (logical)
    UseOvationSMEIon  (logical)

This is for using Dongjies aurora.

    #AEMODEL
    UseAeModel        (logical)

This is for using the FTA Model of the aurora.

    #FTAMODEL
    UseFtaModel        (logical)

This is for using Dongjies aurora.

    #FANGENERGY
    UseFangEnergyDeposition        (logical)

This is for specifying a cusp.

    #USECUSP
    UseCusp        (logical)
    CuspAveE       (real)
    CuspEFlux      (real)

    #AMIEFILES
    cAMIEFileNorth  (string)
    cAMIEFileSouth  (string)

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

The limiter is quite important. It is a value between 1.0 and 2.0, with
1.0 being more diffuse and robust, and 2.0 being less diffuse, but less
robust.

    #LIMITER
    TypeLimiter  (string)
    BetaLimiter  (real between 1.0-minmod and 2.0-mc)

This will turn on all of the NaN checks in the code!

    #NANCHECK
    DoCheckForNans (logical)

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

    #IONLIMITS
    MaxVParallel     (real, default=100 m/s)
    MaxEField        (real, default=0.1 V/m)
    MinIonDensity    (real, default=100 m^-3)

    #PHOTOELECTRON
    PhotoElectronHeatingEfficiency   (real)

    #NEUTRALHEATING
    NeutralHeatingEfficiency   (real)

    #DON4SHACK
    DoN4SHack       (logical)

    #THERMO
    UseSolarHeating   (logical)
    UseJouleHeating   (logical)
    UseAuroralHeating (logical)
    UseNOCooling      (logical)
    UseOCooling       (logical)
    UseConduction     (logical)
    UseTurbulentCond  (logical)
    UseIRHeating      (logical)

    #THERMALDIFFUSION
    KappaTemp0    (thermal conductivity, real)

    #THERMALCONDUCTION
    ThermalConduction_AO2 (Conduction A(O2): 3.6e-4, real)
    ThermalConduction_AO  (Conduction A(O): 5.6e-4, real)
    ThermalConduction_s   (Conduction s: 0.75, real)

    #VERTICALSOURCES
    MaximumVerticalVelocity      (real)

    #AUSMSOLVER
    Use AUSM Solver      (logical)

If you use eddy diffusion, you must specify two pressure levels - under
the first, the eddy diffusion is constant. Between the first and the
second, there is a linear drop-off. Therefore The first pressure must be
larger than the second!

    #DIFFUSION
    UseDiffusion (logical)
    EddyDiffusionCoef (real)
    EddyDiffusionPressure0 (real)
    EddyDiffusionPressure1 (real)

    #FORCING
    UsePressureGradient (logical)
    UseIonDrag          (logical)
    UseNeutralFriction  (logical)
    UseViscosity        (logical)
    UseCoriolis         (logical)
    UseGravity          (logical)

    #MODIFYPLANET
    RotationPeriodInput        (real)
    DaysPerYearInput           (real)
    DaysPerYearInput           (real)

    #USEIMPLICITIONMOMENTUM
    UseImplicitFieldAlignedMomentum      (logical)

    #USEIMPROVEDIONADVECTION
    UseImprovedIonAdvection      (logical)
    UseNighttimeIonBCs           (logical)
    MinTEC                       (real)

    #USETESTVISCOSITY
    TestViscosityFactor      (real)

    #DYNAMO
    UseDynamo              (logical)
    DynamoHighLatBoundary  (real)
    nItersMax              (integer)
    MaxResidual            (V,real)
    IncludeCowling         (logical)
    DynamoLonAverage       (real)

    #IONFORCING
    UseExB                 (logical)
    UseIonPressureGradient (logical)
    UseIonGravity          (logical)
    UseNeutralDrag         (logical)

    #FIXTILT
    IsFixedTilt   (logical)

    #DIPOLE
    MagneticPoleRotation   (real)
    MagneticPoleTilt       (real)
    xDipoleCenter          (real)
    yDipoleCenter          (real)
    zDipoleCenter          (real)

    #APEX
    UseApex (logical)
            Sets whether to use a realistic magnetic
            field (T) or a dipole (F)

For Earth, the AltMin is the only variable used here. The altitudes are
set to 0.3 times the scale height reported by MSIS, at the equator for
the specified F107 and F107a values.

    #ALTITUDE
    AltMin                (real, km)
    AltMax                (real, km)
    UseStretchedAltitude  (logical)

If LatStart and LatEnd are set to \< -90 and \> 90, respectively, then
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

    #NEWSTRETCH
    Poleward Edge of Stretch Region (real, degrees)
    StretchWidth  (real, 1.0-20.0, deg)
    StretchingPercentage  (real, 0-1)
    Example (auroral zone):
    #NEWSTRETCH
    65.0 ! location of minimum grid spacing
    5.0         ! Width of stretched region
    0.6         ! Amount of stretch 0 (none) to 1 (lots)

    #STRETCH
    ConcentrationLatitude (real, degrees)
    StretchingPercentage  (real, 0-1)
    StretchingFactor      (real)
    Example (no stretching):
    #STRETCH
    65.0 ! location of minimum grid spacing
    0.0         ! Amount of stretch 0 (none) to 1 (lots)
    1.0  

    #TOPOGRAPHY
    UseTopography (logical)
    AltMinUniform (real)

    #RESTART
    DoRestart (logical)

This allows you to change the output cadence of the files for a limited
time. If you have an event then you can output much more often during
that event.

    #PLOTTIMECHANGE
    yyyy mm dd hh mm ss ms (start)
    yyyy mm dd hh mm ss ms (end)

For satellite files, you can have one single file per satellite, instead
of one for every output. This makes GITM output significantly less
files. It only works for satellite files now.

    #APPENDFILES
    DoAppendFiles    (logical)

Typicaly file is named (e.g.) 1DALL_yymmdd_hhmmss.bin With this it will
be named 1DALL_GITM_yyyy-mm-ddThh-mm-ss.bin

    #CCMCFILENAME
    UseCCMCFileName    (logical)

The DtRestart variable sets the time in between writing full restart
files to the UA/restartOUT directory.\
This sets the output files. The most common type is 3DALL, which outputs
all primary state variables. Types include : 3DALL, 3DNEU, 3DION, 3DTHM,
3DCHM, 3DUSR, 3DGLO, 2DGEL, 2DMEL, 2DUSR, 1DALL, 1DGLO, 1DTHM, 1DNEW,
1DCHM, 1DCMS, 1DUSR. DtPlot sets the frequency of output

    #SAVEPLOT
    DtRestart (real, seconds)
    nOutputTypes  (integer)
    Outputtype (string, 3D, 2D, ION, NEUTRAL, ...)
    DtPlot    (real, seconds)

    #SAVEHIMEPLOT
    HIMEPlotLonStart (real)
    HIMEPlotLonEnd (real)
    HIMEPlotLatStart (real)
    HIMEPlotLatEnd (real)

    #SATELLITES
    nSats     (integer - max = ', nMaxSats, ')
    SatFile1  (string)
    SatOutputtype1  (string, 0DUSR or 1DUSR or other)
    DtPlot1   (real, seconds)
    etc...

    #RCMR
    Input data type to assimilate (RHO/VTEC)
    Output variable to drive (F107/PHOTOELECTRON)
    Initial output estimate
    Number satellites to assimilate
    1st satellite index to assimilate
    2nd satellite index to assimilate
    etc...

    #DART
    useDART (integer, {default 0=no}, 1=master ensemble member, 2=slave ens.)

Sets the time for updating the high-latitude (and low-latitude)
electrodynamic drivers, such as the potential and the aurora.

    #ELECTRODYNAMICS
    DtPotential (real, seconds)
    DtAurora    (real, seconds)

Sets the time delay for the high latitude drivers and solar EUV inputs.

    #INPUTTIMEDELAY
    TimeDelayHighLat (real, seconds)
    TimeDelayEUV     (real, seconds)

    #LTERadiation
    DtLTERadiation (real)

You can only have an AMIE input file for this now. Make sure you put the
ions in the AMIE file!!!

    #IONPRECIPITATION
    UseIonPrecipitation     (logical)

You really want a log file. They are very important. It is output in
UA/data/log\*.dat. You can output the log file at whatever frequency you
would like, but if you set dt to some very small value, you will get an
output every iteration, which is probably a good thing.

    #LOGFILE
    DtLogFile   (real, seconds)

If EUVAC = true, use EUVAC Model If Tobiska = true, use Tobiska91 Model
If both are true, average them together If UseAboveHigh, then extend the
spectrum to the longer wavelengths If UseBelowLow, then extend the
spectrum to the shorter wavelengths

    #EUV_DATA
    UseEUVData            (logical)
    cEUVFile              (string)

This is for a FISM or some other solar spectrum file.

    #EUV_DATA
    UseEUVData            (logical)
    cEUVFile              (string)

    #ECLIPSE
    start y m d h m s            (int x 6)
    end   y m d h m s            (int x 6)
    start-y-gse                  (real)
    start-z-gse                  (real)
    end-y-gse                    (real)
    end-z-gse                    (real)
    peak                         (real)
    Max Distance                 (real)
    ExpAmp                       (real)
    ExpWidth                     (real)
