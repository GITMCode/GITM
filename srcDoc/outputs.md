# Outputs {#outputs.sec}

Now that you have managed to successfully complete a GITM run you've
found yourself with a bunch of output files. All of the GITM output is
in mks units and this data is contained within several files located in
the `UA/data` directory

After [postprocessing](postprocessing.md), you will find yourself with a
`log*.dat` file,and many `.bin` files in whichever formats you specified in
SAVEPLOT (see [`#SAVEPLOT`](set_inputs.md#saveplot)

The log file provides useful information about the run, such as whether a
restart was performed, which physical processes were used, and a list of some
variables which represent the state of the modeled system and its drivers. This
file can be very useful when sharing runs with other users, when revisiting an
old run, or merely ensuring that GITM performed as expected. Some example log
files can be found in the `srcTest/auto_test/ref_solns/` directory.

## Possible Output Variables

The output binary files can contain the following atmospheric
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


Currently the following output files are available:

## 2DANC

    1 Longitude
    2 Latitude
    3 Altitude
    4 Local Time
    5 Solar Zenith Angle
    6 Vertical TEC
    7 AltIntJouleHeating (W/m2)
    8 AltIntHeatingTransfer (W/m2)
    9 AltIntEuvHeating (W/m2)
    10 AltIntPhotoElectronHeating (W/m2)
    11 AltIntChamicalHeating (W/m2)
    12 AltIntRadCooling (W/m2)
    12 AltIntCO2Cooling (W/m2)
    12 AltIntNOCooling (W/m2)
    12 AltIntOCooling (W/m2)

## 2DGEL

    1 Longitude
    2 Latitude
    3 Altitude
    4 Potential
    5 Pedersen Conductance
    6 Hall Conductance
    7 Electron_Average_Energy
    8 Electron_Energy_Flux
    9 DivJuAlt
    10 Pedersen FL Conductance
    11 Hall FL Conductance
    12 DivJu FL
    13 FL Length

## 2DHME

    1 Longitude
    2 Latitude
    3 Altitude
    4 Local Time
    5 Vertical TEC

## 2DMEL

    1 Longitude
    2 Latitude
    3 Altitude
    4 MLT
    5 GeoLat
    6 GeoLon
    7 Pedersen Conductance
    8 Hall Conductance
    9 DivJuAlt
    10 Field Line Length
    11 Sigma PP
    12 Sigma LL
    13 Sigma H
    14 Sigma C
    15 Sigma PL
    16 Sigma LP
    17 K^D_{m\phi}
    18 K^D_{m\lamda}
    19 Solver A
    20 Solver B
    21 Solver C
    22 Solver D
    23 Solver E
    24 Solver S
    25 DynamoPotential
    26 Ed1new
    27 Ed2new
    28 Kphi
    29 Klamda

## 2DTEC

    1 Longitude
    2 Latitude
    3 Altitude
    4 Solar Zenith Angle
    5 Vertical TEC


## 2DUSR

    1 Longitude
    2 Latitude
    3 Altitude
    4 Potential (kV)
    5 Average Energy (keV)
    6 Total Energy (ergs)
    7 Discrete Average Energy (keV)
    8 Discrete Total Energy (ergs)
    9 Wave Average Energy (keV)
    10 Wave Total Energy (ergs)
    11 Flux@5.000E+05eV (/cm2/s)
    12 Flux@4.009E+05eV (/cm2/s)
    13 Flux@3.215E+05eV (/cm2/s)
    14 Flux@2.578E+05eV (/cm2/s)
    15 Flux@2.067E+05eV (/cm2/s)
    16 Flux@1.658E+05eV (/cm2/s)
    17 Flux@1.329E+05eV (/cm2/s)
    18 Flux@1.066E+05eV (/cm2/s)
    19 Flux@8.547E+04eV (/cm2/s)
    20 Flux@6.853E+04eV (/cm2/s)
    21 Flux@5.495E+04eV (/cm2/s)
    22 Flux@4.407E+04eV (/cm2/s)
    23 Flux@3.533E+04eV (/cm2/s)
    24 Flux@2.833E+04eV (/cm2/s)
    25 Flux@2.272E+04eV (/cm2/s)
    26 Flux@1.822E+04eV (/cm2/s)
    27 Flux@1.461E+04eV (/cm2/s)
    28 Flux@1.171E+04eV (/cm2/s)
    29 Flux@9.393E+03eV (/cm2/s)
    30 Flux@7.532E+03eV (/cm2/s)
    31 Flux@6.040E+03eV (/cm2/s)
    32 Flux@4.843E+03eV (/cm2/s)
    33 Flux@3.884E+03eV (/cm2/s)
    34 Flux@3.114E+03eV (/cm2/s)
    35 Flux@2.497E+03eV (/cm2/s)
    36 Flux@2.002E+03eV (/cm2/s)
    37 Flux@1.606E+03eV (/cm2/s)
    38 Flux@1.287E+03eV (/cm2/s)
    39 Flux@1.032E+03eV (/cm2/s)
    40 Flux@8.278E+02eV (/cm2/s)
    41 Flux@6.638E+02eV (/cm2/s)
    42 Flux@5.323E+02eV (/cm2/s)
    43 Flux@4.268E+02eV (/cm2/s)
    44 Flux@3.423E+02eV (/cm2/s)
    45 Flux@2.744E+02eV (/cm2/s)
    46 Flux@2.201E+02eV (/cm2/s)
    47 Flux@1.765E+02eV (/cm2/s)
    48 Flux@1.415E+02eV (/cm2/s)
    49 Flux@1.135E+02eV (/cm2/s)
    50 Flux@9.099E+01eV (/cm2/s)
    51 Flux@7.296E+01eV (/cm2/s)
    52 Flux@5.850E+01eV (/cm2/s)
    53 Flux@4.691E+01eV (/cm2/s)
    54 Flux@3.762E+01eV (/cm2/s)
    55 Flux@3.016E+01eV (/cm2/s)
    56 Flux@2.419E+01eV (/cm2/s)
    57 Flux@1.940E+01eV (/cm2/s)
    58 Flux@1.555E+01eV (/cm2/s)
    59 Flux@1.247E+01eV (/cm2/s)
    60 Flux@1.000E+01eV (/cm2/s)

## 3DALL

    1 Longitude
    2 Latitude
    3 Altitude
    4 Rho
    5 [O(!U3!NP)           ]
    6 [O!D2!N              ]
    7 [N!D2!N              ]
    8 [N(!U4!NS)           ]
    9 [NO                  ]
    10 [He                  ]
    11 [N(!U2!ND)           ]
    12 [N(!U2!NP)           ]
    13 [H                   ]
    14 [CO!D2!N             ]
    15 [O(!U1!ND)           ]
    16 Temperature
    17 V!Dn!N (east)
    18 V!Dn!N (north)
    19 V!Dn!N (up)
    20 V!Dn!N (up,O(!U3!NP)           )
    21 V!Dn!N (up,O!D2!N              )
    22 V!Dn!N (up,N!D2!N              )
    23 V!Dn!N (up,N(!U4!NS)           )
    24 V!Dn!N (up,NO                  )
    25 V!Dn!N (up,He                  )
    26 [O_4SP_!U+!N         ]
    27 [NO!U+!N             ]
    28 [O!D2!U+!N           ]
    29 [N!D2!U+!N           ]
    30 [N!U+!N              ]
    31 [O(!U2!ND)!U+!N      ]
    32 [O(!U2!NP)!U+!N      ]
    33 [H!U+!N              ]
    34 [He!U+!N             ]
    35 [e-                  ]
    36 eTemperature
    37 iTemperature
    38 V!Di!N (east)
    39 V!Di!N (north)
    40 V!Di!N (up)

## 3DCHM

    1 Longitude
    2 Latitude
    3 Altitude
    4 N!D2!U+!N + e
    5 O!D2!U+!N + e
    6 N!D2!U+!N + O
    7 NO!U+!N + e
    8 N!U+!N + O!D2!N
    9 NO + N
    10 O!U+!N + O!D2!N
    11 N + O!D2!N
    12 O!D2!U+!N + N
    13 O!D2!U+!N + NO
    14 O!D2!U+!N + N2
    15 N!D2!U+!N + O!D2!N
    16 N!U+!N + O
    17 O!+!N + N!D2!N
    18 O(1D) + N!D2!N
    19 O(1D) + O!D2!N
    20 O(1D) + O
    21 O(1D) + e
    22 N(2D) + O!D2!N
    23 O!U+!N(2D)+e
    24 N(2D) + O
    25 N(2D) + e
    26 O!U+!N(2D + N!D2!N
    27 O!U+!N(2P) + e
    28 O!U+!N(2P) + O
    29 O!U+!N(2P) + N!D2!N
    30 Chemical Heating Rate

## 3DGLO

    1 Longitude
    2 Latitude
    3 Altitude
    4 6300 A Emission
    5 PhotoElectronUp
    6 PhotoElectronDown

## 3DHME

    1 Longitude
    2 Latitude
    3 Altitude
    4 Rho
    5 [O(!U3!NP)           ]
    6 [O!D2!N              ]
    7 [N!D2!N              ]
    8 [N(!U4!NS)           ]
    9 [NO                  ]
    10 [He                  ]
    11 [N(!U2!ND)           ]
    12 [N(!U2!NP)           ]
    13 [H                   ]
    14 [CO!D2!N             ]
    15 [O(!U1!ND)           ]
    16 Temperature
    17 V!Dn!N (east)
    18 V!Dn!N (north)
    19 V!Dn!N (up)
    20 V!Dn!N (up,O(!U3!NP)           )
    21 V!Dn!N (up,O!D2!N              )
    22 V!Dn!N (up,N!D2!N              )
    23 V!Dn!N (up,N(!U4!NS)           )
    24 V!Dn!N (up,NO                  )
    25 V!Dn!N (up,He                  )
    26 [O_4SP_!U+!N         ]
    27 [NO!U+!N             ]
    28 [O!D2!U+!N           ]
    29 [N!D2!U+!N           ]
    30 [N!U+!N              ]
    31 [O(!U2!ND)!U+!N      ]
    32 [O(!U2!NP)!U+!N      ]
    33 [H!U+!N              ]
    34 [He!U+!N             ]
    35 [e-                  ]
    36 eTemperature
    37 iTemperature
    38 V!Di!N (east)
    39 V!Di!N (north)
    40 V!Di!N (up)
    41 PhotoElectron Heating
    42 Joule Heating
    43 Auroral Heating
    44 Specific Heat
    45 Magnetic Latitude
    46 Magnetic Longitude
    47 B.F. East
    48 B.F. North
    49 B.F. Vertical
    50 B.F. Magnitude
    51 Potential
    52 PotentialY
    53 E.F. East
    54 E.F. North
    55 E.F. Vertical

## 3DION

    1 Longitude
    2 Latitude
    3 Altitude
    4 [O_4SP_!U+!N         ]
    5 [NO!U+!N             ]
    6 [O!D2!U+!N           ]
    7 [N!D2!U+!N           ]
    8 [N!U+!N              ]
    9 [O(!U2!ND)!U+!N      ]
    10 [O(!U2!NP)!U+!N      ]
    11 [H!U+!N              ]
    12 [He!U+!N             ]
    13 [e-                  ]
    14 eTemperature
    15 iTemperature
    16 V!Di!N (east)
    17 V!Di!N (north)
    18 V!Di!N (up)
    19 Ed1
    20 Ed2
    21 Je1
    22 Je2
    23 Magnetic Latitude
    24 Magnetic Longitude
    25 B.F. East
    26 B.F. North
    27 B.F. Vertical
    28 B.F. Magnitude
    29 Potential
    30 E.F. East
    31 E.F. North
    32 E.F. Vertical
    33 E.F. Magnitude
    34 IN Collision Freq
    35 PressGrad (east)
    36 PressGrad (north)
    37 PressGrad (up)

## 3DLST

    1 Longitude
    2 Latitude
    3 Altitude
    4 Rho (kg/m3)
    5 [O(!U3!NP)           ] (/m3)
    6 [O!D2!N              ] (/m3)
    7 [N!D2!N              ] (/m3)
    8 [NO                  ] (/m3)
    9 Vn (east) (m/s)
    10 Vn (north) (m/s)
    11 Vn (up) (m/s)
    12 [O_4SP_!U+!N         ] (/m3)
    13 [NO!U+!N             ] (/m3)
    14 [O!D2!U+!N           ] (/m3)
    15 [N!D2!U+!N           ] (/m3)
    16 [e-                  ] (/m3)
    17 Vi (east) (m/s)
    18 Vi (north) (m/s)
    19 Vi (up) (m/s)
    20 Neutral Temperature (K)

## 3DMAG

    1 Longitude
    2 Latitude
    3 Altitude
    4 Magnetic Latitude
    5 Magnetic Longitude
    6 B.F. East
    7 B.F. North
    8 B.F. Vertical
    9 B.F. Magnitude

## 3DNEU

    1 Longitude
    2 Latitude
    3 Altitude
    4 Rho
    5 [O(!U3!NP)           ]
    6 [O!D2!N              ]
    7 [N!D2!N              ]
    8 [N(!U4!NS)           ]
    9 [NO                  ]
    10 [He                  ]
    11 [N(!U2!ND)           ]
    12 [N(!U2!NP)           ]
    13 [H                   ]
    14 [CO!D2!N             ]
    15 [O(!U1!ND)           ]
    16 Temperature
    17 V!Dn!N (east)
    18 V!Dn!N (north)
    19 V!Dn!N (up)
    20 V!Dn!N (up,O(!U3!NP)           )
    21 V!Dn!N (up,O!D2!N              )
    22 V!Dn!N (up,N!D2!N              )
    23 V!Dn!N (up,N(!U4!NS)           )
    24 V!Dn!N (up,NO                  )
    25 V!Dn!N (up,He                  )

## 3DTHM

    1 Longitude
    2 Latitude
    3 Altitude
    4 EUV Heating
    5 Conduction
    6 Molecular Conduction
    7 Eddy Conduction
    8 Eddy Adiabatic Conduction
    9 Chemical Heating
    10 Auroral Heating
    11 Joule Heating
    12 NO Cooling
    13 O Cooling
    14 Total Abs EUV
    15 Cp
    16 Rho
    17 E-Field Mag
    18 Sigma Ped

## 3DUSR

    1 Longitude
    2 Latitude
    3 Altitude
    4 Joule Heating
    5 JPara