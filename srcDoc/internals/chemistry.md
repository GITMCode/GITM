# Chemistry

Chemical reactions in GITM for Earth.

The source of truth is `src/calc_chemistry.Earth.f90`;
this page summarizes what is active in that routine. Conventions:

- Two-body rate coefficients $`k`$ are in $`\mathrm{m^3\,s^{-1}}`$ (the code
  works in SI; values quoted in $`\mathrm{cm^3\,s^{-1}}`$ in the literature are
  multiplied by $`10^{-6}`$). Radiative (first-order) rates are in
  $`\mathrm{s^{-1}}`$. The one three-body rate is in $`\mathrm{m^6\,s^{-1}}`$.
- $`T_n`$, $`T_i`$, $`T_e`$ are the neutral, ion and electron temperatures in K.
  Temperatures enter the rates as ratios, for example
  $`(T_e/300)^{-0.5}`$.
- The rate of a reaction $`A + B \rightarrow \ldots`$ is $`k\,n_A\,n_B`$ and
  $`k\,n_A`$ for a first-order loss. Photo-ionization, photo-dissociation and
  auroral rates come from the EUV and aurora modules; the one exception is NO
  photodissociation, which this routine computes itself (see below).
- The "Heating" column is the exothermic energy (eV per reaction) added to the
  chemical heating rate. The code splits each energy into a "neutral" and an
  "ion" part, but both are added to the neutral heating rate: the ion and
  electron chemical heating rates are always zero. A dash means the reaction
  adds no heating.
- Each species is advanced once per model time step, with the sources
  explicit and the losses implicit: $`n^{new} = (n + S\,\Delta t)/(1 + L\,\Delta t/n)`$.
  There is no sub-cycling.
- `#CHEMISTRY` switches: with `UseIonChemistry` false the whole routine is
  skipped, neutral chemistry included. With `UseNeutralChemistry` false the
  neutral sources and losses are zeroed while the ion chemistry still runs.

A temperature-dependent ion-drift "effective temperature" is computed at the
top of the routine, but it is not used by any rate below.

---

## Neutral Chemistry

### Three-body recombination

```math
\begin{align}
O + O + M \rightarrow O_2 + M, \qquad M = O_2,\ O,\ N_2
\end{align}
```

```math
\begin{align}
k = 4.7\times 10^{-45} \left(\frac{300}{T_n}\right)^{2} \ \mathrm{m^6\,s^{-1}},
\qquad R = k\,[O]^2\,\left([O_2] + [O] + [N_2]\right)
\end{align}
```

Each reaction consumes two O and produces one $`O_2`$; heating is 5.12 eV.

### Reactions of N(4S)

| Reaction | $`k`$ (m³ s⁻¹) | Heating (eV) |
|---|---|---|
| $`N(^4S) + O_2 \rightarrow NO + O`$ | $`1.5\times10^{-20}\,T_n\,e^{-3270/T_n}`$ | 1.385 |
| $`N(^4S) + NO \rightarrow N_2 + O`$ | $`3.4\times10^{-17}`$ | 3.25 |

### Reactions of N(2D)

| Reaction | $`k`$ (m³ s⁻¹) | Heating (eV) |
|---|---|---|
| $`N(^2D) + O_2 \rightarrow NO + O`$ | $`9.7\times10^{-18}\,e^{-185/T_n}`$ | see below |
| $`N(^2D) + NO \rightarrow N_2 + O`$ | $`6.7\times10^{-17}`$ | 5.63 |
| $`N(^2D) + O \rightarrow N(^4S) + O`$ | $`6.9\times10^{-19}`$ | see below |
| $`N(^2D) + e^- \rightarrow N(^4S) + e^-`$ | $`3.86\times10^{-16}\,(T_e/300)^{0.81}`$ | 2.38 |
| $`N(^2D) \rightarrow N(^4S) + h\nu\ (5200\,\mathrm{\AA})`$ | $`1.06\times10^{-5}\ \mathrm{s^{-1}}`$ | -- |

- $`N(^2D) + O_2`$: the oxygen product is split 90% $`O(^3P)`$ and 10%
  $`O(^1D)`$, with heating $`0.9\times3.76 + 0.1\times1.80`$ eV. (The two
  channels release 3.76 and 1.80 eV.)
- $`N(^2D) + O`$: $`N(^4S)`$ is produced; $`O`$ is returned as 90%
  $`O(^3P)`$ (2.38 eV) and 10% $`O(^1D)`$ (0.42 eV), so the heating is
  $`0.9\times2.38 + 0.1\times0.42`$ eV.
- Both $`O(^1D)`$ branches apply when $`O(^1D)`$ is a species,
  which it always is in this routine (see the $`O(^1D)`$ section). If
  it were switched off, $`N(^2D)+O_2`$ would use only 3.76 eV and
  $`N(^2D)+O`$ only 2.38 eV.

### Reactions of N(2P)

| Reaction | $`k`$ (m³ s⁻¹) | Heating (eV) |
|---|---|---|
| $`N(^2P) \rightarrow N(^2D) + h\nu\ (10400\,\mathrm{\AA})`$ | $`7.9\times10^{-2}\ \mathrm{s^{-1}}`$ | -- |
| $`N(^2P) + e^- \rightarrow N(^2D) + e^-`$ | $`9.5\times10^{-15}`$ | 1.19 |
| $`N(^2P) + e^- \rightarrow N(^4S) + e^-`$ | $`2.04\times10^{-16}\,(T_e/300)^{0.85}`$ | 3.57 |
| $`N(^2P) + NO \rightarrow N(^4S) + NO`$ | $`1.8\times10^{-16}`$ | 3.44 |
| $`N(^2P) + O \rightarrow N(^2D) + O`$ | $`1.7\times10^{-17}`$ | 1.19 |
| $`N(^2P) + O_2 \rightarrow NO + O`$ | $`3.09\times10^{-18}\,e^{-60/T_n}`$ | 4.95 |

### Reactions of O(1D)

These are active only if $`O(^1D)`$ is used as a species (`UseNeutralConstituent(iO_1D_)`,
which is set true in the routine). All of them return $`O(^3P)`$.

| Reaction | $`k`$ | Heating (eV) |
|---|---|---|
| $`O(^1D) \rightarrow O(^3P) + h\nu\ (6300\,\mathrm{\AA})`$ | $`7.1\times10^{-3}\ \mathrm{s^{-1}}`$ | -- |
| $`O(^1D) \rightarrow O(^3P) + h\nu\ (6364\,\mathrm{\AA})`$ | $`2.2\times10^{-3}\ \mathrm{s^{-1}}`$ | -- |
| $`O(^1D) + e^- \rightarrow O(^3P) + e^-`$ | $`2.87\times10^{-16}\,(T_e/300)^{0.91}`$ | 1.96 |
| $`O(^1D) + N_2 \rightarrow O(^3P) + N_2`$ | $`1.8\times10^{-17}\,e^{107/T_n}`$ | 1.96 |
| $`O(^1D) + O_2 \rightarrow O(^3P) + O_2`$ | $`3.2\times10^{-17}\,e^{67/T_n}`$ | 1.96 |
| $`O(^1D) + O \rightarrow O(^3P) + O`$ | $`6.47\times10^{-18}\,(T_n/300)^{0.14}`$ | 1.96 |

### NO photodissociation

```math
\begin{align}
NO + h\nu \rightarrow N(^4S) + O
\end{align}
```

```math
\begin{align}
J_{NO} = 4.5\times10^{-6}\left[1 + 0.11\,\frac{F_{10.7} - 65}{165}\right]
\exp\!\left[-10^{-8}\left(10^{-4}\,\mathcal{C}_{O_2}\right)^{0.38}\right]\ \mathrm{s^{-1}}
\end{align}
```

Here $`\mathcal{C}_{O_2}`$ is the $`O_2`$ Chapman (column) integral. $`F_{10.7}`$
is the model F10.7, replaced by its B10.7 equivalent when `#B107` is in use. The rate is set
to zero where $`\mathcal{C}_{O_2} \ge 0.5\,\mathcal{C}_{shadow}`$, with
$`\mathcal{C}_{shadow} = 10^{26}`$ (i.e. in the planetary shadow; without this
cutoff there would be no solar zenith angle dependence). No heating is applied.

---

## Photo-dissociation and Photo-ionization / Auroral Sources

These terms use rates computed elsewhere (`EuvDissRateS`, `EuvIonRateS`,
`AuroralIonRateS`, `IonPrecipIonRateS`). Each ionization removes one neutral
and adds one ion.

### EUV dissociation

```math
\begin{align}
O_2 + h\nu \rightarrow 2\,O(^3P)
\end{align}
```

```math
\begin{align}
N_2 + h\nu \rightarrow 0.25\,N(^4S) + 0.60\,N(^2D)
\end{align}
```

The $`N_2`$ dissociation rate is multiplied by 0.25 into $`N(^4S)`$ and 0.60
into $`N(^2D)`$, as coded (one $`N_2`$ is lost per event).

### EUV ionization

| Reaction | Ion source |
|---|---|
| $`N_2 + h\nu \rightarrow N_2^+ + e^-`$ | $`N_2^+`$ |
| $`O_2 + h\nu \rightarrow O_2^+ + e^-`$ | $`O_2^+`$ |
| $`O(^3P) + h\nu \rightarrow O^+(^4S) + e^-`$ | $`O^+(^4S)`$ |
| $`O(^3P) + h\nu \rightarrow O^+(^2D) + e^-`$ | $`O^+(^2D)`$ |
| $`O(^3P) + h\nu \rightarrow O^+(^2P) + e^-`$ | $`O^+(^2P)`$ |
| $`N(^4S) + h\nu \rightarrow N^+ + e^-`$ | $`N^+`$ |
| $`He + h\nu \rightarrow He^+ + e^-`$ | $`He^+`$ |
| $`NO + h\nu \rightarrow NO^+ + e^-`$ | $`NO^+`$ |

### Auroral and ion-precipitation ionization

The auroral ionization rate plus the ion-precipitation rate is applied for
$`N_2`$ and $`O_2`$:

```math
\begin{align}
N_2 + e^*_{aur} \rightarrow N_2^+ + 2e^-, \qquad O_2 + e^*_{aur} \rightarrow O_2^+ + 2e^-
\end{align}
```

For atomic oxygen, the auroral plus ion-precipitation ionization rate of
$`O(^3P)`$ is partitioned 0.4 / 0.4 / 0.2 into $`O^+(^4S)`$, $`O^+(^2D)`$ and
$`O^+(^2P)`$ respectively:

```math
\begin{align}
O(^3P) + e^*_{aur} \rightarrow 0.4\,O^+(^4S) + 0.4\,O^+(^2D) + 0.2\,O^+(^2P) + 2e^-
\end{align}
```

Helium is ionized by the auroral rate only (the ion-precipitation rate is not
included for He):

```math
\begin{align}
He + e^*_{aur} \rightarrow He^+ + 2e^-
\end{align}
```

---

## Ion–Neutral Reactions

Temperature dependences use $`T_i`$ (the ion temperature). "Heating" is the
exothermic energy applied to the chemical heating rate.

### N$`_2^+`$ as reactant

| Reaction | $`k`$ (m³ s⁻¹) | Heating (eV) |
|---|---|---|
| $`N_2^+ + O_2 \rightarrow O_2^+ + N_2`$ | $`5.1\times10^{-17}\,(T_i/300)^{-1.16}`$ for $`T_i \le 1000`$ K<br>$`1.26\times10^{-17}\,(T_i/1000)^{0.67}`$ otherwise | 3.53 |
| $`N_2^+ + O \rightarrow NO^+ + N(^2D)`$ | $`1.33\times10^{-16}\,(T_i/300)^{-0.44}`$ for $`T_i \le 1500`$ K<br>$`6.55\times10^{-17}\,(T_i/1500)^{0.2}`$ otherwise | 0.70 |
| $`N_2^+ + O \rightarrow O^+(^4S) + N_2`$ | $`7.0\times10^{-18}\,(T_i/300)^{-0.23}`$ | 1.96 |
| $`N_2^+ + N(^4S) \rightarrow N_2 + N^+`$ | $`1.0\times10^{-17}`$ | 2.48 |
| $`N_2^+ + NO \rightarrow NO^+ + N_2`$ | $`3.6\times10^{-16}`$ | 6.33 |

The $`N_2^+ + O \rightarrow NO^+ + N`$ channel produces $`N(^2D)`$ only.

### O$`_2^+`$ as reactant

| Reaction | $`k`$ (m³ s⁻¹) | Heating (eV) |
|---|---|---|
| $`O_2^+ + N(^4S) \rightarrow NO^+ + O`$ | $`1.0\times10^{-16}`$ | 4.21 |
| $`O_2^+ + N(^2D) \rightarrow NO^+ + O`$ | $`1.8\times10^{-16}`$ | 6.519 |
| $`O_2^+ + N(^2D) \rightarrow N^+ + O_2`$ | $`8.65\times10^{-17}`$ | 0 |
| $`O_2^+ + N(^2P) \rightarrow O_2^+ + N(^4S)`$ | $`2.2\times10^{-17}`$ | 3.565 |
| $`O_2^+ + NO \rightarrow NO^+ + O_2`$ | $`4.5\times10^{-16}`$ | 2.813 |
| $`O_2^+ + N_2 \rightarrow NO^+ + NO`$ | $`5.0\times10^{-22}`$ | 0.9333 |

In the last reaction the code adds one NO and removes one $`N_2`$ per event.

### O$`^+(^4S)`$ as reactant

| Reaction | $`k`$ (m³ s⁻¹) | Heating (eV) |
|---|---|---|
| $`O^+(^4S) + O_2 \rightarrow O_2^+ + O`$ | $`1.6\times10^{-17}\,(T_i/300)^{-0.52}`$ for $`T_i \le 900`$ K<br>$`9\times10^{-18}\,(T_i/900)^{0.92}`$ otherwise | 1.55 |
| $`O^+(^4S) + N_2 \rightarrow NO^+ + N(^4S)`$ | $`1.2\times10^{-18}\,(T_i/300)^{-0.45}`$ for $`T_i \le 1000`$ K<br>$`7.0\times10^{-19}\,(T_i/1000)^{2.12}`$ otherwise | 1.10 |
| $`O^+(^4S) + NO \rightarrow NO^+ + O`$ | $`7.0\times10^{-19}\,(T_i/300)^{0.87}`$ | 4.36 |
| $`O^+(^4S) + N(^2D) \rightarrow N^+ + O`$ | $`1.3\times10^{-16}`$ | 1.45 |

### O$`^+(^2D)`$ as reactant

| Reaction | $`k`$ (m³ s⁻¹) | Heating (eV) |
|---|---|---|
| $`O^+(^2D) + N_2 \rightarrow N_2^+ + O`$ | $`1.5\times10^{-16}\,(T_i/300)^{0.55}`$ | 1.35 |
| $`O^+(^2D) + N_2 \rightarrow NO^+ + N(^4S)`$ | $`2.5\times10^{-17}`$ | 4.41 |
| $`O^+(^2D) + N_2 \rightarrow O^+(^4S) + N_2`$ | $`8.0\times10^{-16}`$ | 3.31 |
| $`O^+(^2D) + O_2 \rightarrow O_2^+ + O`$ | $`7.0\times10^{-16}`$ | 4.865 |
| $`O^+(^2D) + NO \rightarrow NO^+ + O`$ | $`1.2\times10^{-15}`$ | 4.37 |
| $`O^+(^2D) + N(^4S) \rightarrow N^+ + O`$ | $`1.5\times10^{-16}`$ | 1.0 |
| $`O^+(^2D) + O \rightarrow O^+(^4S) + O`$ | $`1.0\times10^{-17}`$ | see below |
| $`O^+(^2D) + e^- \rightarrow O^+(^4S) + e^-`$ | $`6.03\times10^{-14}\,(T_e/300)^{-0.5}`$ | 3.31 |
| $`O^+(^2D) \rightarrow O^+(^4S) + h\nu\ (3726\,\mathrm{\AA})`$ | $`7.7\times10^{-5}\ \mathrm{s^{-1}}`$ | -- |

For $`O^+(^2D) + O`$, the neutral oxygen is returned 50% as $`O(^3P)`$
(3.31 eV) and 50% as $`O(^1D)`$ (1.35 eV), so the heating is
$`0.5\times3.31 + 0.5\times1.35 = 2.33`$ eV. If $`O(^1D)`$ were not a
species, the heating would be 3.31 eV and no $`O(^1D)`$ would be produced.

### O$`^+(^2P)`$ as reactant

| Reaction | $`k`$ (m³ s⁻¹) | Heating (eV) |
|---|---|---|
| $`O^+(^2P) + N_2 \rightarrow N_2^+ + O`$ | $`2.0\times10^{-16}\,(T_i/300)^{0.55}`$ | 3.05 |
| $`O^+(^2P) + N_2 \rightarrow N^+ + NO`$ | $`1.0\times10^{-16}`$ | 0.70 |
| $`O^+(^2P) + O_2 \rightarrow O_2^+ + O`$ | $`1.3\times10^{-16}`$ | 6.54 |
| $`O^+(^2P) + O_2 \rightarrow O^+(^4S) + O_2`$ | $`1.3\times10^{-16}`$ | 5.016 |
| $`O^+(^2P) + N(^4S) \rightarrow N^+ + O`$ | $`1.0\times10^{-16}`$ | 2.7 |
| $`O^+(^2P) + O \rightarrow O^+(^4S) + O`$ | $`4.0\times10^{-16}`$ | 5.0 |
| $`O^+(^2P) + e^- \rightarrow O^+(^4S) + e^-`$ | $`3.03\times10^{-14}\,(T_e/300)^{-0.5}`$ | 5.0 |
| $`O^+(^2P) + e^- \rightarrow O^+(^2D) + e^-`$ | $`1.84\times10^{-13}\,(T_e/300)^{-0.5}`$ | 1.69 |
| $`O^+(^2P) \rightarrow O^+(^4S) + h\nu\ (2470\,\mathrm{\AA})`$ | $`0.047\ \mathrm{s^{-1}}`$ | -- |
| $`O^+(^2P) \rightarrow O^+(^2D) + h\nu\ (7320\,\mathrm{\AA})`$ | $`0.171\ \mathrm{s^{-1}}`$ | -- |

### N$`^+`$ as reactant

The $`N^+ + O_2`$ reaction has five channels, each with the same temperature
form, $`k = a\,(T_i/300)^{-0.45}`$ for $`T_i \le 1000`$ K and a constant $`b`$
above.

| Reaction | $`a`$ (m³ s⁻¹) | $`b`$ (m³ s⁻¹) | Heating (eV) |
|---|---|---|---|
| $`N^+ + O_2 \rightarrow O_2^+ + N(^4S)`$ | $`1.925\times10^{-16}`$ | $`3.325\times10^{-16}`$ | 2.5 |
| $`N^+ + O_2 \rightarrow O_2^+ + N(^2D)`$ | $`0.825\times10^{-16}`$ | $`1.425\times10^{-16}`$ | 0.1 |
| $`N^+ + O_2 \rightarrow O^+(^4S) + NO`$ | $`0.275\times10^{-16}`$ | $`0.475\times10^{-16}`$ | 2.31 |
| $`N^+ + O_2 \rightarrow NO^+ + O(^3P)`$ | $`0.495\times10^{-16}`$ | $`0.855\times10^{-16}`$ | 6.67 |
| $`N^+ + O_2 \rightarrow NO^+ + O(^1D)`$ | $`1.98\times10^{-16}`$ | $`3.42\times10^{-16}`$ | 4.71 |

In the last channel $`O(^1D)`$ is produced because it is a species; otherwise
$`O(^3P)`$ would be produced.

| Reaction | $`k`$ (m³ s⁻¹) | Heating (eV) |
|---|---|---|
| $`N^+ + NO \rightarrow N_2^+ + O`$ | $`8.33\times10^{-17}\,(T_i/300)^{-0.24}`$ | 2.2 |
| $`N^+ + NO \rightarrow NO^+ + N(^4S)`$ | $`4.72\times10^{-16}\,(T_i/300)^{-0.24}`$ | 3.4 |
| $`N^+ + O \rightarrow O^+(^4S) + N(^4S)`$ | $`2.2\times10^{-18}`$ | 0.93 |

### He$`^+`$ as reactant

| Reaction | $`k`$ (m³ s⁻¹) | Heating (eV) |
|---|---|---|
| $`He^+ + N_2 \rightarrow N^+ + N(^4S) + He`$ | $`7.8\times10^{-16}`$ | 0.28 |
| $`He^+ + N_2 \rightarrow N_2^+ + He`$ | $`5.2\times10^{-16}`$ | -- |
| $`He^+ + O_2 \rightarrow O^+(^4S) + O(^3P) + He`$ | $`9.7\times10^{-16}`$ | -- |

---

## Dissociative and Radiative Recombination

These use the electron temperature $`T_e`$ except where noted.

### N$`_2^+`$ + e

```math
\begin{align}
N_2^+ + e^- \rightarrow 2\,N(^2D) \quad (56\%), \qquad
N_2^+ + e^- \rightarrow 2\,N(^4S) \quad (44\%)
\end{align}
```

```math
\begin{align}
k = 2.2\times10^{-13}\left(\frac{T_e}{300}\right)^{-0.39} \ \mathrm{m^3\,s^{-1}}
\end{align}
```

Heating is $`0.56\times1.04 + 0.44\times5.77`$ eV (1.04 eV for the
$`N(^2D)`$ channel, 5.77 eV for the $`N(^4S)`$ channel).

### O$`_2^+`$ + e

```math
\begin{align}
k = \begin{cases}
1.95\times10^{-13}\,(T_e/300)^{-0.7} & T_i \le 1200\ \mathrm{K} \\
7.39\times10^{-14}\,(T_e/1200)^{-0.56} & T_i > 1200\ \mathrm{K}
\end{cases} \ \mathrm{m^3\,s^{-1}}
\end{align}
```

The switch between the two forms is made on $`T_i`$, not on $`T_e`$. The
product branching (applies when $`O(^1D)`$ is a species):

| Products | Fraction | Heating (eV) |
|---|---|---|
| $`O(^3P) + O(^3P)`$ | 0.22 | 6.99 |
| $`O(^3P) + O(^1D)`$ | 0.42 | 5.02 |
| $`O(^1D) + O(^1D)`$ | 0.31 | 3.06 |
| $`O(^3P) + O(^1D)`$ | 0.05 | -- |

Total heating is $`0.22\times6.99 + 0.42\times5.02 + 0.31\times3.06`$ eV. (The
last 0.05 is a placeholder in the code; it is meant to go to $`O(^1D)`$ and
$`O(^1S)`$.) If $`O(^1D)`$ were not a species, every reaction would give
$`2\,O(^3P)`$ with 5.0 eV.

### NO$`^+`$ + e

```math
\begin{align}
k = 3.4\times10^{-13}\left(\frac{T_e}{300}\right)^{-0.85}: \quad NO^+ + e^- \rightarrow O + N(^2D) \quad (0.38\ \mathrm{eV})
\end{align}
```

```math
\begin{align}
k = 0.6\times10^{-13}\left(\frac{T_e}{300}\right)^{-0.85}: \quad NO^+ + e^- \rightarrow O + N(^4S) \quad (2.77\ \mathrm{eV})
\end{align}
```

Both $`k`$ are in $`\mathrm{m^3\,s^{-1}}`$.

### He$`^+`$ radiative recombination

```math
\begin{align}
He^+ + e^- \rightarrow He + h\nu, \qquad
k = 4.8\times10^{-18}\left(\frac{250}{T_e}\right)^{0.7} \ \mathrm{m^3\,s^{-1}}
\end{align}
```

No heating is applied.
