# Psychrometric Relationships
#
# The purpose is to document a consistent and reliable set 
# of reference equations for psychrometric data validation.

# Bruce Wernick
# 20 November 2026

# Properties
# pa  Barometric Pressure        [kPa]
# db  Dry Bulb Temperature       [°C]
# wb  Wet Bulb Temperature       [°C]
# ha  Specific Enthalpy          [kJ/kg]
# wa  Absolute Humidity Ratio    [kg/kg]
# rh  Relative Humidity          [fraction]
# td  Dew Point Temperature      [°C]
# pv  Vapor Pressure             [kPa]
# ps  Saturated Vapor Pressure   [kPa]
# va  Specific Volume            [m³/kg]
# ro  Moist Air Density          [kg (water + dry air) / m³]
# cp  Specific Heat              [kJ/kg-K]
# za  Sigma Heat                 [kJ/kg]
# ds  Degree of Saturation       [fraction]
# sa  Specific Entropy           [kJ/kg-K]

## Reference conditions
P0 = 101.325          # barometric pressure at sea level [kPa]
K0 = 273.15           # °C to Kelvin offset
PS0 = 0.611213        # Saturation pressure at 0°C [kPa]
TMIN = -100.0         # Minimum valid temperature (°C)
TMAX = 200.0          # Maximum valid temperature (°C)
TRIPLE = 0.01         # Triple point temperature (°C)

# Air and Water Constants

Ma = 28.9645077       # molecular weight of air
Mw = 18.01528         # molecular weight of water
Mwa = Mw / Ma         # Ratio of molar masses (water / dry air)
RU = 8.31446261815324 # universal gas constant
Ra = RU / Ma          # air gas constant
Rw = RU / Mw          # water gas constant
Raw = Ra / Rw         # gas constant ratio

CPa = 1.006           # Dry air specific heat [kJ/kg-K]
CPv = 1.82            # Water vapor specific heat (Gatley pg 228)
CPw = 4.186           # Liquid water specific heat
CPi = 2.093           # Ice specific heat
Hfgo = 2500.77        # Latent heat of vaporization at 0°C (kJ/kg)
Hsgo = 2834.2         # Latent heat of sublimation at 0°C (kJ/kg)
Sgo = Hfgo / K0       # saturated entropy of gas at 0°C

fs = 1.0044           # molecular interaction coeff used by Gatley

# ---------------------------------------------------------

# Absolute Humidity Ratio (wa [kg/kg])
wa = Mwa * f * pv / (pa - f * pv)

# ---------------------------------------------------------

# Relative Humidity (rh [fraction])
rh = pv / ps

# ---------------------------------------------------------

# Specific Enthalpy (ha)
ha = CPa * db + wa * (Hfgo + Cpv * db)

# ---------------------------------------------------------

# Adiabatic Saturation (introducing wet bulb)
# wa_from_dbwb()

wsw = wsat(wb)
wb >= TRIPLE
  # Above Freezing 
  num = wsw * (Hfgo - wb * (Cpw - Cpv)) - Cpa*(db - wb)
  den = Hfgo + Cpv * db - Cpw * wb 
wb < TRIPLE
  # Below Freezing wb < TRIPLE
  num = wsw * (Hsgo - wb * (Cpi - Cpv)) - Cpa*(db - wb)
  den = Hsgo + Cpv * db - Cpi * wb 
wa = num / den

# ---------------------------------------------------------

# Wet bulb temperature
# Inverse of Adiabatic Saturation using root finder.

Given pa, db, wa 
let wb0 be a guess of wb
fsolve(lambda x: wa - wa_from_dbwb(pa, db, x), wb0)

# ---------------------------------------------------------

# Saturation Pressure (ps = f(t))

# Magnus equation (t in °C and p in kPa)
# (simple, fast but not accurate)
# listed in Wikipedia as Tetens equation
# also known as August–Roche–Magnus formula
a, b, c = 0.61078, 17.27, 237.3 # above zero
a, b, c = 0.61078, 21.875, 265.5 # below zero (Murray-1967)
p = a * exp(b * t / (t + c))
x = ln(p / a)
t = c * x / (b - x)

# ASHRAE curve fit (recently added inverse)
# Saturation pressure over ice: 
# 173.15 K to 273.15 K (-100°C to 0°C)
C1 = -5.6745359e3
C2 = 6.3925247
C3 = -9.6778430e-3
C4 = 6.2215701e-7
C5 = 2.0747825e-9
C6 = -9.4840240e-13
C7 = 4.1635019
poly = C2 + T * (C3 + T * (C4 + T * (C5 + T * C6)))
ln_pa = C1 / T + poly + C7 * ln(T)

# Saturation pressure over liquid water: 
# 273.15 K to 473.15 K (0°C to 200°C)
C8 = -5.8002206e3
C9 = 1.3914993
C10 = -4.8640239e-2
C11 = 4.1764768e-5
C12 = -1.4452093e-8
C13 = 6.5459673
poly = C9 + T * (C10 + T * (C11 + T * (C12)))
ln_pa = C8 / T + poly + C13 * ln(T)

# Convert ln(Pa) to ln(kPa) by subtracting ln(1000)
ps = exp(ln_pa - 6.907755278982137)

# ASHRAE inverse function

a = log(p_kpa)

# 0°C to 93°C (liquid water)
C14 = 6.54
C15 = 14.526
C16 = 0.7389
C17 = 0.09486
C18 = 0.4569
T = C14 + a*(C15 + a*(C16 + a*(C17))) + C18*(p_kpa**0.1984)

# -60°C to 0°C (over ice)
T = 6.09 + a*(12.608 + 0.4959*a)


# IAPWS - saturation pressure coefficients (Region 4)
# Note, IAPWS Region 4 is strictly a liquid-vapor saturation 
# equation (273.15 K to 647.096 K). When fed temperatures below 0°C, 
# it extrapolates the saturation pressure of supercooled liquid water.
# Use IAPWS-08—specifically for sublimation pressure over ice.

nps = (1167.0521452767, -724213.16703206,
       -17.073846940092, 12020.82470247,
       -3232555.0322333, 14.91510861353,
       -4823.2657361591, 405113.40542057,
       -0.23855557567849, 650.17534844798)

# P from T
g = T + nps[8] / (T - nps[9])
a = g*(g + nps[0]) + nps[1]
b = g*(g*nps[2] + nps[3]) + nps[4]
c = g*(g*nps[5] + nps[6]) + nps[7]
d = -b + sqrt(b*b - 4*a*c)
p = 1e3*(2*c/d)**4

# T from P
g = (1e-3*p)**0.25
a = (g + nps[2]) * g + nps[5]
b = (nps[0] * g + nps[3]) * g + nps[6]
c = nps[7] + (nps[1] * g + nps[4]) * g
d = -2*c/(b + sqrt(b*b - 4*a*c))
nd = nps[9] + d
T = 0.5*(nd - sqrt(nd*nd - 4*(nps[9]*d + nps[8])))

IAPWS over ice
T = t + K0
theta = T / TRIPLE
a = (-21.2144006, 27.3203819, -6.10598130)
b = (0.00333333333, 1.20666667, 1.70333333)
ln_pi = (1.0 / theta) * sum(ai * (theta**bi) for ai, bi in zip(a, b))
p = 0.611657 * exp(ln_pi)

# ---------------------------------------------------------

## Dew Point Temperature [°C]

# Inverse of ps = psat(t)
# Calculate by using root finder or a correlation
# Use Magnus as a starting guess

# ---------------------------------------------------------

## Specific Volume [m³/kg]
va = Ra * (db + K0) * (1.0 + wa / Mwa) / pa

# ---------------------------------------------------------

## Moist Air Density [kg/m³]
ro = (1.0 + wa) / va

# ---------------------------------------------------------

## Specific Heat [kJ/kg-K]
cp = CPa + wa * CPv

# ---------------------------------------------------------

## Degree of saturation [fraction]
ds = wa / ws

# ---------------------------------------------------------

## Sigma Heat
za = ha - wa * Cpw * wb

# ---------------------------------------------------------

## Entropy [kJ/kg-K]

# Dry air
sa = Cpa * ln(T / K0) - Ra * ln((pa - pv) / P0)

# Moist air
sg = Sgo + Cpv * ln(T / K0) - Raw * ln(pv / PS0)
sa = sd + wa * sg  

# ---------------------------------------------------------

## Calculate all properties

# Given these three properties:
pa = 101.325 # kPa sea level barometric pressure
db = 23.0    # °C typical indoor comfort conditions
wb = 16.0    # °C round figure resulting in about 50% rh

# This, I will declare as my reference condition
# for comparing to all other calculators.
# It is possible to calculate all other directly 
# (no iteration required)

ps = psat(db) # ASHRAE pressure eqn
wa = wa_from_dbwb(pa, db, wb) # adi sat
pv = pv_from_wa(pa, wa) # wa eqn
ha = ha_from_wa(db, wa) # enthalpy
rh = rh_from_pvps(pv, ps)
td = td_from_wa(pa, wa) # inverse psat using pv
va = va_from_dbwa(pa, db, wa)
ro = ro_from_wava(wa, va)
cp = cp_from_wa(wa)
za = za_from_wb(wb)

# In my function name, I leave out the reference to
# barometric pressure because it is never calculated.

# ---------------------------------------------------------

# Reference Point Benchmark
# (pa=101.325 kPa, db=23.0 °C, wb=16.0 °C)

# Property             Sym    Unit      Gatley     Ideal Gas   Wernick
#                                       f=1.0044   f=1
# ---------------------------------------------------------------------
# Sat. Vapor Pressure   ps    kPa       2.810442   2.810442    2.81044
# Humidity Ratio        wa    g/kg      8.514      8.463       8.51405
# Vapor Pressure        pv    kPa       1.362281   1.360242    1.36228
# Specific Enthalpy     ha    kJ/kg    44.786080  44.657250   44.7861
# Relative Humidity     rh    %        48.472125  48.399560   48.4721
# Dew Point Temp        td    °C       11.558561  11.535923   11.5586
# Specific Volume       va    m³/kg     0.850487   0.850419    0.850487
# Moist Air Density     ro    kg/m³     1.185808   1.185843    1.18581
# Specific Heat         cp    kJ/kg-K   1.021496   1.021403    1.0215
# Sigma Heat            za    kJ/kg    44.215843  44.090407   44.2158
# Degree of Saturation  ds    frac      0.477668   0.476974

# The Wernick calculations are done in a Python program 
# entitled psy_wernick.py.  These agree fully with the
# equations posted in this report.

# ---------------------------------------------------------

## Limits stored in air_meta.toml data file

[pa]
label = "Barometric Pressure"
name = "Barom"
value = 101.325
minval = 78
maxval = 119
ndigs = 3
sf = 1.0
units = "kPa"

[db]
label = "Dry Bulb Temperature"
name = "DryBulb"
value = 23.0
minval = -60
maxval = 90
ndigs = 3
sf = 1.0
units = "°C"

[wb]
label = "Wet Bulb Temperature"
name = "WetBulb"
value = 16.0
minval = -60
maxval = 90
ndigs = 3
sf = 1.0
units = "°C"

[ha]
label = "Specific Enthalpy"
name = "Enthalpy"
value = 44.786
minval = -60
maxval = 355
ndigs = 3
sf = 1.0
units = "kJ/kg"

[wa]
label = "Absolute Humidity"
name = "Humidity"
value = 8.514
minval = 0.0
maxval = 100
ndigs = 3
sf = 1.0
units = "g/kg"

[rh]
label = "Relative Humidity"
name = "RelHumid"
value = 48.47
minval = 0.0
maxval = 100.0
ndigs = 3
sf = 1.0
units = "%"

[td]
label = "Dew Point Temperature"
name = "DewPoint"
value = 11.558
minval = -60
maxval = 90
ndigs = 3
sf = 1.0
units = "°C"

[pv]
label = "Vapor Pressure"
name = "VapPress"
value = 1.36228
minval = 0.001
maxval = 12.4
ndigs = 3
sf = 1.0
units = "kPa"

[ps]
label = "Saturation Pressure"
name = "SatPress"
value = 2.81
minval = 0.001
maxval = 12.4
ndigs = 3
sf = 1.0
units = "kPa"

[va]
label = "Specific Volume"
name = "SpecVol"
value = 0.85048
minval = 0.514
maxval = 1.41
ndigs = 4
sf = 1.0
units = "m³/kg"

[ro]
label = "True Moist Air Density"
name = "Density"
value = 1.1858147
minval = 0.792
maxval = 1.95
ndigs = 4
sf = 1.0
units = "kg/m³"

[cp]
label = "Specific Heat"
name = "SpecHeat"
value = 1.02149
minval = 0.5
maxval = 5.0
ndigs = 4
sf = 1.0
units = "kJ/kg-K"

[za]
label = "Sigma Heat"
name = "SigHeat"
value = 44.21585
minval = -60
maxval = 330
ndigs = 3
sf = 1.0
units = "kJ/kg"

[ws]
label = "Saturated Abs Humidity"
name = "SatHumidity"
value = 17.8242
minval = 0
maxval = 117
ndigs = 3
sf = 1.0
units = "g/kg"

[ds]
label = "Degree of Saturation"
name = "DegSat"
value = 0.4776678
minval = 0.0
maxval = 1.0
ndigs = 3
sf = 1.0
units = "-"

[sa]
label = "Specific Entropy"
name = "Entropy"
value = 0.1601728
minval = -0.30
maxval = 1.25
ndigs = 3
sf = 1.0
units = "kJ/kg-K"
