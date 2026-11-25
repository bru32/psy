# Psychrometric Relationships & Moist Air State Calculations (Python)

A Python reference library and validation benchmark for **psychrometry**, **moist air thermodynamics**, and **HVAC air-conditioning** calculations. 

Based on **ASHRAE** formulations, **Donald Gatley's** *Understanding Psychrometrics*, and **PsychroLib**, this repository provides verified equations for **dry bulb temperature**, **wet bulb temperature**, **dew point**, **relative humidity**, **absolute humidity ratio**, **specific enthalpy**, and **psychrometric charts (psychart)**


**Author:** Bruce Wernick  
**Date:** 20 November 2026  

Copyright © (2026) Bruce Wernick,  All rights reserved...

---

## Properties

| Symbol | Property | Unit |
| :--- | :--- | :--- |
| `pa` | Barometric Pressure | $\text{kPa}$ |
| `db` | Dry Bulb Temperature | $^\circ\text{C}$ |
| `wb` | Wet Bulb Temperature | $^\circ\text{C}$ |
| `ha` | Specific Enthalpy | $\text{kJ/kg}$ |
| `wa` | Absolute Humidity Ratio | $\text{kg/kg}$ |
| `rh` | Relative Humidity | $\text{fraction}$ |
| `td` | Dew Point Temperature | $^\circ\text{C}$ |
| `pv` | Vapor Pressure | $\text{kPa}$ |
| `ps` | Saturated Vapor Pressure | $\text{kPa}$ |
| `va` | Specific Volume | $\text{m}^3/\text{kg}$ |
| `ro` | Moist Air Density | $\text{kg (water + dry air) / m}^3$ |
| `cp` | Specific Heat | $\text{kJ/kg}\cdot\text{K}$ |
| `za` | Sigma Heat | $\text{kJ/kg}$ |
| `ds` | Degree of Saturation | $\text{fraction}$ |
| `sa` | Specific Entropy | $\text{kJ/kg}\cdot\text{K}$ |

---

## Reference Conditions

```python
P0 = 101.325          # Barometric pressure at sea level [kPa]
K0 = 273.15           # °C to Kelvin offset
PS0 = 0.611213        # Saturation pressure at 0°C [kPa]
TMIN = -100.0         # Minimum valid temperature [°C]
TMAX = 200.0          # Maximum valid temperature [°C]
TRIPLE = 0.01         # Triple point temperature [°C]
```

## Air and Water Constants

```python
Ma = 28.9645077       # Molecular weight of air
Mw = 18.01528         # Molecular weight of water
Mwa = Mw / Ma         # Ratio of molar masses (water / dry air)
RU = 8.31446261815324 # Universal gas constant
Ra = RU / Ma          # Air gas constant
Rw = RU / Mw          # Water gas constant
Raw = Ra / Rw         # Gas constant ratio

CPa = 1.006           # Dry air specific heat [kJ/kg-K]
CPv = 1.82            # Water vapor specific heat (Gatley pg 228) [kJ/kg-K]
CPw = 4.186           # Liquid water specific heat [kJ/kg-K]
CPi = 2.093           # Ice specific heat [kJ/kg-K]
Hfgo = 2500.77        # Latent heat of vaporization at 0°C [kJ/kg]
Hsgo = 2834.2         # Latent heat of sublimation at 0°C [kJ/kg]
Sgo = Hfgo / K0       # Saturated entropy of gas at 0°C

fs = 1.0044           # Molecular interaction coeff used by Gatley
```

---

## Core Psychrometric Equations

### Absolute Humidity Ratio $w_a$ $\text{kg/kg}$
$$w_a = \frac{M_{wa} \cdot f \cdot p_v}{p_a - f \cdot p_v}$$

### Relative Humidity $rh$ $\text{fraction}$
$$rh = \frac{p_v}{p_s}$$

### Specific Enthalpy $h_a$ $\text{kJ/kg}$
$$h_a = CP_a \cdot db + w_a \cdot (H_{fgo} + CP_v \cdot db)$$

### Adiabatic Saturation (Introducing Wet Bulb)
Let $w_{sw} = w_{\text{sat}}(wb)$:

* **Above Freezing $wb \ge \text{TRIPLE}$:**
  $$w_a = \frac{w_{sw} \cdot \left(H_{fgo} - wb \cdot (CP_w - CP_v)\right) - CP_a \cdot (db - wb)}{H_{fgo} + CP_v \cdot db - CP_w \cdot wb}$$

* **Below Freezing ($wb < \text{TRIPLE}$):**
  $$w_a = \frac{w_{sw} \cdot \left(H_{sgo} - wb \cdot (CP_i - CP_v)\right) - CP_a \cdot (db - wb)}{H_{sgo} + CP_v \cdot db - CP_i \cdot wb}$$

### Wet Bulb Temperature $wb$ $^\circ\text{C}$
Calculated as the inverse of Adiabatic Saturation using a root finder. Given $p_a$, $db$, and $w_a$, let $wb_0$ be an initial guess of $wb$:

```python
fsolve(lambda x: wa - wa_from_dbwb(pa, db, x), wb0)
```

---

## Saturation Pressure $p_s = f(t)$

### 1. Magnus Equation
*(Simple and fast, but less accurate. Listed in Wikipedia as the Tetens equation; also known as the August–Roche–Magnus formula, where $t$ is in $^\circ\text{C}$ and $p$ is in $\text{kPa}$.)*

```python
# Above zero (°C)
a, b, c = 0.61078, 17.27, 237.3

# Below zero (°C) (Murray-1967)
a, b, c = 0.61078, 21.875, 265.5

# Forward and Inverse
p = a * exp(b * t / (t + c))
x = ln(p / a)
t = c * x / (b - x)
```

### 2. ASHRAE Curve Fit

#### Saturation pressure over ice: $173.15\text{ K}$ to $273.15\text{ K}$ ($-100^\circ\text{C}$ to $0^\circ\text{C}$)
```python
C1 = -5.6745359e3
C2 = 6.3925247
C3 = -9.6778430e-3
C4 = 6.2215701e-7
C5 = 2.0747825e-9
C6 = -9.4840240e-13
C7 = 4.1635019

poly = C2 + T * (C3 + T * (C4 + T * (C5 + T * C6)))
ln_pa = C1 / T + poly + C7 * ln(T)
```

#### Saturation pressure over liquid water: $273.15\text{ K}$ to $473.15\text{ K}$ ($0^\circ\text{C}$ to $200^\circ\text{C}$)
```python
C8 = -5.8002206e3
C9 = 1.3914993
C10 = -4.8640239e-2
C11 = 4.1764768e-5
C12 = -1.4452093e-8
C13 = 6.5459673

poly = C9 + T * (C10 + T * (C11 + T * (C12)))
ln_pa = C8 / T + poly + C13 * ln(T)
```

#### Convert $\ln(\text{Pa})$ to $\ln(\text{kPa})$ by subtracting $\ln(1000)$:
```python
ps = exp(ln_pa - 6.907755278982137)
```

#### ASHRAE Inverse Function
Let $a = \ln(p_{\text{kPa}})$:

```python
# 0°C to 93°C (liquid water)
C14 = 6.54
C15 = 14.526
C16 = 0.7389
C17 = 0.09486
C18 = 0.4569
T = C14 + a * (C15 + a * (C16 + a * C17)) + C18 * (p_kpa ** 0.1984)

# -60°C to 0°C (over ice)
T = 6.09 + a * (12.608 + 0.4959 * a)
```

### 3. IAPWS Saturation Pressure Coefficients (Region 4)
> **Note:** IAPWS Region 4 is strictly a liquid-vapor saturation equation ($273.15\text{ K}$ to $647.096\text{ K}$). When fed temperatures below $0^\circ\text{C}$, it extrapolates the saturation pressure of supercooled liquid water. Use **IAPWS-08** specifically for sublimation pressure over ice.

```python
nps = (1167.0521452767, -724213.16703206,
       -17.073846940092, 12020.82470247,
       -3232555.0322333, 14.91510861353,
       -4823.2657361591, 405113.40542057,
       -0.23855557567849, 650.17534844798)

# P from T
g = T + nps[8] / (T - nps[9])
a = g * (g + nps[0]) + nps[1]
b = g * (g * nps[2] + nps[3]) + nps[4]
c = g * (g * nps[5] + nps[6]) + nps[7]
d = -b + sqrt(b * b - 4 * a * c)
p = 1e3 * (2 * c / d) ** 4

# T from P
g = (1e-3 * p) ** 0.25
a = (g + nps[2]) * g + nps[5]
b = (nps[0] * g + nps[3]) * g + nps[6]
c = nps[7] + (nps[1] * g + nps[4]) * g
d = -2 * c / (b + sqrt(b * b - 4 * a * c))
nd = nps[9] + d
T = 0.5 * (nd - sqrt(nd * nd - 4 * (nps[9] * d + nps[8])))

# IAPWS over ice
T = t + K0
theta = T / TRIPLE
a = (-21.2144006, 27.3203819, -6.10598130)
b = (0.00333333333, 1.20666667, 1.70333333)
ln_pi = (1.0 / theta) * sum(ai * (theta ** bi) for ai, bi in zip(a, b))
p = 0.611657 * exp(ln_pi)
```

---

## Additional Derived Properties

* **Dew Point Temperature $t_d$ $^\circ\text{C}$:**  
  Inverse of $p_s = p_{\text{sat}}(t)$. Calculate using a root finder or correlation with the Magnus equation as a starting guess.

* **Specific Volume $v_a$ $\text{m}^3/\text{kg}$:**  
  $$v_a = \frac{R_a \cdot (db + K_0) \cdot \left(1.0 + \frac{w_a}{M_{wa}}\right)}{p_a}$$

* **Moist Air Density $\rho$ $\text{kg/m}^3$:**  
  $$\rho = \frac{1.0 + w_a}{v_a}$$

* **Specific Heat $c_p$ $\text{kJ/kg}\cdot\text{K}$:**  
  $$c_p = CP_a + w_a \cdot CP_v$$

* **Degree of Saturation $ds$ $\text{fraction}$:**  
  $$ds = \frac{w_a}{w_s}$$

* **Sigma Heat $z_a$ $\text{kJ/kg}$:**  
  $$z_a = h_a - w_a \cdot CP_w \cdot wb$$

* **Specific Entropy $s_a$ $\text{kJ/kg}\cdot\text{K}$:**  
  * *Dry air:*
    $$s_d = CP_a \cdot \ln\left(\frac{T}{K_0}\right) - R_a \cdot \ln\left(\frac{p_a - p_v}{P_0}\right)$$
  * *Moist air:*
    $$s_g = S_{go} + CP_v \cdot \ln\left(\frac{T}{K_0}\right) - R_{aw} \cdot \ln\left(\frac{p_v}{PS_0}\right)$$
    $$s_a = s_d + w_a \cdot s_g$$

---

## Calculate All Properties (Direct Workflow)

Given these three reference properties (typical indoor comfort conditions at sea level, resulting in ~50% RH), it is possible to calculate all other properties directly without iteration:

```python
pa = 101.325 # kPa sea level barometric pressure
db = 23.0    # °C typical indoor comfort conditions
wb = 16.0    # °C round figure resulting in about 50% rh

ps = psat(db)                 # ASHRAE pressure eqn
wa = wa_from_dbwb(pa, db, wb) # Adiabatic saturation
pv = pv_from_wa(pa, wa)       # Vapor pressure from wa eqn
ha = ha_from_wa(db, wa)       # Specific enthalpy
rh = rh_from_pvps(pv, ps)     # Relative humidity
td = td_from_wa(pa, wa)       # Dew point (inverse psat using pv)
va = va_from_dbwa(pa, db, wa) # Specific volume
ro = ro_from_wava(wa, va)     # Moist air density
cp = cp_from_wa(wa)           # Specific heat
za = za_from_wb(wb)           # Sigma heat
```
*(Note: Function names omit the reference to barometric pressure because it is never calculated.)*

---

## Reference Point Benchmark
**Conditions:** $p_a = 101.325\text{ kPa}$, $db = 23.0^\circ\text{C}$, $wb = 16.0^\circ\text{C}$

| Property | Sym | Unit | Gatley ($f=1.0044$) | Ideal Gas ($f=1$) | Wernick |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **Sat. Vapor Pressure** | `ps` | $\text{kPa}$ | 2.810442 | 2.810442 | 2.81044 |
| **Humidity Ratio** | `wa` | $\text{g/kg}$ | 8.514 | 8.463 | 8.51405 |
| **Vapor Pressure** | `pv` | $\text{kPa}$ | 1.362281 | 1.360242 | 1.36228 |
| **Specific Enthalpy** | `ha` | $\text{kJ/kg}$ | 44.786080 | 44.657250 | 44.7861 |
| **Relative Humidity** | `rh` | $\%$ | 48.472125 | 48.399560 | 48.4721 |
| **Dew Point Temp** | `td` | $^\circ\text{C}$ | 11.558561 | 11.535923 | 11.5586 |
| **Specific Volume** | `va` | $\text{m}^3/\text{kg}$ | 0.850487 | 0.850419 | 0.850487 |
| **Moist Air Density** | `ro` | $\text{kg/m}^3$ | 1.185808 | 1.185843 | 1.18581 |
| **Specific Heat** | `cp` | $\text{kJ/kg}\cdot\text{K}$ | 1.021496 | 1.021403 | 1.0215 |
| **Sigma Heat** | `za` | $\text{kJ/kg}$ | 44.215843 | 44.090407 | 44.2158 |


## Search Keywords
  **abs humidity**,
  **air**,
  **air-conditioning**,
  **air-water-mixture**,
  **ASHRAE**,
  **dew-point temperature**,
  **drybulb**,
  **enthalpy**,
  **gatley**,
  **humidity**,
  **moist air state**,
  **moisture**,
  **psychart**,
  **psychrometry**,
  **psychrometrics**,
  **python**,
  **relative humidity**,
  **temperature**,
  **wetbulb**
