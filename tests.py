
from psy_wernick import *


# Utils for testing

def dicts_are_close(d1, d2, rel_tol=1e-9, abs_tol=0.0):
  # Ensure both dictionaries have the exact same keys
  if d1.keys() != d2.keys():
    return False

  # Check that every corresponding float pair is close
  return all(
    math.isclose(d1[k], d2[k], rel_tol=rel_tol, abs_tol=abs_tol) for k in d1)

from itertools import combinations

def compare(air1, air2, key):
  print(f"({key}) passed ✔️") if air1.is_close(air2) else print(f"({key}) failed ❌")

def frange(start: float, stop: float, step: float):
  """Yield floating-point numbers from start up to stop (inclusive).
  """
  steps = int(round((stop - start) / step))
  for i in range(steps + 1):
    yield round(start + i * step, 10)

# =========================================================

#####################
# Some general tests
#####################

def test0():
  """ Compare astate and calcstate with same inputs.
  """
  print()
  print("Test 0: Compare astate and calcstate at ref conditions.")
  print()
  pa, db, wb = 101.325, 23.0, 16.0
  air1 = astate(pa, db, wb)

  props = ("wb", "ha", "wa", "rh", "td", "pv", "va", "ro", "za")
  cases = {k: getattr(air1, k) for k in props}
  for param, val in cases.items():
    air2 = calcstate(pa, db=db, **{param: val})
    compare(air1, air2, f"db,{param}")

  print()
  print("-"*60)

def test1():
  """ Try all calcstate combinations.
  """

  print("\nTest 1: Try all calcstate combinations.\n")

  pa, db, wb = 101.325, 23.0, 16.0
  air = astate(pa, db, wb)
  print(air.report("Ref Air State"))
  print()

  keys = "db,wb,ha,wa,rh,td,pv,ps,va,ro,za".split(',')
  for k1, k2 in combinations(keys, 2):
    kwargs = {k1: air[k1], k2: air[k2]}
    label = f"{k1},{k2}"
    try:
      compare(air, calcstate(pa, **kwargs), label)
    except Exception as e:
      print(f"({k1},{k2}) ❌ {e}")
  
  print()
  print("-"*60)

def test2():
  """ Test typical astate and standard air.
  """

  print("\nTest 2: Typical astate and standard air\n")
  
  # Typical aircon conditions
  pa, db, wb = 101.325, 23.0, 16.0
  air = astate(pa, db, wb)
  print(air.report("Typical Air"))
  print()

  # Standard air 
  # (sea level, 20.0 °C 1.200 kg/m³)
  air = calcstate(pa, db=20, ro=1.200)
  print(air.report("Standard Air"))

  print()
  print("-"*60)

def test3():
  """ Calculate moist air states 
      across the full psychrometric chart.
  """

  print("\nTest 3: Property extremes across FULL range.\n")

  pa = 101.325 # kPa
  db0, db1 = -60.0, 90.0 # dry bulb range [°C]
  pc0, pc1 = 0.1, 99.9 # percent saturation
  n, m = 40, 20

  k = 0

  # Dictionary to store {'property': {'min': val, 'max': val}}
  extremes = {}

  for i in range(n):
    
    db = db0 + (db1 - db0) * i / (n - 1)
    w_max = wsat(pa, db) # kg/kg

    for j in range(m):
      
      # Step in percent of saturation, pc [%]
      pc = pc0 + (pc1 - pc0) * j / (m - 1)
      wa = w_max * pc / 100

      wb = wb_from_dbwa(pa, db, wa)
      air = astate(pa, db, wb)
      
      # Convert dataclass instance to a dictionary
      state_dict = asdict(air)

      if not extremes:
        extremes = {prop: {"min": val, "max": val} for prop, val in state_dict.items()}
      else:
        for prop, val in state_dict.items():
          if val < extremes[prop]["min"]:
            extremes[prop]["min"] = val
          elif val > extremes[prop]["max"]:
            extremes[prop]["max"] = val
      k += 1

  for key, value in extremes.items():
    a = round(value['min'], 6)
    b = round(value['max'], 6)
    print(f"{key}: min={a:14}, max={b:14}")

  print()
  print("-"*60)

def test4():
  """ Test astate variants.
  """

  print("\nTest 4: astate variants.\n")

  pa, db, wb = 101.325, 23.0, 16.0
  air = astate(pa, db, wb)
  
  wa = air.wa
  air = astate2(pa, db, wa)
  
  ha = air.ha
  air = astate3(pa, ha, wa)
  print(air.report())
  print()
  
  print(f"Passed ✔️") if math.isclose(db, air.db) and math.isclose(wb, air.wb) else print(f"Failed ❌")
  print()
  print("-"*60)

def test5():
  """ Compare sigma heat calcs
  """

  print("\nTest 5: Test sigma heat equation.\n")

  pa = 101.325 # kPa
  db0, db1 = -60.0, 90.0 # dry bulb range [°C]
  pc0, pc1 = 0.1, 99.9 # percent saturation
  n, m = 40, 20

  all_ok = True
  for i in range(n):
    db = db0 + (db1 - db0) * i / (n - 1)
    w_max = wsat(pa, db) # kg/kg
    for j in range(m):
      # Step in percent of saturation, pc [%]
      pc = pc0 + (pc1 - pc0) * j / (m - 1)
      wa = w_max * pc / 100
      wb = wb_from_dbwa(pa, db, wa)
      za1 = sigma_heat(db, wb, wa)
      za2 = za_from_wb(pa, wb)
      if not math.isclose(za1, za2):
        all_ok = False
        print(f"za1: {za1:10.6f} za2: {za2:10.6f} kJ/kg ❌")
  
  if all_ok:
    print(f"All za tests passed ✔️\n")

  print()  
  print("-"*60)

def test7():

  print("\nTest 7: External call\n")

  # call externally
  pa, db, wb = 101.325, 23.0, 16.0
  #pa, db, wb = map(float, sys.argv[1:4])
  air = astate(pa, db, wb)
  print(air.report())
  print()
  print("-"*60)

def test8():
  """ 1. Testing the Magnus and Magnus inverse.
      2. Testing the ASHRAE inverse correlation.
  """

  print("\nTest8: Pressure Equation Tests\n")
  
  print("tsat by root find using ASHRAE psat")
  n = 20
  tmin, tmax = -60.0, 90.0
  for i in range(n):
    t = float(tmin + (tmax - tmin) * i / (n - 1))
    p = psat(t)
    t2 = tsat(p)
    if math.isclose(t, t2):
      print(f"{100*(t-t2)/t:6.3f}: passed ✔️")
    else:
      print(f"{100*(t-t2)/t:6.3f} failed ❌")
  print()
  print("-"*60)


test0()
test1()
test2()
test3()
test4()
test5()
test7()
test8()
