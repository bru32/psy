# Brent method root finder

import sys

eps = sys.float_info.epsilon

def brent(f, guess, xtol=1e-9, maxits=100):
  """Pure-Python Brent's method root-finder 
  on bracket [a, b].
  """
  a, b = guess
  fa = f(a)
  fb = f(b)
  if fa == 0.0:
    return a
  if fb == 0.0:
    return b
  if fa * fb > 0.0:
    raise ValueError("Root is not bracketed on [a, b].")

  c, fc = a, fa
  d = e = b - a

  for _ in range(maxits):
    
    #print(f"{_:3d}: {b:12.8f}, {fb:12.8f}")

    if fb * fc > 0.0:
      c, fc = a, fa
      d = e = b - a
    if abs(fc) < abs(fb):
      a, b, c = b, c, b
      fa, fb, fc = fb, fc, fb

    tol = 2.0 * eps * abs(b) + 0.5 * xtol
    xm = 0.5 * (c - b)

    if abs(xm) <= tol or fb == 0.0:
      return b

    if abs(e) >= tol and abs(fa) > abs(fb):
      s = fb / fa
      if a == c:
        p_step = 2.0 * xm * s
        q_step = 1.0 - s
      else:
        q_step = fa / fc
        r_step = fb / fc
        p_step = s * (2.0 * xm * q_step * (q_step - r_step) - (b - a) * (r_step - 1.0))
        q_step = (q_step - 1.0) * (r_step - 1.0) * (s - 1.0)
      if p_step > 0.0:
        q_step = -q_step
      p_step = abs(p_step)

      min1 = 3.0 * xm * q_step - abs(tol * q_step)
      min2 = abs(e * q_step)
      if 2.0 * p_step < (min1 if min1 < min2 else min2):
        e = d
        d = p_step / q_step
      else:
        d = xm
        e = d
    else:
      d = xm
      e = d

    a, fa = b, fb
    if abs(d) > tol:
      b += d
    else:
      b += tol if xm > 0.0 else -tol
    fb = f(b)

  return b

# ---------------------------------------------------------------------

if __name__ == "__main__":

  # Test with "narrow_well" function.
  from math import exp
  f = lambda x: 1-2*exp(-2*(x-3)**2)
  print(brent(f,(2.4,2.5)))
  print(brent(f,(3.5,3.6)))
