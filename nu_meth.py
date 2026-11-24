# nu_meth.py
#
# Minimal Numerical Methods Toolbox.
# (No external libraries)
#
# Includes:
#   1D root finders (Bisect, Secant, Brent, Newton)
#   Linear systems (Gaussian Elimination)
#   2D root finders (Newton and Broyden)
#   Scaled Broyden method (for moist air state)
#   1D Bounded Minimization (Brent's Method)
#   Polynomial Eval (low and high order coeff)
#   Polynomial Roots (Aberth-Ehrlich Method)
#   Interpolation (Lagrange)
#   Linear regression
#   Singular Valued Decomposition
#   Multiple Linear Regression
#   ODE Integration: Fixed-Step RK4 & Adaptive RK45
#
# Author: Bruce Wernick
# Date: 24 November 2026
# Version 0.0.0.1
#
# Copyright © (2026) Bruce Wernick,  All rights reserved...
# ---------------------------------------------------------

import math
import sys
  
# IEEE 754 float64 machine epsilon
eps = sys.float_info.epsilon # 2.220446049250313e-16

# ---------------------------------------------------------

# 1D Root Finders 
# (Bisection, Secant, Newton, Brent)

def bisect(f, a, b, tol=1e-12, maxits=100):
  fa, fb = f(a), f(b)
  if fa * fb > 0:
    raise ValueError("Root not bracketed: f(a) and f(b) must have opposite signs.")
  for its in range(1, maxits + 1):
    mid = a + 0.5 * (b - a)
    fmid = f(mid)
    if abs(fmid) <= tol or 0.5 * (b - a) <= tol:
      return mid, its, True
    if fa * fmid < 0:
      b, fb = mid, fmid
    else:
      a, fa = mid, fmid
  return 0.5 * (a + b), maxits, False

def secant(f, x0, x1, tol=1e-12, maxits=100):
  f0, f1 = f(x0), f(x1)
  for its in range(1, maxits + 1):
    if abs(f1 - f0) < 1e-15:
      return x1, its, False  # Prevent division by zero
    x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
    if abs(x2 - x1) <= tol:
      return x2, its, True
    x0, f0, x1, f1 = x1, f1, x2, f(x2)
  return x1, maxits, False

def newton_1d(f, x0, df=None, tol=1e-12, maxits=100):
  """
  Finds a root of f(x) = 0 starting from x0 using Newton-Raphson.
  df: derivative function f'(x). Uses central finite differences if None.
  Returns (root, iterations, converged).
  """
  x = float(x0)
  tiny = 1e-8  # Step size for finite-difference fallback

  for its in range(1, maxits + 1):
    fx = f(x)
    if abs(fx) <= tol:
      return x, its, True

    # Evaluate analytical slope or approximate via central difference
    if df is not None:
      dfx = df(x)
    else:
      h = tiny * max(1.0, abs(x))
      dfx = (f(x + h) - f(x - h)) / (2.0 * h)

    # Guard against stationary points / horizontal tangents
    if abs(dfx) < 1e-15:
      return x, its, False

    dx = fx / dfx
    x_next = x - dx

    if abs(dx) <= tol * max(1.0, abs(x_next)):
      return x_next, its, True

    x = x_next

  return x, maxits, False

def brent(f, a, b, tol=1e-12, maxits=100):
  """
  Finds a root of f(x) = 0 on [a, b] using Brent's method.
  Requires f(a) and f(b) to have opposite signs.
  Returns (root, iterations, converged).
  """
  fa = f(a)
  fb = f(b)
  if fa * fb > 0.0:
    raise ValueError("Root not bracketed: f(a) and f(b) must have opposite signs.")
  if abs(fa) < abs(fb):
    a, b = b, a
    fa, fb = fb, fa

  c = a
  fc = fa
  d = b - a
  e = d

  for its in range(1, maxits + 1):
    if fb == 0.0:
      return b, its, True

    # Ensure |f(b)| <= |f(c)| so b is always the best current approximation
    if (fb > 0.0 and fc > 0.0) or (fb < 0.0 and fc < 0.0):
      c = a
      fc = fa
      d = b - a
      e = d

    if abs(fc) < abs(fb):
      a, b, c = b, c, b
      fa, fb, fc = fb, fc, fb

    # Convergence check (accounts for floating-point precision limits)
    tol1 = 2.0 * eps * abs(b) + 0.5 * tol
    xm = 0.5 * (c - b)
    if abs(xm) <= tol1 or fb == 0.0:
      return b, its, True

    # Decide between interpolation (IQI / Secant) and Bisection
    if abs(e) >= tol1 and abs(fa) > abs(fb):
      s = fb / fa
      if a == c:
        # Secant step
        p = 2.0 * xm * s
        q = 1.0 - s
      else:
        # Inverse Quadratic Interpolation (IQI) step
        q = fa / fc
        r = fb / fc
        p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0))
        q = (q - 1.0) * (r - 1.0) * (s - 1.0)

      if p > 0.0:
        q = -q
      p = abs(p)

      # Accept interpolation only if it stays well inside the bracket
      min1 = 3.0 * xm * q - abs(tol1 * q)
      min2 = abs(e * q)
      if 2.0 * p < (min1 if min1 < min2 else min2):
        e = d
        d = p / q
      else:
        # Fallback to Bisection
        d = xm
        e = d
    else:
      # Fallback to Bisection
      d = xm
      e = d

    # Move previous best guess to a, and advance b
    a = b
    fa = fb
    if abs(d) > tol1:
      b += d
    else:
      b += math.copysign(tol1, xm)
    fb = f(b)

  return b, maxits, False


# ---------------------------------------------------------

# Linear Systems
# Tiny Gaussian Elimination 
# (Engine for 2D/ND Newton & Splines)

def solve_linear(A, b):
  """Solves Ax = b using Gaussian elimination 
     with partial pivoting.
  """
  n = len(b)
  M = [row[:] + [b[i]] for i, row in enumerate(A)]  # Augmented matrix

  for col in range(n):
    # Partial pivoting
    pivot = max(range(col, n), key=lambda r: abs(M[r][col]))
    M[col], M[pivot] = M[pivot], M[col]
    if abs(M[col][col]) < 1e-15:
      raise ValueError("Singular matrix.")
    
    for row in range(col + 1, n):
      factor = M[row][col] / M[col][col]
      for k in range(col, n + 1):
        M[row][k] -= factor * M[col][k]

  # Back substitution
  x = [0.0] * n
  for i in range(n - 1, -1, -1):
    s = sum(M[i][j] * x[j] for j in range(i + 1, n))
    x[i] = (M[i][n] - s) / M[i][i]
  return x

# ---------------------------------------------------------

# Newton 2D (with Finite-Difference Fallback)
# & Broyden 2D with scaled finite differences.

def approx_jac(F, x, tiny=1e-8):
  n = len(x)
  f0 = F(x)
  J = [[0.0] * n for _ in range(n)]
  for j in range(n):
    h = tiny * max(1.0, abs(x[j]))
    xh = x[:]
    xh[j] += h
    fh = F(xh)
    for i in range(n):
      J[i][j] = (fh[i] - f0[i]) / h
  return J

def newton_nd(F, x0, jac=None, tol=1e-10, maxits=50):
  x = list(x0)
  for its in range(1, maxits + 1):
    fx = F(x)
    if max(abs(v) for v in fx) <= tol:
      return x, its, True
    J = jac(x) if jac else approx_jac(F, x)
    dx = solve_linear(J, [-v for v in fx])
    x = [x[i] + dx[i] for i in range(len(x))]
    if max(abs(d) for d in dx) <= tol:
      return x, its, True
  return x, maxits, False

# Scaled Broyden method ...
# For psychrometrics and thermodynamic state solvers.

def broyden_nd(F, x0, tol=1e-10, maxits=100):
  """Damped Good Broyden method 
     with scaled finite differences.
  """
  x = list(x0)
  n = len(x)
  fx = F(x)
  norm_fx = sum(v * v for v in fx)

  # Initialize inverse Jacobian B once
  J = approx_jac(F, x)
  B = [[0.0] * n for _ in range(n)]
  for i in range(n):
    e = [1.0 if j == i else 0.0 for j in range(n)]
    col = solve_linear(J, e)
    for r in range(n):
      B[r][i] = col[r]

  for its in range(1, maxits + 1):
    if max(abs(v) for v in fx) <= tol:
      return x, its, True

    dx = [-sum(B[i][j] * fx[j] for j in range(n)) for i in range(n)]

    # Backtracking line search to prevent non-physical overshoots
    alpha = 1.0
    for _ in range(10):
      step = [alpha * d for d in dx]
      x_new = [x[i] + step[i] for i in range(n)]
      try:
        fx_new = F(x_new)
        norm_new = sum(v * v for v in fx_new)
        if norm_new < norm_fx or alpha < 0.05:
          break
      except (ValueError, OverflowError):
        pass  # Shrink step if correlation went out of physical range
      alpha *= 0.5
    else:
      return x, its, False

    df = [fx_new[i] - fx[i] for i in range(n)]

    # Sherman-Morrison rank-1 update using actual damped step
    B_df = [sum(B[i][j] * df[j] for j in range(n)) for i in range(n)]
    denom = sum(step[i] * B_df[i] for i in range(n))
    if abs(denom) < 1e-15:
      return x_new, its, False

    u = [(step[i] - B_df[i]) / denom for i in range(n)]
    for i in range(n):
      for j in range(n):
        B[i][j] += u[i] * sum(step[k] * B[k][j] for k in range(n))

    x, fx, norm_fx = x_new, fx_new, norm_new

  return x, maxits, False

# ---------------------------------------------------------

# 1D Bounded Minimization (Brent's Method)

def brent_min(f, a, b, tol=1e-8, maxits=100):
  """Finds a local minimum of f(x) on the interval [a, b] 
     using Brent's method
     (Golden Section Search + Parabolic Interpolation).
     Zero derivatives needed.
     Returns (x_min, f_min, iterations, converged).
  """
  if a > b:
    a, b = b, a

  golden = 0.5 * (3.0 - math.sqrt(5.0))  # ~0.38196601125010515
  sqrt_eps = math.sqrt(eps)

  x = w = v = a + golden * (b - a)
  fx = fw = fv = f(x)
  d = 0.0
  e = 0.0

  for its in range(1, maxits + 1):
    xm = 0.5 * (a + b)
    tol1 = sqrt_eps * abs(x) + tol / 3.0
    tol2 = 2.0 * tol1

    # Check if bracket width is within tolerance
    if abs(x - xm) <= (tol2 - 0.5 * (b - a)):
      return x, fx, its, True

    use_golden = True
    if abs(e) > tol1:
      # Fit parabola through (v, fv), (w, fw), (x, fx)
      r = (x - w) * (fx - fv)
      q = (x - v) * (fx - fw)
      p = (x - v) * q - (x - w) * r
      q = 2.0 * (q - r)
      if q > 0.0:
        p = -p
      q = abs(q)
      etemp = e
      e = d

      # Accept parabolic step only if it falls inside (a, b) and shrinks bracket
      if abs(p) < abs(0.5 * q * etemp) and p > q * (a - x) and p < q * (b - x):
        d = p / q
        u = x + d
        if (u - a) < tol2 or (b - u) < tol2:
          d = math.copysign(tol1, xm - x)
        use_golden = False

    if use_golden:
      e = (b - x) if x < xm else (a - x)
      d = golden * e

    # Ensure step is at least tol1 to avoid evaluating identical floats
    u = x + (d if abs(d) >= tol1 else math.copysign(tol1, d))
    fu = f(u)

    # Update bracket [a, b] and three best points (x, w, v)
    if fu <= fx:
      if u < x:
        b = x
      else:
        a = x
      v, fv = w, fw
      w, fw = x, fx
      x, fx = u, fu
    else:
      if u < x:
        a = u
      else:
        b = u
      if fu <= fw or w == x:
        v, fv = w, fw
        w, fw = u, fu
      elif fu <= fv or v == x or v == w:
        v, fv = u, fu

  return x, fx, maxits, False

# ---------------------------------------------------------

# Polynomial Evaluation

def poly_eval(coeffs, x):
  """Evaluates polynomial
     P(x) = coeffs[0]*x^n + ... + coeffs[-1]
     using Horner's method.
  """
  res = 0
  for c in coeffs:
    res = res * x + c
  return res

def poly_eval_deriv(coeffs, x):
  """Returns (P(x), P'(x)) in a single O(n) pass.
  """
  p = 0
  dp = 0
  for c in coeffs:
    dp = dp * x + p
    p = p * x + c
  return p, dp

def poly_eval_asc(coeffs, x):
  """Polynomial (ascending coefficient)
     Evaluates 
     P(x) = coeffs[0] + coeffs[1]*x + ... + coeffs[n]*x^n.
  """
  res = 0
  for c in reversed(coeffs):
    res = res * x + c
  return res

def poly_eval_deriv_asc(coeffs, x):
  """Returns (P(x), P'(x)) 
     for ascending coefficients [c_0, c_1, ..., c_n].
  """
  p = 0
  dp = 0
  for c in reversed(coeffs):
    dp = dp * x + p
    p = p * x + c
  return p, dp

# ---------------------------------------------------------

# Polynomial Roots (Aberth-Ehrlich Method)

def poly_roots(coeffs, tol=1e-12, max_iter=200):
  """Finds all roots of 
     P(x) = coeffs[0]*x^n + ... + coeffs[n]
     using theAberth-Ehrlich method.
  """
  # Strip leading zeros
  while len(coeffs) > 1 and coeffs[0] == 0:
    coeffs = coeffs[1:]
  n = len(coeffs) - 1
  if n <= 0:
    return []
  
  # Normalize monic polynomial
  c = [complex(a / coeffs[0]) for a in coeffs]
  
  # Initial guesses on a circle with slight asymmetry
  radius = 1.0 + max(abs(a) for a in c[1:])
  roots = [radius * complex(math.cos(2 * math.pi * k / n + 0.1), 
    math.sin(2 * math.pi * k / n + 0.1)) for k in range(n)]

  for _ in range(max_iter):
    max_step = 0.0
    for i in range(n):
      z = roots[i]
      # Horner's method for P(z) and P'(z)
      pz, dpz = c[0], 0j
      for a in c[1:]:
        dpz = dpz * z + pz
        pz = pz * z + a
      if abs(pz) <= tol:
        continue
      newton_step = pz / dpz
      sum_inv = sum(1.0 / (z - roots[j]) for j in range(n) if j != i)
      w = newton_step / (1.0 - newton_step * sum_inv)
      roots[i] -= w
      max_step = max(max_step, abs(w))
    if max_step <= tol:
      break
  return roots

# ---------------------------------------------------------

# Interpolation

def lagrange_interp(xs, ys):
  """Constructs a Barycentric Lagrange interpolant 
     through (xs, ys).
     Returns a callable function P(x).
  """
  n = len(xs)
  if n != len(ys) or n == 0:
    raise ValueError("xs and ys must have the same non-zero length.")

  # Precompute barycentric weights w_j = 1 / prod_{k != j} (x_j - x_k)
  weights = [1.0] * n
  for j in range(n):
    xj = xs[j]
    prod = 1.0
    for k in range(n):
      if k != j:
        diff = xj - xs[k]
        if diff == 0.0:
          raise ValueError(f"Duplicate x value detected at index {j}: {xj}")
        prod *= diff
    weights[j] = 1.0 / prod

  # Copy lists into closure to prevent external mutation bugs
  x_nodes = list(xs)
  y_nodes = list(ys)

  def interpolant(x):
    num = 0.0
    den = 0.0
    for j in range(n):
      diff = x - x_nodes[j]
      # Exact node hit: return y_j immediately to avoid 0/0 division
      if diff == 0.0:
        return y_nodes[j]
      term = weights[j] / diff
      num += term * y_nodes[j]
      den += term
    return num / den

  return interpolant

# ---------------------------------------------------------

# 2D Bilinear Interpolation 
# (Catalog / Compressor / Coil Maps)

def interp_2d(xs, ys, grid):
  """Bilinear interpolation over a 2D rectangular grid.
     xs: sorted list of x-coordinates (length nx)
     ys: sorted list of y-coordinates (length ny)
     grid: 2D list of shape [ny][nx] where grid[j][i] = f(xs[i], ys[j])
     Returns a callable f(x, y) with edge clamping.
  """
  nx, ny = len(xs), len(ys)
  if len(grid) != ny or any(len(row) != nx for row in grid):
    raise ValueError("Grid dimensions must match [len(ys)][len(xs)].")

  def _find_cell(arr, val):
    if val <= arr[0]:
      return 0, 0.0
    if val >= arr[-1]:
      return len(arr) - 2, 1.0
    lo, hi = 0, len(arr) - 1
    while hi - lo > 1:
      mid = (lo + hi) // 2
      if arr[mid] <= val:
        lo = mid
      else:
        hi = mid
    t = (val - arr[lo]) / (arr[lo + 1] - arr[lo])
    return lo, t

  def interpolant(x, y):
    i, tx = _find_cell(xs, x)
    j, ty = _find_cell(ys, y)
    z00 = grid[j][i]
    z10 = grid[j][i + 1]
    z01 = grid[j + 1][i]
    z11 = grid[j + 1][i + 1]
    return (1.0 - tx) * (1.0 - ty) * z00 + tx * (1.0 - ty) * z10 + \
           (1.0 - tx) * ty * z01 + tx * ty * z11

  return interpolant

# ---------------------------------------------------------

# Natural Cubic Spline (Smooth Refrigerant / Steam Tables)

def cubic_spline(xs, ys):
  """Constructs a Natural Cubic Spline through sorted (xs, ys).
     Solved in O(n) via the Thomas tridiagonal algorithm.
     Returns a callable S(x).
  """
  n = len(xs)
  if n != len(ys) or n < 2:
    raise ValueError("xs and ys must have matching length >= 2.")

  h = [xs[i + 1] - xs[i] for i in range(n - 1)]
  if any(step <= 0.0 for step in h):
    raise ValueError("xs must be strictly increasing.")

  # Set up tridiagonal system for second derivatives M[0..n-1]
  alpha = [0.0] * n
  for i in range(1, n - 1):
    alpha[i] = (3.0 / h[i]) * (ys[i + 1] - ys[i]) - (3.0 / h[i - 1]) * (ys[i] - ys[i - 1])

  l = [1.0] + [0.0] * (n - 1)
  mu = [0.0] * n
  z = [0.0] * n

  for i in range(1, n - 1):
    l[i] = 2.0 * (xs[i + 1] - xs[i - 1]) - h[i - 1] * mu[i - 1]
    mu[i] = h[i] / l[i]
    z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i]

  l[n - 1] = 1.0
  c = [0.0] * n
  b = [0.0] * (n - 1)
  d = [0.0] * (n - 1)

  for j in range(n - 2, -1, -1):
    c[j] = z[j] - mu[j] * c[j + 1]
    b[j] = (ys[j + 1] - ys[j]) / h[j] - h[j] * (c[j + 1] + 2.0 * c[j]) / 3.0
    d[j] = (c[j + 1] - c[j]) / (3.0 * h[j])

  x_nodes = list(xs)
  y_nodes = list(ys)

  def spline(x):
    if x <= x_nodes[0]:
      i = 0
    elif x >= x_nodes[-1]:
      i = n - 2
    else:
      lo, hi = 0, n - 1
      while hi - lo > 1:
        mid = (lo + hi) // 2
        if x_nodes[mid] <= x:
          lo = mid
        else:
          hi = mid
      i = lo
    dx = x - x_nodes[i]
    return y_nodes[i] + b[i] * dx + c[i] * dx * dx + d[i] * dx * dx * dx

  return spline

# ---------------------------------------------------------

# Numerical Integration (Adaptive Simpson's Quadrature)

def integrate_simpson(f, a, b, tol=1e-10, max_depth=20):
  """Computes definite integral of f(x) 
     from a to b.
     Uses Adaptive Simpson's Rule.
  """
  def _simpson(fa, fm, fb, step):
    return (step / 6.0) * (fa + 4.0 * fm + fb)

  def _recurse(a, b, fa, fm, fb, whole, tol, depth):
    mid = 0.5 * (a + b)
    lm = 0.5 * (a + mid)
    rm = 0.5 * (mid + b)
    flm = f(lm)
    frm = f(rm)
    left = _simpson(fa, flm, fm, mid - a)
    right = _simpson(fm, frm, fb, b - mid)
    delta = left + right - whole
    if depth <= 0 or abs(delta) <= 15.0 * tol:
      return left + right + delta / 15.0
    return (_recurse(a, mid, fa, flm, fm, left, 0.5 * tol, depth - 1) +
            _recurse(mid, b, fm, frm, fb, right, 0.5 * tol, depth - 1))

  if a == b:
    return 0.0
  sign = 1.0
  if b < a:
    a, b = b, a
    sign = -1.0
  fa, fm, fb = f(a), f(0.5 * (a + b)), f(b)
  whole = _simpson(fa, fm, fb, b - a)
  return sign * _recurse(a, b, fa, fm, fb, whole, tol, max_depth)


# ---------------------------------------------------------

# Linear Least Squares Regression

def linear_fit(xs, ys):
  """Fits y = m*x + b via Ordinary Least Squares.
     Returns (predict_fn, (m, b), r_squared).
  """
  n = len(xs)
  if n != len(ys) or n < 2:
    raise ValueError("xs and ys must have equal length >= 2.")

  mean_x = sum(xs) / n
  mean_y = sum(ys) / n

  s_xx = 0.0
  s_yy = 0.0
  s_xy = 0.0
  for x, y in zip(xs, ys):
    dx = x - mean_x
    dy = y - mean_y
    s_xx += dx * dx
    s_yy += dy * dy
    s_xy += dx * dy

  if s_xx == 0.0:
    raise ValueError("All x values are identical; slope is undefined.")

  m = s_xy / s_xx
  b = mean_y - m * mean_x
  r_squared = (s_xy * s_xy) / (s_xx * s_yy) if s_yy > 0.0 else 1.0

  def model(x):
      return m * x + b

  return model, (m, b), r_squared


# ---------------------------------------------------------

def svd_thin(A, tol=1e-12, max_sweeps=100):
  """Computes Thin SVD (A = U * diag(S) * Vt) 
     for an m x n matrix (m >= n)
     using One-Sided Jacobi rotations.
     Zero external dependencies.
  """
  m = len(A)
  n = len(A[0])
  if m < n:
    raise ValueError("Matrix must have rows >= columns (m >= n).")

  # Work on column-major copies for fast dot products and rotations
  U_cols = [[float(A[r][c]) for r in range(m)] for c in range(n)]
  V_cols = [[1.0 if r == c else 0.0 for r in range(n)] for c in range(n)]

  for _ in range(max_sweeps):
    converged = True
    for i in range(n - 1):
      for j in range(i + 1, n):
        ui, uj = U_cols[i], U_cols[j]
        
        # Compute 2x2 Gram submatrix entries
        a = sum(x * x for x in ui)
        b = sum(y * y for y in uj)
        c = sum(x * y for x, y in zip(ui, uj))

        if abs(c) <= tol * math.sqrt(a * b):
          continue
        converged = False

        # Jacobi rotation angle to orthogonalize columns i and j
        zeta = (b - a) / (2.0 * c)
        t = math.copysign(1.0, zeta) / (abs(zeta) + math.sqrt(1.0 + zeta * zeta))
        cos = 1.0 / math.sqrt(1.0 + t * t)
        sin = cos * t

        # Apply rotation to U_cols and V_cols
        for r in range(m):
          u_ir, u_jr = ui[r], uj[r]
          ui[r] = cos * u_ir - sin * u_jr
          uj[r] = sin * u_ir + cos * u_jr

        vi, vj = V_cols[i], V_cols[j]
        for r in range(n):
          v_ir, v_jr = vi[r], vj[r]
          vi[r] = cos * v_ir - sin * v_jr
          vj[r] = sin * v_ir + cos * v_jr

    if converged:
      break

  # Extract singular values (norms of U columns) and normalize U
  S = []
  for j in range(n):
    norm = math.sqrt(sum(x * x for x in U_cols[j]))
    S.append(norm)
    if norm > 0.0:
      U_cols[j] = [x / norm for x in U_cols[j]]

  # Sort singular values descending
  order = sorted(range(n), key=lambda k: S[k], reverse=True)
  S = [S[k] for k in order]
  U_cols = [U_cols[k] for k in order]
  V_cols = [V_cols[k] for k in order]

  # Convert U_cols back to row-major (m x n); V_cols is already Vt (n x n)
  U = [[U_cols[c][r] for c in range(n)] for r in range(m)]
  Vt = V_cols
  return U, S, Vt

def svd_multiple_regression(X, y, rcond=1e-12):
  """Fits y = b0 + b1*x1 + ... + bk*xk
     using SVD pseudoinverse.
     X: list of rows [[x1_1, x1_2, ...], [x2_1, x2_2, ...], ...]
     y: list of target values
     Returns (predict_fn, coeffs, R, R_squared, singular_values).
  """
  m = len(X)
  if m != len(y) or m == 0:
    raise ValueError("X and y must have matching non-zero length.")

  # Prepend column of 1.0s for the intercept b0
  A = [[1.0] + list(row) for row in X]
  n = len(A[0])

  U, S, Vt = svd_thin(A)

  # Threshold small singular values to filter out collinear noise
  cutoff = rcond * S[0] if S else 0.0

  # Compute coeffs = V * diag(1/S) * U^T * y
  # Step 1: ut_y = U^T * y (length n)
  ut_y = [sum(U[r][c] * y[r] for r in range(m)) for c in range(n)]

  # Step 2: w = diag(1/S) * ut_y (length n)
  w = [(ut_y[i] / S[i]) if S[i] > cutoff else 0.0 for i in range(n)]

  # Step 3: coeffs = V * w (since Vt rows are V columns)
  coeffs = [sum(Vt[r][c] * w[r] for r in range(n)) for c in range(n)]

  # Compute R^2 and Multiple Correlation R = sqrt(R^2)
  mean_y = sum(y) / m
  ss_tot = sum((yi - mean_y) ** 2 for yi in y)
  ss_res = sum((y[i] - sum(coeffs[j] * A[i][j] for j in range(n))) ** 2 for i in range(m))
  
  r_squared = max(0.0, 1.0 - (ss_res / ss_tot)) if ss_tot > 0.0 else 1.0
  multiple_r = math.sqrt(r_squared)

  def predict(x_row):
    return coeffs[0] + sum(c * val for c, val in zip(coeffs[1:], x_row))

  return predict, coeffs, multiple_r, r_squared, S

# ---------------------------------------------------------

# ODE Integration: Fixed-Step RK4 
# & Adaptive RK45 (Dormand-Prince)

def rk4(f, t_span, y0, steps=100):
  """Solves dy/dt = f(t, y) on [t0, tf] 
     using fixed-step 4th-order Runge-Kutta.
     y0 can be a float (scalar ODE) 
     or a list of floats (system of ODEs).
     Returns (ts, ys).
  """
  t0, tf = float(t_span[0]), float(t_span[1])
  h = (tf - t0) / steps
  scalar = isinstance(y0, (int, float))
  y = [float(y0)] if scalar else [float(v) for v in y0]
  n = len(y)

  ts = [t0]
  ys = [y[0] if scalar else y[:]]
  t = t0

  def _eval(ti, yi):
    res = f(ti, yi[0] if scalar else yi)
    return [float(res)] if scalar else res

  for _ in range(steps):
    k1 = _eval(t, y)
    y2 = [y[i] + 0.5 * h * k1[i] for i in range(n)]
    k2 = _eval(t + 0.5 * h, y2)
    y3 = [y[i] + 0.5 * h * k2[i] for i in range(n)]
    k3 = _eval(t + 0.5 * h, y3)
    y4 = [y[i] + h * k3[i] for i in range(n)]
    k4 = _eval(t + h, y4)

    y = [y[i] + (h / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) for i in range(n)]
    t += h
    ts.append(t)
    ys.append(y[0] if scalar else y[:])

  return ts, ys

def rk45(f, t_span, y0, rtol=1e-6, atol=1e-9, max_steps=10000):
  """Solves dy/dt = f(t, y) on [t0, tf] 
     using adaptive Dormand-Prince RK45.
     y0 can be a float (scalar ODE) 
     or a list of floats (system of ODEs).
     Returns (ts, ys, converged).
  """
  t0, tf = float(t_span[0]), float(t_span[1])
  direction = 1.0 if tf >= t0 else -1.0
  span = abs(tf - t0)
  if span == 0.0:
    return [t0], [y0], True

  scalar = isinstance(y0, (int, float))
  y = [float(y0)] if scalar else [float(v) for v in y0]
  n = len(y)

  def _eval(ti, yi):
    res = f(ti, yi[0] if scalar else yi)
    return [float(res)] if scalar else res

  # Dormand-Prince RK45 Butcher tableau coefficients
  c2, c3, c4, c5, c6 = 1/5, 3/10, 4/5, 8/9, 1.0
  a21 = 1/5
  a31, a32 = 3/40, 9/40
  a41, a42, a43 = 44/45, -56/15, 32/9
  a51, a52, a53, a54 = 19372/6561, -25360/2187, 64448/6561, -212/729
  a61, a62, a63, a64, a65 = 9017/3168, -355/33, 46732/5247, 49/176, -5103/18656

  # 5th-order solution weights
  b1, b3, b4, b5, b6 = 35/384, 500/1113, 125/192, -2187/6784, 11/84
  # Error coefficients (b_i - b*_i, difference between 5th and 4th order)
  e1 = 35/384 - 5179/57600
  e3 = 500/1113 - 7571/16695
  e4 = 125/192 - 393/640
  e5 = -2187/6784 - (-92097/339200)
  e6 = 11/84 - 187/2100
  e7 = -1/40

  t = t0
  h = direction * (span / 100.0)
  ts = [t]
  ys = [y[0] if scalar else y[:]]
  k1 = _eval(t, y)

  for _ in range(max_steps):
    if direction * (t - tf) >= -1e-15:
      return ts, ys, True

    # Clamp step so we land squarely on tf
    if direction * (t + h - tf) > 0.0:
      h = tf - t

    y2 = [y[i] + h * a21 * k1[i] for i in range(n)]
    k2 = _eval(t + c2 * h, y2)

    y3 = [y[i] + h * (a31 * k1[i] + a32 * k2[i]) for i in range(n)]
    k3 = _eval(t + c3 * h, y3)

    y4 = [y[i] + h * (a41 * k1[i] + a42 * k2[i] + a43 * k3[i]) for i in range(n)]
    k4 = _eval(t + c4 * h, y4)

    y5 = [y[i] + h * (a51 * k1[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i]) for i in range(n)]
    k5 = _eval(t + c5 * h, y5)

    y6 = [y[i] + h * (a61 * k1[i] + a62 * k2[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i]) for i in range(n)]
    k6 = _eval(t + c6 * h, y6)

    y_new = [y[i] + h * (b1 * k1[i] + b3 * k3[i] + b4 * k4[i] + b5 * k5[i] + b6 * k6[i]) for i in range(n)]
    k7 = _eval(t + h, y_new)  # FSAL (First Same As Last)

    # RMS scaled error norm
    err_sq = 0.0
    for i in range(n):
      err_i = h * (e1 * k1[i] + e3 * k3[i] + e4 * k4[i] + e5 * k5[i] + e6 * k6[i] + e7 * k7[i])
      scale = atol + rtol * max(abs(y[i]), abs(y_new[i]))
      err_sq += (err_i / scale) ** 2
    err_norm = math.sqrt(err_sq / n)

    if err_norm <= 1.0:
      # Accept step
      t += h
      y = y_new
      k1 = k7  # Reuse k7 as k1 for the next step!
      ts.append(t)
      ys.append(y[0] if scalar else y[:])

    # Adjust step size for next attempt
    factor = 5.0 if err_norm == 0.0 else 0.9 * (err_norm ** -0.2)
    factor = max(0.2, min(5.0, factor))
    h *= factor

  return ts, ys, False

# ---------------------------------------------------------

# HVAC and Engineering solutions

# Colebrook Friction factor 

def colebrook(Re, k, D):
  """Darcy friction factor via Colebrook-White 
     using 1D Newton-Raphson.
  """
  if Re <= 0.0:
    raise ValueError("Reynolds number must be positive.")
  if Re < 2300.0:
    return 64.0 / Re  # Laminar flow

  C0 = (k / D) / 3.7
  C1 = 2.51 / Re
  ln10 = math.log(10.0)

  g = lambda x: x + 2.0 * math.log10(C0 + x * C1)
  dg = lambda x: 1.0 + (2.0 * C1) / (ln10 * (C0 + x * C1))

  # First guess using Haaland equation
  a0 = -1.8 * math.log10((C0)**1.11 + 6.9 / Re)

  a, _, _ = newton_1d(g, a0, df=dg, tol=1e-9)
  return 1.0 / (a * a)



# ---------------------------------------------------------

if __name__ == "__main__":

  #g = lambda x: poly_eval_asc([400, 4, 0.4, 0.04], x)
  g = lambda x: poly_eval([0.04, 0.4, 4, 400], x)
  xs = [1,2,3,4,5]
  ys = [g(x) for x in xs]
  f = lagrange_interp(xs, ys)
  x = 2.5
  print(f(x), g(x))

  f, (m, b), r2  = linear_fit(xs, ys)
  print(f"R² = {r2:0.4f}")
  for x in xs:
    print(f"{x:2d}, {g(x):6.2f}, {f(x):6.2f}")


  # Data where x2 is nearly 2 * x1 (severe multicollinearity)
  X = [
      [1.0, 2.01],
      [2.0, 3.99],
      [3.0, 6.02],
      [4.0, 7.98],
      [5.0, 10.01],
  ]
  y = [5.1, 8.9, 13.2, 16.8, 21.1]

  predict, coeffs, R, R2, S = svd_multiple_regression(X, y)

  print(f"Intercept & Slopes: {[round(c, 4) for c in coeffs]}")
  print(f"Multiple Correlation R: {R:.6f}, R²: {R2:.6f}")
  print(f"Singular Values: {[round(s, 4) for s in S]}")  # Notice tiny 3rd singular value!
  print(f"Prediction at [6.0, 12.0]: {predict([6.0, 12.0]):.4f}")
  print()

  # Tricky function with flat regions and steep slopes: f(x) = (x + 3) * (x - 1)^2 ... 
  # Let's use a simple odd-multiplicity root: f(x) = x^3 - 2*x - 5
  f = lambda x: x**3 - 2*x - 5
  df = lambda x: 3*x**2 - 2

  root, its, ok = brent(f, 1.0, 3.0)
  print(f"Brent: {root:.15f} in {its} iterations")

  root, its, ok = secant(f, 1.0, 3.0)
  print(f"Secant: {root:.15f} in {its} iterations")

  root, its, ok = bisect(f, 1.0, 3.0)
  print(f"Bisect: {root:.15f} in {its} iterations")

  # 1. With analytical slope (quadratic convergence: ~4-5 iterations)
  root, its, ok = newton_1d(f, x0=1.0, df=df)
  print(f"Analytical slope: {root:.15f} in {its} iterations (converged={ok})")  

  # 2. Without analytical slope (central finite-difference fallback)
  root, its, ok = newton_1d(f, x0=1.0)
  print(f"Approx slope: {root:.15f} in {its} iterations (converged={ok})")
  print()

  def F(pt):
    x, y = pt
    return [x**2 + y**2 - 4.0, x * y - 1.0]

  def J(pt):
    x, y = pt
    return [[2.0 * x, 2.0 * y], [y, x]]

  x0 = [1.8, 0.4]  # Initial guess near one intersection

  sol_a, it_a, ok_a = newton_nd(F, x0, jac=J)
  sol_b, it_b, ok_b = newton_nd(F, x0, jac=None)
  sol_c, it_c, ok_c = broyden_nd(F, x0)
  sol_d, it_d, ok_d = broyden_nd(F, x0)
  exact_x = math.sqrt(2.0 + math.sqrt(3.0))
  exact_y = math.sqrt(2.0 - math.sqrt(3.0))
  print(f"Exact Solution: ({exact_x:.6f}, {exact_y:.6f})")
  print(f"Newton (Jc): ({sol_a[0]:.6f}, {sol_a[1]:.6f}) its={it_a}")
  print(f"Newton (FD): ({sol_b[0]:.6f}, {sol_b[1]:.6f}) its={it_b}")
  print(f"Broydn (2D): ({sol_c[0]:.6f}, {sol_c[1]:.6f}) its={it_c}")
  print()


  # Classic thermodynamic rule of thumb: P_int ≈ sqrt(P_evap * P_cond)
  # Let's minimize total two-stage isentropic compressor work for R-134a-like properties
  p_evap = 200.0   # kPa
  p_cond = 1200.0  # kPa
  gamma = 1.12     # Isentropic exponent

  def total_work(p_int):
    exp = (gamma - 1.0) / gamma
    w_low = (p_int / p_evap) ** exp - 1.0
    w_high = (p_cond / p_int) ** exp - 1.0
    return w_low + w_high

  p_opt, w_min, its, ok = brent_min(total_work, p_evap, p_cond)
  print(f"Optimal P_int: {p_opt:.4f} kPa (exact geometric mean: {math.sqrt(p_evap*p_cond):.4f}) in {its} its")

  # Friction factor
  # Example: Water in a 0.1 m commercial steel pipe  
  Re = 1e5
  e = 4.5e-5  # m
  D = 100e-3  # m
  f = colebrook(Re, e, D)
  print(f"D'Arcy friction factor = {f:0.6f}")


  # Parameters (SI units: Watts, Joules/K, W/K)
  C_air = 1.5e5    # Room air thermal capacitance
  C_mass = 2.0e6   # Walls/furniture thermal capacitance
  UA_env = 120.0   # Envelope conductance to outdoor air (W/K)
  UA_int = 450.0   # Convective conductance between air and thermal mass (W/K)
  T_out = 35.0     # Outdoor temperature (°C)
  Q_cool = 5000.0  # Cooling coil capacity (W)

  def room_ode(t, state):
    T_air, T_mass = state
    dT_air = (UA_env * (T_out - T_air) + UA_int * (T_mass - T_air) - Q_cool) / C_air
    dT_mass = (UA_int * (T_air - T_mass)) / C_mass
    return [dT_air, dT_mass]

  # Simulate 2 hours (7200 seconds) starting at 30°C
  ts, ys, ok = rk45(room_ode, (0.0, 7200.0), [30.0, 30.0])

  print(f"Converged: {ok} in {len(ts) - 1} adaptive steps")
  print(f"After 2 hours -> T_air: {ys[-1][0]:.2f}°C, T_mass: {ys[-1][1]:.2f}°C")
