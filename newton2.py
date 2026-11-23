# =====================================================================
# Forward State Evaluator & 2D Newton Psychrometric Solver
# =====================================================================

def newton2(f, guess, bounds=None, scales=None, tol=1e-10, maxits=40):
  """Solve a 2D system f(x) = [0, 0] starting from guess = (x0, x1)
  using damped Newton-Raphson with finite-difference Jacobian.
  """
  x0, x1 = float(guess[0]), float(guess[1])
  s0, s1 = (1.0, 1.0) if scales is None else (float(scales[0]), float(scales[1]))

  def clamp(v0, v1):
    if bounds is not None:
      (lo0, hi0), (lo1, hi1) = bounds
      if lo0 is not None: v0 = max(lo0, v0)
      if hi0 is not None: v0 = min(hi0, v0)
      if lo1 is not None: v1 = max(lo1, v1)
      if hi1 is not None: v1 = min(hi1, v1)
    return v0, v1

  for _ in range(maxits):
    r0, r1 = f((x0, x1))
    norm = (r0 / s0) ** 2 + (r1 / s1) ** 2

    #for debugging
    #print(_, x0, x1)

    if max(abs(r0) / s0, abs(r1) / s1) < tol:
      return x0, x1

    # Numerical Jacobian via forward finite differences
    h0 = 1e-6 * max(1.0, abs(x0))
    h1 = 1e-6 * max(1.0, abs(x1))

    r_x0_0, r_x0_1 = f((x0 + h0, x1))
    r_x1_0, r_x1_1 = f((x0, x1 + h1))

    j00 = (r_x0_0 - r0) / h0
    j10 = (r_x0_1 - r1) / h0
    j01 = (r_x1_0 - r0) / h1
    j11 = (r_x1_1 - r1) / h1

    det = j00 * j11 - j01 * j10
    if abs(det) < 1e-16:
      raise ValueError("Singular Jacobian: equations are dependent or ill-conditioned!")

    # Newton step: dx = -J^-1 * r
    dx0 = -( j11 * r0 - j01 * r1) / det
    dx1 = -(-j10 * r0 + j00 * r1) / det

    # Backtracking line search within physical bounds
    alpha = 1.0
    while alpha > 1e-4:
      nx0, nx1 = clamp(x0 + alpha * dx0, x1 + alpha * dx1)
      try:
        nr0, nr1 = f((nx0, nx1))
        if (nr0 / s0) ** 2 + (nr1 / s1) ** 2 < norm:
          x0, x1 = nx0, nx1
          break
      except ValueError:
        pass
      alpha *= 0.5
    else:
      x0, x1 = nx0, nx1

  raise ValueError("newton2 failed to converge within maximum iterations!")
