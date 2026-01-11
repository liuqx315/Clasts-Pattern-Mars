"""
Pure-Python implementation of Bartlett's isotropic estimator.

- No numba
- No rpy2
- Fully NumPy 2.0 compatible
- Vectorized for speed
"""

import numpy as np
import scipy.special as sc
from scipy.spatial.distance import pdist

from structure_factor.spatial_windows import BallWindow, UnitBallWindow


# ---------------------------------------------------------------------
#  PURE PYTHON Bessel helper (scalar-safe)
# ---------------------------------------------------------------------
def bessel1(order, x):
    """Return J_order(x). Handles special cases for speed."""
    if order == 0.0:
        return sc.j0(x)
    if order == 1.0:
        return sc.j1(x)
    return sc.jv(order, x)


# ---------------------------------------------------------------------
#  PURE PYTHON vectorized version of Bartlett estimator kernel
# ---------------------------------------------------------------------
def bartlett_kernel(k_norm, norm_xi_xj, order):
    """
    Vectorized Bartlett kernel:

        2 * sum_{i<j} J_{order}(k ||xi-xj||) / (k||xi-xj||)^{order}

    Parameters
    ----------
    k_norm : array of shape (m,)
    norm_xi_xj : array of pairwise distances shape (N(N-1)/2,)
    order : float = d/2 - 1

    Returns
    -------
    result : array shape (m,)
    """

    # shape (m, n_pairs)
    k_r = np.outer(k_norm, norm_xi_xj)

    # J_{order}(k*r)
    J = bessel1(order, k_r)

    if order > 0:
        denom = np.power(k_r, order, where=(k_r > 0), out=np.zeros_like(k_r))
        J = np.divide(J, denom, out=np.zeros_like(J), where=(denom != 0))

    # per-k sum over all pairs
    return 2.0 * np.sum(J, axis=1)


# ---------------------------------------------------------------------
#  MAIN ESTIMATOR
# ---------------------------------------------------------------------
def bartlett_estimator(k_norm, point_pattern):
    """
    Bartlett's isotropic structure-factor estimator.

    Pure Python / NumPy-2.0 version (no numba).
    """
    if not isinstance(point_pattern.window, BallWindow):
        raise TypeError(
            "Window must be a BallWindow. Use point_pattern.restrict_to_window() first."
        )

    window = point_pattern.window
    d = window.dimension
    X = np.atleast_2d(point_pattern.points)

    # Pairwise distances ||xi - xj||
    norm_xi_xj = pdist(X, metric="euclidean")

    k_norm = np.asarray(k_norm, dtype=float)
    order = float(d / 2 - 1)

    # main computation
    sf_estimated = bartlett_kernel(k_norm, norm_xi_xj, order)

    # prefactor
    surface = UnitBallWindow(np.zeros(d)).surface
    volume = window.volume
    rho = float(point_pattern.intensity)

    sf_estimated *= (2.0 * np.pi) ** (d / 2) / (surface * volume * rho)
    sf_estimated += 1.0

    return sf_estimated


# ---------------------------------------------------------------------
#  ALLOWED k VALUES
# ---------------------------------------------------------------------
def allowed_k_norm_bartlett_isotropic(dimension, radius, nb_values=60):
    """
    Allowed wavenumbers for even-dimensional d.

    k_n = (zeros of J_{d/2}) / R
    """
    d2, mod = divmod(dimension, 2)
    if mod != 0:
        raise ValueError("Dimension must be even to use allowed wavenumbers.")

    zeros = sc.jn_zeros(d2, nb_values)
    return zeros / radius
