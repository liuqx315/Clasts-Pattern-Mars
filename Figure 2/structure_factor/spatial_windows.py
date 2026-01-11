"""
Collection of classes representing observation windows (box, ball, etc).

- BallWindow: D-dimensional Euclidean ball
- BoxWindow:  D-dimensional hyperrectangle

R interface and spatstat dependencies removed.
Compatible with NumPy >= 2.0.
"""

from abc import ABCMeta, abstractmethod
import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
from matplotlib.patches import Circle, Rectangle

from structure_factor.utils import get_random_number_generator


# -------------------------------------------------------------------------
#  ABSTRACT BASE CLASS
# -------------------------------------------------------------------------
class AbstractSpatialWindow(metaclass=ABCMeta):
    """Encapsulate the notion of spatial window in R^d."""

    @property
    @abstractmethod
    def dimension(self):
        """Ambient dimension."""

    @property
    @abstractmethod
    def volume(self):
        """Volume of window."""

    @abstractmethod
    def __contains__(self, point):
        """True if point in window."""

    def indicator_function(self, points):
        """
        Given a single point (d-vector) or array (n×d), return indicator mask.
        """
        points = np.asarray(points)

        if points.ndim == 1 and points.size == self.dimension:
            return bool(self.__contains__(points))

        return np.apply_along_axis(self.__contains__, axis=1, arr=points)

    @abstractmethod
    def rand(self, n=1, seed=None):
        """Generate `n` uniform random points inside window."""

    # ---- NEW: pure-Python export helper ----
    def to_dict(self):
        """Return a Python representation (no R dependency)."""
        return {
            "dimension": self.dimension,
            "volume": float(self.volume),
            "type": self.__class__.__name__,
        }


# -------------------------------------------------------------------------
#  BALL WINDOW
# -------------------------------------------------------------------------
class BallWindow(AbstractSpatialWindow):
    """
    d-dimensional closed ball B(center, radius).
    """

    def __init__(self, center, radius=1.0):
        center = np.asarray(center)
        if center.ndim != 1:
            raise ValueError("center must be a 1D array")
        if radius <= 0:
            raise ValueError("radius must be positive")

        self.center = center
        self.radius = float(radius)

    @property
    def dimension(self):
        return self.center.size

    @property
    def surface(self):
        d, r = self.dimension, self.radius
        if d == 1:
            return 0.0
        if d == 2:
            return 2 * np.pi * r
        if d == 3:
            return 4 * np.pi * r**2

        return 2 * np.pi ** (d / 2) * r ** (d - 1) / sp.special.gamma(d / 2)

    @property
    def volume(self):
        d, r = self.dimension, self.radius
        if d == 1:
            return 2 * r
        if d == 2:
            return np.pi * r**2
        if d == 3:
            return 4 / 3 * np.pi * r**3

        return np.pi ** (d / 2) * r**d / sp.special.gamma(d / 2 + 1)

    def __contains__(self, point):
        point = np.asarray(point)
        if point.ndim != 1 or point.size != self.dimension:
            raise ValueError("point must be 1D coordinate vector")
        return np.linalg.norm(point - self.center) <= self.radius

    def indicator_function(self, points):
        points = np.asarray(points)
        return np.linalg.norm(points - self.center, axis=-1) <= self.radius

    def rand(self, n=1, seed=None):
        """
        Sample uniformly from the ball using 'dropped coordinates' method.
        """
        rng = get_random_number_generator(seed)
        d = self.dimension
        # sample on sphere in (d+2) dims → project down
        pts = rng.standard_normal((n, d + 2))
        pts /= np.linalg.norm(pts, axis=-1, keepdims=True)
        pts = pts[:, :d]  # drop extra coords
        return self.center + self.radius * pts

    # Removed: R conversion (to_spatstat_owin)

    def plot(self, axis=None, **kwargs):
        """Plot circle in 2D."""
        if self.dimension != 2:
            raise NotImplementedError("Plotting only implemented for 2D")

        if axis is None:
            _, axis = plt.subplots(figsize=(5, 5))

        kwargs.setdefault("fill", False)
        circ = Circle(self.center, self.radius, **kwargs)
        axis.add_patch(circ)
        axis.set_aspect("equal")
        return axis


class UnitBallWindow(BallWindow):
    """Ball with radius = 1."""

    def __init__(self, center):
        super().__init__(center, radius=1.0)


# -------------------------------------------------------------------------
#  BOX WINDOW
# -------------------------------------------------------------------------
class BoxWindow(AbstractSpatialWindow):
    """d-dimensional axis-aligned box ∏ᵢ [aᵢ, bᵢ]."""

    def __init__(self, bounds):
        bounds = np.asarray(bounds)
        bounds = np.atleast_2d(bounds)

        if bounds.ndim != 2 or bounds.shape[1] != 2:
            raise ValueError("bounds must be (d × 2) array")

        if np.any(np.diff(bounds, axis=1) <= 0):
            raise ValueError("Each row must satisfy a_i < b_i")

        # Store as 2×d for easier use (a and b rows)
        self._bounds = bounds.T  # shape (2, d)

    @property
    def bounds(self):
        return self._bounds.T  # back to d×2 format

    @property
    def dimension(self):
        return self._bounds.shape[1]

    @property
    def volume(self):
        return np.prod(np.diff(self._bounds, axis=0))

    def __contains__(self, point):
        point = np.asarray(point)
        if point.ndim != 1 or point.size != self.dimension:
            raise ValueError("point must be a 1D coordinate vector")

        a, b = self._bounds
        return np.all(a <= point) and np.all(point <= b)

    def indicator_function(self, points):
        points = np.asarray(points)
        a, b = self._bounds
        return np.logical_and(
            np.all(a <= points, axis=-1),
            np.all(points <= b, axis=-1),
        )

    def rand(self, n=1, seed=None):
        rng = get_random_number_generator(seed)
        a, b = self._bounds
        d = self.dimension

        if n == 1:
            return rng.uniform(a, b)
        return rng.uniform(a, b, size=(n, d))

    # Removed: to_spatstat_owin (R dependency)

    def plot(self, axis=None, **kwargs):
        """Plot rectangle in 2D."""
        if self.dimension != 2:
            raise NotImplementedError("Plotting only implemented for 2D boxes")

        if axis is None:
            _, axis = plt.subplots(figsize=(5, 5))

        kwargs.setdefault("fill", False)
        a = self._bounds[0]
        w, h = np.diff(self._bounds, axis=0).ravel()
        rect = Rectangle(a, w, h, **kwargs)
        axis.add_patch(rect)
        axis.set_aspect("equal")
        return axis


class UnitBoxWindow(BoxWindow):
    """Centered box with side length 1."""

    def __init__(self, center):
        center = np.asarray(center)
        if center.ndim != 1:
            raise ValueError("center must be a 1D array")
        bounds = np.add.outer(center, [-0.5, 0.5])
        super().__init__(bounds)


# -------------------------------------------------------------------------
#  UTILITY CHECKS
# -------------------------------------------------------------------------
def check_cubic_window(window):
    if not isinstance(window, BoxWindow):
        raise TypeError("window must be BoxWindow")
    lengths = np.diff(window.bounds, axis=1)
    if not np.allclose(lengths, lengths[0]):
        raise ValueError("window must be cubic")
    return None


def check_centered_window(window):
    if isinstance(window, BoxWindow):
        if not np.allclose(window.bounds.sum(), 0):
            raise ValueError("BoxWindow not centered at origin.")
    elif isinstance(window, BallWindow):
        if not np.allclose(window.center, 0):
            raise ValueError("BallWindow not centered at origin.")
    return None


def subwindow_parameter_max(window, subwindow_type="BoxWindow"):
    """
    Find largest sub-window contained in `window`.
    """
    if subwindow_type not in ["BoxWindow", "BallWindow"]:
        raise ValueError("subwindow_type must be 'BoxWindow' or 'BallWindow'")

    check_centered_window(window)

    if isinstance(window, BallWindow):
        if subwindow_type == "BallWindow":
            return window.radius
        else:
            return window.radius * 2 / np.sqrt(2)

    if isinstance(window, BoxWindow):
        min_side = np.min(np.diff(window.bounds, axis=1))
        if subwindow_type == "BallWindow":
            return min_side / 2
        return min_side
