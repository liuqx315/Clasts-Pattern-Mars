import matplotlib.pyplot as plt
import numpy as np

from structure_factor.spatial_windows import AbstractSpatialWindow


class PointPattern:
    """
    Encapsulate one realization of a point process, the observation window,
    and the intensity of the underlying point process.

    * NumPy-2.0 compatible
    * Pure Python (no rpy2, no R, no spatstat)
    """

    def __init__(self, points, window, intensity=None, **params):
        points = np.asarray(points)
        assert points.ndim == 2, "Points should be an N×d array"
        self.points = points

        assert isinstance(window, AbstractSpatialWindow)
        self.window = window

        # Estimate intensity if not provided
        if intensity is None:
            intensity = self.points.shape[0] / window.volume

        if intensity <= 0:
            raise ValueError("Intensity must be positive.")

        self.intensity = float(intensity)
        self.params = params

    # ------------------------------------------------------------------
    @property
    def dimension(self):
        """Ambient dimension of the point space."""
        return self.points.shape[1]

    # ------------------------------------------------------------------
    def restrict_to_window(self, window):
        """
        Returns a new PointPattern restricted to the given window.
        """
        assert isinstance(window, AbstractSpatialWindow)
        mask = window.indicator_function(self.points)
        return PointPattern(self.points[mask], window, self.intensity)

    # ------------------------------------------------------------------
    # OLD (removed):
    # def convert_to_spatstat_ppp(...)
    #
    # NEW:
    def to_dict(self):
        """
        Return a pure-Python dictionary representation of the point pattern.

        Useful for:
        * JSON export
        * writing to .csv / .npz
        * interoperability without needing R
        """
        return {
            "points": self.points.copy(),
            "window_bounds": self.window.bounds if hasattr(self.window, "bounds") else None,
            "intensity": self.intensity,
            "params": self.params,
        }

    def to_npz(self, filename):
        """Save point pattern and metadata to a .npz file."""
        data = self.to_dict()
        np.savez(filename, **data)

    # ------------------------------------------------------------------
    def plot(self, axis=None, window=None, show_window=False, file_name="", **kwargs):
        """
        Scatter plot of the point pattern.
        """
        if axis is None:
            fig, axis = plt.subplots(figsize=(5, 5))

        if window is None:
            window = self.window
            points = self.points
        else:
            assert isinstance(window, AbstractSpatialWindow)
            mask = window.indicator_function(self.points)
            points = self.points[mask]

        if show_window:
            window.plot(axis=axis)

        kwargs.setdefault("c", "k")
        kwargs.setdefault("s", 1.0)
        axis.scatter(points[:, 0], points[:, 1], **kwargs)
        axis.set_aspect("equal")

        if file_name:
            axis.get_figure().savefig(file_name, bbox_inches="tight")

        return axis
