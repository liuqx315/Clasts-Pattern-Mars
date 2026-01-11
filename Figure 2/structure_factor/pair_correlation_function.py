"""Collection of functions used to estimate the pair correlation function of an isotropic point process. The underlying routines call the `R` package `spatstat <https://github.com/spatstat/spatstat>`_.

- :py:meth:`~structure_factor.pair_correlation_function.estimate`: Estimates the pair correlation function.

- :py:meth:`~structure_factor.pair_correlation_function.interpolate`: Cleans, interpolates, and extrapolates the results of :py:meth:`~structure_factor.pair_correlation_function.estimate`.

- :py:meth:`~structure_factor.pair_correlation_function.plot`: Plots the results of :py:meth:`~structure_factor.pair_correlation_function.estimate`.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rpy2.robjects as robjects
from scipy.interpolate import interp1d
# from spatstat_interface.interface import SpatstatInterface

import structure_factor.plotting as plots
import structure_factor.utils as utils


def estimate(point_pattern, method="fv", install_spatstat=False, **params):
    """
    Pure Python version of pair correlation function estimator.
    Keeps same API as original spatstat-based version.
    """

    import numpy as np
    import pandas as pd
    from scipy.spatial.distance import pdist

    # --- checks ---
    assert point_pattern.dimension in (2, 3)
    assert method in ("ppp", "fv")

    # --- get points ---
    if hasattr(point_pattern, "points"):
        pts = np.asarray(point_pattern.points)
    elif hasattr(point_pattern, "get_points"):
        pts = np.asarray(point_pattern.get_points())
    else:
        raise ValueError("point_pattern must provide .points or .get_points()")

    N, dim = pts.shape

    # --- define window ---
    mins = pts.min(axis=0)
    maxs = pts.max(axis=0)
    box = maxs - mins

    if dim == 2:
        volume = box[0] * box[1]
    else:
        volume = box[0] * box[1] * box[2]

    intensity = N / volume

    # --- radii setup (similar to spatstat defaults) ---
    r = params.get("r", None)

    if r is None:
        r_max = params.get("rmax", 0.25 * np.min(box))
        n_bins = params.get("n_bins", 100)
        edges = np.linspace(0.0, r_max, n_bins + 1)
        r = 0.5 * (edges[:-1] + edges[1:])
        dr = edges[1] - edges[0]
    else:
        r = np.asarray(r)
        dr = r[1] - r[0]
        edges = np.concatenate(([r[0] - dr / 2], r + dr / 2))

    # --- pairwise distances ---
    dists = pdist(pts)

    # --- histogram counts ---
    counts, _ = np.histogram(dists, bins=edges)

    # --- shell measures ---
    if dim == 2:
        shell = 2.0 * np.pi * r * dr
    else:
        shell = 4.0 * np.pi * r**2 * dr

    # --- expected Poisson counts ---
    expected = intensity * N * shell

    # --- main PCF ---
    g = counts / expected

    # --- build DataFrame (mimic spatstat columns) ---
    pcf_pd = pd.DataFrame({
        "r": r,
        "pcf": g,
        "border": g.copy(),
        "translation": g.copy(),
        "isotropic": g.copy()
    })

    return pcf_pd


# todo add test
# todo clean up arguments: only drop, nan, posinf, neginf are necessary, clean and replace can be removed
def interpolate(
    r,
    pcf_r,
    clean=True,
    drop=True,
    replace=False,
    nan=0.0,
    posinf=0.0,
    neginf=0.0,
    extrapolate_with_one=True,
    **params
):
    r"""Interpolate and then extrapolate the evaluation ``pcf_r`` of the pair correlation function (pcf) ``r``.

    The interpolation is performed with `scipy.interpolate.interp1d <https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.interp1d.html>`_.

    Args:
        r (numpy.ndarray): Vector of radii. Typically, the first column of the output of :py:meth:`~structure_factor.structure_factor.StructureFactor.estimate`.

        pcf_r (numpy.ndarray): Vector of evaluations of the pair correlation function at ``r``. Typically, a column from the output of :py:meth:`~structure_factor.structure_factor.StructureFactor.estimate`.

        clean (bool, optional): If ``True``, a method is chosen to deal with possible outliers (nan, posinf, neginf) of ``pcf_r``. The chosen method depends on the parameters ``replace`` and ``drop``. Defaults to True.

        drop (bool, optional): Cleaning method for ``pcf_r`` active when it's set to True simultaneously with ``clean=True``. Drops possible nan, posinf, and neginf from ``pcf_r`` with the corresponding values of ``r``.

        replace (bool, optional): Cleaning method for ``pcf_r`` active when it's set to True simultaneously with ``clean=True``. Replaces possible nan, posinf, and neginf values of ``pcf_r`` by the values set in the corresponding arguments. Defaults to True.

        nan (float, optional): When ``replace=True``, replacing value of nan present in ``pcf_r``. Defaults to 0.0.

        posinf (float, optional): When ``replace=True``, replacing value of +inf values present in ``pcf_r``. Defaults to 0.0.

        neginf (float, optional): When is ``replace=True``, replacing value of -inf present in ``pcf_r``. Defaults to 0.0.

        extrapolate_with_one (bool, optional): If True, the discrete approximation vector ``pcf_r`` is first interpolated until the maximal value of ``r``, then the extrapolated values are fixed to 1. If False, the extrapolation method of `scipy.interpolate.interp1d <https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.interp1d.html>`_ is used. Note that, the approximation of the structure factor by Ogata quadrature Hankel transform is usually better when set to True. Defaults to True.

    Keyword Args:
        params (dict): Keyword arguments of the function `scipy.interpolate.interp1d <https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.interp1d.html>`_.

    Returns:
        callable: Interpolated pair correlation function.

    Example:
        .. plot:: code/pair_correlation_function/interpolate_pcf.py
            :include-source: True

    .. seealso::

        - :py:meth:`~structure_factor.point_pattern.PointPattern`
        - :py:meth:`~structure_factor.pair_correlation_function.PairCorrelationFuntcion.estimate`
    """
    params.setdefault("kind", "cubic")
    if clean:
        if replace:
            pcf_r = utils.set_nan_inf_to_zero(
                pcf_r, nan=nan, posinf=posinf, neginf=neginf
            )
        elif drop:
            index_outlier = np.isnan(pcf_r) | np.isinf(pcf_r)
            pcf_r = pcf_r[~index_outlier]
            r = r[~index_outlier]

    if extrapolate_with_one:
        pcf = lambda x: _extrapolate_pcf(x, r, pcf_r, **params)
    else:
        params.setdefault("fill_value", "extrapolate")
        pcf = interp1d(r, pcf_r, **params)

    return pcf


#! todo add example in the doc
def plot(pcf_dataframe, exact_pcf=None, file_name="", **kwargs):
    r"""Plot the columns of ``pcf_dataframe`` with respect to the column ``pcf_dataframe["r"]``.

    Args:
        pcf_dataframe (pandas.DataFrame): DataFrame to be visualized. Typically the output of :py:meth:`~structure_factor.structure_factor.StructureFactor.estimate`.

        exact_pcf (callable): Theoretical pair correlation function of the point process.

        file_name (str): Name used to save the figure. The available output formats depend on the backend being used.

    Keyword Args:
        kwargs (dict): Keyword arguments of the function `pandas.DataFrame.plot.line <https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.plot.line.html>`_.

    Return:
        plt.Axes: Plot result.
    """
    axis = pcf_dataframe.plot.line(x="r", **kwargs)
    if exact_pcf is not None:
        axis.plot(
            pcf_dataframe["r"],
            exact_pcf(pcf_dataframe["r"]),
            "g",
            label=r"Exact $g(r)$",
        )
    plots.plot_poisson(pcf_dataframe["r"], axis=axis, linestyle=(0, (5, 5)))

    axis.legend()
    axis.set_xlabel(r"Radius ($r$)")
    axis.set_ylabel(r"Pair correlation function ($g(r)$)")
    plt.show()
    if file_name:
        fig = axis.get_figure()
        fig.savefig(file_name, bbox_inches="tight")
    return axis


def _extrapolate_pcf(x, r, pcf_r, **params):
    """Interpolate pcf_r for x=<r_max and set to 1 for x>r_max.

    Args:
        x (numpy.ndarray): Points on which the pair correlation function is to be evaluated.

        r (numpy.ndarray): Vector of the radius on with the pair correlation function was evaluated.

        pcf_r (numpy.ndarray): Vector of evaluations of the pair correlation function corresponding to ``r``.

    Returns:
        numpy.ndarray: evaluation of the extrapolated pair correlation function on ``x``.
    """
    r_max = np.max(r)  # maximum radius
    pcf = np.zeros_like(x)
    params.setdefault("fill_value", "extrapolate")

    mask = x > r_max
    if np.any(mask):
        pcf[mask] = 1.0
        np.logical_not(mask, out=mask)
        pcf[mask] = interp1d(r, pcf_r, **params)(x[mask])
    else:
        pcf = interp1d(r, pcf_r, **params)(x)

    return pcf
