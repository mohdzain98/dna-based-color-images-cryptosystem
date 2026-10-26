"""Diagnostics for the 2D Multiple Collapse Chaotic Map (2D-MCCM).

Contains the map definitions plus bifurcation, phase trajectory, value
distribution and Lyapunov exponent routines. Plots are saved to the given
path (if any) and shown only when `show` is True.
"""

from collections.abc import Callable
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray

from src.logger import AppLogger

logger = AppLogger.get_logger(__name__)

MapFunction = Callable[[float, float, float, float], tuple[float, float]]

A_LABELS = ["[0.0,0.2)", "[0.2,0.4)", "[0.4,0.6)", "[0.6,0.8)", "[0.8,1.0]"]
A_COLORS = ["#cce5ff", "#99ccff", "#66b2ff", "#3399ff", "#007acc"]  # light to dark


def twoDMCCM(a: float, b: float, x: float, y: float) -> tuple[float, float]:
    """Apply one iteration of the 2D-MCCM and return the next (x, y)."""
    x1 = np.arctan((b / (10 * a * y)) + np.tan(a * np.pi * x))
    y1 = np.arctan((b / (a * x)) + np.tan(10 * a * np.pi * y))
    return x1, y1


def icm(a: float, b: float, x_i: float, y_i: float) -> tuple[float, float]:
    """Apply one iteration of the iterative chaotic map with sine feedback."""
    x_next = np.sin(a / y_i) * np.sin(b / x_i)
    y_next = np.sin(a / x_i) * np.sin(b / y_i)
    return x_next, y_next


def licm(x: float, y: float, a: float, b: float) -> tuple[float, float]:
    """Apply one iteration of the logistic-style coupled map.

    Clamps x and y away from zero; returns (nan, nan) on numerical failure.
    """
    try:
        epsilon = 1e-10  # avoid division by zero
        if abs(x) < epsilon:
            x = epsilon
        if abs(y) < epsilon:
            y = epsilon

        denominator1 = b * (y + 3) * a * x * (1 - a * x)
        xi_next = np.sin(21 / denominator1)

        denominator2 = b * (a * xi_next + 3) * y * (1 - y)
        yi_next = np.sin(21 / denominator2)
        return xi_next, yi_next
    except Exception as e:  # noqa: BLE001
        logger.error("Error in licm: %s", e)
        return np.nan, np.nan


def _finish(save_path: str | Path | None, show: bool) -> None:
    """Save the current figure if a path is given, then show or close it."""
    if save_path:
        Path(save_path).parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(save_path)
    if show:
        plt.show()
    else:
        plt.close()


def bifurcation_diagram(
    a_min: float = 0.1,
    a_max: float = 1,
    a_steps: int = 500,
    iterations: int = 2000,
    last: int = 100,
    b: float = 1,
    save_path: str | Path | None = None,
    show: bool = True,
) -> None:
    """Plot x against parameter `a` after discarding `iterations` transients."""
    a_values = np.linspace(a_min, a_max, a_steps)
    x1_vals = []
    y1_vals = []

    for a in a_values:
        x, y = 0.4, -0.6  # initial conditions
        for _ in range(iterations):
            x, y = twoDMCCM(x, y, a, b)  # arg order kept from the notebook
        for _ in range(last):
            x, y = twoDMCCM(x, y, a, b)
            x1_vals.append(a)
            y1_vals.append(x)

    plt.figure(figsize=(8, 6))
    plt.scatter(x1_vals, y1_vals, s=0.05, color="blue")
    plt.title("Bifurcation Diagram of 2D MCCM")
    plt.xlabel("a")
    plt.ylabel("x")
    _finish(save_path, show)


def generate_bifurcation_data(
    map_function: MapFunction,
    a_min: float = 0.0,
    a_max: float = 1.0,
    b: float = 21,
    a_steps: int = 500,
    discard: int = 500,
    keep: int = 1000,
    x0: float = 0.4,
    y0: float = 0.6,
) -> tuple[NDArray, NDArray]:
    """Iterate `map_function(x, y, a, b)` per `a`; return kept (xs, as) samples."""
    xs = []
    as_ = []
    for a in np.linspace(a_min, a_max, a_steps):
        x, y = x0, y0
        for _ in range(discard):
            x, y = map_function(x, y, a, b)
        for _ in range(keep):
            x, y = map_function(x, y, a, b)
            xs.append(x)
            as_.append(a)
    return np.array(xs), np.array(as_)


def compute_distribution(
    xs: NDArray, as_: NDArray, a_bins: int = 5, x_bins: int = 10
) -> tuple[NDArray, NDArray]:
    """Return x bin centers and the per-`a`-bin proportion of x values.

    Rows of the distribution sum to 1, or stay 0 when an `a` bin is empty.
    """
    a_bin_edges = np.linspace(0, 1, a_bins + 1)
    x_bin_edges = np.linspace(np.min(xs), np.max(xs), x_bins + 1)
    distribution = np.zeros((a_bins, x_bins))

    for i in range(a_bins):
        mask_a = (as_ >= a_bin_edges[i]) & (as_ < a_bin_edges[i + 1])
        xs_in_bin = xs[mask_a]
        if len(xs_in_bin) > 0:
            counts, _ = np.histogram(xs_in_bin, bins=x_bin_edges)
            distribution[i, :] = counts / np.sum(counts)

    x_bin_centers = (x_bin_edges[:-1] + x_bin_edges[1:]) / 2
    return x_bin_centers, distribution


def plot_stacked_bar(
    x_centers: NDArray,
    distribution: NDArray,
    save_path: str | Path | None = None,
    show: bool = True,
) -> None:
    """Draw the distribution as stacked bars, one color per `a` range."""
    _, ax = plt.subplots(figsize=(10, 6))
    bottoms = np.zeros(len(x_centers))
    for i in range(distribution.shape[0]):
        ax.bar(
            x_centers,
            distribution[i],
            bottom=bottoms,
            color=A_COLORS[i],
            width=(x_centers[1] - x_centers[0]) * 0.8,
            label=f"a in {A_LABELS[i]}",
        )
        bottoms += distribution[i]

    ax.set_xlabel("x(i)")
    ax.set_ylabel("Proportion")
    ax.set_title("Bifurcation Value Distribution")
    ax.legend()
    plt.tight_layout()
    _finish(save_path, show)


def phase_trajectory(
    a: float = 1,
    b: float = 21,
    total_steps: int = 20000,
    save_path: str | Path | None = None,
    show: bool = True,
) -> None:
    """Scatter-plot the (x, y) orbit of the map over `total_steps` iterations."""
    x, y = 0.4, 0.6
    xs, ys = [], []
    for _ in range(total_steps):
        x, y = twoDMCCM(x, y, a, b)  # arg order kept from the notebook
        xs.append(x)
        ys.append(y)

    plt.figure(figsize=(8, 6))
    plt.scatter(xs, ys, s=0.1, color="blue")
    plt.title("Phase Trajectory of 2D-MCCM", fontsize=14)
    plt.xlabel("xi", fontsize=12)
    plt.ylabel("yi", fontsize=12)
    _finish(save_path, show)


def lyapunov_exponents_2d_mccm(
    a: float, b: float = 21, steps: int = 20000, delta: float = 1e-6
) -> NDArray:
    """Estimate both Lyapunov exponents via QR on finite-difference Jacobians."""
    x, y = 0.4, 0.6
    q = np.identity(2)
    lyap = np.zeros(2)

    for _ in range(steps):
        fx, fy = twoDMCCM(x, y, a, b)  # arg order kept from the notebook
        x_dx, y_dx = twoDMCCM(x + delta, y, a, b)
        x_dy, y_dy = twoDMCCM(x, y + delta, a, b)

        jacobian = np.array(
            [
                [(x_dx - fx) / delta, (y_dy - fy) / delta],
                [(y_dx - fy) / delta, (x_dy - fx) / delta],
            ]
        )  # same layout as the notebook
        q, r = np.linalg.qr(jacobian @ q)
        lyap += np.log(np.abs(np.diag(r)))
        x, y = fx, fy

    return lyap / steps


def plot_lyapunov_sweep(
    a_values: NDArray | None = None,
    b: float = 21,
    save_path: str | Path | None = None,
    show: bool = True,
) -> tuple[list[float], list[float]]:
    """Sweep `a`, plot both exponents and return them as (lambda1, lambda2)."""
    if a_values is None:
        a_values = np.linspace(0.1, 1.0, 50)
    lambda1, lambda2 = [], []
    logger.info("Calculating Lyapunov exponents...")
    for a in a_values:
        l1, l2 = lyapunov_exponents_2d_mccm(a, b=b)
        lambda1.append(l1)
        lambda2.append(l2)

    plt.figure(figsize=(10, 6))
    plt.plot(a_values, lambda1, label="λ1", color="red")
    plt.plot(a_values, lambda2, label="λ2", color="blue")
    plt.axhline(0, color="black", linestyle="--", lw=1)
    plt.title("Lyapunov Exponents 2D-MCCM", fontsize=14)
    plt.xlabel("a", fontsize=13)
    plt.ylabel("Lyapunov Exponent", fontsize=13)
    plt.legend()
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.tight_layout()
    _finish(save_path, show)
    return lambda1, lambda2
