"""Tests for the analysis helpers: differential, metrics, performance, chaos, nist."""

import numpy as np
import pytest

from src.analysis.chaos import compute_distribution, twoDMCCM
from src.analysis.differential import npcr, uaci
from src.analysis.metrics import experimentalAnalysis
from src.analysis.nist import array_to_binary
from src.analysis.performance import calculate_et_and_cpb


@pytest.fixture
def metrics(tmp_path) -> experimentalAnalysis:
    """Return a metrics helper writing into a temp dir."""
    return experimentalAnalysis(np.zeros((64, 64), np.uint8), str(tmp_path))


def test_npcr_identical_and_different() -> None:
    """NPCR is 0 for equal images and 100 when every pixel differs."""
    a = np.zeros((8, 8), np.uint8)
    assert npcr(a, a) == 0
    assert npcr(a, a + 1) == 100


def test_npcr_shape_mismatch_raises() -> None:
    """Mismatched shapes are rejected."""
    with pytest.raises(ValueError):
        npcr(np.zeros((2, 2)), np.zeros((3, 3)))


def test_uaci_does_not_wrap_uint8() -> None:
    """0 vs 200 is a 200/255 difference, not a wrapped uint8 value."""
    a = np.zeros((4, 4), np.uint8)
    assert uaci(a, np.full((4, 4), 200, np.uint8)) == pytest.approx(200 / 255 * 100)
    assert uaci(a, a) == 0


def test_entropy_bounds(metrics: experimentalAnalysis) -> None:
    """Random data is near 8 bits, constant data is 0."""
    rng = np.random.default_rng(1)
    assert metrics.entropy(rng.integers(0, 256, (256, 256))) > 7.9
    assert metrics.entropy(np.zeros((64, 64), np.uint8)) == pytest.approx(0, abs=1e-6)


def test_variance_matches_numpy(metrics: experimentalAnalysis) -> None:
    """The numba variance agrees with numpy's population variance."""
    img = np.random.default_rng(2).integers(0, 256, (64, 64)).astype(np.float64)
    assert metrics.variance_numba(img) == pytest.approx(np.var(img), abs=1e-2)


def test_psnr_uses_image_shape(metrics: experimentalAnalysis) -> None:
    """PSNR works for non-256 sizes and is inf for identical images."""
    a = np.zeros((16, 32), np.uint8)
    b = np.full((16, 32), 10, np.uint8)
    assert metrics.calculate_PSNR(a, b) == pytest.approx(10 * np.log10(255**2 / 100))
    assert metrics.calculate_PSNR(a, a) == float("inf")


def test_et_and_cpb() -> None:
    """Throughput is bytes over seconds and CPB is clock over throughput."""
    (row,) = calculate_et_and_cpb([(10, 10)], [1.0], cpu_speed_ghz=1.0)
    assert row["ET (bytes/sec)"] == 300.0
    assert row["CPB (cycles/byte)"] == pytest.approx(1e9 / 300, abs=0.01)
    assert row["Image Size"] == "10x10"


def test_compute_distribution_rows_sum_to_one() -> None:
    """Each non-empty `a` bin is a proportion vector."""
    rng = np.random.default_rng(3)
    as_ = rng.uniform(0, 0.999, 2000)
    xs = rng.normal(size=2000)
    centers, dist = compute_distribution(xs, as_)
    assert len(centers) == 10
    assert np.allclose(dist.sum(axis=1), 1)


def test_twodmccm_is_deterministic_and_bounded() -> None:
    """The map is arctan-valued, so outputs stay within (-pi/2, pi/2)."""
    first = twoDMCCM(0.5, 21, 0.4, 0.6)
    assert first == twoDMCCM(0.5, 21, 0.4, 0.6)
    assert all(abs(v) < np.pi / 2 for v in first)


def test_array_to_binary_length() -> None:
    """Every uint8 expands into 8 bits."""
    bits = array_to_binary(np.array([0, 255, 1], np.uint8))
    assert len(bits) == 24
    assert bits[:8].sum() == 0 and bits[8:16].sum() == 8
