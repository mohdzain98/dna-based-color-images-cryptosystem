"""Shared pytest fixtures: headless matplotlib, a fixed key and a small image."""

import os

import numpy as np
import pytest

os.environ.setdefault("MPLBACKEND", "Agg")  # matrices() calls plt.show()

KEY = (
    "10101111110011010100000101001011010110110101001010010101010001101100111011011001"
    "01011111100101100001000101010000110001011000110110111111100110011111101011111111"
    "01001001111100111010101110000010000011101100100110000100000111101000100010101101"
)
SIZE = 32


@pytest.fixture(scope="session")
def image() -> np.ndarray:
    """Return a small deterministic random BGR image."""
    rng = np.random.default_rng(0)
    return rng.integers(0, 256, (SIZE, SIZE, 3), dtype=np.uint8)
