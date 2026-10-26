"""Shared pytest fixtures: headless matplotlib and a key-derived cipher setup."""

import os
import sys
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")  # matrices() calls plt.show()
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402
import pytest  # noqa: E402

from src.cryptosystem import Cryptosystem  # noqa: E402

KEY = (
    "10101111110011010100000101001011010110110101001010010101010001101100111011011001"
    "01011111100101100001000101010000110001011000110110111111100110011111101011111111"
    "01001001111100111010101110000010000011101100100110000100000111101000100010101101"
)
SIZE = 32


def build_params(image: np.ndarray, key: str = KEY) -> dict:
    """Derive the IV, rule maps and DSM from `key`, as the notebooks do."""
    x = Cryptosystem(image)
    a, b, x0, y0 = x.generateInitials(key)
    maps = [x.getChaoticMap(a, b, x0, y0)]
    for _ in range(5):
        maps.append(x.getChaoticMap(a, b, maps[-1][-1][0], maps[-1][-1][1]))
    rm2 = x.getRuleMap(maps[2], 8)
    return {
        "IV": x.getfloorIV(8, maps[0]),
        "rm1": x.getRuleMap(maps[1], 8),
        "rm3": x.getRuleMap(maps[4], 8),
        "rm4": x.getRuleMap(maps[5], 8),
        "dsm": x.getDSMMap(maps[3], rm2),
    }


@pytest.fixture(scope="session")
def image() -> np.ndarray:
    """Return a small deterministic random BGR image."""
    rng = np.random.default_rng(0)
    return rng.integers(0, 256, (SIZE, SIZE, 3), dtype=np.uint8)


@pytest.fixture(scope="session")
def params(image: np.ndarray) -> dict:
    """Return cipher parameters derived from the fixed test key."""
    return build_params(image)
