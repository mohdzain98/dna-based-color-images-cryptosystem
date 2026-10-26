"""Differential attack measures between two cipher images.

Implements NPCR (number of pixels change rate) and UACI (unified average
changing intensity). Used to evaluate plaintext sensitivity of the cipher.
"""

import numpy as np
from numpy.typing import NDArray


def npcr(image1: NDArray, image2: NDArray) -> float:
    """Return the percentage of pixel positions that differ between two images.

    Raises ValueError when the image shapes do not match.
    """
    if image1.shape != image2.shape:
        raise ValueError("Images must be the same size.")
    return float(np.sum(image1 != image2) / image1.size * 100)


def uaci(img1: NDArray, img2: NDArray) -> float:
    """Return the mean absolute intensity difference as a percentage of 255.

    Casts to int first so uint8 subtraction cannot wrap around.
    """
    if img1.shape != img2.shape:
        raise ValueError("Images must be the same size.")
    diff = np.abs(img1.astype(np.int64) - img2.astype(np.int64))
    return float(diff.sum() / (img1.size * 255) * 100)
