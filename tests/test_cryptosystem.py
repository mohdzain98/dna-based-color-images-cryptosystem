"""Round-trip and key-dependence tests for the DNA cryptosystem."""

import numpy as np
import pytest

from src.pipeline import decrypt_image, encrypt_image
from tests.conftest import KEY


def test_encrypt_decrypt_roundtrip(image: np.ndarray) -> None:
    """Decrypting the cipher image returns the original pixels."""
    cipher = encrypt_image(image, KEY)
    assert np.array_equal(decrypt_image(cipher, KEY), image)


def test_cipher_differs_from_plain(image: np.ndarray) -> None:
    """The cipher image shares the shape but not the content of the plain image."""
    cipher = encrypt_image(image, KEY)
    assert cipher.shape == image.shape
    assert not np.array_equal(cipher, image)


def test_encryption_is_deterministic(image: np.ndarray) -> None:
    """The same key and image always give the same cipher."""
    assert np.array_equal(encrypt_image(image, KEY), encrypt_image(image, KEY))


def test_different_key_gives_different_cipher(image: np.ndarray) -> None:
    """Flipping one key bit changes the cipher image."""
    flipped = ("0" if KEY[0] == "1" else "1") + KEY[1:]
    assert not np.array_equal(encrypt_image(image, KEY), encrypt_image(image, flipped))


@pytest.mark.parametrize("bad_key", ["", "1010", "2" * 256])
def test_invalid_key_rejected(image: np.ndarray, bad_key: str) -> None:
    """Short or non-binary keys raise ValueError."""
    with pytest.raises(ValueError):
        encrypt_image(image, bad_key)
