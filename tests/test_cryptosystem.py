"""Round-trip and key-dependence tests for the DNA cryptosystem."""

import numpy as np

from src.cryptosystem import Decryption, Encryption
from tests.conftest import KEY, build_params


def encrypt(image: np.ndarray, p: dict) -> np.ndarray:
    """Encrypt `image` with the given parameter dict."""
    return Encryption(image).applyEncryption(
        image, p["IV"], p["rm1"], p["rm3"], p["rm4"], p["dsm"]
    )


def test_encrypt_decrypt_roundtrip(image: np.ndarray, params: dict) -> None:
    """Decrypting the cipher image returns the original pixels."""
    cipher = encrypt(image, params)
    p = params
    plain = Decryption(cipher).applyDecryption(
        cipher, p["IV"], p["rm1"], p["rm3"], p["rm4"], p["dsm"]
    )
    assert np.array_equal(plain, image)


def test_cipher_differs_from_plain(image: np.ndarray, params: dict) -> None:
    """The cipher image shares the shape but not the content of the plain image."""
    cipher = encrypt(image, params)
    assert cipher.shape == image.shape
    assert not np.array_equal(cipher, image)


def test_encryption_is_deterministic(image: np.ndarray, params: dict) -> None:
    """The same key and image always give the same cipher."""
    assert np.array_equal(encrypt(image, params), encrypt(image, params))


def test_different_key_gives_different_cipher(image: np.ndarray, params: dict) -> None:
    """Flipping one key bit changes the cipher image."""
    flipped = ("0" if KEY[0] == "1" else "1") + KEY[1:]
    other = build_params(image, flipped)
    assert not np.array_equal(encrypt(image, params), encrypt(image, other))
