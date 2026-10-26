"""High-level encrypt/decrypt entry points built on the cipher classes.

Derives the IV, rule maps and substitution map from a binary key, then runs
the full encryption or decryption algorithm. Shared by the CLI, tests and
notebooks so the key schedule lives in one place.
"""

import secrets
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from src.cryptosystem import Cryptosystem, Decryption, Encryption

KEY_BITS = 256
MIN_KEY_BITS = 240  # generateInitials reads the first 240 bits
MAP_COUNT = 6


@dataclass
class CipherParams:
    """Key-derived inputs of the algorithm: IV, rule maps and DSM."""

    IV: list
    rm1: list
    rm3: list
    rm4: list
    dsm: list


def generate_key() -> str:
    """Return a random binary key string of KEY_BITS bits."""
    return format(secrets.randbits(KEY_BITS), f"0{KEY_BITS}b")


def build_params(image: NDArray[np.uint8], key: str) -> CipherParams:
    """Derive the cipher parameters for `image` from a binary `key`.

    Raises ValueError when the key is shorter than MIN_KEY_BITS or not binary.
    """
    if len(key) < MIN_KEY_BITS or set(key) - {"0", "1"}:
        raise ValueError(f"Key must be a binary string of >= {MIN_KEY_BITS} bits.")
    system = Cryptosystem(image)
    a, b, x0, y0 = system.generateInitials(key)
    maps = [system.getChaoticMap(a, b, x0, y0)]
    for _ in range(
        MAP_COUNT - 1
    ):  # each map starts from the last point of the previous
        last_x, last_y = maps[-1][-1]
        maps.append(system.getChaoticMap(a, b, last_x, last_y))
    return CipherParams(
        IV=system.getfloorIV(8, maps[0]),
        rm1=system.getRuleMap(maps[1], 8),
        rm3=system.getRuleMap(maps[4], 8),
        rm4=system.getRuleMap(maps[5], 8),
        dsm=system.getDSMMap(maps[3], system.getRuleMap(maps[2], 8)),
    )


def encrypt_image(image: NDArray[np.uint8], key: str) -> NDArray[np.uint8]:
    """Encrypt a BGR image with the binary `key`."""
    p = build_params(image, key)
    return Encryption(image).applyEncryption(image, p.IV, p.rm1, p.rm3, p.rm4, p.dsm)


def decrypt_image(cipher: NDArray[np.uint8], key: str) -> NDArray[np.uint8]:
    """Decrypt a BGR cipher image with the binary `key` used to encrypt it."""
    p = build_params(cipher, key)
    return Decryption(cipher).applyDecryption(cipher, p.IV, p.rm1, p.rm3, p.rm4, p.dsm)
