"""End-to-end test of the command line interface."""

from pathlib import Path

import cv2
import numpy as np

from src.__main__ import main


def test_cli_encrypt_then_decrypt(tmp_path: Path, image: np.ndarray) -> None:
    """encrypt writes a cipher and key; decrypt restores the original image."""
    src = tmp_path / "plain.png"
    cv2.imwrite(str(src), image)
    cipher = tmp_path / "plain_cipher.png"

    assert main(["encrypt", str(src), "-o", str(cipher)]) == 0
    key_file = cipher.with_suffix(".key")
    assert key_file.exists()

    restored = tmp_path / "restored.png"
    assert main(["decrypt", str(cipher), "-k", str(key_file), "-o", str(restored)]) == 0
    assert np.array_equal(cv2.imread(str(restored)), image)


def test_cli_decrypt_requires_key(tmp_path: Path, image: np.ndarray) -> None:
    """decrypt without --key-file fails."""
    src = tmp_path / "c.png"
    cv2.imwrite(str(src), image)
    assert main(["decrypt", str(src)]) == 1


def test_cli_missing_image(tmp_path: Path) -> None:
    """A nonexistent image path returns a non-zero exit code."""
    assert main(["encrypt", str(tmp_path / "nope.png")]) == 1
