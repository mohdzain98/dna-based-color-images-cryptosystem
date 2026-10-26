"""Command line interface: python -m src encrypt|decrypt IMAGE.

Encrypt writes the cipher image and a key file (generated unless --key-file
is given); decrypt needs the same key file. Runs headless (Agg backend).
"""

import argparse
import os
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")  # cipher classes call plt.show()

import cv2  # noqa: E402

from src.logger import AppLogger  # noqa: E402
from src.pipeline import decrypt_image, encrypt_image, generate_key  # noqa: E402

logger = AppLogger.get_logger(__name__)


def build_parser() -> argparse.ArgumentParser:
    """Create the argument parser with encrypt and decrypt subcommands."""
    parser = argparse.ArgumentParser(prog="python -m src", description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ("encrypt", "decrypt"):
        cmd = sub.add_parser(name, help=f"{name} an image")
        cmd.add_argument("image", type=Path, help="input image path")
        cmd.add_argument("-o", "--output", type=Path, help="output image path")
        cmd.add_argument("-k", "--key-file", type=Path, help="binary key text file")
    return parser


def main(argv: list[str] | None = None) -> int:
    """Run the CLI and return a process exit code."""
    args = build_parser().parse_args(argv)
    image = cv2.imread(str(args.image))
    if image is None:
        logger.error("Cannot read image: %s", args.image)
        return 1

    suffix = "cipher" if args.command == "encrypt" else "decrypted"
    output = args.output or args.image.with_name(f"{args.image.stem}_{suffix}.png")

    if args.key_file:
        key = args.key_file.read_text().strip()
    elif args.command == "encrypt":
        key = generate_key()
        key_path = output.with_suffix(".key")
        key_path.write_text(key)
        logger.info("Generated key saved to %s", key_path)
    else:
        logger.error("decrypt requires --key-file")
        return 1

    run = encrypt_image if args.command == "encrypt" else decrypt_image
    cv2.imwrite(str(output), run(image, key))
    logger.info("Saved %s", output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
