"""Logging helpers for consistent console and file output.

Provides colored console formatting, optional file logging, and a
`logger.print()` pretty-printer toggled by the DNA_PRINT_ENABLED env var.
Call AppLogger.init() once; modules then use AppLogger.get_logger(__name__).
"""

import json
import logging
import os
from pathlib import Path

from rich.console import Console
from rich.pretty import pretty_repr

PROJECT_ROOT = Path(__file__).resolve().parents[1]

LEVEL_COLORS = {
    "DEBUG": "\033[36m",
    "INFO": "\033[32m",
    "WARNING": "\033[33m",
    "ERROR": "\033[31m",
    "CRITICAL": "\033[35m",
}
RESET_COLOR = "\033[0m"

_STANDARD_LOG_RECORD_FIELDS = set(logging.makeLogRecord({}).__dict__)
_INTERNAL_FORMATTER_FIELDS = {"colored_levelname", "error_suffix", "extra_suffix"}


def _format_extra_value(value: object) -> str:
    """Render an extra log field as a compact string."""
    if isinstance(value, str):
        return value
    try:
        return json.dumps(value, ensure_ascii=False, default=str)
    except TypeError:
        return str(value)


def _env_flag(name: str) -> bool:
    """Return True when the environment variable is set to a truthy value."""
    return os.environ.get(name, "").strip().lower() in {"1", "true", "yes", "on"}


class SafeExtraFormatter(logging.Formatter):
    """Appends any `extra=` fields on a record as `key=value` pairs."""

    def format(self, record: logging.LogRecord) -> str:
        extra_items = []
        for key, value in sorted(record.__dict__.items()):
            if key in _STANDARD_LOG_RECORD_FIELDS or key in _INTERNAL_FORMATTER_FIELDS:
                continue
            if key.startswith("_"):
                continue
            extra_items.append(f"{key}={_format_extra_value(value)}")
        record.extra_suffix = f" | {' | '.join(extra_items)}" if extra_items else ""
        return super().format(record)


class ColorFormatter(SafeExtraFormatter):
    """Adds an ANSI-colored, padded level name to console output."""

    def format(self, record: logging.LogRecord) -> str:
        padded_level = f"{record.levelname + ':':<9}"
        color = LEVEL_COLORS.get(record.levelname, "")
        record.colored_levelname = (
            f"{color}{padded_level}{RESET_COLOR}" if color else padded_level
        )
        return super().format(record)


_console = Console()


class DebugAwareLogger(logging.Logger):
    """Logger subclass adding `.print()`, a no-op unless DNA_PRINT_ENABLED is true.

    The flag is read on every call, so notebooks can toggle it via os.environ.
    Replaces scattered `print()` calls with one switchable pretty-printer.
    """

    def print(
        self, label: str, data: object = None, *, raw_keys: set[str] | None = None
    ) -> None:
        """Pretty-print `data` under a `label` rule.

        `raw_keys` names dict keys whose string values print as literal text
        (real newlines) instead of escaped reprs.
        """
        if not _env_flag("DNA_PRINT_ENABLED"):
            return

        _console.rule(f"[bold cyan]{label}[/bold cyan]", style="cyan")
        if data is None:
            return

        if isinstance(data, dict) and raw_keys:
            remaining = {}
            for key, value in data.items():
                if key in raw_keys and isinstance(value, str):
                    _console.print(f"[bold]{key}:[/bold]")
                    _console.print(value)
                else:
                    remaining[key] = value
            if remaining:
                _console.print(pretty_repr(remaining))
            return

        _console.print(pretty_repr(data))


# Must run before any getLogger() call so loggers get the print() method
logging.setLoggerClass(DebugAwareLogger)


class AppLogger:
    """Central logging setup and logger factory.

    init() is idempotent and configures the root logger once. File logging
    writes WARNING and above to logs/ under the project root.
    """

    _configured: bool = False

    @classmethod
    def init(
        cls,
        level: int | None = None,
        log_to_file: bool = False,
        filename: str = "dna_cryptosystem.log",
    ) -> None:
        """Configure the root logger once; later calls are no-ops."""

        if cls._configured:
            return
        cls._configured = True
        if level is None:
            level = logging.DEBUG if _env_flag("DNA_DEBUG") else logging.INFO

        root_logger = logging.getLogger()
        root_logger.setLevel(level)

        for handler in root_logger.handlers[:]:
            root_logger.removeHandler(handler)

        console_handler = logging.StreamHandler()
        console_handler.setLevel(level)
        console_format = "%(colored_levelname)s %(name)s | %(message)s%(extra_suffix)s"
        console_handler.setFormatter(ColorFormatter(console_format))
        root_logger.addHandler(console_handler)

        if log_to_file:
            logs_dir = PROJECT_ROOT / "logs"
            logs_dir.mkdir(parents=True, exist_ok=True)
            file_path = logs_dir / filename

            file_handler = logging.FileHandler(file_path, encoding="utf-8")
            file_handler.setLevel(logging.WARNING)
            file_format = (
                "%(asctime)s | %(levelname)-8s | %(name)s | "
                "%(filename)s:%(lineno)d | %(message)s%(extra_suffix)s"
            )
            file_handler.setFormatter(
                SafeExtraFormatter(file_format, datefmt="%Y-%m-%d %H:%M:%S")
            )
            root_logger.addHandler(file_handler)

    @staticmethod
    def get_logger(name: str | None = None) -> "DebugAwareLogger":
        """Return a logger named after the last dotted segment of `name`."""

        logger_name = name if name is not None else __name__
        short_name = logger_name.split(".")[-1]
        return logging.getLogger(short_name)
