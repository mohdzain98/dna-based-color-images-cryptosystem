"""NIST SP 800-22 randomness testing for cipher images.

Runs the nistrng battery per color channel and keeps the scores and
PASSED/FAILED statuses in two dataframes that can be loaded from and saved
to CSV. Used by the quality checks notebook.
"""

from pathlib import Path

import numpy as np
import pandas as pd
from nistrng import (
    SP800_22R1A_BATTERY,
    check_eligibility_all_battery,
    run_all_battery,
    run_by_name_battery,
)
from numpy.typing import NDArray

from src.logger import AppLogger

logger = AppLogger.get_logger(__name__)

RESULTS_FILE = "nist_test_results.csv"
STATUS_FILE = "nist_test_results_status.csv"
CHANNEL_NAMES = ("R", "G", "B")


def array_to_binary(flat_array: NDArray[np.uint8]) -> NDArray[np.uint8]:
    """Unpack a flat uint8 array into its bit sequence (8 bits per value)."""
    return np.unpackbits(flat_array)


class NistRunner:
    """Runs the NIST battery on images and accumulates results.

    Holds a score table and a pass/fail table keyed by (image, channel);
    re-running an existing pair overwrites its row instead of duplicating.
    """

    def __init__(self, csv_dir: str | Path | None = None) -> None:
        """Load existing result CSVs from `csv_dir` if present, else start empty."""
        self.csv_dir = Path(csv_dir) if csv_dir else None
        self.results_df = pd.DataFrame()
        self.pf_df = pd.DataFrame()
        if self.csv_dir and self.csv_dir.exists():
            self.results_df = pd.read_csv(self.csv_dir / RESULTS_FILE)
            self.pf_df = pd.read_csv(self.csv_dir / STATUS_FILE)

    def run_and_save(
        self, binary_seq: NDArray, image_name: str, channel_name: str, size: int
    ) -> None:
        """Run all eligible tests on `binary_seq` and store the row in memory."""
        eligible = check_eligibility_all_battery(binary_seq, SP800_22R1A_BATTERY)
        results = run_all_battery(binary_seq, eligible, False)

        scores = {"image": image_name, "channel": channel_name, "size": size}
        statuses = dict(scores)
        for result, _elapsed in results:
            scores[result.name] = np.round(result.score, 3)
            statuses[result.name] = "PASSED" if result.passed else "FAILED"
        self._upsert(scores, statuses)

    def _upsert(self, scores: dict, statuses: dict) -> None:
        """Overwrite the matching (image, channel) row, or append a new one."""
        if self.results_df.empty:
            mask = pd.Series(dtype=bool)
        else:
            mask = (self.results_df["image"] == scores["image"]) & (
                self.results_df["channel"] == scores["channel"]
            )
        if mask.any():
            logger.info("already exist, overwriting")
            for col in scores:
                self.results_df.loc[mask, col] = scores[col]
                self.pf_df.loc[mask, col] = statuses[col]
        else:
            logger.info("appending new result")
            self.results_df = pd.concat(
                [self.results_df, pd.DataFrame([scores])], ignore_index=True
            )
            self.pf_df = pd.concat(
                [self.pf_df, pd.DataFrame([statuses])], ignore_index=True
            )

    def process_image(self, img: NDArray, image_name: str) -> None:
        """Run the battery on channels 0, 1, 2 of the image, labeled R, G, B."""
        if img.ndim != 3 or img.shape[2] != 3:
            logger.error("Image %s is not a valid RGB image!", image_name)
            return
        size = img.shape[0]
        # index order kept from the notebook (cv2 images are BGR, see plans)
        for idx, channel in enumerate(CHANNEL_NAMES):
            logger.info("Running NIST tests for %s - %s channel", image_name, channel)
            bits = array_to_binary(img[:, :, idx].flatten())
            self.run_and_save(bits, image_name, channel, size)

    def save(self, directory: str | Path) -> None:
        """Write both result tables as CSV into `directory`."""
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
        self.results_df.to_csv(directory / RESULTS_FILE, index=False)
        self.pf_df.to_csv(directory / STATUS_FILE, index=False)
        logger.info("NIST test results saved to %s", directory)


def run_test(test_name: str, bit_seq: NDArray) -> None:
    """Run one named NIST test and log its score; logs if not eligible."""
    outcome = run_by_name_battery(
        test_name=test_name, bits=bit_seq, battery=SP800_22R1A_BATTERY
    )
    if outcome is None:
        logger.warning("%s test NOT ELIGIBLE", test_name)
        return
    result, elapsed_time = outcome
    status = "PASSED" if result.passed else "FAILED"
    logger.info(
        "%s - score: %s - %s - elapsed time: %s ms",
        status,
        np.round(result.score, 3),
        result.name,
        elapsed_time,
    )
