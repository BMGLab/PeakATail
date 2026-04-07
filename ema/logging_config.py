import logging
import sys
import time
from typing import Optional


def setup_logging(level: str = "INFO", log_file: Optional[str] = None) -> logging.Logger:
    """Configure structured logging for PeakATail.

    Args:
        level: Logging level (DEBUG, INFO, WARNING, ERROR)
        log_file: Optional path to log file (in addition to stderr)

    Returns:
        Root logger for peakatail
    """
    logger = logging.getLogger("peakatail")
    logger.setLevel(getattr(logging, level.upper(), logging.INFO))

    # Clear existing handlers
    logger.handlers.clear()

    formatter = logging.Formatter(
        "%(asctime)s | %(levelname)-5s | %(name)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )

    # stderr handler
    stderr_handler = logging.StreamHandler(sys.stderr)
    stderr_handler.setFormatter(formatter)
    logger.addHandler(stderr_handler)

    # optional file handler
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    return logger


class PeakCallingLogger:
    """Structured logger for peak calling progress and decisions."""

    def __init__(self, strategy_name: str, direction: bool):
        self.logger = logging.getLogger("peakatail.peakcalling")
        self.strategy = strategy_name
        self.direction = "negative" if direction else "positive"
        self.start_time = time.time()
        self.reads_processed = 0
        self.peaks_found = 0
        self.pas_found = 0
        self.current_chrom = None

    def log_start(self):
        self.logger.info(
            f"peak_calling started | strategy={self.strategy} | strand={self.direction}"
        )

    def log_chromosome(self, chrom: str, peaks: int):
        if self.current_chrom and self.current_chrom != chrom:
            self.logger.info(
                f"chromosome={self.current_chrom} | peaks_found={peaks}"
            )
        self.current_chrom = chrom

    def log_peak_decision(self, peak_id: int, height: int, lambda_local: float,
                          p_value: float, significant: bool, pas_count: int = 0):
        level = logging.DEBUG
        status = "significant" if significant else "filtered"
        self.logger.log(level,
            f"peak_id={peak_id} | height={height} | lambda={lambda_local:.1f} | "
            f"p_value={p_value:.2e} | {status} | pas_count={pas_count}"
        )

    def log_read_progress(self, count: int):
        self.reads_processed = count
        if count % 1_000_000 == 0:
            elapsed = time.time() - self.start_time
            self.logger.info(
                f"reads_processed={count:,} | elapsed={elapsed:.1f}s | "
                f"strand={self.direction}"
            )

    def log_complete(self, total_peaks: int, total_pas: int):
        elapsed = time.time() - self.start_time
        self.logger.info(
            f"run_complete | strategy={self.strategy} | strand={self.direction} | "
            f"total_peaks={total_peaks} | total_pas={total_pas} | "
            f"runtime={elapsed:.1f}s"
        )
