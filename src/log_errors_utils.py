"""Handling errors"""

import logging
import os
import sys
from pathlib import Path


class CustomException(Exception):
    """Exceptions handling."""

    def __init__(self, msg):
        # Call the base class constructor with the custom message
        super().__init__(msg)


def _check_file(in_file, log=False, msg="FILE"):
    try:
        if not os.path.isfile(in_file):
            raise CustomException("File is missing")
        if os.path.getsize(in_file) == 0:
            process_warning(f"{msg}\t{in_file}: is empty")
    except Exception as e:
        process_exception(f"{msg}\t{in_file}: {e}")
    else:
        if log:
            logging.info(f"{msg}\t{in_file}")


def check_file(in_file):
    _check_file(in_file, log=False)


def log_file(in_file):
    _check_file(in_file, log=True)


def process_exception(msg):
    logging.exception(msg)
    print(f"EXCEPTION\t{msg}", file=sys.stderr)
    sys.exit(1)


def process_error(msg):
    logging.error(msg)
    print(f"ERROR\t{msg}", file=sys.stderr)
    sys.exit(1)


def process_warning(msg):
    logging.warning(msg)
    print(f"WARNING\t{msg}", file=sys.stderr)


# ---------------------------------------------------------------------------- #
#                        Files And Directories Function                        #
# ---------------------------------------------------------------------------- #
def create_directory(in_dir_list: list[Path]) -> None:
    """Create directories."""
    for in_dir in in_dir_list:
        if in_dir:
            in_dir.mkdir(parents=True, exist_ok=True)
