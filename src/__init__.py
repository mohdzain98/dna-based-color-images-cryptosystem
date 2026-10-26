"""DNA based color image cryptosystem package.

Holds the core cipher (cryptosystem), analysis metrics and shared logging.
Importing the package configures logging once so notebooks and scripts get
consistent console output.
"""

from src.logger import AppLogger

AppLogger.init()
