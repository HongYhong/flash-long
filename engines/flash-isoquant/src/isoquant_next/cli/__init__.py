"""Command-line entry points for IsoQuant 4.0 refactored (`isoquant-next`)."""

from isoquant_next.cli.detect_barcodes import main as detect_barcodes_main
from isoquant_next.cli.main import main as main
from isoquant_next.cli.visualize import main as visualize_main

__all__ = ["detect_barcodes_main", "main", "visualize_main"]
