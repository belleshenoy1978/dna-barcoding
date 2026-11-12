#!/usr/bin/env python
"""Thin wrapper so the pipeline can be run from a clone without installing.

    python run_barcode.py --input traces/ --email you@lab.org

Once the package is installed (`pip install -e .`) the same interface is
available as the `barcode-pipeline` command.
"""
from barcode_pipeline.cli import main

if __name__ == "__main__":
    main()
