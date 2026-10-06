# TOI contamination analysis

This package performs a TESS ON/OFF contamination analysis for a target star or TOI.
It queries Gaia and TIC catalogs around the target, evaluates nearby contaminant sources,
and estimates the dilution and contamination levels affecting the photometric aperture.

The tool is packaged as a small Python CLI and produces summary tables and diagnostic plots
that help assess whether a detected signal is likely on-source or contaminated by nearby stars.

## Requirements

- Python 3.10
- A dedicated conda environment is recommended
- The astronomy dependencies needed by the analysis workflow

## Install

Create and activate a dedicated conda environment before installing the package:

```bash
conda create -y -n toi-contam python=3.10
conda activate toi-contam
python -m pip install -e .
```

## Quick start

Run the analysis for a TOI target:

```bash
toi-contam --target TOI1011.01 --output-dir on-off
```

## Optional arguments

```bash
toi-contam --target TOI1011.01 --output-dir on-off --custom-aperture 3 --dark
```

Available options include:

- `--target`: target name or folder path
- `--output-dir`: output directory for generated files and plots
- `--custom-aperture`: use a custom aperture radius in pixels
- `--dark`: enable dark plotting style

## Outputs

The script writes contamination summaries and diagnostic plots into the selected output directory,
including the main catalog-based contamination estimates and comparison figures.
