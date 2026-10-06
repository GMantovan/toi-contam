# TOI contamination analysis

This package performs a contamination analysis for a TESS target star (specifying a TOI).
It queries the Gaia DR3 catalogue around the target, evaluates nearby contaminant sources,
and estimates the *dilution* and *contamination factors* affecting the TESS photometric aperture 
("default" from the TPF file or a "custom" one, e.g. a circular aperture).

The tool produces summary tables and diagnostic plots that help assess whether a detected signal 
is likely on-source or potentially contaminated by nearby stars.

## Requirements

- Python 3.10
- A dedicated conda environment is recommended

## Install

Create and activate a dedicated conda environment before installing the package:

```bash
conda create -y -n toi-contam python=3.10
conda activate toi-contam
python -m pip install -e .
```

## Quick start

Run the analysis for a specific TOI target:

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
- `--custom-aperture`: to use a custom circular aperture, radius in pixels
- `--dark`: enable dark plotting style

## Outputs

The script writes contamination summaries and diagnostic plots into the selected output directory (see example in *example_folder*),
including the main catalogue-based contamination estimates and comparison figures.

## Reference

Please cite [Mantovan et al., 2022](https://academic.oup.com/mnras/article/516/3/4432/6692879?login=true) if you use this tool.
