# TOI contamination analysis

[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

This package performs a contamination analysis for a TESS target star (specifying a TOI).
It queries the Gaia DR3 catalogue around the target, evaluates nearby contaminant sources,
and estimates the *dilution factor* (DF) and *contamination ratio* (CR) affecting the TESS photometric aperture 
("default" from the TPF file or a "custom" one, e.g. a circular aperture).

The tool produces summary tables and diagnostic plots that help assess whether a detected signal 
is likely on-source or potentially contaminated by nearby stars.

Possible uses:

- To check if nearby stars could also generate the detected signal. If not, then the signal is 100% on-source. See 
  an example of a false positive (TOI-154.01) in *example_folder/toi-contam:example.ipynb*;
- To perform your own flux dilution corrections on TESS Simple Aperture Photometry (SAP) or your favourite TESS
  diluted photometry (e.g. PATHOS). This is useful if you want to preserve stellar variability, which is often flattened
  out in the TESS PDCSAP flux; 
- To perform proper dilution corrections in the case of binary stars. This is particularly helpful, for example,
  when you want to correct for flux contamination from the primary star and accurately estimate the planetary radius
  in cases where the secondary star is the stellar host.

### Definitions used in the code

For a target star with flux $F_{\star}$ and contaminant sources with fluxes $F_i$, the quantities used here are:

- Dilution factor:

  $DF = \frac{\sum_i F_i}{F_{\star} + \sum_i F_i}$

  where $F_{\star} + \sum_i F_i$ is the total flux in the aperture, and measures the fraction of flux that 
  comes from nearby contaminating stars.

- Contamination ratio:

  $CR = \frac{\sum_i F_i}{F_{\star}}$

  which measures how much contaminant flux is present relative to the target star itself.

## Release notes

### v1.0.0

- First public release of the TOI contamination analysis tool.
- Added support for Gaia DR3 based contamination estimates and TESS aperture diagnostics. Gaia to TESS magnitude 
  conversion is done following Consolini et al. (in prep).
- Included a command-line interface for quick TOI analysis runs.
- Improved handling of output folders, sector and custom aperture choices.

## Requirements

- Python 3.10+
- A conda environment is recommended
- Internet access for Gaia and MAST queries

## Installation

Create and activate a dedicated environment, then install the package in editable mode:

```bash
conda create -y -n toi-contam python=3.10
conda activate toi-contam
python -m pip install -e .
```

## Quick start

Run the analysis for a TOI:

```bash
toi-contam --target TOI1011.01 --output-dir TOI1011_contam
```

## Optional arguments

```bash
toi-contam --target TOI1011.01 --output-dir TOI1011_contam --custom-aperture 3 --sector 35 --dark
```

## Command-line options

- `--target`: target name or folder path
- `--output-dir`: output directory for generated files and plots
- `--custom-aperture`: to use a custom circular aperture, radius in pixels (default = tpf aperture)
- `--sector`: to perform the analysis for a specific TESS sector (default = first observed sector)
- `--dark`: enable dark plotting style

## Outputs

The script writes contamination summaries and diagnostic plots into the selected output directory (see example in *example_folder*),
including the main catalogue-based contamination estimates and comparison figures.

## Citation

Please cite [Mantovan et al., 2022](https://academic.oup.com/mnras/article/516/3/4432/6692879?login=true) if you use this tool.
