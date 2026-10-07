# TOI contamination analysis

This package performs a contamination analysis for a TESS target star (specifying a TOI).
It queries the Gaia DR3 catalogue around the target, evaluates nearby contaminant sources,
and estimates the *dilution* and *contamination factors* affecting the TESS photometric aperture 
("default" from the TPF file or a "custom" one, e.g. a circular aperture).

The tool produces summary tables and diagnostic plots that help assess whether a detected signal 
is likely on-source or potentially contaminated by nearby stars.

Possible uses:

- To check if nearby stars could also generate the detected signal. If not, then the signal is 100% on-source;
- To perform your own flux dilution corrections on TESS Simple Aperture Photometry (SAP) or your favourite TESS
  diluted photometry (e.g. PATHOS). This is useful if you want to preserve stellar variability, which is often flattened
  out in the TESS PDCSAP flux; 
- To perform proper dilution corrections in the case of binary stars. This is particularly helpful, for example,
  when you want to correct for flux contamination from the primary star and accurately estimate the planetary radius
  in cases where the secondary star is the stellar host.

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
- `--custom-aperture`: to use a custom circular aperture, radius in pixels (default = tpf aperture)
- `--sector`: to perform the analysis for a specific TESS sector (default = first observed sector)
- `--dark`: enable dark plotting style

## Outputs

The script writes contamination summaries and diagnostic plots into the selected output directory (see example in *example_folder*),
including the main catalogue-based contamination estimates and comparison figures.

## Reference

Please cite [Mantovan et al., 2022](https://academic.oup.com/mnras/article/516/3/4432/6692879?login=true) if you use this tool.
