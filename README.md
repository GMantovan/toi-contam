# TOI contamination analysis

This package performs a TESS contamination analysis for a target star or TOI. It queries the Gaia catalog around the target, evaluates nearby contaminant sources, and computes the expected contamination contribution using the target aperture and nearby stellar sources. The code produces plots and analysis outputs that can be used to assess whether a detected signal is likely on source. Moreover, it estimates the contamination and dilution factors.

## Install locally

Before installing the package, create and activate a dedicated conda environment:

```bash
conda create -y -n toi-contam python=3.10
conda activate toi-contam
python -m pip install -e .
```

## Run

```bash
toi-contam --target TOI1011.01 --output-dir on-off
```

You can also pass the other supported options, such as `--custom-aperture` and `--dark`.
