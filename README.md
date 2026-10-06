# TOI contamination analysis

This project packages the TESS contamination script as an installable Python command.

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
