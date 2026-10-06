# TOI contamination analysis

This project packages the TESS ON/OFF contamination script as an installable Python command.

## Install locally

```bash
python3 -m pip install -e .
```

## Run

```bash
toi-contam --target TOI1011.01 --output-dir on-off
```

You can also pass the other supported options, such as `--custom-aperture` and `--dark`.
