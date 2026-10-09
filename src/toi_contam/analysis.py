from __future__ import annotations

import runpy
from pathlib import Path

from .config import AnalysisConfig


def _legacy_script_path() -> Path:
    return Path(__file__).resolve().parents[2] / "toi-contam_dr3.py"


def run_analysis(target=None, custom_aperture=0, sector=None, dark=False, output_dir="on-off"):
    """Run the scientific analysis through the package entry point.

    The long-term goal is to move the actual analysis code into the package itself,
    but this keeps the project installable and centralizes the configuration in a
    dedicated dataclass while the legacy script remains as a compatibility layer.
    """
    config = AnalysisConfig.from_arguments(
        target=target,
        custom_aperture=custom_aperture,
        sector=sector,
        dark=dark,
        output_dir=output_dir,
    )
    script_path = _legacy_script_path()
    namespace = runpy.run_path(str(script_path))

    main = namespace.get("main")
    if main is None:
        raise RuntimeError(f"Could not find main() in legacy script: {script_path}")

    return main(config.as_argv())
