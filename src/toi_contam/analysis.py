from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path


def _load_script_module():
    script_path = Path(__file__).resolve().parents[2] / "toi-contam_dr3.py"
    spec = spec_from_file_location("toi_contam_script", script_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load script from {script_path}")
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run_analysis(target=None, custom_aperture=0, dark=False, output_dir="on-off"):
    module = _load_script_module()
    argv = []

    if target is not None:
        argv.extend(["--target", str(target)])
    if custom_aperture:
        argv.extend(["--custom-aperture", str(custom_aperture)])
    if dark:
        argv.append("--dark")
    if output_dir != "on-off":
        argv.extend(["--output-dir", str(output_dir)])

    return module.main(argv)
