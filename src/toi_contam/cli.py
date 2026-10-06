import argparse

from .analysis import run_analysis


def parse_args():
    parser = argparse.ArgumentParser(description="Run the TESS contamination analysis.")
    parser.add_argument("--target", type=str, default=None,
                        help="Target folder or TOI name, e.g. TOI1011.01 or a path.")
    parser.add_argument("--custom-aperture", type=int, default=0,
                        help="Use a custom circular aperture radius in pixels. 0 means default behavior.")
    parser.add_argument("--dark", action="store_true",
                        help="Use the dark plotting style.")
    parser.add_argument("--output-dir", type=str, default="on-off",
                        help="Directory where output plots and files are saved.")
    return parser.parse_args()


def main():
    args = parse_args()
    return run_analysis(
        target=args.target,
        custom_aperture=args.custom_aperture,
        dark=args.dark,
        output_dir=args.output_dir,
    )
