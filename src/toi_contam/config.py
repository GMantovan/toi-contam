from __future__ import annotations

from dataclasses import dataclass


@dataclass
class AnalysisConfig:
    """Configuration used by the package entry point and the legacy script wrapper."""

    target: str | None = None
    output_dir: str = "on-off"
    custom_aperture: int = 0
    dark: bool = False
    search_radius_pixels: int = 10
    pixel_size_arcsec: int = 21
    prf_col_offset: int = 44
    prf_image_size: int = 13

    @classmethod
    def from_arguments(
        cls,
        *,
        target: str | None = None,
        custom_aperture: int = 0,
        dark: bool = False,
        output_dir: str = "on-off",
    ) -> "AnalysisConfig":
        return cls(
            target=target,
            output_dir=output_dir,
            custom_aperture=int(custom_aperture),
            dark=bool(dark),
        )

    def as_argv(self) -> list[str]:
        argv: list[str] = []
        if self.target is not None:
            argv.extend(["--target", str(self.target)])
        if self.custom_aperture:
            argv.extend(["--custom-aperture", str(self.custom_aperture)])
        if self.dark:
            argv.append("--dark")
        if self.output_dir != "on-off":
            argv.extend(["--output-dir", str(self.output_dir)])
        return argv
