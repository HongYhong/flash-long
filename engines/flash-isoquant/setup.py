from __future__ import annotations

import os
from pathlib import Path

import numpy as np
from Cython.Build import cythonize
from setuptools import Extension, find_packages, setup

ROOT = Path(__file__).resolve().parent
SRC_DIR = ROOT / "src"


def build_extensions() -> list[Extension]:
    extensions: list[Extension] = []
    compiler_args = ["-O3", "-fno-strict-aliasing"]
    include_dirs = [np.get_include()]

    for pyx_path in sorted((SRC_DIR / "isoquant_next").rglob("*.pyx")):
        rel = pyx_path.relative_to(SRC_DIR).with_suffix("")
        mod_name = ".".join(rel.parts)
        extensions.append(
            Extension(
                name=mod_name,
                sources=[str(pyx_path.relative_to(ROOT))],
                include_dirs=include_dirs,
                extra_compile_args=compiler_args,
            )
        )
    return extensions


os.environ.setdefault("CC", "/usr/bin/gcc")

setup(
    name="isoquant-next",
    version="0.1.0",
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    package_data={"isoquant_next": ["data/*.json"]},
    ext_modules=cythonize(
        build_extensions(),
        compiler_directives={
            "language_level": "3",
            "boundscheck": False,
            "wraparound": False,
            "cdivision": True,
            "initializedcheck": False,
        },
        nthreads=max(1, min(4, os.cpu_count() or 1)),
    ),
)
