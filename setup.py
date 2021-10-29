import sys

try:
    from skbuild import setup
except ImportError:
    print(
        "Please update pip, you need pip 10 or greater,\n"
        " or you need to install the PEP 518 requirements in pyproject.toml yourself",
        file=sys.stderr,
    )
    raise

from setuptools import find_packages

setup(
    name="pylimer_tools",
    version="0.1.0",
    description="A collection of utility python functions for handling LAMMPS output and polymers in Python ",
    author="Tim Bernhard",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    cmake_install_dir="src/pylimer_tools",
    cmake_args=[
        
    ],
    include_package_data=True,
    extras_require={"test": ["unittest"]},
    python_requires=">=3.7",
)
