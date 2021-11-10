import os
import sys
import warnings

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

cmake_args = []

if (os.getenv('VCPKG_ROOT')):
    toolchainFile = os.path.join(
        os.getenv('VCPKG_ROOT'), "scripts", "buildsystems", "vcpkg.cmake")
    if (os.path.isfile(toolchainFile)):
        cmake_args.append(
            "-DCMAKE_TOOLCHAIN_FILE={}".format(toolchainFile.replace("\\", "/")))
        # cmake_args.append("-DVCPKG_TARGET_TRIPLET=x86-windows-static")
        print("Using toolchain \"{}\"".format(toolchainFile))
    else:
        warnings.warn("Detected VCPKG_ROOT. Did not find toolchain file {} though.".format(toolchainFile))


setup(
    name="pylimer_tools",
    version="0.1.0",
    description="A collection of utility python functions for handling LAMMPS output and polymers in Python ",
    author="Tim Bernhard",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    cmake_install_dir="src/pylimer_tools_cpp",
    cmake_args=cmake_args,
    include_package_data=True,
    extras_require={"test": ["unittest"]},
    python_requires=">=3.8"
)
