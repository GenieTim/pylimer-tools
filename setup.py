import multiprocessing
import os
import platform
import re
import shutil
import subprocess
import sys
import warnings
from pathlib import Path

from setuptools import Extension, find_namespace_packages, setup
from setuptools.command.build_ext import build_ext


def get_version():
    """Get version using version management system."""
    try:
        version_script = Path(__file__).parent / "bin" / "get_version.py"
        
        # Ensure version file exists for packaging
        result = subprocess.run(
            [sys.executable, str(version_script), "--ensure-file"],
            cwd=Path(__file__).parent,
            capture_output=True,
            text=True,
            timeout=10
        )
        
        if result.returncode == 0 and result.stdout.strip():
            version = result.stdout.strip()
            print(f"Using version: {version}")
            return version
            
    except (subprocess.SubprocessError, FileNotFoundError, OSError) as e:
        warnings.warn(f"Failed to get version from script: {e}")
    
    # Final fallback
    return "0.0.0"


VERSION = get_version()

# ==============================================================================
# CMake Build Configuration
# ==============================================================================

# Default to Release build for better performance
# Use RelWithDebInfo for debugging while maintaining reasonable performance
default_build_type = "Release"
build_type = os.environ.get("CMAKE_BUILD_TYPE", default_build_type)

cmake_args = [
    f"-DCMAKE_BUILD_TYPE={build_type}",
    f"-Dvendor_suffix=-skbuild-{platform.system()}",
    "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON",  # For IDE support
]

# Add any user-specified CMake arguments
if os.environ.get("CMAKE_ARGS"):
    cmake_args.extend(os.environ.get("CMAKE_ARGS", "").split())

print(f"Building with CMake build type: {build_type}")
print(f"Using {multiprocessing.cpu_count()} parallel jobs")

if os.getenv("VCPKG_ROOT"):
    toolchain_file = os.path.join(
        os.getenv("VCPKG_ROOT", ""), "scripts", "buildsystems", "vcpkg.cmake"
    )
    if os.path.isfile(toolchain_file):
        cmake_args.append(
            "-DCMAKE_TOOLCHAIN_FILE={}".format(toolchain_file.replace("\\", "/"))
        )
        # cmake_args.append("-DVCPKG_TARGET_TRIPLET=x86-windows-static")
        print('Using toolchain "{}"'.format(toolchain_file))
    else:
        warnings.warn(
            "Detected VCPKG_ROOT. Did not find toolchain file {} though.".format(
                toolchain_file
            )
        )
else:
    print("VCPKG_ROOT not set. Not using vcpk dependencies.")

# check env to decide whether we should add high performance flags
if os.getenv("HIGH_PERFORMANCE", False):
    cmake_args.append("-DHIGH_PERFORMANCE=ON")
else:
    cmake_args.append("-DHIGH_PERFORMANCE=OFF")
if os.getenv("NO_SERIALIZATION", False):
    cmake_args.append("-DCEREALIZABLE=OFF")

# delete vendor caches — this is useful if you compile
# this project using CMake (e.g. for tests) as well as skbuild,
# as the two build directories of vendor do not interact well.
vendor_files_to_delete = [
    os.path.abspath(
        os.path.join(
            os.path.dirname(__file__),
            "vendor/igraph-skbuild-{}/src/igraphLib-build".format(platform.system()),
        )
    ),
    os.path.abspath(
        os.path.join(
            os.path.dirname(__file__),
            "vendor/nlopt-skbuild-{}/src/nloptLib-build".format(platform.system()),
        )
    ),
    os.path.abspath(
        os.path.join(
            os.path.dirname(__file__),
            "vendor/cereal-skbuild-{}/src/cerealLib-build".format(platform.system()),
        )
    ),
]
for vendor_file in vendor_files_to_delete:
    if os.path.exists(vendor_file):
        try:
            shutil.rmtree(vendor_file)
        except Exception:
            warnings.warn(
                "Could not delete directory {}. Errors incoming.".format(vendor_file)
            )
    else:
        print("No need to delete {}".format(vendor_file))

# skbuildCaches = os.path.abspath(os.path.join(
#     os.path.dirname(__file__), '_skbuild'))
# if (os.path.exists(skbuildCaches)):
#     try:
#         shutil.rmtree(skbuildCaches)
#     except:
#         warnings.warn(
#             "Could not delete directory {}. Errors incoming.".format(skbuildCaches))

with open("README.md", "r", encoding="utf-8") as file:
    readme_content = file.read()

# Convert distutils Windows platform specifiers to CMake -A arguments
PLAT_TO_CMAKE = {
    "win32": "Win32",
    "win-amd64": "x64",
    "win-arm32": "ARM",
    "win-arm64": "ARM64",
}


# A CMakeExtension needs a sourcedir instead of a file list.
# The name must be the _single_ output extension from the CMake build.
# If you need multiple extensions, see scikit-build.
class CMakeExtension(Extension):
    def __init__(self, name: str, sourcedir: str = "") -> None:
        super().__init__(name, sources=[])
        self.sourcedir = os.fspath(Path(sourcedir).resolve())


class CMakeBuild(build_ext):
    def _convert_excel_to_json(self):
        """Convert Excel data to JSON if needed during build."""
        project_root = Path(__file__).parent
        excel_path = project_root / "src" / "pylimer_tools" / "data" / "everaers_et_al_unit_properties.xlsx"
        json_path = project_root / "src" / "pylimer_tools" / "data" / "everaers_et_al_unit_properties.json"
        convert_script = project_root / "bin" / "convert-excel-to-json.py"
        
        # Check if conversion is needed (Excel newer than JSON or JSON doesn't exist)
        if (not json_path.exists() or 
            (excel_path.exists() and excel_path.stat().st_mtime > json_path.stat().st_mtime)):
            
            print("Converting Excel data to JSON...")
            try:
                subprocess.run([
                    sys.executable, str(convert_script)
                ], cwd=project_root, check=True, capture_output=True, text=True)
                print("Excel to JSON conversion completed successfully.")
            except subprocess.CalledProcessError as e:
                print(f"Warning: Excel to JSON conversion failed: {e}")
                print(f"stdout: {e.stdout}")
                print(f"stderr: {e.stderr}")
                # Don't fail the build if conversion fails and JSON already exists
                if not json_path.exists():
                    raise RuntimeError("JSON file doesn't exist and conversion failed")
            except FileNotFoundError:
                print("Warning: Conversion script not found, skipping Excel to JSON conversion")
                if not json_path.exists():
                    raise RuntimeError("JSON file doesn't exist and conversion script not found")

    def build_extension(self, ext: CMakeExtension) -> None:
        # Convert Excel to JSON if needed before building
        self._convert_excel_to_json()
        
        # Must be in this form due to bug in .resolve() only fixed in Python
        # 3.10+
        ext_fullpath = Path.cwd() / self.get_ext_fullpath(ext.name)
        extdir = ext_fullpath.parent.resolve()

        # Using this requires trailing slash for auto-detection & inclusion of
        # auxiliary "native" libs

        debug = int(os.environ.get("DEBUG", 0)) if self.debug is None else self.debug
        cfg = "Debug" if debug else "RelWithDebInfo"  # "Release"

        # The CMake build process is done in a temporary directory,
        # so we have to create it here.
        build_temp = Path(self.build_temp) / ext.name
        if not build_temp.exists():
            build_temp.mkdir(parents=True)
        print("Building in {}".format(build_temp))

        # CMake lets you override the generator - we need to check this.
        # Can be set with Conda-Build, for example.
        cmake_generator = os.environ.get("CMAKE_GENERATOR", "")

        # Set Python_EXECUTABLE instead if you use PYBIND11_FINDPYTHON
        # EXAMPLE_VERSION_INFO shows you how to pass a value into the C++ code
        # from Python.
        global cmake_args
        cmake_args.extend(
            [
                f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}{os.sep}",
                f"-DPYTHON_EXECUTABLE={sys.executable}",
                "-DCODE_COVERAGE=OFF",
                "-DLEAK_ANALYSIS=OFF",
                f"-DCMAKE_BUILD_TYPE={cfg}",  # not used on MSVC, but no harm
            ]
        )
        build_args = []
        # Adding CMake arguments set as environment variable
        # (needed e.g. to build for ARM OSx on conda-forge)
        if "CMAKE_ARGS" in os.environ:
            cmake_args += [item for item in os.environ["CMAKE_ARGS"].split(" ") if item]

        # Detect Pyodide / Emscripten build environment
        is_emscripten = (
            os.environ.get("PYODIDE_BUILD")
            or sys.platform.startswith("emscripten")
            or "emscripten" in str(sys.platform)
        )
        if is_emscripten:
            cmake_args.append("-DPYODIDE_BUILD=ON")
            print("Detected Pyodide/Emscripten build environment, setting PYODIDE_BUILD=ON")

        # In this example, we pass in the version to C++. You might not need
        # to.
        cmake_args += [f"-DVERSION_NR={self.distribution.get_version()}"]

        if self.compiler.compiler_type != "msvc":
            # Using Ninja-build since it a) is available as a wheel and b)
            # multithreads automatically. MSVC would require all variables be
            # exported for Ninja to pick it up, which is a little tricky to do.
            # Users can override the generator with CMAKE_GENERATOR in CMake
            # 3.15+.
            if (
                not cmake_generator
                or cmake_generator == "Ninja"
                or cmake_generator == ""
            ):
                # import ninja
                # ninja_executable_path = Path(ninja.BIN_DIR) / "ninja"
                import shutil

                ninja_executable_path = shutil.which("ninja")

                if ninja_executable_path:
                    try:
                        print(
                            "Checking whether ninja can be run at {}".format(
                                ninja_executable_path
                            )
                        )
                        subprocess.run(
                            [ninja_executable_path, "--version"],
                            check=True,
                            timeout=5.0,
                            cwd=build_temp,
                        )
                        cmake_args += [
                            "-GNinja",
                            f"-DCMAKE_MAKE_PROGRAM:FILEPATH={ninja_executable_path}",
                        ]
                        using_ninja = True
                        print(
                            "Using Ninja generator (parallelization handled automatically)"
                        )
                    except (
                        ImportError,
                        subprocess.CalledProcessError,
                        PermissionError,
                        subprocess.TimeoutExpired,
                    ):
                        warnings.warn(
                            "Ninja check did not pass, using default generator."
                        )
                        pass
                else:
                    warnings.warn("Ninja is not available, using default generator.")

        else:
            # Single config generators are handled "normally"
            single_config = any(x in cmake_generator for x in {"NMake", "Ninja"})

            # CMake allows an arch-in-generator style for backward
            # compatibility
            contains_arch = any(x in cmake_generator for x in {"ARM", "Win64"})

            # Specify the arch if using MSVC generator, but only if it doesn't
            # contain a backward-compatibility arch spec already in the
            # generator name.
            if not single_config and not contains_arch:
                cmake_args += ["-A", PLAT_TO_CMAKE[self.plat_name]]

            # Multi-config generators have a different way to specify configs
            if not single_config:
                cmake_args += [
                    f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{cfg.upper()}={extdir}"
                ]
                build_args += ["--config", cfg]

        if sys.platform.startswith("darwin"):
            # Cross-compile support for macOS - respect ARCHFLAGS if set
            archs = re.findall(r"-arch (\S+)", os.environ.get("ARCHFLAGS", ""))
            if archs:
                cmake_args += ["-DCMAKE_OSX_ARCHITECTURES={}".format(";".join(archs))]

        # Handle parallelization settings with clear priority to avoid conflicts:
        # 1. CMAKE_BUILD_PARALLEL_LEVEL (as an environment variable)
        # 2. User's self.parallel setting
        # 3. The number of CPU cores available
        if not os.environ.get("CMAKE_BUILD_PARALLEL_LEVEL"):
            if hasattr(self, "parallel") and self.parallel:
                parallel_level = self.parallel
            else:
                parallel_level = min(1, multiprocessing.cpu_count()-1)
            if parallel_level > 1:
                # Set the CMake variable for parallel builds
                cmake_args.append(f"-DCMAKE_BUILD_PARALLEL_LEVEL={parallel_level}")
                os.environ["CMAKE_BUILD_PARALLEL_LEVEL"] = str(parallel_level)
                print(f"Using CMAKE_BUILD_PARALLEL_LEVEL={parallel_level}")

        # Actually run CMake
        configure_cmd = ["cmake", ext.sourcedir, *cmake_args]
        if is_emscripten:
            print("Using emcmake for Emscripten configuration")
            configure_cmd = ["emcmake"] + configure_cmd
        subprocess.run(configure_cmd, cwd=build_temp, check=True)
        subprocess.run(
            ["cmake", "--build", ".", *build_args], cwd=build_temp, check=True
        )


setup(
    name="pylimer_tools",
    version=VERSION,
    description="A toolkit for handling bead-spring polymers and LAMMPS output in Python",
    long_description_content_type="text/markdown",
    long_description=readme_content,
    keywords=["Polymer", "Chemistry", "Network", "LAMMPS", "Science"],
    author="Tim Bernhard",
    author_email="tim@bernhard.dev",
    url="https://github.com/GenieTim/pylimer-tools",
    license="GPL-3.0-or-later",
    packages=find_namespace_packages(where="src", exclude=("tests",)),
    package_dir={"": "src"},
    include_package_data=True,
    extras_require={"test": ["unittest"]},
    python_requires=">=3.9",
    ext_modules=[CMakeExtension("pylimer_tools_cpp")],
    cmdclass={"build_ext": CMakeBuild},
    entry_points={
        "console_scripts": [
            "pylimer-generate-network=pylimer_tools.generate_network:cli",
            "pylimer-analyse-networks=pylimer_tools.analyse_networks:cli",
            "pylimer-basic-lammps-stats=pylimer_tools.basic_lammps_structure_stats:cli",
            "pylimer-displace-randomly=pylimer_tools.displace_randomly:cli",
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
        "Operating System :: OS Independent",
    ],
)
