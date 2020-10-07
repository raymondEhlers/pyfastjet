"""
Adapted from https://github.com/pybind/cmake_example
"""
import os
import platform
import re
import shutil
import subprocess
import sys
import sysconfig
from distutils.command.build_ext import build_ext
from distutils.core import Distribution, Extension
from distutils.version import LooseVersion
from typing import Any, Dict

#from numpy import get_include as get_numpy_include
#from setuptools.command.build_ext import build_ext
#from setuptools.extension import Extension
#from setuptools import find_packages


class CMakeExtension(Extension):
    name: str  # exists, even though IDE doesn't find it

    def __init__(self, name: str, sourcedir: str="") -> None:
        super().__init__(name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class ExtensionBuilder(build_ext):
    def run(self) -> None:
        self.validate_cmake()
        super().run()

    def build_extension(self, ext: Extension) -> None:
        if isinstance(ext, CMakeExtension):
            self.build_cmake_extension(ext)
        else:
            super().build_extension(ext)

    def validate_cmake(self) -> None:
        print("About to run validation")
        cmake_extensions = [x for x in self.extensions if isinstance(x, CMakeExtension)]
        if len(cmake_extensions) > 0:
            try:
                out = subprocess.run(["cmake", "--version"], check=True, stdout=subprocess.PIPE).stdout
            except OSError:
                raise RuntimeError(
                    "CMake must be installed to build the following extensions: "
                    + ", ".join(e.name for e in cmake_extensions)
                )
            if platform.system() == "Windows":
                cmake_version = LooseVersion(re.search(r"version\s*([\d.]+)", out.decode()).group(1))  # type: ignore
                if cmake_version < "3.1.0":
                    raise RuntimeError("CMake >= 3.1.0 is required on Windows")

    def build_cmake_extension(self, ext: CMakeExtension) -> None:
        print("About to build")
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        cmake_args = ["-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=" + extdir, "-DPYTHON_EXECUTABLE=" + sys.executable]

        cfg = "Debug" if self.debug else "Release"
        # cfg = 'Debug'
        build_args = ["--config", cfg]

        if platform.system() == "Windows":
            cmake_args += ["-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}".format(cfg.upper(), extdir)]
            if sys.maxsize > 2 ** 32:
                cmake_args += ["-A", "x64"]
            build_args += ["--", "/m"]
        else:
            cmake_args += ["-DCMAKE_BUILD_TYPE=" + cfg]
            build_args += ["--", "-j4"]
        cmake_args += ["-DPYTHON_INCLUDE_DIR={}".format(sysconfig.get_path("include"))]

        env = os.environ.copy()
        env["CXXFLAGS"] = '{} -DVERSION_INFO=\\"{}\\"'.format(env.get("CXXFLAGS", ""), self.distribution.get_version())

        # Fastjet
        # Grab the fastjet libraries.
        # FASTJET, CGAL, and GMP are needed. It's up to the user to set them (for now).
        # Yes, this is a huge pain...
        # alibuild updated the fastjet path...
        #cmake_args += [f'-DFASTJET={env["FASTJET_ROOT"]}']
        cmake_args += [f'-DFASTJET={env["FASTJET"]}']
        cmake_args += [f'-DCGAL={env["CGAL_ROOT"]}']
        cmake_args += [f'-DGMP={env["GMP_ROOT"]}']
        #cmake_args += [f'-DGMP=/Users/re239/alice/sw/osx_x86-64/GMP/v6.0.0-1/']
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        print("Configurtion...")
        print(f"Source dir: {ext.sourcedir} with: {os.listdir(ext.sourcedir)}")
        print(f"self.build_temp: {self.build_temp} with {os.listdir(self.build_temp)}")
        print(f"cwd: {os.getcwd()}")
        # TODO: Maybe glob this directory?
        # Somehow, I need to run the cwd from the build dir, but scan the package directory.
        #subprocess.run(["cmake", ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env, check=True)
        # Source directory = -S, build directory = -B
        # cmake -S . -B build
        # cmake --build build
        #subprocess.run(["cmake", "-S", ext.sourcedir] + cmake_args, cwd=ext.sourcedir, env=env, check=True)
        subprocess.run(["cmake", "-S", ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env, check=True)
        print("Building...")
        subprocess.run(["cmake", "--build", "."] + build_args, cwd=self.build_temp, check=True)


def build(setup_kwargs: Dict[str, Any]) -> None:
    #cmake_modules = [CMakeExtension("project.package.pybind11_extension", sourcedir="project/package/pybind11_extension")]
    # 29 September 2020: Apparently setting the sourcedir, even if it's in a subdirectory, will breaks things. So don't set it!
    # 2 October 2020 TODO: Update fully to new poetry build system with:
    #   - https://github.com/sdispater/pendulum/pull/488/files
    #   - https://github.com/python-poetry/poetry/issues/2740
    cmake_modules = [CMakeExtension("pyfastjet._src")]
    ext_modules = cmake_modules
    #setup_kwargs.update({
    #    "packages": find_packages(),
    #    "ext_modules": ext_modules,
    #    "cmdclass": dict(build_ext=ExtensionBuilder),
    #    "zip_safe": False,
    #})

    distribution = Distribution({"name": "pyfastjet", "ext_modules": ext_modules})
    distribution.package_dir = "pyfastjet"

    # More or less directly from https://github.com/python-poetry/poetry/issues/2740
    cmd = ExtensionBuilder(distribution)
    cmd.ensure_finalized()
    cmd.run()

    # Copy built extensions back to the project
    for output in cmd.get_outputs():
        relative_extension = os.path.relpath(output, cmd.build_lib)
        print(f"Copying {output} to {relative_extension}")
        shutil.copyfile(output, relative_extension)
        mode = os.stat(relative_extension).st_mode
        mode |= (mode & 0o444) >> 2
        os.chmod(relative_extension, mode)

    #sys.exit(1)

    return setup_kwargs

if __name__ == "__main__":
    build({})
