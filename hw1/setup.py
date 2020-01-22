import sys
from setuptools import setup, Extension, find_packages
from shutil import rmtree, copyfile
import glob
import os


class getPybindInclude(object):
    """Helper class to determine the pybind11 include path
    The purpose of this class is to postpone importing pybind11
    until it is actually installed, so that the ``get_include()``
    method can be invoked.
    https://github.com/pybind/python_example/blob/master/setup.py
    """

    def __init__(self, user=False):
        self.user = user

    def __str__(self):
        import pybind11
        return pybind11.get_include(self.user)


def getIncludes():
    return [
        getPybindInclude(),
        getPybindInclude(user=True)
    ]

sources = [
    'cParticle.cpp',
    'particles.cpp',
]

extra_compile_args = ["--std=c++11", "-fPIC", "-v", "-O3"]
extra_link_args = None
if sys.platform == 'darwin':
    extra_compile_args += ['-stdlib=libc++', '-mmacosx-version-min=10.9']
    extra_link_args = ["-v", '-mmacosx-version-min=10.9']

    # extra_compile_args += ['-mmacosx-version-min=10.9', '-stdlib=libc++']

module = Extension(
    'cParticle',
    # include_dirs=getIncludes(),
    extra_compile_args=extra_compile_args,
    extra_link_args = extra_link_args,
    sources=sources
)

setup(
    ext_modules=[module],
)

if sys.argv[-1] == "build":
    # copy the build module to the local directory
    sopath = glob.glob("build/lib*/cParticle*")[0]
    sofile = os.path.split(sopath)[-1]
    dest = os.path.join(os.getcwd(), sofile)
    copyfile(sopath, dest)
