
from setuptools import setup
from setuptools.extension import Extension

import codecs
import os.path
import sys
import re
import shlex
try:
    from shlex import quote as shlexquote
except ImportError:
    from pipes import quote as shlexquote


try:
    import numpy # prerequisite, numpy.get_include()
except ImportError:
    raise RuntimeError("ERROR: Please install `numpy` first.")

try:
    import pybind11 # prerequisite, pybind11.get_include()
except ImportError:
    raise RuntimeError("ERROR: Please install `pybind` first.")

try:
    import tomographer # prerequisite
except ImportError:
    raise RuntimeError(
        "ERROR: Please install `tomographer` first. "
        "(See https://tomographer.github.io/tomographer/get-started/#python-version )"
    )

# version check
import tomographer.version
if tomographer.version.version_info < (5, 3):
    raise RuntimeError("ERROR: QPtomographer requires `tomographer` >= 5.3. "
                       "Please upgrade `tomographer`.")

import tomographer.include # tomographer.include.get_include()
from tomographer.include import find_include_dir, find_lib, Vars, ensure_str

import tomographer.version # compile_info

vv = Vars([
    'SCS_ROOT',
    'LAPACK_LDFLAGS',
    'CXX_FLAGS'
])

#
# Defaults for LAPACK_LDFLAGS, CXX_FLAGS
#

if sys.platform == 'darwin':
    vv.setDefault('LAPACK_LDFLAGS', '-framework Accelerate')
else:
    vv.setDefault('LAPACK_LDFLAGS', '-llapack -lblas')

vv.setDefault('CXX_FLAGS', " ".join((shlexquote(x) for x in tomographer.version.compile_info['cflags'])))


if not vv.get('SCS_ROOT'):
    raise RuntimeError("Please specify where to look for the SCS library and headers "
                       "using the environment variable SCS_ROOT. Check out the README "
                       "for more info at https://github.com/Tomographer/QPtomographer .")

if not os.path.isdir(vv.get('SCS_ROOT')):
    raise RuntimeError("$SCS_ROOT='" + vv.get('SCS_ROOT') + "' does not point to a "
                       "valid directory. Please specify "
                       "the path to the SCS library and headers using the environment "
                       "variable SCS_ROOT. Check out the README file for more info.")


print(vv.message())



#
# Include directories and other compilation flags -- use same as those used to
# compile tomographer itself
#
include_dirs = [
    numpy.get_include(),
    pybind11.get_include(),
] + tomographer.include.get_include() + [
    os.path.join(vv.get('SCS_ROOT'), 'include'),
    vv.get('SCS_ROOT'),
]
compiler_path = tomographer.version.compile_info.get('compiler',{}).get('path', None)
cflags = shlex.split(vv.get('CXX_FLAGS'))
ldflags = [
    # link against libscsdir.a
    os.path.join(vv.get('SCS_ROOT'), 'out', 'libscsdir.a'),
    # add flags for blas/lapack here, needed for SCS:
    ] + vv.get('LAPACK_LDFLAGS').split() + [
]


#
# Header dependencies
#
common_headers = [
    'cxx/channelspace.h',
    'cxx/utils.h'
    'cxx/diamond_norm_scs.h'
    'cxx/diamond_norm_figofmerit.h'
]


#
# Read version info
#
from QPtomographer._version import version


with open('README.md', 'r') as f:
    README = f.read()

#
# Set up the python package
#

setup(name="QPtomographer",
      version=version,
      author='Philippe Faist',
      author_email='phfaist@caltech.edu',
      description='Reliable Diamond Norm Estimation for Quantum Process Tomography',
      long_description=README,
      packages=[
          'QPtomographer'
      ],
      ext_modules=[
          Extension(name="QPtomographer.channelspace",
                    sources=['cxx/pydnormchannelspace.cxx'],
                    library_dirs=[],
                    libraries=[],
                    include_dirs=include_dirs,
                    extra_compile_args=cflags,
                    extra_link_args=ldflags,
                    language='c++',
                    depends=common_headers),
          Extension(name="QPtomographer.bistates",
                    sources=['cxx/pydnormbistates.cxx'],
                    library_dirs=[],
                    libraries=[],
                    include_dirs=include_dirs,
                    extra_compile_args=cflags,
                    extra_link_args=ldflags,
                    language='c++',
                    depends=common_headers),
          Extension(name="QPtomographer.figofmerit",
                    sources=['cxx/pyfigofmerit.cxx'],
                    library_dirs=[],
                    libraries=[],
                    include_dirs=include_dirs,
                    extra_compile_args=cflags,
                    extra_link_args=ldflags,
                    language='c++',
                    depends=common_headers),
      ],

    # runtime requirements: qutip
    install_requires=[
        'qutip',
    ],

    # Not safe for keeping in a ZIP file (because of our extension)
    zip_safe=False,
)
