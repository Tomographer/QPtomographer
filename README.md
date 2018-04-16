# Reliable Diamond Norm Estimation for Quantum Process Tomography

Derive quantum error bars for quantum processes in terms of the diamond norm to
a reference quantum channel.

Theory: see arXiv:XXXX.XXXXX


## Installation

### Build SCS

Download (or clone) [SCS ≥ 2.0](https://github.com/cvxgrp/scs), say to
`$HOME/Downloads/scs`, and compile it as follows:

First, edit the file `scs.mk` as follows:

    > make CTRLC=0 USE_OPENMP=0 USE_LAPACK=1

On some systems such as on *Mac OS X*, you might have to adjust or specify the
necessary flags for linking to `BLAS/LAPACK`.  Use the variable
`BLASLDFLAGS="..."` for that.  For instance, on *Mac OS X*, you should probably
use:

    > make CTRLC=0 USE_OPENMP=0 USE_LAPACK=1 BLASLDFLAGS="-framework Accelerate"
    

## Install `tomographer` and its prerequisites

Make sure you have already installed the `tomographer` package, [as described
here][tomographer_py_inst].

You need `tomographer` version ≥ 5.4.

[tomographer_py_inst]: https://tomographer.github.io/tomographer/get-started/#python-version


## Building and installing `dnormtomo`

You should only have to run

    > SCS_ROOT=$HOME/Downloads/scs python setup.py install

specifying the path where you compiled SCS using the environment variable
`SCS_ROOT`.

Or, to install as administrator,

    > sudo -H SCS_ROOT=$HOME/Downloads/scs python setup.py install

The good news is that `dnormtomo`'s setup script automatically picks up all the
C++ flags set for `tomographer` itself, and uses those same flags. Thus, if
`tomographer` compiled, `dnormtomo` should compile as well (just make sure you
use the same compiler).

### NOTE: Installing using `pip`
You can also install `dnormtomo` using `pip`, so that it is properly seen as a
package and can be uninstalled easily.  For that, build a source package and
install it by running the following commands:

    > SCS_ROOT=$HOME/Downloads/scs python setup.py sdist  # creates dist/dnormtomo-1.0.tar.gz
    > SCS_ROOT=$HOME/Downloads/scs pip install dist/dnormtomo-1.0.tar.gz  # might need sudo -H as above

Then `dnormtomo` is seen as an installed wheel package, which has a certain
number of advantages.  For instance, you can uninstall easily with `pip
uninstall dnormtomo`.


# License

`dnormtomo` is released under the MIT License (see `LICENSE.txt`).
