# Reliable Diamond Norm Estimation for Quantum Process Tomography

Derive quantum error bars for quantum processes in terms of the diamond norm to
a reference quantum channel.

Theory: see arXiv:XXXX.XXXXX


## Installation

### Build SCS

Download (or clone) [SCS](https://github.com/cvxgrp/scs), say to
`$HOME/Downloads/scs`, and compile it according to the following instructions.

First, edit the file `scs.mk` as follows:

 - Search for a line starting with `CTRLC =`, and make sure it reads `CTRLC =
   0` (we will catch Ctrl+C keystrokes ourselves)
   
 - Search for a line starting with `USE_OPENMP =`, and make sure that it reads
   `USE_OPENMP = 0` (our invokation of SCS is already within parallel tasks, so
   avoid double-parallelism which won't speed up anything)
   
 - Search for a line starting with `USE_LAPACK =`, and make sure that it reads
   `USE_LAPACK = 1` (needed for solving SDPs).
   
   You might have to adjust or specify the necessary flags for linking to
   `BLAS/LAPACK` (the variable `BLASLDFLAGS`). [On Mac OS X, use `-framework
   Accelerate`; on Ubuntu, install for example `openblas` and use `-llapack
   -lblas`.]
   
Finally, compile `scs` by running

    > make


## Install `tomographer` and its prerequisites

Make sure you have already installed the `tomographer`
package, [as described here][tomographer_py_inst].

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


# License

XXXXX

