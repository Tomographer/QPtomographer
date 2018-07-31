# Reliable Quantum Process Tomography

![QPtomographer](QPtomographer.svg)

Derive quantum error bars for quantum processes in terms of the diamond norm to
a reference quantum channel.

Theory: see arXiv:XXXX.XXXXX


## Installation

See detailed installation instructions at
https://QPtomographer.readthedocs.org/latest/install/ .

Basically:

- you need to download & compile [SCS ≥ 2.0](https://github.com/cvxgrp/scs);

- you need to have the `tomographer` package
  installed [as described here][tomographer_py_inst];
  
- you can install `QPtomographer` via `pip`.


Here are quick-and-easy steps for some common setups:

### Ubuntu/Linux

Download & compile [SCS ≥ 2.0](https://github.com/cvxgrp/scs):

    > cd $HOME/Downloads
    > curl -L https://github.com/cvxgrp/scs/archive/v2.0.2.tar.gz -o scs-2.0.2.tar.gz
    > tar xvfz scs-2.0.2.tar.gz && cd scs-2.0.2
    > make CTRLC=0 USE_OPENMP=0 USE_LAPACK=1

Download & compile `tomographer ≥ 5.4` and other prerequisites:

    # If you are using *pip* (or if you're not sure):
    > pip install numpy pybind11 qutip
    > pip install tomographer

    # If you are using *anacoda/conda*:
    > conda install numpy gcc libgcc
    > conda install -c conda-forge pybind11 qutip
    > pip install tomographer

Install `QPtomographer` via `pip`:

    > SCS_ROOT=$HOME/Downloads/scs-2.0.2 pip install QPtomographer


### Apple Mac OS X

Download & compile [SCS ≥ 2.0](https://github.com/cvxgrp/scs):

    > cd $HOME/Downloads
    > curl -L https://github.com/cvxgrp/scs/archive/v2.0.2.tar.gz -o scs-2.0.2.tar.gz
    > tar xvfz scs-2.0.2.tar.gz && cd scs-2.0.2
    > make CTRLC=0 USE_OPENMP=0 USE_LAPACK=1 BLASLDFLAGS="-framework Accelerate"

Download & compile `tomographer` and other prerequisites:

    # If you are using *pip* (or if you're not sure):
    > pip install numpy pybind11 qutip
    > pip install tomographer

    # If you are using *anacoda/conda*:
    > conda install numpy gcc libgcc
    > conda install -c conda-forge pybind11 qutip
    > pip install tomographer

Install `QPtomographer` via `pip`:

    > SCS_ROOT=$HOME/Downloads/scs-2.0.2 pip install QPtomographer



## Alternatively, you can build and install `QPtomographer` from source

You should only have to run

    > SCS_ROOT=$HOME/Downloads/scs python setup.py install

specifying the path where you compiled SCS using the environment variable
`SCS_ROOT`.  Or, to install as administrator,

    > sudo -H SCS_ROOT=$HOME/Downloads/scs python setup.py install

The good news is that `QPtomographer`'s setup script automatically picks up all
the C++ flags set for `tomographer` itself, and uses those same flags. Thus, if
`tomographer` compiled, `QPtomographer` should compile as well (just make sure
you use the same compiler).



# License

`QPtomographer` is released under the MIT License (see `LICENSE.txt`).
