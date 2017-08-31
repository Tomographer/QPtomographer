
Installation instructions
-------------------------

The `dnormtomo` package is uses the standard `setuptools` Python infrastructure,
providing a `setup.py` script like most other Python packages.

There are some little things to set up first however.

Compiling SCS
~~~~~~~~~~~~~

First, you need `SCS <https://github.com/cvxgrp/scs>`_. We make use of this
great library to calculate the diamond norm distance between two channels.
(Unfortunately, it is not enough to install the `scs` Python package, because we
need SCS' C interface.)

Download (or clone) `SCS <https://github.com/cvxgrp/scs>`_, say to
``$HOME/Downloads/scs``, and compile it according to the following instructions.

First, edit the file ``scs.mk`` (in the SCS sources) as follows:

 - Search for a line starting with ``CTRLC =``, and make sure it reads ``CTRLC
   = 0`` (we will catch *Ctrl+C* keystrokes ourselves)
   
 - Search for a line starting with ``USE_OPENMP =``, and make sure that it
   reads ``USE_OPENMP = 0`` (our invokation of *SCS* is already within parallel
   tasks, so avoid double-parallelism which won't speed up anything)
   
 - Search for a line starting with ``USE_LAPACK =``, and make sure that it reads
   ``USE_LAPACK = 1`` (needed for solving SDPs).
   
   You might have to adjust or specify the necessary flags for linking to
   `BLAS/LAPACK` (the variable `BLASLDFLAGS`). [On Mac OS X, use ``-framework
   Accelerate``; on Ubuntu, install for example `openblas` and use ``-llapack
   -lblas``.]
   
Finally, compile `scs` by running::

    > make


Install `tomographer` and its prerequisites
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Make sure you have already installed the `tomographer` package, `as described
here <https://tomographer.github.io/tomographer/get-started/#python-version>`_.



Building and installing `dnormtomo`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Simply run::

  > SCS_ROOT=$HOME/Downloads/scs python setup.py install

specifying the path where you compiled SCS using the environment variable ``SCS_ROOT``.

Or, to install as administrator::

  > sudo -H SCS_ROOT=$HOME/Downloads/scs python setup.py install

.. important:: You should use the same compiler as the one you used to compile
               the `tomographer` package.

               If you leave the default this shouldn't be a problem.
