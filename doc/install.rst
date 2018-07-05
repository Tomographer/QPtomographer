
Installation instructions
-------------------------

The `QPtomographer` package is uses the standard `setuptools` Python
infrastructure, providing a `setup.py` script like most other Python packages.

There are some little things to set up first however.

Compiling SCS
~~~~~~~~~~~~~

First, you need to compile `SCS ≥ 2.0 <https://github.com/cvxgrp/scs>`_. We make
use of this great library to calculate the diamond norm distance between two
channels.  (Unfortunately, it is not enough to install the `scs` Python package,
because we need SCS' C interface.)

Make sure you have downloaded SCS version 2.0.0 or later.

Download (or clone) `SCS ≥ 2.0 <https://github.com/cvxgrp/scs>`_, say to
``$HOME/Downloads/scs``, and compile it with specific options set.  This can be
done with the following steps using the command-line.

First, clone SCS into your `Downloads` directory, and enter that directory.
(You may choose a different directory; the following instructions however assume
you clone SCS in the `Downloads` directory.)

.. code-block:: bash

    > cd Downloads
    > git clone https://github.com/cvxgrp/scs.git
    > cd scs

If you are on *Mac OS X*, compile SCS with the following command::

    > make CTRLC=0 USE_OPENMP=0 USE_LAPACK=1 BLASLDFLAGS="-framework Accelerate"

Otherwise (e.g. on *Linux*), compile SCS with the following command::

    > make CTRLC=0 USE_OPENMP=0 USE_LAPACK=1



*Notes:*

 - ``CTRLC=0`` is needed because we will catch *Ctrl+C* keystrokes ourselves,
   and SCS's keystrokes captures interferes with our own mechanicsm.
   
 - We neeed ``USE_OPENMP=0`` because our invokation of *SCS* is already within
   parallel tasks, so we want to avoid double-parallelism which won't speed up
   anything.
   
 - Lapack is needed for solving SDPs, so it must be enabled with
   ``USE_LAPACK=1``.
   
   You might have to adjust or specify the necessary flags for linking to
   `BLAS/LAPACK` (the variable ``BLASLDFLAGS=...``).  On *Mac OS X*, use
   ``"-framework Accelerate"``; on Ubuntu, install for example `openblas` and
   use ``"-llapack -lblas"``.  The defaults might be satisfactory.  Try and see
   what works.




Install `tomographer` and its prerequisites
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Make sure you have already installed the `tomographer` package (version ≥ 5.4),
`as described here
<https://tomographer.github.io/tomographer/get-started/#python-version>`_.



Building and installing `QPtomographer`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Simply run::

  > SCS_ROOT=$HOME/Downloads/scs python setup.py install

specifying the path where you compiled SCS using the environment variable ``SCS_ROOT``.

Or, to install as administrator::

  > sudo -H SCS_ROOT=$HOME/Downloads/scs python setup.py install

.. important:: You should use the same compiler as the one you used to compile
               the `tomographer` package.

               If you leave the default this shouldn't be a problem.
