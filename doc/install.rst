
Installation instructions
-------------------------

The `QPtomographer` package uses the standard `setuptools` Python
infrastructure, providing a `setup.py` script like most other Python packages.

It may also be installed, as most other Python packages, using `pip`.

There are some little things to set up first however.


Prerequisite: Install BLAS/Lapack
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You need a working installation of BLAS/Lapack.  This is provided by default on
`Mac OS X`, in which case you don't have to do anything.  On Ubuntu, run the
command::

  sudo apt-get install liblapack-dev


Prerequisite: Compile SCS
~~~~~~~~~~~~~~~~~~~~~~~~~

You need to compile `SCS ≥ 2.0 <https://github.com/cvxgrp/scs>`_. We make use of
this great library to calculate the diamond norm distance between two channels.
(Unfortunately, it is not enough to install the `scs` Python package, because we
need SCS's C interface.)

Make sure you have downloaded SCS version 2.0.0 or later.

Download and unpack `SCS <https://github.com/cvxgrp/scs/releases>`_, say in
``$HOME/Downloads/``, and compile it with specific options set.  This can be
done as follows.

First, download `the latest version of SCS
<https://github.com/cvxgrp/scs/releases>`_ into your `Downloads` directory,
unpack that archive, and enter that directory.  (You may work in a different
directory; the following instructions assume you unpacked SCS in the `Downloads`
directory.)

.. code-block:: bash

    > cd Downloads
    > tar xvfz scs-2.0.2.tar.gz
    > cd scs-2.0.2

If you are on *Mac OS X*, compile SCS with the following command::

    > make CTRLC=0 USE_OPENMP=0 USE_LAPACK=1 BLASLDFLAGS="-framework Accelerate"

Otherwise (e.g. on *Linux*), compile SCS with the following command::

    > make CTRLC=0 USE_OPENMP=0 USE_LAPACK=1



*Notes:*

 - ``CTRLC=0`` is needed because we will catch *Ctrl+C* keystrokes ourselves,
   and we don't want SCS to interfere with this.
   
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




Prerequisite: Install `tomographer` and other Python packages
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Make sure you have already installed the `tomographer` package (version ≥ 5.4)
along with its prerequisites, `as described here
<https://tomographer.github.io/tomographer/get-started/#python-version>`_.
Basically, you should run the following commands depending on whether you are
using *anaconda* or *pip*::

  # If you are using *pip* (or if you're not sure):
  pip install numpy pybind11
  pip install tomographer

  # If you are using *conda*:
  conda install numpy gcc libgcc
  conda install -c conda-forge pybind11
  pip install tomographer

All ``pip install`` commands (e.g. ``pip install SOMEPACKAGE``) might have to be
prefixed by ``sudo -H`` (e.g. ``sudo -H pip install SOMEPACKAGE``) if you need
administrator priviledges, or you can use the option ``--user`` to install into
a user directory (e.g. ``pip install SOMEPACKAGE --user``).


The `QPtomographer` package also requires `qutip`.  Install it before
proceeding::

  # if you are using *pip* (or if you're not sure)
  pip install qutip

  # if you are using *conda*
  conda install -c conda-forge qutip



Download and install `QPtomographer`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Simplest solution using pip:**

The simplest solution to install `QPtomographer` is using `pip`.  However, you
need to set the environment variable ``SCS_ROOT`` to the root directory where
you downloaded and compiled `SCS`::

  # Works with both *conda* and *pip* users:
  SCS_ROOT=$HOME/Downloads/scs-2.0.2 pip install QPtomographer


**Alternative: Building and installing QPtomographer from source:**

You can also install `QPtomographer` from source.  First, obtain the most recent
version of `QPtomographer` from the `releases page
<https://github.com/Tomographer/QPtomographer/releases>`_.

Then, after unpacking the archive, run the python setup script::

  > SCS_ROOT=$HOME/Downloads/scs-2.0.2 python setup.py install

specifying the path where you compiled SCS using the environment variable ``SCS_ROOT``.

Or, to install as administrator::

  > sudo -H SCS_ROOT=$HOME/Downloads/scs-2.0.2 python setup.py install

 
.. important:: You should use the same compiler as the one you used to compile
               the `tomographer` package.

               If you leave the default this shouldn't be a problem.
