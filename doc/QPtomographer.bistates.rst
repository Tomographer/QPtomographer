
Bipartite sampling method: `QPtomographer.bistates`
===================================================

.. py:module:: QPtomographer.bistates

This module implements the bipartite state sampling method for the reliable
diamond norm estimation.


.. py:function:: run(...)

   The main execution routine for the bipartite state sampling method.  This
   function behaves analogously to :py:func:`tomographer.tomorun.tomorun()`, but
   in the setting described in |paper_arxiv_ref| and calculating the diamond
   norm to a reference channel as figure of merit.

   As described in the paper, the log-likelihood function is calculated from the
   matrices `Emn` as follows:

   .. math::
      \ln\mathcal{L}(\rho_{AB}\mid E)
      = \sum_k \texttt{Nm[k]} \cdot \ln \mbox{tr}\bigl(\rho_{AB}\, \texttt{Emn[k]}\bigr)\ ,

   Hence, ``Emn[k]`` must be the POVM effect which was observed on the joint
   :math:`A\otimes B` systems (see paper, and see the function argument `Emn=`
   below).

   All arguments should be specified as keyword arguments. Accepted arguments
   are:

      - `dimX`: dimension of input system (required)

      - `dimY`: dimension of output system (required)

      - `Emn`: a list of POVM effects on :math:`X\otimes Y`.  Should be a python
        list where each item is a matrix provided as a complex :math:`d_X
        d_Y\times d_X d_Y` `NumPy` array.

      - `Nm`: the corresponding list of frequency counts. This is a list of
        integers, of the same length as `Emn`, where `Nm[k]` is the number of
        times the outcome described by the POVM effect `Emn[k]` was observed.

      - `fig_of_merit`: which figure of merit to use. May be one of
        'diamond-distance' (the default), 'entanglement-fidelity',
        'worst-entanglement-fidelity', or a custom callable Python object.
        Refer to section :ref:`figures-of-merit`.

        The callable will be invoked with as argument a square complex matrix of
        shape `(d_X*dY, d_X*d_Y)` given as a NumPy array, and where the process
        matrix, a normalized state, can be simply computed as::

          rho_YR = np.dot(T, T.conj().T)   # T * T'

      - `ref_channel_XY`: the reference channel to which to calculate the
        diamond norm distance to.  This argument is not used if `fig_of_merit`
        is not 'diamond-distance'.  Specify the channel by its Choi matrix on
        :math:`X\otimes Y`, as a :math:`d_X d_Y \times d_X d_Y` `NumPy`
        array. The normalization should be such that the reduced state on
        :math:`X` is the identity operator, ie., the trace of the full matrix
        should be equal to :math:`d_X`.

        If `dimX==dimY`, you may set `ref_channel_XY=None`, in which case the
        identity process (with respect to the canonical basis) is used as
        reference channel.

      - `hist_params`: the parameters of the histogram to collect, i.e., range
        and number of bins.  Must be a :py:class:`tomographer.HistogramParams`
        object.

      - `mhrw_params`: the parameters of the random walk, i.e., the step size,
        the sweep size, number of thermalization sweeps and number of live run
        sweeps.  Specify as a :py:class:`tomographer.MHRWParams` object.

      - `jump_mode`: one of ``"full"`` or ``"light"``, depending on the
        requested method of random walk step.  This argument has the same effect
        as the `jumps_method=` argument of
        :py:func:`tomographer.tomorun.tomorun()`.

      - `dnorm_epsilon`: the precision at which to calculate the diamond norm
        (which is calculated by numerically solving the corresponding
        semidefinite program using `SCS <https://github.com/cvxgrp/scs>`_). The
        default is `1e-3`.

      - `num_repeats`: the total number of random walks to run.  By default,
        this is set to the number of available cores.

      - `binning_num_levels`: number of levels in the binning analysis. By
        default, or if the value `-1` is specified, an appropriate number of
        levels is determined automatically.

      - `progress_fn`, `progress_interval_ms`, `ctrl_step_size_params`,
        `ctrl_converged_params`: these parameters are treated the same as for
        :py:func:`tomographer.tomorun.tomorun()`.



.. py:exception:: DNormBiStatesInvalidInputError

   Exception is thrown whenever invalid input to the
   :py:func:`~QPtomographer.bistates.run()` function is encountered.


