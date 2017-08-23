

Channel-space sampling method: `dnormtomo.channelspace`
=======================================================

.. py:module:: dnormtomo.channelspace

This module implements the channel-space sampling method for the reliable
diamond norm estimation.

.. py:data:: dnormtomo.channelspace.RandHermExp

   Numerical constant which signifies to carry out the random walk using the
   ":math:`e^{iH}` jumps mode" (see paper).  This value can be specified to the
   `channel_walker_jump_mode=` argument of :py:func:`~dnormtomo.channelspace.run()`.

.. py:data:: dnormtomo.channelspace.ElemRotations

   Numerical constant which signifies to carry out the random walk using the "elementary
   rotations jumps mode" (see paper).  This value can be specified to the
   `channel_walker_jump_mode=` argument of :py:func:`~dnormtomo.channelspace.run()`.


.. py:function:: dnormtomo.channelspace.run(...)

   The main execution routine for the bipartite state sampling method.

   Documentation here ...............................................

   All arguments should be specified as keyword arguments. Accepted arguments are:

      - `dimX`: dimension of input system

      - `dimY`: dimension of output system

      - `Emn`: a list of POVM effects on :math:`X\otimes Y`.  Should be a python list
        where each item is a POVM effect provided as a complex :math:`d_X d_Y\times d_X
        d_Y` matrix, of type `NumPy` array.

      - `Nm`: the corresponding list of frequency counts. This is a list of integers, of
        the same length as `Emn`, where `Nm[k]` is the number of times the outcome
        described by the POVM effect `Emn[k]` was observed.

      - `ref_channel_XY`: the reference channel to which to calculate the diamond norm
        distance to. Specify the channel by its Choi matrix on :math:`X\otimes Y`, as a
        :math:`d_X d_Y \times d_X d_Y` `NumPy` array. The normalization should be such
        that the reduced state on :math:`X` is the identity operator, ie., the trace of
        the full matrix should be equal to :math:`d_X`.

        If `dimX==dimY`, you may set to `ref_channel_XY=None`, in which case the identity
        process (with respect to the canonical basis) is used as reference channel.

      - `hist_params`: the parameters of the histogram to collect, i.e., range and number
        of bins.  Must be a :py:class:`tomographer.HistogramParams` object.

      - `mhrw_params`: the parameters of the random walk, i.e., the step size, the sweep
        size, number of thermalization sweeps and number of live run sweeps.  Specify as a
        :py:class:`tomographer.MHRWParams` object.

      - `channel_walker_jump_mode`: one of :py:data:`~dnormtomo.channelspace.RandHermExp`
        or :py:data:`~dnormtomo.channelspace.ElemRotations`, depending on the requested
        method of random walk step.

      - `dnorm_epsilon`: the precision at which to calculate the diamond norm (which is
        calculated by numerically solving the corresponding semidefinite program using
        `SCS <https://github.com/cvxgrp/scs>`_). The default is `1e-3`.

      - `num_repeats`: the total number of random walks to run.

      - `binning_num_levels`: number of levels in the binning analysis. By default, or if
        the value `-1` is specified, an appropriate number of levels is determined
        automatically.

      - `progress_fn`, `progress_interval_ms`, `ctrl_step_size_params`,
        `ctrl_converged_params`: these parameters are treated the same as for
        :py:func:`tomographer.tomorun.tomorun()`.



.. py:exception:: dnormtomo.channelspace.DNormChannelSpaceInvalidInputError

   Exception is thrown whenever invalid input to the
   :py:func:`~dnormtomo.channelspace.run()` function is encountered.


