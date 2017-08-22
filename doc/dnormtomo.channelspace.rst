

Channel-space sampling method: `dnormtomo.channelspace`
=======================================================

.. py:module:: dnormtomo.channelspace

This module implements the channel-space sampling method for the reliable
diamond norm estimation.

.. py:function:: dnormtomo.channelspace.run(...)

    The main execution routine for the bipartite state sampling method.

    Documentation here ...............................................

    Arguments:

      - `dimX`: dimension of input system

      - `dimY`: dimension of output system

      - `Emn`: a list of POVM effects on :math:`X\otimes Y`.  Should be a python list
        where each item is a POVM effect provided as a complex :math:`dX dY\times dX dY`
        matrix, of type `NumPy` array.

      - `Nm`: the corresponding list of frequency counts. This is a list of integers, of
        the same length as `Emn`, where `Nm[k]` is the number of times the outcome
        described by the POVM effect `Emn[k]` was observed.

      - `hist_params`: the parameters of the histogram to collect, i.e., range and number
        of bins.  Must be a :py:class:`tomographer.HistogramParams` object.

      - `mhrw_params`: the parameters of the random walk, i.e., the step size, the sweep
        size, number of thermalization sweeps and number of live run sweeps.  Specify as a
        :py:class:`tomographer.MHRWParams` object.

      - `binning_num_levels`: ......

      py::arg("binning_num_levels") = -1,
      py::arg("num_repeats") = std::thread::hardware_concurrency(),
      py::arg("ref_channel_XY") = py::none(),
      py::arg("progress_fn") = py::none(),
      py::arg("progress_interval_ms") = (int)500,
      py::arg("channel_walker_jump_mode") = (int)RandHermExp,
      py::arg("dnorm_epsilon") = (double)1e-3,
      py::arg("ctrl_step_size_params") = py::dict(),
      py::arg("ctrl_converged_params") = py::dict(),



.. py:exception:: dnormtomo.channelspace.DNormChannelSpaceInvalidInputError

   Exception is thrown whenever invalid input to the
   :py:func:`~dnormtomo.channelspace.run()` function is encountered.


