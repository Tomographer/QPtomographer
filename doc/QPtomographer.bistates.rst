
Bipartite sampling method: `QPtomographer.bistates`
===================================================

.. py:module:: QPtomographer.bistates

This module implements the bipartite state sampling method for the reliable
diamond norm estimation.


.. py:function:: run(...)

   The main execution routine for the bipartite state sampling method.  This
   function behaves analogously to :py:func:`tomographer.tomorun.tomorun()`, but
   in the setting described in |paper_arxiv_ref|_ and calculating the diamond
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

        The callable will be invoked with as argument a square complex matrix
        :math:`T` of shape `(d_X*dY, d_X*d_Y)` given as a NumPy array, and where
        the process matrix, a normalized bipartite state, can be simply computed
        as :math:`T T^\dagger`.  If you conjugate by the inverse square root of
        the reduced state on the first system, you get the Choi matrix of the
        process.  The following code computes this::

          def custfigofmerit(T):
              # compute average entanglement fidelity manually with custom
              # callable, assuming a qubit->qubit quantum channel
              #
              # T is a numpy object of size (dX*dY, dX*dY) representing a
              # square root of a bipartite state (the current point of the
              # random walk).
              #print(T)
          
              # To compute the Choi matrix of the process encoded in T:
              #
              # Compute the bipartite state (Choi matrix weighted by input state)
              rho_RB = np.dot(T, T.conj().T)
              rho_q_RB = qutip.Qobj(rho_RB, dims=[[2,2],[2,2]])
              # conjugate by inverse square root of input state to get Choi matrix
              invsqrtmrho_R = qutip.tensor(invsqrtm(rho_q_RB.ptrace([0])), qutip.qeye(2))
              E_RB = invsqrtmrho_R * rho_q_RB * invsqrtmrho_R
          
              #print(E_RB)
              #print(E_RB.ptrace([0]))
              assert npl.norm( E_RB.ptrace([0]).data.toarray() - np.eye(2), 'nuc' ) <= 1e-6
          
              # Maximally entangled ket between R and A:  Phi_RA = (\sum_k |k>_R |k>_A)
              Phi_RA = qutip.Qobj(np.eye(2).reshape([4,1]), dims=[[2,2],[1]])
          
              # <\Phi|_RA E_RB |\Phi>_RA / (2*2) :  divide by (2*2) because \Phi is not normalized
              val = qutip.expect(E_RB, Phi_RA) / (2*2)
              return val

          # helper function: compute matrix inverse square root of a Hermitian `qutip.Qobj` matrix
          def invsqrtm(Aq):
              d, V = npl.eigh(Aq.data.toarray())
              return qutip.Qobj(np.dot(np.dot(V, np.diag(np.reciprocal(np.sqrt(d)))), V.conj().T),
                                dims=Aq.dims)
          

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

      - `jumps_method`: one of ``"full"`` or ``"light"``, depending on the
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


