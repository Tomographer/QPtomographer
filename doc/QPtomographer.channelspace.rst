

Channel-space sampling method: `QPtomographer.channelspace`
===========================================================

.. py:module:: QPtomographer.channelspace

This module implements the channel-space sampling method for the reliable
diamond norm estimation.

.. py:function:: run(...)

   The main execution routine for the channel-space state sampling method.  This
   function behaves analogously to :py:func:`tomographer.tomorun.tomorun()`, but
   in the setting described in |paper_arxiv_ref| and calculating a :ref:`figure
   of merit for quantum processes <figures-of-merit>`.

   Note: The `Emn` matrices must be weighted by input state as follows. The
   log-likelihood function is calculated from the matrices `Emn` as:

   .. math::
      \ln\mathcal{L}(\Lambda_{A\to B}\mid E)
      = \sum_k \texttt{Nm[k]} \cdot \ln \mbox{tr}\bigl(\Lambda_{A\to B}(\Phi_{AP})\,
        \texttt{Emn[k]}\bigr)\ ,

   where :math:`\left|\Phi\right\rangle_{AP} = \sum_k \left|k\right\rangle_A
   \left|k\right\rangle_P`.  Hence, ``Emn[k]`` must be equal to
   :math:`\sigma_{P,j}^{1/2} E \sigma_{P,j}^{1/2}` where :math:`\sigma_{P,j}` is
   the corresponding transposed input state (see paper, and see the function
   argument `Emn=` below).

   All arguments should be specified as keyword arguments. Accepted arguments are:

      - `dimX`: dimension of input system (required)

      - `dimY`: dimension of output system (required)

      - `Emn`: a list of POVM effects *weighted by input state* on
        :math:`X\otimes Y`.  Should be a python list where each item is a matrix
        provided as a complex :math:`d_X d_Y\times d_X d_Y` `NumPy` array.

      - `Nm`: the corresponding list of frequency counts. This is a list of
        integers, of the same length as `Emn`, where `Nm[k]` is the number of
        times the outcome described by the POVM effect `Emn[k]` was observed.

      - `fig_of_merit`: which figure of merit to use. May be one of
        'diamond-distance' (the default), 'entanglement-fidelity',
        'worst-entanglement-fidelity', or a custom callable Python object.
        Refer to section :ref:`figures-of-merit`.

        The callable will be invoked with as argument an isometry of shape
        `(dY*dX*dY, dX)` given as a NumPy array, and where the output system is
        the "outer" system of the output of the isometry.  More precisely,
        `V[j*dX*dY+k*dY+l,i]` is the matrix element of the isometry

        .. math::
           V = \sum_{i,j,k,\ell} V_{(jk\ell),i} |{j}\rangle_Y |{k}\rangle_{\tilde{X}}
              |{\ell}\rangle_{\tilde{Y}} \langle{i}|_X

        which is to be interpreted as a Stinespring dilation of a process
        :math:`X\to Y` with an environment system
        :math:`\tilde{X}\otimes\tilde{Y}` (all indices start at zero).

        The Choi matrix of the process `E_RY` is recovered from `V` as follows::

          T_RY = V.T.reshape((dX*dY*dX*dY,)).reshape((dX*dY,dX*dY,), order='F').T
          E_RY = np.dot(T_RY, T_RY.conj().T)  # T * T'

        where `E_RB` now follows the conventions described in
        :ref:`conventions-channels-states`.  See also below for a more in-depth
        example.

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

      - `channel_walker_jump_mode`: one of
        :py:data:`~QPtomographer.channelspace.RandHermExp` or
        :py:data:`~QPtomographer.channelspace.ElemRotations`, depending on the
        requested method of random walk step.

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

.. rubric:: Custom figure of merit

Here is some sample code to understand how a custom figure of merit can be
implemented.  It computes the average entanglement fidelity manually (you can
compare this to the built-in version of the entanglement fidelity)::

     import numpy as np
     import numpy.linalg as npl

     import qutip # only needed for Method #1 below
     
     def custom_figofmerit(V):
         # compute average entanglement fidelity manually with custom callable
         # (assuming a qubit->qubit channel)
         #
         # This code is a purely illustrative example, since the average entanglement
         # fidelity is already implemented as a built-in figure of merit (and it is
         # written in C++, making it much faster).
         #
         # V is a numpy object of size (dY*dX*dY, dX) representing an isometry
         # X->YXY, where in the output, the first 'Y' is the actual output system of
         # the corresponding process of which this is a purification, and the
         # remaining 'XY' are the environment ancilla.  The actual process can be
         # recovered as follows.

         # Method #1: more explicit using qutip (but potentially slower)
         #
         #print("V=%r"%(V))
         Vq = qutip.Qobj(V, dims=[[2,2,2],[2]])
         # Phi_RA is the maximally entangled ket between R and A:  Phi_RA = (\sum_k |k>_R |k>_A)
         Phi_RA = qutip.Qobj(np.eye(2).reshape([4,1]), dims=[[2,2],[1]])
         # E_RB is the Choi matrix of the current process, E_RB = (id_R\otimes\mathcal{E}_{A\to B})(|\Phi><\Phi|_RA)
         E_RB = qutip.ptrace(
             qutip.tensor(qutip.qeye(2),Vq) * Phi_RA * Phi_RA.dag() * qutip.tensor(qutip.qeye(2),Vq.dag()),
             [0,1]
         )
         #print(E_RB)
         #
         # we can check that the reduced state on the input is indeed the identity operator:
         #print(E_RB.ptrace([0]))
         assert npl.norm(E_RB.ptrace([0]).data.toarray() - np.eye(2), 'nuc') <= 1e-6
         
         # Method #2: use reshape directly (faster, but more obscure)
         #
         # we have: V((jkl),i) |i>_A |jkl>_{BA'B'}
         # we want: T((ij),(kl)) from V((jkl),i)
         # (ij) = i*dJ+j
         # (jkl) = j*dK*dL+k*dL+l
         #
         # V.T is (V.T)(i,(jkl))
         # V.T.reshape((16,),) is (V...)((ijkl),) with (ijkl) = i*dJ*dK*dL+(jkl)
         # V.T.reshape((16,),).reshape((4,4,), order='F') is (V...)((kl),(ij))
         # V.T.reshape((16,),).reshape((4,4,), order='F').T is (V...)((ij),(kl))
         T_RB = qutip.Qobj(V.T.reshape((16,)).reshape((4,4,), order='F').T, dims=[[2,2],[2,2]])
         E2_RB = T_RB*T_RB.dag()
         #print(E2_RB)
     
         # Compare method #1 and method #2 and make sure they give the same thing
         assert npl.norm( (E_RB - E2_RB).data.toarray(), 'nuc') <= 1.e-6 # nuclear norm is Schatten-1 norm

         # Compute the entanglement fidelity using E_RB:
         #
         # <\Phi|_RA E_RB |\Phi>_RA / (2*2) :  divide by (2*2) because \Phi is not normalized
         val = qutip.expect(E_RB, Phi_RA) / (2*2)
         return val


.. py:data:: RandHermExp

   Numerical constant which signifies to carry out the random walk using the
   ":math:`e^{iH}` jumps mode" (see paper).  This value can be specified to the
   `channel_walker_jump_mode=` argument of
   :py:func:`~QPtomographer.channelspace.run()`.

.. py:data:: ElemRotations

   Numerical constant which signifies to carry out the random walk using the
   "elementary rotations jumps mode" (see paper).  This value can be specified
   to the `channel_walker_jump_mode=` argument of
   :py:func:`~QPtomographer.channelspace.run()`.


.. py:exception:: DNormChannelSpaceInvalidInputError

   Exception is thrown whenever invalid input to the
   :py:func:`~QPtomographer.channelspace.run()` function is encountered.


