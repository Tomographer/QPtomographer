
Evaluating the figures of merit: `QPtomographer.figofmerit`
===========================================================

.. py:module:: QPtomographer.figofmerit

This module provides some functions to calculate the diamond norm distance, the
average entanglement fidelity as well as the worst-case entanglement fidelity
(see :ref:`figures-of-merit`).

.. py:function:: diamond_norm_dist(Delta_XY, dim_X, [epsilon=1e-6])

   Compute the diamond norm of `Delta_XY`.  This computes :math:`\frac{1}{2}
   \lVert \Delta_{X\to Y} \rVert_{\diamond}` by solving a semidefinite program,
   as described in :ref:`figures-of-merit`.  Specify the dimension of the input
   system :math:`X` in the argument `dim_X` (the dimension of the output system
   is deduced from the size of `Delta_XY`).

   The Hermiticity-preserving map :math:`\Delta_{X\to Y}(\cdot)` is specified in
   the argument `Delta_XY` as a complex NumPy array, via the Choi matrix, of
   dimension :math:`d_X d_Y \times d_X d_Y`, as specified in
   :ref:`conventions-channels-states`.

   The semidefinite problem is solved using `SCS` via QPtomographer's C++
   wrappers.  This function is in fact a wrapper to the same C++ subroutine code
   that computes the diamond norm when calling
   :py:func:`QPtomographer.channelspace.run()` or
   :py:func:`QPtomographer.bistates.run()`.

   The argument `epsilon` specifies the precision to which the SDP
   characterizing the diamond norm is solved.


.. py:function:: avg_entgl_fidelity2(E_XY, dim_X)

   Compute the average entanglement fidelity (squared) of the process specified
   by `E_XY` to the identity process.  See :ref:`figures-of-merit` for details.
   Specify the dimension of the input system :math:`X` in the argument `dim_X`
   (the dimension of the output system is deduced from the size of `Delta_XY`).

   The completely positive, trace-preserving map :math:`\mathcal{E}_{X\to
   Y}(\cdot)` is specified in the argument `E_XY` as a complex NumPy array, via
   the Choi matrix, of dimension :math:`d_X d_Y \times d_X d_Y`, as specified in
   :ref:`conventions-channels-states`.

   The average entanglement fidelity is calculated by direct computation, via
   the same C++ code that computes the average entanglement fidelity when
   calling :py:func:`QPtomographer.channelspace.run()` or
   :py:func:`QPtomographer.bistates.run()`.


.. py:function:: worst_entgl_fidelity2(E_XY, dim_X, [epsilon=1e-6])

   Compute the worst-case entanglement fidelity (squared) of the process
   specified by `E_XY` to the identity process.  See :ref:`figures-of-merit` for
   details.  Specify the dimension of the input system :math:`X` in the argument
   `dim_X` (the dimension of the output system is deduced from the size of
   `Delta_XY`).

   The completely positive, trace-preserving map :math:`\mathcal{E}_{X\to
   Y}(\cdot)` is specified in the argument `E_XY` as a complex NumPy array, via
   the Choi matrix, of dimension :math:`d_X d_Y \times d_X d_Y`, as specified in
   :ref:`conventions-channels-states`.

   The semidefinite problem is solved using `SCS` via QPtomographer's C++
   wrappers.  This function is in fact a wrapper to the same C++ subroutine code
   that computes the worst-case entanglement fidelity when calling
   :py:func:`QPtomographer.channelspace.run()` or
   :py:func:`QPtomographer.bistates.run()`.

   The argument `epsilon` specifies the precision to which the SDP
   characterizing the worst-case entanglement fidelity is solved.

