
.. _conventions-channels-states:

Conventions for specifying channels and states
==============================================

Channels and states are specified as NumPy :py:class:`array <numpy.ndarray>`
objects, or as QuTip :py:class:`~qutip.Qobj` objects.  A state on a single
system :math:`X` is represented as a complex NumPy array of dimensions
:math:`d_X \times d_X`.  A state on several systems :math:`X, Y, Z` is specified
as a NumPy array of dimension :math:`d_X d_Y d_Z \times d_X d_Y d_Z`.  For
instance, the state :math:`|0\rangle\langle0|_X \otimes |{+}\rangle\langle{+}|_Y
\otimes |1\rangle\langle1|_Z` be formed for instance as follows::

  import numpy as np
  rho_XYZ = np.kron(np.kron(np.array([[1,0],[0,0]]), np.array([[.5,.5],[.5,.5]])), np.array([[0,0],[0,1]]))

  import qutip
  rho2_XYZ = qutip.Qobj(
      np.kron(np.kron(np.array([[1,0],[0,0]]), np.array([[.5,.5],[.5,.5]])), np.array([[0,0],[0,1]])),
      dims=[[2,2,2],[2,2,2]]
      )


Channels are specified by their non-normalized Choi matrix.  That is, let
:math:`\{ |k\rangle_X \}, \{ |\ell\rangle_Y \}` denote the standard bases of the
input and the output systems, let :math:`R\simeq X` be a copy of :math:`X`,
and define

.. math::

   | \tilde{\Phi} \rangle_{XR} = \sum |k\rangle_X \otimes |k\rangle_{R}\ .


Then the Channel :math:`\mathcal{E}_{X\to Y}` is represented by its unnormalized
Choi matrix

.. math::

   E_{R Y} = (\mathcal{E}_{X\to Y}\otimes\operatorname{id}_{R})(
     | \tilde{\Phi} \rangle \langle \tilde{\Phi} |_{XR}
   )\ .

Using this convention, the application of a channel :math:`\mathcal{E}_{X\to Y}`
onto a state :math:`\sigma` is computed as

.. math::

   \rho_Y = \operatorname{tr}_{R}[ E_{R Y} \, \sigma_X^T ]\ ,

where :math:`\sigma_X^T` is the partial transpose of :math:`\sigma_X` during
which the system :math:`X` is relabled as :math:`R`.

The unnormalized Choi matrix can be represented as a NumPy :py:class:`array
<numpy.ndarray>` or as a :py:class:`qutip.Qobj` in the same way as states.
