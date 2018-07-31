
.. _figures-of-merit:

Figures of Merit for Channels
-----------------------------

The current package has built-in capability for computing the following figures
of merit.  The choice of the figure of merit can be specified using the
`fig_of_merit=` parameter to :py:func:`QPtomographer.bistates.run()` or
:py:func:`QPtomographer.channelspace.run()`.


The diamond norm distance
~~~~~~~~~~~~~~~~~~~~~~~~~

The diamond norm distance between two channels :math:`\mathcal{E}_{X\to Y}` and
:math:`\mathcal{F}_{X\to Y}` is defined as

.. math::
   \frac{1}{2} \left\lVert \mathcal{E} - \mathcal{F} \right\rVert_{\diamond}
   = \sup_{\sigma_{XR}} \frac{1}{2} \left\lVert \mathcal{E}(\sigma_{XR})
   - \mathcal{F}(\sigma_{XR}) \right\lVert_{1}\ ,

where the supremum is taken over all bipartite states on :math:`X` and a
reference system :math:`R\simeq X`.  In the right hand side, the channels act on
the system :math:`X` and act as the identity channel on :math:`R`.  The
optimization may be restricted to pure states without loss of generality.

The diamond norm distance can be formulated as the following semidefinite
program, in terms of a positive semidefinite variable :math:`Z_{YR} \geqslant 0`
and a real variable :math:`\alpha`:[#WatrousDiamondSDP]_

.. math::
   \begin{array}{rc}
   \frac{1}{2} \left\lVert \Delta_{X\to Y} \right\rVert_{\diamond} \quad
   = \quad \mbox{minimize:} \quad
   &   \alpha \qquad\ , \\
   \mbox{subject to:}\quad
   & \mbox{tr}_Y(Z_{YR}) \leqslant \alpha\mathbb{I}_{R} \\
   & Z_{YR} \geqslant \Delta_{X\to Y}(\tilde{\Phi}_{X:R}) \end{array}

where we have defined

.. math::
   {|\tilde\Phi\rangle}_{X:R} = \sum_k {|k\rangle}_X\otimes{|k\rangle}_R\ .


The entanglement fidelity for channels
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An alternative quantity which is relevant to quantify the distance of a channel
:math:`\mathcal{E}_{X\to X'}` (where :math:`X'\simeq X`) to the identity channel
:math:`\mbox{id}_{X\to X'}` is the so-called *entanglement fidelity*

.. math::
    F_{e}(\mathcal{E}_{X\to X'}) = F^2(\mathcal{E}_{X\to X'}(\varphi_{XR}), \varphi_{XR})
     = {\langle\varphi|}_{XR} \mathcal{E}_{X\to X'}({|\varphi\rangle}{\langle\varphi|}_{XR}) {|\varphi\rangle}_{XR}\ ,

where we have defined the maximally entangled state :math:`|\varphi\rangle_{XR}`
(using a reference system :math:`R\simeq X`) and the (square) fidelity between
states :math:`F^2(\rho,\sigma)`, as:

.. math::
   |\varphi\rangle_{XR} = \frac{1}{\sqrt{\mbox{dim}(X)}} \,
      \sum_k {|k\rangle}_X\otimes{|k\rangle}_{R}\ ;

.. math::
   F^2(\rho,\sigma) = \left\lVert \sqrt{\rho} \sqrt{\sigma} \right\rVert_{1}^{2} \ .


The worst-case entanglement fidelity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Finally, one may optimize the entanglement fidelity over all possible input
states and seek out the worst case, instead of computing the entanglement
fidelity only for the maximally entangled state.  Hence, the *worst-case
entanglement fidelity* of a channel :math:`\mathcal{E}_{X\to X'}` (where
:math:`X'\simeq X`) is defined as

.. math::
    F_{\mbox{worst}}(\mathcal{E}_{X\to X'}) =
     \inf_{\sigma_{XR}} F^2(\mathcal{E}_{X\to X'}(\sigma_{XR}), \sigma_{XR})\ ,

where the optimization ranges over all bipartite quantum states defined over the
input and a reference system :math:`R\simeq X`, and where :math:`F(\rho,\sigma)`
is defined as above.

The optimization may be restricted to pure states without loss of generality,
because the minimum of a concave function is over a convex set is reached at the
boundary of the set.

The worst-case entanglement fidelity may be computed by evaluating the following
semidefinite program (as described in |paper_arxiv_ref|), in terms of the real
variable :math:`\mu` and positive semidefinite variable :math:`\rho_X \geqslant
0`:

.. math::
   \begin{array}{rcl}
   F_{\mbox{worst}}(\mathcal{E}_{X\to X'})\quad =
   \quad\mbox{minimize:}\quad
   & \mu \qquad            &\ , \\
   \mbox{subject to:}\quad
   & \mbox{tr}(\rho_X) = 1 &\\
   & 
     \left[\begin{array}{c|c}
      \mathbb{I}\vphantom{\Bigg[]} & M_{X'R}^\dagger \mbox{vec}(\rho_X) \\ \hline
      \mbox{vec}(\rho_X)^\dagger M_{X'R} & \mu
     \end{array}\right] \geqslant 0 &
   \end{array}

where :math:`M_{X'R}` is a factorization of the Choi matrix of the process, satisfying:

.. math::
   M_{X'R} \, M_{X'R}^\dagger = \mathcal{E}_{X\to X'}(\tilde\Phi_{XR})\ ,

with :math:`{|\tilde\Phi\rangle}_{XR}` as defined above.  The factorization can be
obtained using a Cholseky or LDLT factorization, for instance; or more generally
by computing any matrix square root.  The unitary freedom of the matrix square
root is irrelevant here.

In the above, the operation :math:`\mbox{vec}(\rho)` "vectorizes" the given
operator by stacking its coefficients into one large column; it is more
precisely defined as

.. math::
   \mbox{vec}(\rho) = (\mathbb{I}\otimes\rho) {|\tilde{\Phi}\rangle}_{XR}\ .



References
~~~~~~~~~~

.. [#WatrousDiamondSDP] Watrous, J. (2009). Semidefinite Programs for Completely
                        Bounded Norms. `Theory of Computing, 5(11), 217–238
                        <https://doi.org/10.4086/toc.2009.v005a011>`_;
                        `arXiv:0901.4709 <https://arxiv.org/abs/0901.4709>`_
