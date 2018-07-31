
# -*- coding: utf-8 -*-

"""
Utility routines for reliable process tomography.
"""

try:
    import numpy as np
    import qutip
except ImportError:
    # Allow readthedocs to parse this file without importing unknown module names
    import os
    if os.environ.get('READTHEDOCS') != 'True':
        raise


# A namespace to hold attributes
#
# Do NOT change this class name or relocate it.  Pickled output may use this
# class and it needs to be reliably there to be able to load pickled files.
#
class _Store(object):
    def __init__(self, **kwargs):
        super(_Store, self).__init__()
        self.__dict__.update(kwargs)

    def _store_items(self):
        return [(k, v) for k, v in self.__dict__.items() if k[:1] != '_']

    def __repr__(self):
        return '_Store({})'.format([','.join(['{}={}'.format(k,repr(v))
                                              for k,v in self._store_items()])])

    def __str__(self):
        return ', '.join(['{}={}'.format(k, str(v)) for k,v in self._store_items()])



# general utilities
def projpaulimat(i, v):
    """
    Projector onto the eigenspace corresponding to the eigenvalue `v` of the `i`-th
    Pauli matrix counted from 1 (as a NumPy :py:class:`array <numpy.ndarray>`)
    """
    if i == 1 and v == 1:
        return np.array([[.5, .5],[.5, .5]])
    elif i == 1 and v == -1:
        return np.array([[.5, -.5],[-.5, .5]])
    elif i == 2 and v == 1:
        return np.array([[.5, -.5j],[.5j, .5]])
    elif i == 2 and v == -1:
        return np.array([[.5, .5j],[-.5j, .5]])
    elif (i == 3 or i == 0 or i == 4):
        if v == 1:
            return np.array([[1,0],[0,0]])
        elif v == -1:
            return np.array([[0,0],[0,1]])
    raise ValueError("Invalid Pauli number or eigenvalue: i="+repr(i)+", v="+repr(v))

def projpauli(i, v):
    """
    Projector onto the eigenspace corresponding to the eigenvalue `v` of the
    `i`-th Pauli matrix counted from 1 (as a :py:class:`qutip.Qobj` object)
    """
    return qutip.Qobj(projpaulimat(i,v))

def process_matrix(sigma_A, E_AB):
    r"""
    Returns the process matrix corresponding to the process `E_AB` (represented
    by its Choi matrix) applied onto the input state `sigma_A`.  Both arguments
    must be :py:class:`qutip.Qobj` objects.

    This is simply :math:`\sigma_A^{1/2}\,E_{AB}\,\sigma_A^{1/2}`.
    """
    if (sigma_A.dims[0][0] != E_AB.dims[0][0] or len(E_AB.dims[0])!=2):
        raise ValueError("Incompatible dimensions, or E_AB has more than two systems")
    if not sigma_A.isherm:
        raise ValueError("Input sigma_A not hermitian")
    
    s = qutip.tensor(sigma_A.sqrtm(), qutip.qeye(E_AB.dims[0][1]))
    return s*E_AB*s   # s is hermitian


def simulate_process_measurements(sigma_A, E_AB, Mk_in, Mk_out, num_samples_per_setting):
    """
    Simulate measurements for process tomography. `Mk_in` and `Mk_out` describe
    the POVM effects applied on the reference system (or prepared as inputs) and
    the output system, respectively.

    Both `sigma_A` and `E_AB` should be :py:class:`qutip.Qobj`
    objects. `sigma_A` is the input state, which determines the probabilities of
    preparing each `Mk_in` effect, or alternatively, the reduced state on the
    reference system if the channel is applied onto half of a pure state.
    `E_AB` is the Choi matrix of the channel being applied, as a
    :py:class:`qutip.Qobj` density matrix (NOT as a superoperator).

    `Mk_in` and `Mk_out` must be lists of lists of POVM effects.  Each POVM
    effect is a positive semidefinite matrix of norm ≤ 1.  `Mk[k][i]` is the
    POVM effect corresponding to outcome `i` of measurement setting `k`; this
    applies both to the input and output effects.
    
    `num_samples_per_setting` specifies the number of repetitions of the same
    pair of measurement settings. (I.e., for each pair of an input and an output
    measurement setting, we collect `num_samples_per_setting` measurement
    results.)

    Returns: an object `d` with properties `d.Emn` and `d.Nm`, representing the
    POVM effects and simulated frequency counts.  They are in a format suitable
    for direct input to the C++ code.
    """

    num_in_settings = len(Mk_in)
    num_out_settings = len(Mk_out)

    num_total_in_effects = np.sum([ len(Mk_in[k]) for k in range(len(Mk_in)) ])
    num_total_out_effects = np.sum([ len(Mk_out[k]) for k in range(len(Mk_out)) ])
    num_total_effects = num_total_in_effects*num_total_out_effects

    rho_AB = process_matrix(sigma_A, E_AB)

    dimA = E_AB.dims[0][0]
    dimB = E_AB.dims[0][1]
    dimAB = dimA*dimB


    # prepare the list of POVM effects, and simulate the measurements
    Emn = []
    Nm = []

    for i in range(num_in_settings):
        for j in range(num_out_settings):

            # get some random numbers for this measurement setting
            x = np.random.rand(num_samples_per_setting)
            proboffset = 0

            # We sample the measurement outcomes as follows: we split the interval [0,1]
            # in number sections equal to the number of possible outcomes, each of length
            # = probability of that outcome.  Then, for each possible outcome, we count
            # the number of random numbers in `x` that fall into the corresponding
            # section.

            for s in range(len(Mk_in[i])):
                for t in range(len(Mk_out[j])):

                    Ek = qutip.tensor(Mk_in[i][s], Mk_out[j][t])

                    p = qutip.expect(Ek, rho_AB) # = trace(Ek * rho)
                    Emn.append( Ek.data.toarray().astype(dtype=complex) )
                    Nm.append( np.count_nonzero( (proboffset <= x) & (x < proboffset+p) ) )

                    proboffset += p

            # sanity check
            assert np.abs(proboffset-1.0) < 1e-6

    d = _Store()
    d.Emn = Emn
    d.Nm = np.array(Nm)
    return d



def simulate_process_prep_measure(E_AB, prep_meas_settings):
    """
    Simulate measurements for process tomography using the prepare-and-measure
    scheme, given a "true" quantum process `E_AB`, a list of input states and
    measurement settings specified by `prep_meas_settings`.

    The process `E_AB` should be a :py:class:`qutip.Qobj` containing the Choi
    matrix of the channel being applied, as a :py:class:`qutip.Qobj` matrix (NOT
    as a superoperator).  It should not be normalized, i.e., we expect
    :math:`\mathrm{tr}_B(E_{AB}) = \mathbb{I}_A`.

    The argument `prep_meas_settings` must be a list of tuples `(sigma_in,
    Mk_out, num_repeats)`, where `sigma_in` is an input state specified as a
    :py:class:`qutip.Qobj` object, where `Mk_out` is a POVM specified as a list
    of POVM effects (each POVM effect is a positive semidefinite matrix of norm
    ≤ 1, and is specified as a :py:class:`qutip.Qobj` matrix), and finally
    `num_repeats` is the number of times to repeat this setting.

    Returns: an object `d` with properties `d.Emn_ch`, `d.Nm`, representing the
    arguments suitable for :py:meth:`QPtomographer.channelspace.run()` and the
    simulated frequency counts.
    """

    Emn_ch = []
    Nm = []

    dimA = E_AB.dims[0][0]
    dimB = E_AB.dims[0][1]
    dimAB = dimA*dimB

    for (sigma_in, Mk_out, num_repeats) in prep_meas_settings:

        # calculate the output state after sending sigma_in through the process
        rho_out = (E_AB*qutip.tensor(qutip.partial_transpose(sigma_in, [1]), qutip.qeye(dimB))).ptrace(1)

        # output POVM
        num_out_effects = len(Mk_out)

        # get some random numbers for this measurement setting
        x = np.random.rand(num_repeats)
        proboffset = 0

        # We sample the measurement outcomes as follows: we split the interval
        # [0,1] in number sections equal to the number of possible outcomes,
        # each of length = probability of that outcome.  Then, for each possible
        # outcome, we count the number of random numbers in `x` that fall into
        # the corresponding section.

        for Mk in Mk_out:
            p = qutip.expect(Mk, rho_out) # = trace(Mk * rho_out)
            Emn_ch.append( qutip.tensor(qutip.partial_transpose(sigma_in, [1]), Mk)
                           .data.toarray().astype(dtype=complex) )
            Nm.append( np.count_nonzero( (proboffset <= x) & (x < proboffset+p) ) )
            proboffset += p

        # sanity check
        assert np.abs(proboffset-1.0) < 1e-6

    d = _Store()
    d.Emn_ch = Emn_ch
    d.Nm = np.array(Nm)
    return d
