
import os.path
import sys
import datetime

import numpy as np
import matplotlib.pyplot as plt

import QPtomographer.channelspace
import QPtomographer.bistates
import QPtomographer.util

import qutip

import tomographer
import tomographer.jpyutil
import tomographer.querrorbars


#
# number of outcomes per Pauli pairs
#
NumSamplesPerSetting = 500
# Pauli measurement settings on one system
PauliMeasSettings = [
    [
        QPtomographer.util.projpauli(i, s)
        for s in [1, -1]
    ]
    for i in [1, 2, 3]
]
# simulate the measurements
sigmareal_X = qutip.Qobj(np.array([[ 0.6, 0.01], [0.01, 0.4]]))
Ereal_XY = qutip.Qobj(np.array([[ 0.95, 0, 0, 0.95],
                                [ 0, 0.05, 0, 0],
                                [ 0, 0,    0, 0],
                                [ 0.95, 0, 0, 1],
                               ]), dims=[[2,2],[2,2]])

d = QPtomographer.util.simulate_process_measurements(sigmareal_X, Ereal_XY, PauliMeasSettings, PauliMeasSettings,
                                                 NumSamplesPerSetting)


# Naive, bipartite sampling method

r_naiveopt = None
with tomographer.jpyutil.RandWalkProgressBar() as prg:
    r_naiveopt = QPtomographer.bistates.run(
        dimX=2, dimY=2, Emn=d.Emn, Nm=np.array(d.Nm),
        hist_params=tomographer.HistogramParams(0, 0.2, 50),
        mhrw_params=tomographer.MHRWParams(0.008, 512, 32768, 32768), # thermalize a lot
        progress_fn=prg.progress_fn,
        jumps_method='light' # use optimized random walk
        )
    prg.displayFinalInfo(r_naiveopt['final_report_runs'])
print("Calculation ran for {!s} seconds".format(datetime.timedelta(seconds=r_naiveopt['elapsed_seconds'])))
print(r_naiveopt['final_report_runs'])


analysis_naiveopt = tomographer.querrorbars.HistogramAnalysis(r_naiveopt['final_histogram'],
                                                              threshold_fraction=1e-3)

analysis_naiveopt.plot()


# Channel-space method

sigmareal_X_sqrtm_eyeY = np.kron(sigmareal_X.sqrtm().data.toarray(), np.eye(2))
Emn_for_channelspace = [
    np.dot(np.dot(sigmareal_X_sqrtm_eyeY, E), sigmareal_X_sqrtm_eyeY)
    for E in d.Emn
]

r_elr = None
with tomographer.jpyutil.RandWalkProgressBar() as prg:
    r_elr = QPtomographer.channelspace.run(
        dimX=2, dimY=2, Emn=Emn_for_channelspace, Nm=np.array(d.Nm),
        hist_params=tomographer.UniformBinsHistogramParams(0, 0.2, 50),
        channel_walker_jump_mode=QPtomographer.channelspace.ElemRotations,
        mhrw_params=tomographer.MHRWParams(0.02, 50, 2048, 32768),
        progress_fn=prg.progress_fn
        )
    prg.displayFinalInfo(r_elr['final_report_runs'])
print("Calculation ran for {!s} seconds".format(datetime.timedelta(seconds=r_elr['elapsed_seconds'])))
print(r_elr['final_report_runs'])


analysis_elr = tomographer.querrorbars.HistogramAnalysis(r_elr['final_histogram'],
                                                         threshold_fraction=1e-3)

analysis_elr.plot()


# Compare the two methods, make a figure in log scale

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlabel('Diamond Norm distance to the identity channel')
ax.set_ylabel('probability density')
# use log scale
ax.set_yscale('log')

clist = 'rgbcmyk'

# these are tomographer.Histogram objects with the normalized histogram
# data. Use h1.bins to access the histogram data as a numpy array.
h1 = r_naiveopt['final_histogram'].normalized()
h2 = r_elr['final_histogram'].normalized()

# bin locations in figure of merit
f1 = h1.values_center
f2 = h2.values_center


flist = np.linspace(np.min(f1), np.max(f1), 100)

ax.errorbar(x=f1, y=h1.bins, yerr=h1.delta, fmt='.', c='b')
ax.plot(flist, np.exp(analysis_naiveopt.fit_fn(analysis_naiveopt.ftox(flist), *analysis_naiveopt.fit_params)), c='b')

ax.errorbar(x=f2, y=h2.bins, yerr=h2.delta, fmt='.', c='r')
ax.plot(flist, np.exp(analysis_elr.fit_fn(analysis_elr.ftox(flist), *analysis_elr.fit_params)), c='r')

ax.set_ylim([1e-4, 1e2])

plt.show()
