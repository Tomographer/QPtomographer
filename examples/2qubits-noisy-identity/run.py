
# coding: utf-8

# # Example: noisy identity process on two qubits
# 
# This example is featured in the paper.
# 
# 
# **Contents:**
# * [Python initializations](#Python-initializations)
# * [Data set-up](#Data-set-up)
# * [Bipartite sampling method](#Bipartite-sampling-method)
# * [Bipartite sampling method, optimized](#Bipartite-sampling-method,-optimized)
# * [Fit with our empirical model \#2](#Fit-with-our-empirical-model-#2)
# * [Prepare data for channel-space methods](#Prepare-data-for-channel-space-methods)
# * [Channel-space method, "$e^{iH}$" variant](#Channel-space-method,-%22$e^{iH}$%22-variant)
# * [Channel-space method, "elementary rotations" variant](#Channel-space-method,-%22elementary-rotation%22-variant)
# * [Grand comparison plots](#Grand-comparison-plots)
# 
# 
# ## Python initializations

# In[1]:


from __future__ import print_function
import os.path
import sys
import datetime

import numpy as np
import numpy.linalg as npl

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.font_manager
matplotlib.rc('font', family='Arial')

from IPython.display import display, Markdown, Latex

import dnormtomo.channelspace
import dnormtomo.bistates
import dnormtomo.util

# use this to load pickle files saved with our old python code
sys.modules['pydnorm_util'] = dnormtomo.util
dnormtomo.util._Ns = dnormtomo.util._Store

import qutip

import tomographer
import tomographer.jpyutil
import tomographer.querrorbars
display(Markdown("Using `tomographer` **{}** and `dnormtomo` **{}**"
                 .format(tomographer.__version__, dnormtomo.__version__)))

# to save & load computation results
try:
    import cPickle as pickle
except:
    import pickle

# interact with plots in notebook
get_ipython().magic('matplotlib notebook')


# In[2]:


# utilities for storing & loading calculation results
def save_to_cache(cachesuffix, x):
    with open('_CACHE_'+cachesuffix+'.dat', 'wb') as f:
        pickle.dump(x, f, 2)
def load_from_cache(cachesuffix):
    cachefile = '_CACHE_'+cachesuffix+'.dat'
    if not os.path.exists(cachefile):
        return None
    with open(cachefile, 'rb') as f:
        display(Markdown("Loading `{}` from cache".format(cachesuffix)))
        return pickle.load(f)


# ## Data set-up

# In[3]:


#
# number of outcomes per Pauli pairs
#
NumSamplesPerSetting = 500
# Pauli measurement settings on one system
PauliMeasSettings = [
    [
        qutip.Qobj(qutip.tensor(dnormtomo.util.projpauli(i1, s1),
                                dnormtomo.util.projpauli(i2, s2)).data.toarray(),
                   dims=[[4],[4]])
        for s1 in [1, -1]
        for s2 in [1, -1]
    ]
    for i1 in [1, 2, 3]
    for i2 in [1, 2, 3]
]
#display(PauliMeasSettings)

#
# the "real" states & process from which we simulate outcomes
#
sigmareal_X = qutip.Qobj(np.array([[ 0.35, 0.00, 0.04, 0.1j],
                                   [ 0.00, 0.15, 0.05, 0.00],
                                   [ 0.04, 0.05, 0.32, 0.00],
                                   [-0.1j, 0.00, 0.00, 0.18]]), dims=[[4],[4]])
display(Markdown("Eigenvalues of sigmareal_X = $[" + 
                 ",".join("{:.4g}".format(x) for x in npl.eigvalsh(sigmareal_X.data.toarray()))
                 + "]$"))
#display(sigmareal_X.eigenstates())
MaxEntgl_XY = qutip.Qobj(np.array([ [ 1. if i==j else 0. ] for i in range(4) for j in range(4)]),
                         dims=[[4,4],[1,1]])
#display(MaxEntgl_XY.data.toarray())
Ereal_XY = 0.9*MaxEntgl_XY*MaxEntgl_XY.dag() + 0.1*qutip.Qobj(np.identity(16), dims=[[4,4],[4,4]])/4
#display(Ereal_XY.tr())#.data.toarray().diagonal())
#rho_AB = pydnorm_util.process_matrix(sigmareal_X, Ereal_XY)
#display(rho_AB)
#display(rho_AB.data.diagonal())
#display(rho_AB.tr())

def simulate_measurements():
    #
    # Simulate the outcomes
    #
    d = dnormtomo.util.simulate_process_measurements(sigmareal_X, Ereal_XY,
                                                     PauliMeasSettings,
                                                     PauliMeasSettings,
                                                     NumSamplesPerSetting)
    
    return d

#
# Only simulate the measurements once. After that, use the same data when comparing methods!!
#
d = load_from_cache('meas_data')
if d is None:
    d = simulate_measurements()
    save_to_cache('meas_data', d)

#print(d.__dict__) # prints Emn, Nm ... long outputs!!


# In[4]:


#
# Stuff for the analysis, later
#

def print_report(r):
    display(Markdown("Calculation ran for **{!s} seconds**".format(datetime.timedelta(seconds=r['elapsed_seconds']))))
    display(Markdown("```\n{}\n```".format(r['final_report_runs'])))

def do_analysis(r, name, plots=False):
    final_histogram = r['final_histogram']
    analysis = tomographer.querrorbars.HistogramAnalysis(final_histogram)
    fitparams = analysis.printFitParameters()
    analysis.printQuantumErrorBars()
    
    p1 = None
    p2 = None
    
    if plots:
        p1 = analysis.plot(show_plot=False) # linear scale
        p1.ax.set_title("Distribution of the diamond norm: " + name)
        p1.ax.set_xlabel('Diamond Norm distance to the identity channel')
        p1.show()

        p2 = analysis.plot(log_scale=True, show_plot=False) # log scale
        p2.ax.set_title("Distribution of the diamond norm: " + name)
        p2.ax.set_xlabel('Diamond Norm distance to the identity channel')
        p2.show()
    
    return {'r': r, 'name': name, 'analysis': analysis, 'fitparams': fitparams, 'p1': p1, 'p2': p2}


# ## Bipartite sampling method

# In[5]:


r_naive = load_from_cache('r_naive')
if r_naive is None:
    # perform calculation 
    with tomographer.jpyutil.RandWalkProgressBar() as prg:
        r_naive = pydnormbistates.run(
            dimX=4, dimY=4, Emn=d.Emn, Nm=np.array(d.Nm),
            hist_params=tomographer.UniformBinsHistogramParams(0.1, 0.3, 100),
            mhrw_params=tomographer.MHRWParams(0.001, 1000, 2048, 32768),
            progress_fn=prg.progress_fn,
            progress_interval_ms=2000,
            )
        prg.displayFinalInfo(r_naive['final_report_runs'])
    save_to_cache('r_naive', r_naive)

print_report(r_naive)


# In[6]:


a_naive = do_analysis(r_naive, 'bip.', plots=True)


# ## Bipartite sampling method, optimized

# In[7]:


r_naiveopt = load_from_cache('r_naiveopt')
if r_naiveopt is None:
    # perform calculation 
    with tomographer.jpyutil.RandWalkProgressBar() as prg:
        r_naiveopt = dnormtomo.bistates.run(
            dimX=4, dimY=4, Emn=d.Emn, Nm=np.array(d.Nm),
            hist_params=tomographer.UniformBinsHistogramParams(0.1, 0.3, 100),
            mhrw_params=tomographer.MHRWParams(0.001, 4000, 8192, 32768),
            progress_fn=prg.progress_fn,
            progress_interval_ms=2000,
            jumps_method='light',
            )
        prg.displayFinalInfo(r_naiveopt['final_report_runs'])
    save_to_cache('r_naiveopt', r_naiveopt)

print_report(r_naiveopt)


# In[8]:


a_naiveopt = do_analysis(r_naiveopt, 'bip., opt.', plots=True)


# ## Fit with our empirical model \#2

# In[9]:


def fit_fn_test(x, a2, a1, m, p, c):
    return -a2* np.square(x) - a1*x - m*np.power(-np.log(x),p) + c
    
a = tomographer.querrorbars.HistogramAnalysis(
    r_naiveopt['final_histogram'], fit_fn=fit_fn_test, bounds=((0,-np.inf,0,0,-np.inf), np.inf)
)
a.printFitParameters()
a.plot(plot_deskewed_gaussian=False, log_scale=True)


# ## Prepare data for channel-space methods

# In[10]:


# we need to encode the input state in the POVM effects

sigmareal_X_sqrtm_eyeY = np.kron(sigmareal_X.sqrtm().data.toarray(), np.eye(4))
Emn_for_channelspace = [
    np.dot(np.dot(sigmareal_X_sqrtm_eyeY, E), sigmareal_X_sqrtm_eyeY)
    for E in d.Emn
]


# ## Channel-space method, "$e^{iH}$" variant

# In[11]:


r_eiH = load_from_cache('r_eiH')
if r_eiH is None:
    # no stored result, perform computation
    with tomographer.jpyutil.RandWalkProgressBar() as prg:
        r_eiH = pydnormchannelspace.run(
            dimX=4, dimY=4, Emn=Emn_for_channelspace, Nm=np.array(d.Nm),
            hist_params=tomographer.UniformBinsHistogramParams(0.1, 0.2, 100),
            channel_walker_jump_mode=pydnormchannelspace.RandHermExp,
            mhrw_params=tomographer.MHRWParams(0.001, 1000, 4096, 32768),
            progress_fn=prg.progress_fn,
            progress_interval_ms=2000,
            ctrl_converged_params={'enabled':False},
            )
        prg.displayFinalInfo(r_eiH['final_report_runs'])
    save_to_cache('r_eiH', r_eiH)
print_report(r_eiH)


# In[12]:


a_eiH = do_analysis(r_eiH, '"e^iH"')


# ## Channel-space method, "elementary-rotation" variant

# In[13]:


r_elr = load_from_cache('r_elr')
if r_elr is None:
    # no stored result, perform computation
    with tomographer.jpyutil.RandWalkProgressBar() as prg:
        r_elr = pydnormchannelspace.run(
            dimX=4, dimY=4, Emn=Emn_for_channelspace, Nm=np.array(d.Nm),
            hist_params=tomographer.UniformBinsHistogramParams(0.1, 0.2, 100),
            channel_walker_jump_mode=pydnormchannelspace.ElemRotations,
            mhrw_params=tomographer.MHRWParams(0.005, 500, 4096, 32768),
            progress_fn=prg.progress_fn,
            progress_interval_ms=2000,
            )
        prg.displayFinalInfo(r_elr['final_report_runs'])
    save_to_cache('r_elr', r_elr)
print_report(r_elr)


# In[14]:


a_elr = do_analysis(r_elr, 'el. rot.')


# # Grand comparison plots

# In[15]:


def fit_fn_q_lnxp(x, a, xq, m, p, c):
    return -a*np.square(x-xq) - m*np.power(-np.log(x)/-np.log(0.16), p) + c


def do_comparison_plot(alist, log_scale=False):

    fig = plt.figure(figsize=(8,3.5))
    ax1 = fig.add_subplot(121)
    ax1.set_xlabel('Diamond Norm distance to the identity channel')
    ax1.set_ylabel('probability density')
    
    ax2 = fig.add_subplot(122)
    ax2.set_xlabel('Diamond Norm distance to the identity channel')
    ax2.set_ylabel('probability density')
    ax2.set_yscale('log')
    
    #fig.subplots_adjust(left=0.08, right=0.99, bottom=0.1, top=.98)
    
    clist = 'crbgmyk'
    
    class _Ns: pass
    thelabels = _Ns()
    thelabels.d = {}
    
    def add_label(s):
        k = len(thelabels.d)
        thelabels.d[k] = s
        return str(k)
    def last_label():
        return thelabels.d.get(len(thelabels.d)-1)
    
    qeb = dict()
    
    for i in range(len(alist)):
        a = alist[i]
        
        print("Taking care of plots for {}".format(a['name']))
        
        r = a['r']
        h = r['final_histogram'].normalized()
        c = clist[i%len(clist)]
        f = h.values_center
        #analysis = a['analysis']
        
        analysis_dflt = tomographer.querrorbars.HistogramAnalysis(h)
        analysis_dflt.printFitParameters()
        qeb[a['name']] = analysis_dflt.printQuantumErrorBars()
        
        analysis = tomographer.querrorbars.HistogramAnalysis(h, fit_fn=fit_fn_q_lnxp, maxfev=10000,
                                                             p0=(12, -13, 0.057, 8.8, 2100),
                                                             bounds=((0,-np.inf,0,0,-np.inf),np.inf))
        analysis.printFitParameters()
        
        ax1.errorbar(x=f, y=h.bins, yerr=h.delta, c=c, fmt='.', label=add_label('{}, data'.format(a['name'])))
        ax2.errorbar(x=f, y=h.bins, yerr=h.delta, c=c, fmt='.', label=last_label())
        
        flist = np.linspace(np.min(f), np.max(f), 100)
        
        plist_dflt = np.exp(analysis_dflt.fit_fn(analysis_dflt.ftox(flist), *analysis_dflt.fit_params))
        ax1.plot(flist, plist_dflt, c=c, ls=':', label=add_label('{}, fit1'.format(a['name'])))
        ax2.plot(flist, plist_dflt, c=c, ls=':', label=last_label())
        
        plist = np.exp(analysis.fit_fn(analysis.ftox(flist), *analysis.fit_params))
        ax1.plot(flist, plist, c=c, label=add_label('{}, fit2'.format(a['name'])))
        ax2.plot(flist, plist, c=c, label=last_label())
    
    ax1.set_xlim([0.1, 0.25])
    ax2.set_xlim([0.1, 0.25])
    ax2.set_ylim([1e-5, 1e2])
    
    handles, labels = ax1.get_legend_handles_labels()
    iorder = sorted(range(len(handles)), key=lambda i: int(labels[i]))
    
    ax1.legend([handles[i] for i in iorder], [thelabels.d[int(labels[i])] for i in iorder])
    
    q1 = qeb['el. rot.']
    ax2.text(0.143, 3e-4, "$v_0=${:.3f}\n$\\Delta=${:.3f}\n$\\gamma=${:.1e}".format(q1.f0, q1.Delta, q1.gamma),
             color='#000080')
    q2 = qeb['bip.']
    ax2.text(0.211, .1, "$v_0=${:.3f}\n$\\Delta=${:.3f}\n$\\gamma=${:.1e}".format(q2.f0, q2.Delta, q2.gamma),
             color='#a00000')
    
    abfont = matplotlib.font_manager.FontProperties()
    abfont.set_weight('bold')
    abfont.set_size(12)
    ax1.text(-0.16, 1, "a.", transform=ax1.transAxes, fontproperties=abfont)
    ax2.text(-0.19, 1, "b.", transform=ax2.transAxes, fontproperties=abfont)
    fig.tight_layout()
    
    fig.savefig('comparison_plot.pdf', format='pdf')
    plt.show()
    
do_comparison_plot([a_naive, a_naiveopt, a_eiH, a_elr])


# In[ ]:





# In[ ]:





# In[ ]:




