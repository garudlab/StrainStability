import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.special as special
from scipy.optimize import curve_fit
import miscellaneous_utils as misc

pval_thresh = 0.05

def powlaw(x, a, b) :
    return a * np.power(x, b)

def linlaw(x, a, b) :
    return a + x * b

def curve_fit_log(xdata, ydata) :
    """Fit data to a power law with weights according to a log scale"""
    # Weights according to a log scale
    # Apply fscalex
    xdata_log = np.log10(xdata)
    # Apply fscaley
    ydata_log = np.log10(ydata)
    # Fit linear
    popt_log, pcov_log = curve_fit(linlaw, xdata_log, ydata_log)
    #print(popt_log, pcov_log)
    # Apply fscaley^-1 to fitted data
    ydatafit_log = np.power(10, linlaw(xdata_log, *popt_log))
    # There is no need to apply fscalex^-1 as original data is already available
    
    return (popt_log, pcov_log, ydatafit_log)

## Adapted from code/derivation by Jacopo Grilli
def get_gamma_prediction(x_range, k):
    x_range = x_range/np.sqrt(2)
    k_digamma = special.digamma(k)
    k_trigamma = special.polygamma(1,k)

    gammalog = k*k_trigamma*x_range - np.exp(np.sqrt(k_trigamma)*x_range + k_digamma) - np.log(special.gamma(k)) + k*k_digamma + np.log10(np.exp(1))

    return 10**(gammalog)


df_am = pd.read_csv("strains_am.csv",index_col=0)
df_ao = pd.read_csv("strains_ao.csv",index_col=0)
df_an = pd.read_csv("strains_an.csv",index_col=0)
df_ae = pd.read_csv("strains_ae.csv",index_col=0)

df_am_chisq = pd.read_csv("/u/home/r/rwolff/strain_stability_revisions/strainstability/analysis/chisq/Poyet/am_strain_chisq_test_cross.txt",index_col=0,header=None)
df_ao_chisq = pd.read_csv("/u/home/r/rwolff/strain_stability_revisions/strainstability/analysis/chisq/Poyet/ao_strain_chisq_test_cross.txt",index_col=0,header=None)
df_an_chisq = pd.read_csv("/u/home/r/rwolff/strain_stability_revisions/strainstability/analysis/chisq/Poyet/an_strain_chisq_test_cross.txt",index_col=0,header=None)
df_ae_chisq = pd.read_csv("/u/home/r/rwolff/strain_stability_revisions/strainstability/analysis/chisq/Poyet/ae_strain_chisq_test_cross.txt",index_col=0,header=None)

df_am = df_am.loc[df_am_chisq.index]
df_ao = df_ao.loc[df_ao_chisq.index]
df_an = df_an.loc[df_an_chisq.index]
df_ae = df_ae.loc[df_ae_chisq.index]

df_am = df_am/df_am.sum()
df_ao = df_ao/df_ao.sum()
df_an = df_an/df_an.sum()
df_ae = df_ae/df_ae.sum()

df_am = df_am.dropna(axis=1,how="all")
df_an = df_an.dropna(axis=1,how="all")
df_ao = df_ao.dropna(axis=1,how="all")
df_ae = df_ae.dropna(axis=1,how="all")

df = df_am.append(df_ao).append(df_an).append(df_ae)
df_am_fail = df_am.loc[(df_am_chisq[df_am_chisq < pval_thresh].dropna()).index]
df_ao_fail = df_ao.loc[(df_ao_chisq[df_ao_chisq < pval_thresh].dropna()).index]
df_an_fail = df_an.loc[(df_an_chisq[df_an_chisq < pval_thresh].dropna()).index]
df_ae_fail = df_ae.loc[(df_ae_chisq[df_ae_chisq < pval_thresh].dropna()).index]

df_am_pass = df_am.loc[(df_am_chisq[df_am_chisq > pval_thresh].dropna()).index]
df_ao_pass = df_ao.loc[(df_ao_chisq[df_ao_chisq > pval_thresh].dropna()).index]
df_an_pass = df_an.loc[(df_an_chisq[df_an_chisq > pval_thresh].dropna()).index]
df_ae_pass = df_ae.loc[(df_ae_chisq[df_ae_chisq > pval_thresh].dropna()).index]

df_fail = df_am_fail.append(df_ao_fail).append(df_an_fail).append(df_ae_fail)
df_pass = df_am_pass.append(df_ao_pass).append(df_an_pass).append(df_ae_pass)

## make taylor's law figure
meansh = df.T.mean()
varh = df.T.var()
meansh_pass = df_pass.T.mean()
varh_pass = df_pass.T.var()
meansh_fail = df_fail.T.mean()
varh_fail = df_fail.T.var()
C=curve_fit_log(meansh.values, varh.values)
popt_log, pcov_log, ydatafit_log = C
popt_log[0] = 10**popt_log[0]
xdata_log = np.log10(meansh.values)
xdata_log_space = np.logspace(np.log10(.3*min(meansh.values)),0,100000)

fig_TL,TL_axs = plt.subplots(figsize=(10,10))
TL_axs.tick_params(axis='both', which='major', labelsize=15)
TL_axs.tick_params(axis='both', which='major', labelsize=15)
TL_axs.spines['top'].set_visible(False)
TL_axs.spines['right'].set_visible(False)
TL_axs.spines['left'].set_linewidth(3)
TL_axs.spines['bottom'].set_linewidth(3)
TL_axs.tick_params('both', length=7, width=2, which='major')
TL_axs.tick_params('both', length=3, width=2, which='minor')
fig_TL.suptitle("",size=30)
TL_axs.scatter(meansh.values, varh.values,color="k",s=75)
TL_axs.scatter(meansh_fail.values,varh_fail.values,color= "red",s=75)
TL_axs.plot(xdata_log_space,powlaw(xdata_log_space,*popt_log),color="green",lw = 12,alpha=.5,label= r"$\sigma_{x_i}^2 \propto \langle x_i \rangle^{%s}$" % np.around(popt_log[1],2),zorder=1000)
TL_axs.loglog()
TL_axs.set_xlim([1e-3,1.5])
TL_axs.legend(loc='lower right', prop={'size': 25})
TL_axs.set_xlabel(r"$\langle x_i \rangle$",size=25)
TL_axs.set_ylabel(r"$\sigma_{x_i}^2$",size=25,rotation=90)
fig_TL.savefig("TaylorsLaw_fig")

## make gamma AFD figure
numbins = 20
xbins = np.linspace(-3,3,numbins)
a = []
for strain in df.index:
    
    s = df.loc[strain].values
    s = s[s>0]
    sfreq = list(np.log10(s))
    a.append(sfreq)

SFREQ = []
A = []
for elem in a:
    elem = np.array(elem)
    rescaled = list((elem - np.mean(elem))/(np.std(elem)))
    A.extend(rescaled)
    h = list(plt.hist(rescaled,density=True,bins=xbins))
    H = h[0]
    SFREQ.append((xbins[1:][H>0],H[H>0]))

numbins = 20
xbins = np.log10(np.logspace(-3,3,numbins))
h = list(plt.hist(A,density=True,bins=xbins))
plt.close()

gm_fitlog,opt=curve_fit(get_gamma_prediction, h[1][1:], h[0])

fig_RG,f2_ax_RG = plt.subplots(figsize=(12,8))
f2_ax_RG.tick_params(axis='both', which='major', labelsize=15)
f2_ax_RG.spines['top'].set_visible(False)
f2_ax_RG.spines['right'].set_visible(False)
f2_ax_RG.spines['left'].set_linewidth(3)
f2_ax_RG.spines['bottom'].set_linewidth(3)
f2_ax_RG.tick_params('both', length=5, width=2, which='major')

x_range = np.logspace(-3,3,100000)

x_lin_space = np.linspace(-3 , 3, 10000)
f2_ax_RG.plot(x_lin_space, get_gamma_prediction(x_lin_space,*gm_fitlog), zorder=200,lw=8,alpha=.6,label="Gamma AFD")

f2_ax_RG.plot(h[1][1:],h[0],'.',markersize=20,c="k")

for e in SFREQ:
    xx,yy = e[0],e[1]
    f2_ax_RG.plot(xx,yy,color="grey",alpha=.1,zorder=1)                                                                                                                            
                                                                                                                              
f2_ax_RG.set_ylim([1e-3,1.3])
f2_ax_RG.semilogy()

f2_ax_RG.tick_params('y', length=3, width=2, which='minor')

f2_ax_RG.set_xlabel("Rescaled log relative abundance",size=20)
f2_ax_RG.set_ylabel("Probability density",size=20);
f2_ax_RG.legend(loc="upper left",bbox_to_anchor=(.8, 1.05),prop={'size': 20})
fig_RG.savefig("GammaAFD")