import pandas as pd
import config
import slm_utils
import sys
import figure_utils
import matplotlib.gridspec as gridspec
import miscellaneous_utils as mu
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rc('text', usetex=True)
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}') 
import seaborn as sns
import numpy as np
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition,
                                                  mark_inset)
import os
from statsmodels.tsa.stattools import adfuller

color_mood = ["#a63a33","#f1c132","#3c414d","#d7d7de","#aaaaaa"]
font = {'family': 'sans-serif',
        'color':  'k',
        'size': 40
        }


def plot_fst(ts_all,host,species=None,fst_fig=None,spec=None,ij=0):
 
    species_full = figure_utils.get_pretty_species_name(species)

    ts_plot = ts_all[host]
    dates = mu.return_dates(host)
    
    if fst_fig == None:
        fst_fig = plt.figure(constrained_layout=True,figsize=(24,6))
    
        gs = fst_fig.add_gridspec(1, 10)

        axs = fst_fig.add_subplot(gs[0, :-1])
        axs_null = fst_fig.add_subplot(gs[0, -1:])
     
    else:
        axs = fst_fig.add_subplot(spec[ij, :-1])
        axs_null = fst_fig.add_subplot(spec[ij, -1:])
              
    axs_null.axis('off')
    axs.set_xlim([-4,dates.iloc[-1] + 4])

    if species == None:
                
        for spc in ts_plot.keys():
            d,spc_frq = ts_plot[spc]
            
            axs.plot(d,spc_frq,label=spc,lw=4)
            axs.scatter(d,spc_frq,color="k",s=20,zorder=100)
            
        axs.legend(loc='upper left', bbox_to_anchor=(1.0, 1.0),ncol=1, fancybox=True, shadow=True,prop={'size': 20,'style':'italic'})
        axs.set_xlabel("Time (days)",fontdict=font)
        axs.set_title(host,fontdict=font,fontstyle="italic",fontweight="bold")
        #axs.semilogy()
        axs.text(1, .7, r"Inter-host ${F_{ST}}'$", fontsize = 25,color="red")
        
    if species != None:

        x,y =  ts_plot[species_full]
        x = np.array(x)
        y = np.array(y)
        
        axs.plot(x,y,color="#378293",label=figure_utils.get_abbreviated_species_name(species),zorder=10,lw=10);
        axs.plot(x,y,color="k",zorder=9,lw=12);

        axs.set_xlabel("",fontdict=font)
        axs.set_title(r"$\textit{%s}$" % figure_utils.get_abbreviated_species_name(species),size=60,fontstyle="italic",fontweight="bold")      
        axs.set_xticks([])
        axs.set_xticklabels([])
        axs.text(1, 1.05, r"Inter-host $\boldsymbol{{F_{ST}}'}$", fontsize = 40,color="red")
      
    axs.axhline(1,color="red",lw=2)
    
    axs.set_ylabel(r"$\boldsymbol{{F_{ST}}'}$",fontdict=font,labelpad=25)
    axs.spines['left'].set_linewidth(1)
    axs.spines['bottom'].set_linewidth(1)
    axs.spines['right'].set_linewidth(1)
    axs.spines['top'].set_linewidth(1)
  
    axs.tick_params('both', length=4, width=2, which='major',labelsize=25)
    axs.tick_params('both', length=2, width=2, which='minor')
    
    axs.set_ylim([-.25,1.5])
        
    return axs


def plot_snv_fig(chosen_strain,host,species,snv_fig=None,spec=None,ij=1):   
  
    if snv_fig == None:
        snv_fig = plt.figure(constrained_layout=True,figsize=(16,12))
    
        gs = snv_fig.add_gridspec(0, 10)

        axs= snv_fig.add_subplot(gs[0, :-1])
        axs_null = snv_fig.add_subplot(gs[0, -1:])
 
    else:
        axs = snv_fig.add_subplot(spec[ij, :-1])
        axs_null = snv_fig.add_subplot(spec[ij, -1:])        

    axs_null.axis('off')
    
    axs.set_ylabel("SNV\nfrequency",fontdict=font,labelpad=25,multialignment='center') 
    axs.spines['left'].set_linewidth(1)
    axs.spines['bottom'].set_linewidth(1)
    axs.spines['right'].set_linewidth(1)
    axs.spines['top'].set_linewidth(1)
  
    axs.tick_params('both', length=3, width=2, which='major',labelsize=25)
    axs.tick_params('both', length=3, width=2, which='minor')

    axs.set_xticklabels([])

    analysis_dir = config.analysis_directory
    cluster_dir = "%s/clusters/Poyet/%s/%s" % (analysis_dir,host,species)
    strain_snv_dic = {}
    strain_centroid_dic = {}
    
    i = 1
    for strain in os.listdir(cluster_dir):

        df = pd.read_csv("%s/%s" % (cluster_dir,strain),index_col=0)
        strain_snv_dic[i] = df
        strain_centroid_dic[i] = df.median()

        i+=1

    if len(strain_centroid_dic) == 1:
        strain_snv_dic[2] = 1 - strain_snv_dic[1]
        strain_centroid_dic[2] = 1 - strain_centroid_dic[1]
    
    letter_list = ["A","B","C"]

    dates = mu.return_dates(host)
    
    dates = dates.loc[strain_snv_dic[1].columns].sort_values()

    axs.set_xlim([-4,dates.iloc[-1] + 4])
    
    strain_lines = []

    i = 0
    for strain in strain_centroid_dic:

        if letter_list[i] == chosen_strain:
            xx = 5000
            yy = .01
            lww = 2
            C = strain_snv_dic[strain].sample(min(strain_snv_dic[strain].shape[0],xx)).T.values
            axs.plot(dates.values,strain_centroid_dic[strain],color="k",label=letter_list[i],lw=4,alpha=.75,zorder=10);
            axs.plot(dates.values,C,lw=lww,color="#378293",alpha=yy);
            strain_line = Line2D([0], [0], label=letter_list[i], color=color_mood[i])
            strain_lines.append(strain_line)
            
        else:
            xx = 2500
            yy = .01
            lww = 2
            C = strain_snv_dic[strain].sample(min(strain_snv_dic[strain].shape[0],xx)).T.values
            axs.plot(dates.values,C,lw=lww,color="grey",alpha=yy);
            axs.plot(dates.values,strain_centroid_dic[strain],color="k",label=letter_list[i],lw=2,alpha=.75,zorder=10);
            strain_line = Line2D([0], [0], label=letter_list[i], color=color_mood[i])
            strain_lines.append(strain_line)
            
        i+=1
        
    axs.set_xticks([])
    plt.subplots_adjust(wspace=0.01, hspace=0);
    
def plot_slm_fig(host,species,chosen_strain,fig3=None,spec=None,ij=2):
        
    df = pd.read_csv("strains_%s.csv" % host,index_col=0) 
    dates = mu.return_dates(host)
    
    strain = "%s_%s" % (species,chosen_strain)
    obs_data = df.loc[strain]
    obs_data = obs_data[obs_data.notna()]

    dates = dates.loc[obs_data.index].sort_values()

    obs_data = obs_data.loc[dates.index]
    
    train_num = len(obs_data)//3

    params = slm_utils.fit_SLM_params(obs_data,n=train_num)

    ## initialize SLM 
    S = slm_utils.slm(sigma=params["sigma"],K=params["K"],tau=1,delta_t=1.0/100,init_val=obs_data[0])
    S.run_sim(num_iters=int(dates.iloc[-1]/S.delta_t),num_reps=4000,record_steps=True);
    
    if fig3 == None:
        fig3 = plt.figure(constrained_layout=True,figsize=(16,12))
        
        gs = fig3.add_gridspec(1, 10)

        f3_ax1 = fig3.add_subplot(gs[0, :-1])
        f3_ax1_H = fig3.add_subplot(gs[0, -1:])
    
    else:
        f3_ax1 = fig3.add_subplot(spec[ij, :-1])
        f3_ax1_H = fig3.add_subplot(spec[ij, -1:])                
    
    f3_ax1.set_xticks([])
    f3_ax1.set_ylabel("Strain\nfrequency",fontdict=font,labelpad=25,multialignment='center') 

    f3_ax1.plot(dates.values,obs_data.values,"-",color="#378293",lw=6,zorder=100)
    f3_ax1.plot(dates.values,obs_data.values,"-",color="k",lw=8,zorder=99)

    f3_ax1.set_xlim([-4,dates.iloc[-1] + 4])
    sim_data = np.array(S.trajectory)[[int(d/S.delta_t) for d in dates.values]]
    
    plot_min = min(min(obs_data),min([min(s) for s in sim_data]))*.75
    plot_max = 1.1*max([max(s) for s in sim_data])
    
    plot_max = max((max(1.1*obs_data),plot_max))
        f3_ax1.plot(dates.values,sim_data,color="red",alpha=.01);

    f3_ax1_H.tick_params(labelleft=False, labelbottom=False)

    f3_ax1_H.set_xticks([])
    f3_ax1_H.set_yticks([])

    f3_ax1.spines['right'].set_linewidth(1)
    f3_ax1_H.spines['left'].set_linewidth(1)
    f3_ax1.tick_params('y', length=4, width=2, which='major',labelsize=25)
    f3_ax1.tick_params('y', length=2, width=2, which='minor')

    f3_ax1.spines['right'].set_linewidth(1)
    f3_ax1.spines['top'].set_linewidth(1)
    f3_ax1.spines['bottom'].set_linewidth(1)
    f3_ax1.spines['left'].set_linewidth(1)

    f3_ax1_H.hist(obs_data, bins=15, density=True, color="#378293", orientation='horizontal',edgecolor="k",lw=.3)

    f3_ax1.set_ylim([plot_min,plot_max])
    f3_ax1_H.set_ylim([plot_min,plot_max])

    xx = np.linspace(plot_min,plot_max,100000)

    f3_ax1_H.plot(S.afd.pdf(xx),xx,color="red",alpha=.5,lw=5)
    f3_ax1_H.fill_betweenx(xx,S.afd.pdf(xx),color="red",alpha=.3)
    f3_ax1_H.axis('off')

    plt.subplots_adjust(wspace=0.01, hspace=0);
    
def make_sampling_ax(host,species,samp_fig=None,spec=None,ij=3):
    
    df = pd.read_csv("strains_%s.csv" % host,index_col=0) 

    dates = mu.return_dates(host)
    
    strain = "%s_%s" % (species,chosen_strain)
    obs_data = df.loc[strain]
    obs_data = obs_data[obs_data.notna()]

    dates = dates.loc[obs_data.index].sort_values()    
    
    sampling_ax = samp_fig.add_subplot(spec[ij, :-1])
      
    sampling_ax.set_xlabel("Time (days)",fontdict=font)
 
    sampling_ax.set_yticks([])
    
    sampling_ax.set_ylim([.25,.75])
            
    sampling_ax.scatter(dates.values,.5*np.ones(len(dates.values)),s=1500,facecolor="tab:blue",marker="|",zorder=20)

    sampling_ax.set_xlim([-4,dates.iloc[-1] + 4])
 
    sampling_ax.tick_params(axis='x', which='major',length=4, width=3, labelsize=35)

    
### Collect strain abundance and Fst time series for each species in each host   
ts_all = {}
F_st_all = {}
for host in ["am","an","ae","ao"]:
    
    dates = mu.return_dates(host)
    
    pi_dir = "%s/pi/Poyet/%s" % (config.analysis_directory,host)

    good_species = [s[:-7] for s in os.listdir(pi_dir)]

    all_samples = dates.index
    
    Fst_prime = pd.read_csv(f"fst_{host}.csv",index_col=0)
    Fst_prime.columns = np.array(Fst_prime.columns).astype(int)
    Fst_prime = Fst_prime.astype(float)
    
    H = {}
    ts = {}
    for s in Fst_prime.index:
        st = s
        y = Fst_prime.loc[s].dropna()
        ts[st] = (y.index,y.values)
    
    H_all[host] = H
    ts_all[host] = ts
    F_st_all[host] = Fst_prime    
    
    
##### For a chosen species/strain enter the host, species name, and chosen strain letter

host = sys.argv[1]
species = sys.argv[2]
chosen_strain = sys.argv[3]  
      
fig = plt.figure(figsize=(24,18))
spec = gridspec.GridSpec(ncols=10, nrows=4, height_ratios = (7,7,7,1), figure=fig,hspace=.025)
plot_fst(ts_all,host,species=species,fst_fig=fig,spec=spec)
plot_snv_fig(chosen_strain,host,species,snv_fig=fig,spec=spec)
plot_slm_fig(host,species,chosen_strain,fig3=fig,spec=spec)
make_sampling_ax(host,species,samp_fig=fig,spec=spec)
fig.savefig(f"{species}_{chosen_strain}_{host}_fig1")