import config

import pandas as pd
import matplotlib.pyplot as plt

import seaborn as sns
import numpy as np
import miscellaneous_utils as mu
import figure_utils
from statsmodels.tsa.stattools import adfuller

cohort = "Poyet"
from matplotlib import rc,rcParams

def make_heatmap(host,ax=None,fig=None):
    
    df = pd.read_csv("%s/chisq/%s/%s_strain_chisq_test_cross.txt" % (config.analysis_directory,cohort,host),header=None,index_col=0)
    df_species = pd.read_csv("%s/chisq/%s/%s_species_chisq_test.txt" % (config.analysis_directory,cohort,host),header=None,index_col=0)

    df_species = df_species.squeeze()
    df_species.index = [[figure_utils.get_pretty_species_name(d) for d in df_species.index]]

    species_names = [figure_utils.get_pretty_species_name(d[:-2]) for d in df.index]
    df["species"] = species_names

    df_gb = df.groupby("species")

    df_strain_species = {}
    for species in df_gb.groups:
        df_strain_species[species] = df_gb.get_group(species)[1].values
        a=df_strain_species[species]
        df_strain_species[species] = np.pad(a, (0,3 - len(df_strain_species[species])), mode='constant',constant_values=(np.nan,))

    cmap = plt.get_cmap('RdBu', 10)
    cmap.set_under("red")
    cmap.set_over("blue")

    df_strain_species = pd.DataFrame(df_strain_species).T

    df_strain_species.columns = ["Strain A","Strain B", "Strain C"]

    df_strain_species["Species"] = np.nan

    for ind in df_strain_species.index:
        df_strain_species.loc[ind,"Species"] = df_species.loc[ind][0]
    
    if ax == None:
        fig,ax = plt.subplots(figsize=(12,8))        
    
    
    sns.heatmap(df_strain_species, ax=ax,cmap=cmap,vmin=0.05,vmax=.050001,
                linewidths=3,linecolor="white",cbar=None,annot=True,annot_kws={"fontsize":20,"color":"white"});

    ax.set_title(host,size=30, fontstyle='italic',pad=25)

    ax.tick_params(axis='x', which='major', labelsize=20,rotation=0,pad=20)
    ax.tick_params(axis='y', which='major', labelsize=20)

    ax.set_yticklabels(ax.get_yticklabels(), fontstyle='italic');
    ax.set_xticklabels(ax.get_xticklabels(), fontweight="bold");
    
pass_strain = 0
all_strain = 0
for host in ["am","ao","an","ae"]:
    df = pd.read_csv("%s/chisq/%s/%s_strain_chisq_test_cross.txt" % (config.analysis_directory,cohort,host),header=None,index_col=0)
    pass_strain +=(df[1] > .05).sum()
    all_strain += df[1].shape[0]

## 79% of strains pass
pass_strain_percentage=pass_strain/all_strain
print(f"{np.around(pass_strain_percentage,2)*100} percent of strains pass SLM \n")

pass_species = 0
all_species = 0
for host in ["am","ao","an","ae"]:
    df = pd.read_csv("%s/chisq/%s/%s_species_chisq_test.txt" % (config.analysis_directory,cohort,host),header=None,index_col=0)
    pass_species += (df[1] > .05).sum()
    all_species += df[1].shape[0]

## 86% of species pass 
pass_species_percentage=pass_species/all_species
print(f"{np.around(pass_species_percentage,2)*100} percent of species pass SLM \n")


## Now, plot figure 2 (p-value heatmap)
fig2,axs = plt.subplots(2,2,figsize=(29,17.6))

axs = axs.ravel()
make_heatmap("am",ax=axs[0],fig=fig2)
make_heatmap("ao",ax=axs[1],fig=fig2)
make_heatmap("an",ax=axs[2],fig=fig2)
make_heatmap("ae",ax=axs[3],fig=fig2)

fig2.subplots_adjust(hspace=.3,wspace=.8)
fig2.savefig("Figure_2.pdf",bbox_inches='tight')
