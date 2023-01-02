import state_utils
import pandas as pd
import numpy as np
import config
import os

letter_list = ["A","B","C"]
analysis_dir = config.analysis_directory
all_species = os.listdir("%s/clusters/Poyet/%s" % (analysis_dir,config.host))
data_directory = config.data_directory
spec_df = pd.read_csv("%sspecies/relative_abundance.txt.bz2" % data_directory,index_col=0,sep="\t")
tot_strain_df = pd.DataFrame(columns=config.host_samples)

for species in all_species:

    cluster_dir = "%s/clusters/Poyet/%s/%s" % (analysis_dir,config.host,species)

    strain_snv_dic = {}
    strain_centroid_dic = {}
    i = 0

    for strain in os.listdir(cluster_dir):

        strain_label = "%s_%s" % (species,letter_list[i])
        df = pd.read_csv("%s/%s" % (cluster_dir,strain),index_col=0)
        strain_snv_dic[strain_label] = df
        strain_centroid_dic[strain_label] = df.median()

        print(strain_centroid_dic[strain_label].min() == 1)

        i+=1

    strain_centroid_dic = pd.DataFrame(strain_centroid_dic)

    if strain_centroid_dic.shape[1] < 2:
        if strain_centroid_dic["%s_%s" % (species,letter_list[0])].min() != 1.0:
            strain_centroid_dic["%s_%s" % (species,letter_list[1])] = 1 - strain_centroid_dic["%s_%s" % (species,letter_list[0])]

    spec_strain_df = spec_df.loc[species]*strain_centroid_dic.T

    for idx in spec_strain_df.index:
        tot_strain_df.loc[idx] = spec_strain_df.loc[idx]
    
print(tot_strain_df)

#tot_strain_df.to_csv("strains_%s.csv" % config.host)