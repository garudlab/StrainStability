import pandas as pd
import config
import figure_utils
import numpy as np
import os
import miscellaneous_utils as mu

## mean between host Fst calculated from HMP data
def calculate_BH(spcs):
    
    pi_df_BH = pd.read_csv("/u/home/r/rwolff/strain_stability_revisions/strainstability/analysis/pi/HMP/HMP/%s_pi.txt" % spcs,index_col=0)
    
    Fst = {}
    for ind1 in pi_df_BH.index:
        for ind2  in pi_df_BH.index:
            if ind1 != ind2:
                pi_w1 = pi_df_BH.loc[ind1,ind1]
                pi_w2 = pi_df_BH.loc[ind2,ind2]
            
                Fst[(ind1,ind2)] = 1 - np.mean([pi_w1,pi_w2])/pi_df_BH.loc[ind1,ind2]
            
    return Fst

## return fst_prime trajectory
def return_fst_prime(host,spc):
    dates = mu.return_dates(host)
    all_samples = mu.return_host_samples(host)
    pi_dir = "%s/pi/Poyet/%s" % (config.analysis_directory,host)
    pi_df = pd.read_csv("%s/%s_pi.txt" % (pi_dir,spc),index_col=0)
    pi_df_T = pd.DataFrame(columns=all_samples,index=all_samples)
    pi_df_T.loc[pi_df.index,pi_df.columns] = pi_df
    species_dates = dates.loc[pi_df.index].sort_values()
    first_date = species_dates.index[0]
    pi_df_BT = pi_df_T.loc[first_date]
    pi_df_W = pd.Series(np.diag(pi_df_T),index=pi_df_BT.index)
    
    Fst = pd.Series(dtype='float64' )
    for ind in pi_df_BT.index:
        Fst.loc[ind] = 1 - ((pi_df_W.loc[first_date] + pi_df_W.loc[ind])/2)/pi_df_BT.loc[ind]
        
    return Fst

for host in ["am","an","ae","ao"]:
    
    dates = mu.return_dates(host)
    
    pi_dir = "%s/pi/Poyet/%s" % (config.analysis_directory,host)

    good_species = [s[:-7] for s in os.listdir(pi_dir)]
    
    all_samples = dates.index
    
    Fst_prime = pd.DataFrame(index=good_species,columns=dates.index)
    
    
    for species in good_species:
        Fst_prime.loc[species] = return_fst_prime(host,species)
    
    Fst_prime.columns = np.array(dates.values).astype(int)
        
    Fst_prime.to_csv("fst_%s.csv" % host)
    
    