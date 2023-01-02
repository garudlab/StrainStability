import pandas as pd
import config
import slm_utils
import figure_utils
import os
from statsmodels.tsa.stattools import adfuller

def return_test():

    H_all = {}
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
            H[st] = adfuller(y.values[1:])[1]

        H_all[host] = H
 
    return(H_all)

### create species p-value files w/ results of adf test and slm test
H_all = return_test
df_species_all = {}
for host in ["am","ao","an","ae"]:
    
    df_species = pd.read_csv("%s/chisq/%s/%s_species_chisq_test.txt" % (config.analysis_directory,cohort,host),
                             header=None,index_col=0)
    df_species.index = [figure_utils.get_pretty_species_name(f) for f in df_species.index]
    df_species.columns = ["slm"]
    df_species["adf"] = H_all[host]
    df_species_all[host] = df_species
    df_species.to_csv(f"{host}_pvalues.csv")