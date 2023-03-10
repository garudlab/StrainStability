from numba import njit 

import matplotlib.pyplot as plt
import cPickle
import pandas as pd
import config
import numpy
import random as rand

from random import randint,sample
from math import log

import sys
import os 

### use of numba/njit below *substantially* increases computational efficiency

## calculate distances for forward polarization 
@njit
def D_mat_fun1(num,F,D,D_mat):   

    for k in xrange(num - 1):
        
        O = numpy.zeros(num)
        
        di = D[k]
        fi = F[k]
        
        for i in xrange(num - k - 1):

            j = i + k + 1

            fj = F[j]
            dj = D[j]

            O[j] = 2*numpy.nanmean((di + dj)*((fi - fj)**2)/((fi + fj)*(1 - fi + 1 - fj)))        
        
        D_mat[k] = O
    
    return D_mat

## calculate distances for reverse polarization
@njit
def D_mat_fun2(num,F,D,D_mat_in):   

    for k in xrange(num - 1):
        
        O = numpy.zeros(num)       
        di = D[k]
        
        fi = 1-F[k]
        
        for i in xrange(num - k - 1):

            j = i + k + 1
            
            fj = F[j]
            dj = D[j]

            O[j] = 2*numpy.nanmean((di + dj)*((fi - fj)**2)/((fi + fj)*(1 - fi + 1 - fj)))        
        
        D_mat_in[k] = O
    
    return D_mat_in

def return_clus(D_mat_close,Fs_sub):
    D_mat_close_sorted_sum = D_mat_close.sum().sort_values()
    desired_idx = D_mat_close_sorted_sum.index[-1]
    clus_idxs = D_mat_close.loc[D_mat_close[desired_idx]].index
    
    ### only return indices which are co-clustered w/ at least .25 of other points
    idxtrue = (D_mat_close.loc[clus_idxs,clus_idxs].T.mean() > .25)
    idxtrue = idxtrue[idxtrue].index
    clus_idxs = idxtrue
    clus = Fs_sub.loc[clus_idxs]
     
    return clus,clus_idxs

def drop_clus_idxs(D_mat_close,clus_idxs):
    D_mat_close_out = D_mat_close.drop(clus_idxs).drop(clus_idxs,axis=1)
    return D_mat_close_out

def polarize_clus(clus,clus_idxs,D_mat_1,D_mat_2):
    
    ## polarize whole cluster based on polarization of first cluster element
    clus_to_pol = clus_idxs[numpy.where(D_mat_1.loc[clus_idxs[:1],clus_idxs] > D_mat_2.loc[clus_idxs[:1],clus_idxs])[1]]
    pol_clus = 1 - clus.loc[clus_to_pol]
    clus_non_pol = clus_idxs[numpy.where(D_mat_1.loc[clus_idxs[:1],clus_idxs] < D_mat_2.loc[clus_idxs[:1],clus_idxs])[1]]
    non_pol_plus = clus.loc[clus_non_pol]
    clus_pol = pd.concat([pol_clus,non_pol_plus],ignore_index=True)
    
    return(clus_pol)

@njit
def symmetrize(D_mat):
    for i in range(D_mat.shape[0]-1):
        for j in range(i,D_mat.shape[0]):
            D_mat[j][i] = D_mat[i][j]
    return(D_mat)


species = sys.argv[1]

cluster_min_SNV_size = config.cluster_min_SNV_size
strainfinder_dir = config.strainfinder_directory 
dates = config.dates

min_coverage = config.cluster_min_coverage
max_d = config.max_d
min_cluster_size = config.cluster_min_SNV_size
## optionally, threshold number of SNVs considered for computational efficiency
#max_num_snvs = 5000

sys.stderr.write("Processing %s \n\n" % species)

snp_alignment = pd.read_pickle("%s/%s.strainfinder.p" %  (strainfinder_dir ,species))
samples = pd.read_pickle("%s/%s.strainfinder.samples.p" % (strainfinder_dir ,species))
samples = [s.decode("utf-8") for s in samples]

snp_locations = pd.read_pickle("%s/%s.strainfinder.locations.p" % (strainfinder_dir,species))

cluster_As = []
cluster_Ds = []
for snp_idx in range(0,snp_alignment.shape[1]):
    
    Ds = snp_alignment[:,snp_idx,:].sum(axis=1)
    As = snp_alignment[:,snp_idx,0]
    As = numpy.reshape(As, (1,len(As)))
    Ds = numpy.reshape(Ds, (1,len(Ds)))

    cluster_As.append(As[0])
    cluster_Ds.append(Ds[0])

cluster_As = numpy.array(cluster_As)
cluster_Ds = numpy.array(cluster_Ds)

As = pd.DataFrame(cluster_As,columns=samples,index=snp_locations)
Ds = pd.DataFrame(cluster_Ds,columns=samples,index=snp_locations)

samples_sorted=list(dates.loc[samples].sort_values().index)

As = As.loc[:,samples_sorted]
Ds = Ds.loc[:,samples_sorted]

F = As/Ds

Ass = As
Dss = Ds.loc[Ass.index]
min_coverage = 10
Ass = Ass.mask(Dss < min_coverage)

Ass = Ass.mask(Ass == 0)
Ass = Ass.mask(Ass == Dss)

Fs = Ass/Dss

Fs = Fs.loc[Fs.notna().T.sum() > int(F.shape[1]/4)]

Ass = Ass.loc[Fs.index]
Dss = Dss.loc[Fs.index]

fss = Ass.values/(Dss.values + (Dss.values == 0))

cluster_As = Ass.values
cluster_Ds = Dss.values
cluster_fs = cluster_As/(cluster_Ds + (cluster_Ds == 0))

## for compatibility in case of threshold number of SNVs
#num = min(max_num_snvs,Fs.shape[0])
num = Fs.shape[0]

sys.stderr.write("Processing %s SNVs" % num)

## simply shuffles indices if no threshold is specified
i_list = sample(range(Fs.shape[0]),num)
i_list_idx = Fs.iloc[i_list].index

Ass_sub = Ass.loc[i_list_idx]
Dss_sub = Dss.loc[i_list_idx]
Fs_sub = Fs.loc[i_list_idx]

fss_sub = Ass_sub.values/(Dss_sub.values + (Dss_sub.values == 0))

cluster_As_sub = Ass_sub.values
cluster_Ds_sub = Dss_sub.values
cluster_fs_sub = cluster_As_sub/(cluster_Ds_sub + (cluster_Ds_sub == 0))

D_mat = numpy.zeros([num,num])
D_mat_1 = D_mat_fun1(num,fss_sub,cluster_Ds_sub,D_mat)
D_mat = numpy.zeros([num,num]) 
D_mat_2 = D_mat_fun2(num,fss_sub,cluster_Ds_sub,D_mat)

D_mat = numpy.fmin(D_mat_1,D_mat_2)
D_mat = symmetrize(D_mat)

D_mat_1 = pd.DataFrame(D_mat_1,index=Fs_sub.index,columns=Fs_sub.index)
D_mat_2 = pd.DataFrame(D_mat_2,index=Fs_sub.index,columns=Fs_sub.index)

D_mat_close = pd.DataFrame(D_mat < max_d) 

D_mat_close.index = Fs_sub.index
D_mat_close.columns = Fs_sub.index

### extract up to six clusters. in practice, no species had more than three clusters of sufficient size to 
### to be considered a strain

clus1,clus_idxs_1 = return_clus(D_mat_close,Fs_sub)
clus1_pol = polarize_clus(clus1,clus_idxs_1,D_mat_1,D_mat_2)
D_mat_close_1 = drop_clus_idxs(D_mat_close,clus_idxs_1)

clus2,clus_idxs_2 = return_clus(D_mat_close_1,Fs_sub)
clus2_pol = polarize_clus(clus2,clus_idxs_2,D_mat_1,D_mat_2)
D_mat_close_2 = drop_clus_idxs(D_mat_close_1,clus_idxs_2)

clus3,clus_idxs_3 = return_clus(D_mat_close_2,Fs_sub)
clus3_pol = polarize_clus(clus3,clus_idxs_3,D_mat_1,D_mat_2)
D_mat_close_3 = drop_clus_idxs(D_mat_close_2,clus_idxs_3)

# clus4,clus_idxs_4 = return_clus(D_mat_close_3,Fs_sub)
# clus4_pol = polarize_clus(clus4,clus_idxs_4,D_mat_1,D_mat_2) 
# D_mat_close_4 = drop_clus_idxs(D_mat_close_3,clus_idxs_4)

# clus5,clus_idxs_5 = return_clus(D_mat_close_4,Fs_sub)
# clus5_pol = polarize_clus(clus5,clus_idxs_5,D_mat_1,D_mat_2)
# D_mat_close_4 = drop_clus_idxs(D_mat_close_4,clus_idxs_5)

# clus6,clus_idxs_6 = return_clus(D_mat_close_5,Fs_sub)
# clus6_pol = polarize_clus(clus5,clus_idxs_6,D_mat_1,D_mat_2)
# D_mat_close_5 = drop_clus_idxs(D_mat_close_5,clus_idxs_4)


if not os.path.isdir("/u/home/r/rwolff/strain_stability_revisions/strainstability/analysis/clusters/Poyet/%s/%s" % (config.host,species)):    
    os.mkdir("/u/home/r/rwolff/strain_stability_revisions/strainstability/analysis/clusters/Poyet/%s/%s" % (config.host,species))

i = 1
for clus in [clus1_pol, clus2_pol, clus3_pol]:    

#for clus in [clus1_pol, clus2_pol, clus3_pol,clus4_pol,clus5_pol,clus6_pol]:    
    print(clus)
    if clus.shape[0] > cluster_min_SNV_size:
        clus.to_csv("/u/home/r/rwolff/strain_stability_revisions/strainstability/analysis/clusters/Poyet/%s/%s/%s_cluster_%s.csv.gz" % 
                   (config.host,species,species,str(i)),compression='gzip')
        i+=1

sys.stderr.write("%s strains detected \n" % i)

### if no clusters of sufficient size detected, only a single strain is present whose relative frequency in the species is 1
if i == 1:
    clus = pd.DataFrame(1.0, columns=samples_sorted,index=[species])
    clus.to_csv("/u/home/r/rwolff/strain_stability_revisions/strainstability/analysis/clusters/Poyet/%s/%s/%s_cluster_%s.csv.gz" % 
                   (config.host,species,species,str(i)),compression='gzip')
    
sys.stderr.write("Finished")