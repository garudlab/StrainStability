import pandas as pd
import sys

host = sys.argv[1]

base_dir = "/u/scratch/r/rwolff/Poyet_midas_output/merged_midas_output"

df = pd.read_csv("%s/%s/species/coverage.txt" % (base_dir,host),index_col=0,sep='\t')

pd.Series((df.loc[(df > 3).T.sum() > 1].index)).to_csv("/%s/%s/species_union.txt" % (base_dir,host),index=None,header=None)
