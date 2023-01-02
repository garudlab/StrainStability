#!/usr/bin/env python 
### This script runs the necessary post-processing of the MIDAS output so that we can start analyzing
import os
import sys
import parse_midas_data

########################################################################################
#
# Standard header to read in argument information
#
########################################################################################
if len(sys.argv)>1:
    if len(sys.argv)>2:
        debug=True # debug does nothing in this script
        species_name=sys.argv[2]
    else:
        debug=False
        species_name=sys.argv[1]
else:
    sys.stderr.write("Usage: python postprocess_midas_data.py [debug] species_name")
########################################################################################


sys.stderr.write('Postprocessing species: %s\n' % species_name)

# the following creates this file: marker_coverage.txt.bz2
# It consists of a line recapitulating MIDAS output in terms of coverage for the species of interest
# It also outputs a line summing over the coverage across all species for each sample. 
sys.stderr.write('Calculating species-specific marker gene coverages...\n')
os.system('python %scalculate_marker_gene_coverage.py %s' % (parse_midas_data.scripts_directory, species_name))   
sys.stderr.write('Done calculating species-specific marker gene coverages!\n')

# Calculate error pvalues
# this produces the file annotated_snps.txt.bz2, which contains SNPs that fall between 0.3*median and 3*median, where median=median coverage of a SNP in a sample. The output is in the form of Alt, Ref, where Ref=consensus allele across samples (so, the output is polarized relative to the major allele in the sample). 
sys.stderr.write('Calculating error pvalues...\n')
os.system('python %scalculate_error_pvalues.py %s' % (parse_midas_data.scripts_directory, species_name))
sys.stderr.write('Done calculating error pvalues!\n')

sys.stderr.write("Done postprocessing %s!\n\n" % species_name)
