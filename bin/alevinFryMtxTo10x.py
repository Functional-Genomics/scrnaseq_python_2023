#!/usr/bin/env python

# Read the results of Alevin-fry and write outputs to a .mtx file readable by tools
# expecting 10X outputs. Adapted from 
# https://github.com/ebi-gene-expression-group/scxa-droplet-quantification-workflow/blob/develop/bin/alevinMtxTo10x.py
# which is adapted from 
# https://github.com/k3yavi/vpolo/blob/master/vpolo/alevin/parser.py 

from __future__ import print_function
from collections import defaultdict
from struct import Struct
import pandas as pd
import gzip
import sys
import os
from scipy.io import mmread,mmwrite
from scipy.sparse import *
from shutil import copyfile
import pathlib
import numpy as np
import argparse
import pyroe

parser = argparse.ArgumentParser(description='Convert Alevin outputs to 10X .mtx.')
parser.add_argument('alevin_fry_quant', help = 'Alevin output directory')
parser.add_argument('mtx_out', help = 'Output directory for converted results')
parser.add_argument('--cell_prefix', dest='cell_prefix', default='', help = 'Prefix to apply to cell barcodes')
parser.add_argument('exp_type', help='AEExperimentType')
args = parser.parse_args() 

alevin_fry_quant=args.alevin_fry_quant
mtx_out=args.mtx_out
cell_prefix=args.cell_prefix
exp_type=args.exp_type

# Run some checks in the Alevin-fry output

if not os.path.isdir(alevin_fry_quant):
    print("{} is not a directory".format( alevin_fry_quant ))
    sys.exit(1)
    
# set exp_type to format that pyroe.load_fry can understand
if exp_type == "single_cell":
    format = "scRNA"
elif exp_type == "single_nucleus":
    format = "snRNA"
else:
    print("{} is not an allowed AEExperimentType for droplet scxa".format( exp_type ))
    sys.exit(1)
 

# Read mtx from alevin_fry_quant 
# format is either snRNA or scRNA and decides if reads mapping to unspliced transcripts are quantified
# returns anndata 
ad = pyroe.load_fry(alevin_fry_quant, output_format=format)
# Trying below based on recommendations here: https://www.sc-best-practices.org/introduction/raw_data_processing.html#
# ad = pyroe.load_fry(alevin_fry_quant, output_format={'X': ['U','S','A']})

cb_names = [cell_prefix + s for s in ad.obs_names]
gene_names = ad.var_names
umi_counts = ad.X
    
# Write outputs to a .mtx file readable by tools expecting 10X outputs.
# Barcodes file works as-is, genes need to be two-column, duplicating the
# identifiers. Matrix itself needs to have genes by row, so we transpose. 

pathlib.Path(mtx_out).mkdir(parents=True, exist_ok=True)
mmwrite('%s/matrix.mtx' % mtx_out, umi_counts.transpose(), field='real') 

genes_frame = pd.DataFrame(zip(gene_names, gene_names))
genes_frame.to_csv(path_or_buf='%s/genes.tsv' % mtx_out, index=False, sep="\t", header = False)

with open('%s/barcodes.tsv' % mtx_out, 'w') as f:
    f.write("\n".join(cb_names))    
    f.write("\n")    
