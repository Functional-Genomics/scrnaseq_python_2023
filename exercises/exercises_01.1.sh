#!/usr/bin/env bash

# Set working directory path and subdirs
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

data_dir="../../training_data"
fastq_dir="${data_dir}/ex_read_fastq"
ref_dir="${data_dir}/ex_mouse_ref"

# Add scripts path to $PATH
export PATH=$SCRIPT_DIR/../bin:$PATH

# The conda envs relevant in this are exercise are
# pyroe, af, dropletutils, and sc_py_training
conda_base=$CONDA_PREFIX

# Make splici reference
source $conda_base/bin/activate pyroe

pyroe make-splici \
${ref_dir}/mus_musculus_genome.fa \
${ref_dir}/mus_musculus_genes.gtf \
98 \
splici_rl98_ref &&	

source $conda_base/bin/deactivate

# Index the reference
source $conda_base/bin/activate af

salmon index \
-t $(ls splici_rl98_ref/*\.fa) \
-i salmon_index \
-p 8 &&

# Mapping
salmon alevin \
-i salmon_index \
-l ISR \
-1 ${fastq_dir}/SLX-7632.TAAGGCGA.N701.s_1.r_1.fq-400k.fastq \
-2 ${fastq_dir}/SLX-7632.TAAGGCGA.N701.s_1.r_2.fq-400k.fastq \
-p 8 \
-o salmon_alevin \
--dropseq \
--sketch &&

# Cell barcode correction with force
# "INFO total number of distinct corrected barcodes : 61,328"
alevin-fry generate-permit-list \
-f 100000 \
-d fw \
-i salmon_alevin \
-o alevin_fry_gpl_f &&

# # Alternative: Cell barcode correction with knee distance
# # "INFO total number of distinct corrected barcodes : 23,247"
# alevin-fry generate-permit-list \
# --knee-distance \
# -d fw \
# -i salmon_alevin \
# -o alevin_fry_gpl_k &&

# Filter mapping information
alevin-fry collate \
-i alevin_fry_gpl_f \
-r salmon_alevin \
-t 8 &&

# UMI resolution + quantification
alevin-fry quant \
-r cr-like \
-m $(ls splici_rl98_ref/*3col.tsv) \
-i alevin_fry_gpl_f \
-o alevin_fry_quant \
-t 8 &&

source $conda_base/bin/deactivate


# # (Optional) Convert matrix to 10x
# source $conda_base/bin/activate pyroe

# alevinFryMtxTo10x.py alevin_fry_quant alevin_fry_parsed single_cell

# source $conda_base/bin/deactivate


# # (Optional) Remove empty droplets
# source $conda_base/bin/activate dropletutils

# dropletutils-read-10x-counts.R -s alevin_fry_parsed -c TRUE -o alevin_fry_parsed/matrix.rds

# mkdir -p empty_drops

# dropletutils-empty-drops.R -i alevin_fry_parsed/matrix.rds --lower 5 --niters 1000 --filter-empty TRUE --filter-fdr 0.01 -o empty_drops/nonempty.rds -t empty_drops/nonempty.txt

# source $conda_base/bin/deactivate


# # (Optional) Convert 10X data to anndata format
# source $conda_base/bin/activate sc_py_training

# Rscript -e 'library(sceasy)' \
# -e 'sce <- readRDS("empty_drops/nonempty.rds")' \
# -e 'sceasy::convertFormat(sce, outFile="empty_drops/nonempty.h5ad", from="sce", to="anndata", main_layer="counts")' \
# -e 'print(sce)'
