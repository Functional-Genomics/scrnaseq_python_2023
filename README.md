# Single-cell RNA-seq analysis using Python 2023 
Materials for the EBI course ["Single-cell RNA-seq analysis using Python 2023"](https://www.ebi.ac.uk/training/events/single-cell-rna-seq-analysis-2023/)

Materials adapted from: [Single-cell best practices](www.sc-best-practices.org)

## Table of contents
 * [Overview](#overview)
 * [Repo contents](#repocontents)
 * [From raw reads to feature selection](#Fromrawreadstofeatureselection)
 * [Dimensionality reduction, clustering, and annotation](#Dimensionalityreduction,clustering,andannotation)
 * [Batch correction and data integration](#Batchcorrectionanddataintegration)
 * [Group Projects](#GroupProjects)

### Overview

This repositoty contains all the hands-on materials taught in the course "Single-cell RNA-seq analysis using Python 2023". Here, we learn how to analyse single-cell data starting **from raw reads until the cell-type annotation of our data**. Moreover, we explore and perform **batch correction and data integration methods** for our data. The hands-on materials include practicals(demos), exercises with answer keys, and project exercises for further hands-on practice.

### Repo contents 

- bin  : scripts necessary for preprocessing analysis
- envs : yml files, environments used for the practicals
- exercises : jupyter notebooks with exercises based on the practicals and their answers
- practicals : jupyter notebooks for each practical( demo, follow along tutorial for each session)
- projects : Coming soon

### From raw reads to feature selection

practical_1 notebook is a follow along tutorial for the preprocessing analysis of single-cell data. 
It includes the following steps : 

1. Raw Data Processing 
   * Build the reference index (pyroe and salmon)
   * Perform mapping and quantification (alevin and alevin-fry)
   * remove empty drops (bioconductor-dropletutils) 
   * Alternative method for raw data preprocessing (simpleleaf method)
2. Quality Control (QC)
   * Filtering low quality barcodes(scanpy)
   * Correction of ambient RNA(SoupX)
   * Doublet Detection(scDblFinder)
4. Normalisation methods
   * library size 10e4 and log-transform(scanpy)
   * scran-normalisation(scran)
6. Feature Selection
   * Deviance (scry)
   * Variance(scanpy)

### Dimensionality reduction, clustering, and annotation

### Batch correction and data integration

### Group Projects 














