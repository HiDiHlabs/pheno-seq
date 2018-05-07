# pheno-seq

This repository contains KNIME workflows and R-code required for ‘pheno-seq’. pheno-seq is a method for combined imaging and gene expression profiling of single cell derived spheroids from three-dimensional cell culture systems to understand functional tumor cell heterogeneity. The pheno-seq manuscript is available as pre-print: 
*	https://www.biorxiv.org/content/early/2018/05/01/311472

The methodology is based on the TakaraBio iCELL8™ technology that uses an automated dispenser and barcoded nanowells and was initially developed for single cell genomics applications. 

This repository contains:
*	Detailed pheno-seq protocol for microscopy and image analysis
*	KNIME image analysis workflow for automated processing of nanowell microscopy images (PhenoSeq_preprocessing_SP8.knwf)
*	Shiny web app R code for interactive analysis and selection of spheroids for sequencing (PhenoSelect.R)
*	Additional KNIME image analysis workflows for functional validation assays
*	R code for pheno-seq RNA-seq analysis


#### Example dataset

An example image dataset for testing the KNIME pre-processing workflow (raw) and a folder with already processed images & additional files for the PhenoSelect shiny web app (processed) can be downloaded [here](https://drive.google.com/drive/folders/1opFAoC7HZ2u5pSblKgjz3bDiD1n4SJ0M?usp=sharing).
