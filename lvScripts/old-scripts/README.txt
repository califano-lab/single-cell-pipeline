### INTRODUCTION ###

These scripts are modularized verisons of the CZI pipeline, aiming at a better engineered pipeline. Each subfunction of the pipeline should be abstracted to its own script that can be run independently of all others.

As a general rule, this pipeline aims to shift the file format pipelien from .rda files to .rds. The former are still necessary for use with ARACNe (later versions will hopefully rectify this), but are inefficient methods of storing single matrices, as they often are.

### PACKAGE INSTALLATION ###

The following packages are necessary at some point thorughout the pipeline:
 - stringr
 - optparse
 - viper

What follows is an installation guide for these packages, if you're missing any of them:
install.packages('optparse')
install.packages('stringr')
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("viper", version = "3.8")

### PRE-PROCESSING ###



### QUALITY CONTROL ###


bash /ifs/scratch/c2b2/ac_lab/atg2142/Fall_2018/QSUB_METAVIPER_SCRIPTS/atg_metaviper_qsamples=largeTissue_raw_sampled.rda --ref_samples=largeTissue_raw_sampled.rda --regulon_list=mrc_gts.rda --output_dir=lv_full --core_num=4 --per_num=1000 --cleanup=TRUE