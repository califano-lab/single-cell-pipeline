#!/bin/bash
#$ -l mem=20G,time=12::

Rscript /ifs/scratch/c2b2/ac_lab/CZI/lv_pipeline/code/regProcess.R $1 $2 $3 $4
