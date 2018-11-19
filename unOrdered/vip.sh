#!/bin/bash
#$ -l mem=20G,time=3::

Rscript /ifs/scratch/c2b2/ac_lab/CZI/lv_pipeline/code/vip.R $1 $2 $3
