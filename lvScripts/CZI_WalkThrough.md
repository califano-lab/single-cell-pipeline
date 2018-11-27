### PACKAGES ###

To get started, make sure you have the following R packages installed:
stringr
optparse
viper

You'll also need to install the following python libraries (using python 3):
numpy
pandas
scanpy
igraph


### INTRO ###

This walkthrough will use the following two files from Peter Sims's CZI data as a starting point:
PP001swap.filtered.matrix.txt.bz2
PP002swap.filtered.matrix.txt.bz2

Note that when these two files are specified as arguments to command line functions, you'll need to provide the full path to the files. The same goes for all files you specify as inputs. 

To start, make a directory for the outputs of the analysis (which will be referred to as 'YOUR-OUT-DIRECTORY'). You'll also need a directory with the GTEx bulk interactomes (referred to as 'GTEX-DIRECTORY') as well as a conversion dictionary for different gene name conventions (referred to as CONVERT-DICT).


### STEP 0: cziMerge.R ###

The CZI data arrives in two files; one for acting and one for resting. We'll merge these two files, giving the cells appropriate names, and create a unified data frame for all donstream analysis. This is a 'Step 0' because it is unique to the CZI data set. 

Run the following command:

```Rscript cziMerge.R --rest_file=PP001.swap.filtered.matrix.txt.bz2 --act_file=PP002.swap.filtered.matrix.txt.bz2 --out_name=d1-lung --out_dir=YOUR-OUT-DIRECTORY```


### STEP 1: preProcess.R ###

Next, we will take the raw data and perform some pre-processing steps. We'll remove cells with either too few or too many counts as well as any genes with no reads. We'll also save a CPM normalized copy of the data.

Run the following command:

```Rscript preProcess.R --raw_file=YOUR-OUT-DIRECTORY/d1-lung_mergedRaw.rds --out_name=d1-lung --out_dir=YOUR-OUT-DIRECTORY --min_count=1000 --max_count=100000```


### STEP 2: metaViper.R ###

We'll now compute protein activity of the data using metaVIPER and the GTEx bulk interactomes (see Intro). This will give us protein activity that we can perform our initial clustering with.

Run the following command:

```Rscript gtex-metaViper.R --cpm_file=d1-lung_mergedCPM.rds --net_dir=GTEX-DIRECTORY --convert_dict=CONVERT-DICT --out_name=d1-lung --out_dir=YOUR-OUT-DIRECTORY```


### STEP 2.5: scanpy-prep.R ###

Scanpy needs an anontation vector (indication source file of each cell) and cannot take RDS. Thus, a preprocessing step to create a scanpy friendly text file and an annotation vector is necessary. NOTE: I know this is ineffeicient and terrible, we're working on making it unnecessary.

Run the following command:

```Rscript scanpy-prep.R --in_file=d1-lung_GTExActivity.rds --out_name=d1-lung_GTExActivity --out_dir=YOUR-OUT-DIRECTORY```


### STEP 3: scanpy-clust.py ###

Now, we'll cluster the protein activity using scanpy. This will use seurat's dispersion metrics to select the most variable proteins, then perform Louvain clustering on the cells and generate a text file with cluster labels for each cell.

Run the following command:

```python3 scanpy-clust.py d1-lung_GTExActivity_scanpy-inDat.txt d1-lung_GTExActivity_scanpy-inDat.txt d1-lung_GTEx YOUR-OUT-DIRECTORY```


### STEP 4: clustSubset.R ###

With our cluster labels in hand from scanpy, we will generate raw count files for each cluster that passes a specified size threshold.

Run the following command:

```Rscript clustSubset.R --input_file=d1-lung_mergedFiltered.rds --cluster_labels=d1-lung_GTEx_pAct-louvainClust.txt --sample_index=1 --cluster_index=3 --out_name=d1-lung_GTEx-clust --out_dir=YOUR-OUT-DIRECTORY```


## STEP 5: makeMetaCells.R ###

Now, we generate meta cells for each cluster subset. For each subset from the previous step (where the cluster number replaces #), we create meta cells that can be run with ARACNe. In addition to the raw count matrix (from which meta cells will be computed and normalized), this script requires a matrix of protein activity, which is used to calculate inter-cell distance useing viper similarity.

Run the followign command (for each cluster):

```Rscript makeMetaCells.R --input_file=d1-lung_GTEx-clust_cluster-#.rds --out_dir=YOUR-OUT-DIRECTORY```


## STEP 6: run ARACNe ###

Finally, we run ARACNe on the meta cells from each cluster. Use the subClustAracne.sh script on each of the meta cell files from the previous file