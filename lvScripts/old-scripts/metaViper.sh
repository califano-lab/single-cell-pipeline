## Author = Aaron T. Griffin
## Date = 12/03/18
## Group = Califano Lab

## Objective = perform metaVIPER analysis on the single cell RNA-Seq data set from Lukas and Alec using the GTEx interactomes from Alessandro

# set running parameters
current_test_samples=$1
current_reference_samples=$2
current_regulon_list=$3
current_permutation_number=0
current_output_directory=$4
current_normalization_method=rank
current_combination_method=Stouffer
current_keep_settings=vpmat

# run qsub metaVIPER

bash /ifs/scratch/c2b2/ac_lab/atg2142/Fall_2018/QSUB_METAVIPER_SCRIPTS/atg_qsub_metaviper_wrapper.sh --test_samples=${current_test_samples} --ref_samples=${current_reference_samples} --regulon_list=${current_regulon_list} --per_num=${current_permutation_number} --norm_method=${current_normalization_method} --combine_method=${current_combination_method} --output_dir=${current_output_directory} --keep=${current_keep_settings} >> ${current_output_directory}/qsub_metaviper_logfile.txt

