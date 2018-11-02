#!/bin/bash
#$ -wd /ifs/scratch/c2b2/ac_lab/CZI/lv_pipeline/logs/
#$ -l mem=20G,time=4::

# README: 
# this script will take in two files (one activated and one resting) and a name
# the data will be merged and preprocessed
# ARACNe networks are computed for all four sets of regulons (tfs, cotfs, signaling, surface)
# those networks are used to perform VIPER, and the data is concatenated into one matrix, then subset for lineage markers
# the data is then clustered on those lineage markers
# users should look at the clustering results before choosing an optimal one and proceeding to the phase 2 script

# parameters
dir='/ifs/scratch/c2b2/ac_lab/CZI/lv_pipeline/results'
regDir='/ifs/scratch/c2b2/ac_lab/CZI/lv_pipeline/regSets'
logDir='/ifs/scratch/c2b2/ac_lab/CZI/lv_pipeline/logs/'
# code paths
aracne='/ifs/scratch/c2b2/ac_lab/CZI/lv_pipeline/code/aracneScript.sh'
regConvert='/ifs/scratch/c2b2/ac_lab/CZI/lv_pipeline/code/regProcess.sh'
# arguments
infile=$1
name=$2

mkdir ${dir}/${name}

# ARACNe
mkdir ${dir}/${name}/'fullSampleARACNe/'
regs='cotf sig surface tf'
expFile=$infile
jobNames=()
for i in $regs; do
	echo $i' ARACNe run'
	workDir=${dir}/${name}/'fullSampleARACNe/'$i'ARACNe/'
	mkdir $workDir
	regFile=${regDir}/$i'-ensembl.txt'
	sh $aracne $workDir $regFile $expFile $name $i
	jName='cl_'${name}'-'$i
	jName=${jName:0:10}
	jobNames+=($jName) 	
done
echo ${jobNames[*]}
while qstat -u $USER | grep -q -e ${jobNames[0]} -e ${jobNames[1]} -e ${jobNames[2]} -e ${jobNames[3]} ; do
	sleep 1
done
echo Completed ARACNe

# merge tsvs and make unique
cat ${dir}/${name}/'fullSampleARACNe/cotfARACNe/'${name}'-cotf-network.tsv' ${dir}/${name}/'fullSampleARACNe/sigARACNe/'${name}'-sig-network.tsv' ${dir}/${name}/'fullSampleARACNe/surfaceARACNe/'${name}'-surface-network.tsv' ${dir}/${name}/'fullSampleARACNe/tfARACNe/'${name}'-tf-network.tsv' > ${dir}/${name}/'fullSampleARACNe/merged-network.tsv'
sort -u ${dir}/${name}/'fullSampleARACNe/merged-network.tsv' > ${dir}/${name}/'fullSampleARACNe/merged-uniq-network.tsv'

# aracne2regulon
netFile=${dir}/${name}/'fullSampleARACNe/merged-uniq-network.tsv'
outFile=${dir}/${name}/'fullSampleARACNe/'${name}'_mergedNet.rds'
qsub -wd $logDir $regConvert $netFile $expFile $outFile
