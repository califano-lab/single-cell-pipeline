#!/bin/bash
#$ -wd /ifs/scratch/c2b2/ac_lab/CZI/lv_pipeline/logs/
#$ -l mem=30G,time=8::

# README: 
# this script will take in two files (one activated and one resting) and a name
# the data will be merged and preprocessed
# ARACNe networks are computed for all four sets of regulons (tfs, cotfs, signaling, surface)
# those networks are used to perform VIPER, and the data is concatenated into one matrix, then subset for lineage markers
# the data is then clustered on those lineage markers
# users should look at the clustering results before choosing an optimal one and proceeding to the phase 2 script

# parameters
linMarkers='/ifs/scratch/c2b2/ac_lab/CZI/lv_pipeline/misc/lin-markers-symbol.txt'
convertDict='/ifs/scratch/c2b2/ac_lab/CZI/lv_pipeline/misc/convertDict.txt'
dir='/ifs/scratch/c2b2/ac_lab/CZI/lv_pipeline/results'
regDir='/ifs/scratch/c2b2/ac_lab/CZI/lv_pipeline/regSets'
logDir='/ifs/scratch/c2b2/ac_lab/CZI/lv_pipeline/logs/'
# code paths
preProcess='/ifs/scratch/c2b2/ac_lab/CZI/lv_pipeline/code/preProcess.R'
aracne='/ifs/scratch/c2b2/ac_lab/CZI/lv_pipeline/code/aracneScript.sh'
regConvert='/ifs/scratch/c2b2/ac_lab/CZI/lv_pipeline/code/regProcess.sh'
viper='/ifs/scratch/c2b2/ac_lab/CZI/lv_pipeline/code/vip.sh'
concat='/ifs/scratch/c2b2/ac_lab/CZI/lv_pipeline/code/vipConcat-linMark.R'
clust='/ifs/scratch/c2b2/ac_lab/CZI/lv_pipeline/code/hClust.R'
# arguments
restingFile=$1
activeFile=$2
name=$3

# preprocess
mkdir ${dir}/${name} 
mkdir ${dir}/${name}/'preProcess/'
Rscript $preProcess $activeFile $restingFile ${dir}/${name}/'preProcess/' $name

# ARACNe
mkdir ${dir}/${name}/'fullSampleARACNe/'
regs='cotf sig surface tf'
expFile=${dir}/${name}/'preProcess/'${name}'_ARACNeSample.rda'
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

# aracne2regulon
mkdir ${dir}/${name}/'fullSampleNetworks/'
jobIDs=()
pruneSize=50
for i in $regs; do
	netFile=${dir}/${name}/'fullSampleARACNe/'$i'ARACNe/'${name}'-'$i'-network.tsv'
	outFile=${dir}/${name}/'fullSampleNetworks/'${name}'_'$i'Net.rds'
	jobSub=$(qsub -wd $logDir $regConvert $netFile $expFile $outFile $pruneSize)
	arrJob=(${jobSub// / })
	jobIDs+=(${arrJob[2]})
done	
while qstat -u $USER | grep -q -e ${jobIDs[0]} -e ${jobIDs[1]} -e ${jobNames[2]} -e ${jobNames[3]} ; do
	sleep 1
done
echo Converted Regulons

# viper
mkdir ${dir}/${name}/'fullSampleViper/'
expFile=${dir}/${name}/'preProcess/'${name}'_mergedCPM.rds'
jobIDs=()
for i in $regs; do
	regFile=${dir}/${name}/'fullSampleNetworks/'${name}'_'$i'Net.rds'
	outFile=${dir}/${name}/'fullSampleViper/'${name}'_'$i'_pAct.rds'
	jobSub=$(qsub -wd $logDir $viper $expFile $regFile $outFile)
	arrJob=(${jobSub// / })
	jobIDs+=(${arrJob[2]})
done
echo ${jobIDs[*]}
while qstat -u $USER | grep -q -e ${jobIDs[0]} -e ${jobIDs[1]} -e ${jobNames[2]} -e ${jobNames[3]} ; do
	sleep 1
done
echo Done with Viper

# concatenate
pActDir=${dir}/${name}/'fullSampleViper/'
outFile=${dir}/${name}/'fullSampleViper/'${name}'_linMarkers_pAct.rds'
Rscript $concat $convertDict $pActDir $linMarkers $outFile
echo Done Concatenating

# cluster
mkdir ${dir}/${name}/'linMarkerClust/'
outPath=${dir}/${name}/'linMarkerClust/'
Rscript $clust $outFile $outPath $name






