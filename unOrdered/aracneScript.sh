aracne='/ifs/scratch/c2b2/ac_lab/CZI/single-cell-pipeline/Modules/ARACNe/aracne.sh'
workDir=$1
regFile=$2
expFile=$3
jobName=$4
regName=$5

sh $aracne $workDir $jobName $expFile $regName $regFile
