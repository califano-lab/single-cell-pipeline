##### Arguments
wdir=$1
#working directory
acro=$2
#acronym
ifn=$3
#expression matrix
reg=$4
#regulator name
rfn=$5
#regulator file

##### Pre-defined values
bst=200
adir="/ifs/scratch/c2b2/ac_lab/CZI/single-cell-pipeline/Modules/ARACNe/"
rscript="/nfs/apps/R/3.3.1/bin/Rscript"
java="/nfs/apps/java/1.7.0_25/bin/java"

##### Expression file preparation
command="$rscript $adir/expmat.r ${wdir} ${ifn} ${acro} 0.1 30000"
echo $command | qsub -l mem=32G,time=1:: -N dat_${acro}-${reg} -j yes -cwd -o ${wdir}/dat_${acro}-${reg}.log

##### Threshold
command="$java -Xmx16000M -jar $adir/aracne.jar -e ${wdir}/${acro}-expmat.dat -o ${wdir}/${acro}-${reg}/ --tfs ${rfn} --pvalue 0.00000001 --seed 1 --calculateThreshold"
echo $command | qsub -l mem=20G,time=1:: -N thr_${acro}-${reg} -j yes -o ${wdir}/thr_${acro}-${reg}.log -cwd -hold_jid dat_${acro}-${reg}

##### Run nobootstrap and noDPI
command="mkdir ${wdir}/${acro}-${reg}_all\n
cp ${wdir}/${acro}-${reg}/miThreshold*.txt ${wdir}/${acro}-${reg}_all/\n
rename p1E-8 p1E0 ${wdir}/${acro}-${reg}_all/miThreshold*.txt\n
echo "0" > ${wdir}/${acro}-${reg}_all/miThreshold*.txt\n
$java -Xmx36000M -jar $adir/aracne.jar -e ${wdir}/${acro}-expmat.dat -o ${wdir}/${acro}-${reg}_all/ --tfs ${rfn} --pvalue 1 --threads 4 --nodpi --nobootstrap"
echo -e $command | qsub -l mem=40G,time=2:: -N ar_${acro}-${reg} -j yes -o ${wdir}/ar_${acro}-${reg}_all.log -cwd -hold_jid thr_${acro}-${reg}

##### Run bootstraps
for i in $(seq 1 $bst); do
command="$java -Xmx16000M -jar $adir/aracne.jar -e ${wdir}/${acro}-expmat.dat -o ${wdir}/${acro}-${reg}/ --tfs ${rfn} --pvalue 0.00000001 --threads 8 --seed ${i}"
echo $command | qsub -l mem=20G,time=1:: -N ar_${acro}-${reg} -j yes -o ${wdir}/ar_${acro}-${reg}_${i}.log -cwd -hold_jid thr_${acro}-${reg} 
done

##### Check all jobs finished
command="$rscript $adir/checkjobs.r ${wdir} ${acro} ${reg} ${rfn} $adir $rscript $java"
echo $command | qsub -l mem=1G,time=1:: -N chk_${acro}-${reg} -j yes -o ${wdir}/chk_${acro}-${reg}.log -cwd -hold_jid ar_${acro}-${reg}

##### Consolidate
command="$rscript $adir/consolidate.r ${wdir}/${acro}-${reg} ${wdir}/${acro}-expmat.dat ${rfn} none 0.0001 $bst"
echo $command | qsub -l mem=32G,time=4:: -N net_${acro}-${reg} -j yes -o ${wdir}/net_${acro}-${reg}.log -cwd -hold_jid ar_${acro}-${reg}

##### Clean-up
command="cp ${wdir}/${acro}-${reg}/finalNetwork_4col.tsv ${wdir}/${acro}-${reg}_4col.tsv\n
cp ${wdir}/${acro}-${reg}_all/network.txt ${wdir}/${acro}-${reg}-network.tsv\n
$rscript $adir/regulon.r ${wdir} ${acro} ${reg} ${wdir}/${acro}-${reg}-regulon.rda"
echo -e $command | qsub -l mem=20G,time=8:: -N cl_${acro}-${reg} -j yes -o ${wdir}/cl_${acro}-${reg}.log -cwd -hold_jid net_${acro}-${reg}
