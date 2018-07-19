args <- commandArgs(TRUE)
wd <- args[1]#working directory
acro <- args[2]#acronym
reg <- args[3]#regulator name
rfn <- args[4]#regulator file
#arguments

ad <- "/ifs/scratch/c2b2/ac_lab/CZI/Pipeline/Modules/ARACNe/"
rscript="/nfs/apps/R/3.3.1/bin/Rscript --vanilla"
java="/nfs/apps/java/1.7.0_25/bin/java"
#pre-defined values

fn <- list.files(pattern=paste("ar_", acro, "-", reg, sep=""))
tmp <- sapply(fn, function(x) {
    tmp <- readLines(x)
    tmp[length(tmp)]
})
names(tmp) <- sapply(strsplit(fn, "_"), function(x) sub(".log", "", x[3]))
tmp <- as.numeric(names(tmp))[-grep("Total time elapsed", tmp)]
if (length(tmp)>0) {
    for (i in tmp) {
        qlog <- paste(wd, "/ar_", acro, "-", reg, "_", i, ".log", sep="")
        command <- paste(java, " -Xmx16000M -jar ", ad, "/aracne.jar", " -e ./", acro, "-expmat.dat -o ./", acro, "-", reg, "/ --tfs ", rfn, " --pvalue 0.00000001 --threads 4 --seed ", i, sep="")
        system(paste("echo \"", command, "\" | qsub -l mem=20G,time=1:: -pe smp 4 -N ar_", acro, "-", reg, " -j yes -o ", qlog, " -cwd", sep=""))
    }
    qlog <- paste(wd, "/chk_", acro, "-", reg, ".log", sep="")
    command <- paste(rscript, " checkjobs.r", wd, acro, reg, rfn, sep=" ")
    system(paste("echo \"", command, "\" | qsub -l mem=1G,time=1:: -N chk_", acro, "-", reg, " -j yes -o ", qlog, " -cwd -hold_jid ar_", acro, "-", reg, sep=""))
}
#re-submit aborted bootstraps
