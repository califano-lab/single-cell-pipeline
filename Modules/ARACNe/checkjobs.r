args <- commandArgs(TRUE)
wd <- args[1]#working directory
acro <- args[2]#acronym
reg <- args[3]#regulator name
rfn <- args[4]#regulator file
ad <- args[5]#aracne directory
rscript <- args[6]#rscript
java <- args[7]#java
#arguments

fn <- list.files(pattern=paste("ar_", acro, "-", reg, sep=""))
tmp <- sapply(fn, function(x) {
    tmp <- readLines(x)
    tmp[length(tmp)]
})
names(tmp) <- sapply(strsplit(fn, "_"), function(x) sub(".log", "", x[length(x)]))
names(tmp)[which(names(tmp) == "all")] <- 0#ar_acro-reg_all.log
tmp <- as.numeric(names(tmp))[-grep("Total time elapsed", tmp)]
if (length(tmp)>0) {
    for (i in tmp) {
        if (i == 0){
            message("Re-submitting no-bootstrap no-dpi network")
            qlog <- paste(wd, "/ar_", acro, "-", reg, "_all.log", sep="")
            command <- paste(java, " -Xmx36000M -jar ", ad, "/aracne.jar", " -e ", wd, "/", acro, "-expmat.dat -o ", wd, "/", acro, "-", reg, "_all/ --tfs ", rfn, " --pvalue 1 --threads 4 --nodpi --nobootstrap", sep="")
            system(paste("echo \"", command, "\" | qsub -l mem=40G,time=2:: -N ar_", acro, "-", reg, " -j yes -o ", qlog, " -cwd", sep=""))
            
        }else{
            message(paste("Re-submitting bootstrap network #", i, sep = ""))
            qlog <- paste(wd, "/ar_", acro, "-", reg, "_", i, ".log", sep="")
            command <- paste(java, " -Xmx16000M -jar ", ad, "/aracne.jar", " -e ", wd, "/", acro, "-expmat.dat -o ", wd, "/", acro, "-", reg, "/ --tfs ", rfn, " --pvalue 0.00000001 --threads 4 --seed ", i, sep="")
            system(paste("echo \"", command, "\" | qsub -l mem=20G,time=1:: -N ar_", acro, "-", reg, " -j yes -o ", qlog, " -cwd", sep=""))
        }
    }
    qlog <- paste(wd, "/chk_", acro, "-", reg, ".log", sep="")
    command <- paste(rscript, " ", ad, "checkjobs.r ", wd, " ", acro, " ", reg, " ", rfn, " ", ad, " ", rscript, " ", java, " ", sep="")
    system(paste("echo \"", command, "\" | qsub -l mem=1G,time=1:: -N chk_", acro, "-", reg, " -j yes -o ", qlog, " -cwd -hold_jid ar_", acro, "-", reg, sep=""))
}
message("Finished!")
#re-submit aborted bootstraps
