#args[1] : workdir
#args[2] : prefix


args = commandArgs(trailingOnly=TRUE)

if (length(args)!=2) {
  stop("[workdir] [prefix]", call.=FALSE)
}

# Read in file
valid <- read.table(paste(args[1],paste(args[2],'valid.sample',sep='.'), sep='/'), header=T, colClasses = c("character", "character"),sep="\t")
dat <- read.table(paste(args[1],paste(args[2],'QC.sexcheck',sep='.'), sep='/'), header=T, colClasses = c("character", "character", "integer","integer","character","numeric"), sep="")
valid <- subset(dat, STATUS=="OK" & FID %in% valid$FID)
write.table(valid[,c("FID", "IID")], paste(args[1],paste(args[2],'QC.valid',sep='.'), sep='/'), row.names=F, col.names=F, sep="\t", quote=F) 