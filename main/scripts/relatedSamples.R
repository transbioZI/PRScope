rm(list=ls())
gc()

packages = c("dplyr", "rlang","vctrs","plinkQC")

## Now load or install&load all
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

args = commandArgs(trailingOnly=TRUE)
indir <- args[1]
name <- args[2]
fail <- args[3]
plink <- args[4]
genome_build <- args[5]
ibd <- as.numeric(args[6])
sub=sub(pattern = "(.*)\\..*$", replacement = "\\1", name)
output <- paste0(indir,"/",sub,".",fail)

res = check_relatedness(indir = indir,name,filter_high_ldregion=F,genomebuild=genome_build, imissTh = 1, path2plink = plink, run.check_relatedness = T, highIBDTh = ibd)

write.table(res$failIDs, output, quote = FALSE, row.names = FALSE, col.names = FALSE )
