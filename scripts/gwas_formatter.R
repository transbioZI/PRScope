rm(list=ls())

format_numeric <- function(df, column) {
  df[,column] = gsub("[^0-9| ^e | ^E | ^.| ^-]", "", df[,column])
  df = df[df[,column] != "[e | E]",] 
  df = df[df[,column] != ".",] 
  df = df[df[,column] != "-",]
  df = df[df[,column] != "",]
  df = df[!is.na(df[,column]),]
  df = df[!is.na(as.numeric(df[,column])),]
  return(df)
}

args = commandArgs(trailingOnly=TRUE)
gwas = "/data/PGC_GWASs/raw/CUD_2020/CUD_EUR_casecontrol_public_11.14.2020.tsv" #args[1]
output = "/data/PGC_GWASs/formatted/CUD_2020/CUD_EUR_casecontrol_public_11.14.2020_formatted_37.tsv"

##Format
ss = read.table(gwas,header = T, sep = " ")
ss$effect_allele_frequency = NA
dim(ss)
ss[,"chromosome"] = gsub("[chr | CHR | Chr | chromosome | chrom]","",ss[,"chromosome"])
ss = format_numeric(ss,"chromosome")
dim(ss)
ss = format_numeric(ss,"base_pair_location")
dim(ss)
ss = format_numeric(ss,"info")
ss[ss$info > 1.0,]$info = 1.0
dim(ss)
ss$info = NA
if("odds_ratio" %in% colnames(ss)) {
  ss <- format_numeric(ss,"odds_ratio")
}
if("beta" %in% colnames(ss)) {
  ss <- format_numeric(ss,"beta")
}

ss = ss[grep('^rs', ss$rsid),]

ss = ss[with(ss, order(as.numeric(chromosome), as.numeric(base_pair_location))), ]

columns = c('chromosome', 'base_pair_location', 'effect_allele', 'other_allele', 'beta', 'standard_error', 
            'effect_allele_frequency', 'p_value','info','rsid')

write.table(ss[,columns], output, col.names = T, row.names = F, quote = F, sep = "\t")
