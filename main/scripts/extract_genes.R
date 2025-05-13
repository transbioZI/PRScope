library(TxDb.Hsapiens.UCSC.hg38.knownGene)

dbGENES = genes(TxDb.Hsapiens.UCSC.hg38.knownGene, single.strand.genes.only = F) 
tableGenes = as.data.frame(dbGENES)

library('org.Hs.eg.db')
geneSymbols <- mapIds(org.Hs.eg.db, tableGenes$group_name, 'SYMBOL', 'ENTREZID')

tableGenes = (cbind(tableGenes, geneSymbols))
tableGenes = tableGenes[,c(2,3,4,5,7,8)]
tableGenes$seqnames = gsub("chr", "", tableGenes$seqnames)
a = tableGenes[!duplicated(tableGenes$group_name),]
head(a)
head(tableGenes)

write.table(a, "/data/references/MAGMA/UCSC_TxDb.gene.loc", quote = F, row.names = F, col.names = F, sep = "\t")

b = fread("/zi/flstorage/HITKIP/share/Ersoy/magma/UCSC_TxDb.gene.loc")
