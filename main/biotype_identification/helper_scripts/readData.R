DiscoveryData = read.table("./results/PRSs/discovery_prs_500.tsv", sep=" ", header=TRUE)  # 500 higher gives good results
DiscoveryMeta = read.table("./results/meta_data/discovery.fam", sep=" ", header=F)
DiscoveryPCs  = read.table("./results/meta_data/discovery.POP_STRATIFICATION.eigenvec", sep=" ", header=FALSE)
DiscoveryPCs  = DiscoveryPCs[,-1]
DiscoveryPCs = DiscoveryPCs[,c(1:21)]
print("Outlier removing using pop stratification PCs")
homogenous = log(dist_ogk(as.matrix(DiscoveryPCs[,c(-1)])))
hist(homogenous)
indices = which(homogenous >= 4.5)
if(length(indices) != 0) {
  print(paste0("outlier removed : ", length(indices)))
  DiscoveryPCs  = DiscoveryPCs[-indices,]
}

ValidationData = read.table("./results/PRSs/validation_prs_500.tsv", sep=" ", header=TRUE)
ValidationMeta = read.table("./results/meta_data/validation.fam", sep=" ", header=F)
ValidationPCs  = read.table("./results/meta_data/validation.POP_STRATIFICATION.eigenvec", sep=" ", header=FALSE)
ValidationMeta = ValidationMeta[,-1]
ValidationPCs  = ValidationPCs[,-1]
ValidationPCs = ValidationPCs[,c(1:21)]

homogenous = log(dist_ogk(as.matrix(ValidationPCs[,c(-1)])))
hist(homogenous)
indices = which(homogenous >= 4.5 )
if(length(indices) != 0) {
  print(paste0("outlier removed : ", length(indices)))
  ValidationPCs  = ValidationPCs[-indices,]
}

DiscoveryInter = intersect(intersect(DiscoveryData[,1],DiscoveryMeta[,2]),DiscoveryPCs[,1])
ValidationInter = intersect(intersect(ValidationData[,1],ValidationMeta[,1]),ValidationPCs[,1])

DiscoveryData = DiscoveryData[match(DiscoveryInter,DiscoveryData[,1]),]
DiscoveryMeta = DiscoveryMeta[match(DiscoveryInter,DiscoveryMeta[,2]),]
DiscoveryPCs = DiscoveryPCs[match(DiscoveryInter,DiscoveryPCs[,1]),]

stopifnot(identical(DiscoveryData[,1],DiscoveryMeta[,2]))
stopifnot(identical(DiscoveryPCs[,1],DiscoveryMeta[,2]))

ValidationData = ValidationData[match(ValidationInter,ValidationData[,1]),]
ValidationMeta = ValidationMeta[match(ValidationInter,ValidationMeta[,1]),]
ValidationPCs = ValidationPCs[match(ValidationInter,ValidationPCs[,1]),]

stopifnot(identical(ValidationData[,1],ValidationMeta[,1]))
stopifnot(identical(ValidationPCs[,1],ValidationMeta[,1]))
## match datasets

DiscoveryGroup = DiscoveryMeta[match(DiscoveryData[,1],DiscoveryMeta[,2]),6]
ValidationGroup  = ValidationMeta[match(ValidationData[,1],ValidationMeta[,1]),5]

DiscoverySex = DiscoveryMeta[match(DiscoveryData[,1],DiscoveryMeta[,2]),5]
ValidationSex = ValidationMeta[match(ValidationData[,1],ValidationMeta[,1]),4]

table(DiscoveryGroup)
table(ValidationGroup)

## residualization

DiscoveryRES = apply(DiscoveryData[,-1],2,function(x){ lm( x ~ as.matrix(DiscoveryPCs[,-1]))$res })
ValidationRES  = apply(ValidationData[,-1],2,function(x){ lm( x ~ as.matrix(ValidationPCs[,-1]))$res })

inter = intersect(colnames(DiscoveryRES),colnames(ValidationRES))

DiscoveryRES = DiscoveryRES[,inter]
ValidationRES = ValidationRES[,inter]

######## scale data

DiscoveryRES = scale(DiscoveryRES)
ValidationRES  = scale(ValidationRES)

print(dim(DiscoveryRES))
print(dim(ValidationRES))

print(table(DiscoveryGroup))
print(table(ValidationGroup))
