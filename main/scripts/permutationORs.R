rm(list=ls())
gc()
##################
##### load libraries

library(MASS)
library(igraph)
library(dplyr)
library("FactoMineR")
library("factoextra")
library(liver)
library(stringr)
suppressMessages(library(dplyr))
suppressMessages(library(tidyverse))
suppressMessages(library(gwasrapidd))
suppressMessages(library(rtracklayer))
suppressMessages(library(magrittr))
suppressMessages(library(purrr))
library(stringr)

path = "/data/Output/"
output_name = "permutationORs_51_0_60_cov"
repet = 51
threads = 26
q3 = 2
corr_value = 0.60
prs_path <- "/data/Output/IMAGEN_RANDOM_STUDIES"
files = list.files(prs_path, pattern="*.all_score")
files = sapply(files, function(x) {
  str_split(x,pattern = "\\.")[[1]][1]
})
ll <- c()
studies_to_remove = c()
for(f in files[1:length(files)]) {
  ff= sprintf("%s/%s.prsice", prs_path, f)
  if(file.exists(ff)) {
    df2 <- read.table(ff, header = T)
    
    tt <- c()
    tt <- c(tt, df2[1,4])
    tt <- c(tt, df2[2,4])
    tt <- c(tt, df2[3,4])
    tt <- c(tt, df2[4,4])
    tt <- c(tt, df2[5,4])
    tt <- c(tt, df2[6,4])
    tt <- c(tt, df2[7,4])
    tt <- c(tt, df2[8,4])
    tt <- c(tt, df2[9,4])
    tt <- c(tt, df2[10,4])
    ll<- c(ll, max(tt, na.rm = T))
    if(max(tt, na.rm = T) < 1000) {
      studies_to_remove = c(studies_to_remove,str_split(f,"\\.")[[1]][1])
    }
  } else {
    print(f)
    studies_to_remove = c(studies_to_remove,f)
  }
}

files_exist = c()
for(f in files) {
  if(!(f %in% studies_to_remove)) {
    files_exist = c(files_exist, f)
  } else {
    file.remove(paste0("/data/transformed_random/",f,".transformed.h.tsv"))
  }
}

writeLines(files_exist,"/zi/home/ersoy.kocak/Desktop/Projects/PRSCalculator/resources/final_random_study_list.txt")
for(ze in c(1:100)) {
files_exist_r = files_exist[sample(c(1:length(files_exist)), 276)]

df <- read.table(sprintf("%s/%s.all_score", prs_path, files_exist_r[1]), header = T)
colnames(df) <- c("FID","IID",paste0(colnames(df)[3:ncol(df)], sprintf(".%s",strsplit(files_exist_r[1], '[.]')[[1]][1])))
for(f in files_exist_r[2:length(files_exist_r)]) {
  df2 <- read.table(sprintf("%s/%s.all_score", prs_path, f), header = T)
  colnames(df2) <- c("FID","IID",paste0(colnames(df2)[3:ncol(df2)], sprintf(".%s",strsplit(f, '[.]')[[1]][1])))
  df <- merge(df,df2, by = c("FID","IID"))
}

Sys.setenv(OMP_NUM_THREADS=1)
start_time <- Sys.time()

##################
##### read data
PRS_data = df[,-1]
PC_data = read.table("/zi/flstorage/group_transbio/projects/imagen_analysis/Results/PCA/all.QCED.POP_STRATIFICATION.eigenvec", sep=" ")
meta_data = read.table("/zi/flstorage/group_transbio/projects/imagen_analysis/Results/meta_data/IMAGEN_demographics.csv", sep=",", header=TRUE)
AUDIT = read.table("/zi/flstorage/group_transbio/projects/imagen_analysis/Results/meta_data/IMAGEN-IMGN_AUDIT_FU3_copyES.csv", sep=",", header=TRUE)
ESPAD = read.csv("/zi/flstorage/group_transbio/projects/imagen_analysis/Results/meta_data/IMAGEN-IMGN_ESPAD_FU3.csv", header=TRUE)

##################
##### align data, metadata, PRS scores, AUD rating scales

IDinter   = intersect(intersect(intersect(PRS_data[,1], meta_data[,1]), AUDIT[,1]),ESPAD[,1])
PRS_data  = PRS_data[match(IDinter, PRS_data[,1]),]
meta_data = meta_data[match(IDinter, meta_data[,1]),]
PC_data   = PC_data[match(IDinter, PC_data[,2]),]
AUDIT     = AUDIT[match(IDinter, AUDIT[,1]),]
ESPAD     = ESPAD[match(IDinter, ESPAD[,1]),]

stopifnot(identical(PRS_data[,1],meta_data[,1]))
stopifnot(identical(PRS_data[,1],AUDIT[,1]))

sexvec = rep(0, nrow(PRS_data))
sexvec[which(meta_data$sex == "F")] = 1

PRS_data_RES = apply(PRS_data[,-1],2,function(x){ lm(x ~
                                                       meta_data$recruitment.centre
                                                     + PC_data[,3]
                                                     + PC_data[,4]
                                                     + PC_data[,5]
                                                     + PC_data[,6]
                                                     + PC_data[,7]
                                                     + PC_data[,8]
                                                     + PC_data[,9]
                                                     + PC_data[,10]
                                                     + PC_data[,11]
                                                     + PC_data[,12]
                                                     + PC_data[,13]
                                                     + PC_data[,14]
                                                     + PC_data[,15]
                                                     + PC_data[,16]
                                                     + PC_data[,17]
                                                     + PC_data[,18]
                                                     + PC_data[,19]
                                                     + PC_data[,20]
                                                     + PC_data[,21]
                                                     + PC_data[,22]
                                                     
)$res })

##################
##### Combine correlated PRS features into first principal component
##### The code below extract clusters of correlated genes (not only gene pairs)

dataCOR = cor(PRS_data_RES)		## linear correlation

dataCOR[dataCOR<corr_value] = 0		## use cut-off of 0.7
dataCOR[lower.tri(dataCOR)] = 0
diag(dataCOR) = 0
dataCOR[dataCOR!=0] = 1

data_graph = graph_from_adjacency_matrix(dataCOR, mode="undirected")

components <- decompose(data_graph, min.vertices=5)

max(clusters(data_graph)$csize)
hist(clusters(data_graph)$csize, breaks = 50)
data_clusters = clusters(data_graph)$membership

newDATA = numeric()			## combine correlated clusters into first principal component

for(i in unique(data_clusters)){
  curclu = which(data_clusters == i)
  if(length(curclu) > 1){
    pcmod = prcomp(PRS_data_RES[,curclu], scale=TRUE)$x
    px = pcmod[,1:min(c(1,ncol(pcmod) ))]
    if( mean(cor(px,PRS_data_RES[,curclu]) ) < 0) px = px * (-1)
    newDATA = cbind(newDATA, px)
  }else{
    newDATA = cbind(newDATA, PRS_data_RES[,curclu])
  }
}

colnames(newDATA) =c(1:ncol(newDATA))
dim(newDATA)


##################
##### define outcome

AUDIText = cbind(AUDIT$audit1,AUDIT$audit2,AUDIT$audit3,AUDIT$audit4,AUDIT$audit5,
                 AUDIT$audit6,AUDIT$audit7,AUDIT$audit8,AUDIT$audit9,AUDIT$audit10)
AUDIText[is.na(AUDIText)] = 0

curflag = rep(0, nrow(AUDIT));
curflag[intersect(which(AUDIText[,3] >= q3), which(rowSums(AUDIText)>7))] = 1
#severe_1 = intersect(which(AUDIText[,3] >= 2), which(rowSums(AUDIText)>7))
#severe_2 = intersect(which(AUDIText[,3] >= 3), which(rowSums(AUDIText)>7))
#curflag[intersect(severe_1,severe_2)] = 1
#diff = setdiff(severe_1,severe_2)
#curflag[diff] = NA

sites = meta_data$recruitment.centre
site_flag = table(sites, curflag)
site_flag

usesites = sites

####### Multivariate analysis

sexvec_conf = apply(PRS_data[,-1],2,function(x){ summary(lm(sexvec ~ x))$coefficients[,4] })
min(p.adjust(sexvec_conf, method="BH"))

sites_conf = apply(PRS_data[,-1],2,function(x){ summary(lm(as.numeric(as.factor(sites)) ~ x))$coefficients[,4] })
min(p.adjust(sites_conf, method="BH"))

print(anova(lm(curflag ~ sexvec + sites)))

##################
##### Binarize PRS data using upper 80% percentile as cut-off

## newDATA/PRS_data_RES

TRAINDATA = apply(newDATA,2,function(x){
  out = rep(0,length(x));
  out[x>=quantile(x,0.7)]=1;
  out
})
#saveRDS(TRAINDATA,paste0(path,"/TRAINDATA_",output_name,".rds"))
#TRAINDATA = apply(newDATA,2,rank)

##################
##### Define feature selection function - linear correlation

#OR2 = function(data,out2,train){
#	curcor = apply(data,2,function(x){ fisher.test(x[train], out2[train])$estimate })
#	curcor[which(curcor<1)] = 1/curcor[which(curcor<1)]
#	return(curcor)
#}

OR2 = function(data,out2,train){
  curcor = cor(data[train,], out2[train], use="pairwise.complete.obs")
  return(curcor)
}

##################
##### Define function for stratified cross-validation


runCV = function(usesites, min=TRUE){
  CV_strat = list()
  for(si in 1:length(unique(usesites)) ){
    cursub = which(is.element(usesites, unique(usesites)[si]) ==TRUE)
    checksex = table(curflag[cursub], sexvec[cursub])
    
    if(min == FALSE){
      maxpanel = max(checksex)
      s1 = sample( intersect( which(curflag[cursub] == 0), which(sexvec[cursub] == 0)), maxpanel, replace=TRUE )
      s2 = sample( intersect( which(curflag[cursub] == 1), which(sexvec[cursub] == 0)), maxpanel, replace=TRUE )
      s3 = sample( intersect( which(curflag[cursub] == 0), which(sexvec[cursub] == 1)), maxpanel, replace=TRUE )
      s4 = sample( intersect( which(curflag[cursub] == 1), which(sexvec[cursub] == 1)), maxpanel, replace=TRUE )
    }else{
      maxpanel = min(checksex)
      s1 = sample( intersect( which(curflag[cursub] == 0), which(sexvec[cursub] == 0)), maxpanel, replace=FALSE )
      s2 = sample( intersect( which(curflag[cursub] == 1), which(sexvec[cursub] == 0)), maxpanel, replace=FALSE )
      s3 = sample( intersect( which(curflag[cursub] == 0), which(sexvec[cursub] == 1)), maxpanel, replace=FALSE )
      s4 = sample( intersect( which(curflag[cursub] == 1), which(sexvec[cursub] == 1)), maxpanel, replace=FALSE )
    }
    
    cursel = cursub[c(s1,s2,s3,s4)]
    
    CV_strat[[length(CV_strat) + 1]] = cursel
  }
  return(CV_strat)
}


##################
##### Define function that runs AI pipeline

AImodel = function(TRAIN, TEST, train){
  m1 = lda( TRAIN , as.factor(curflag[train]))
  pred = predict(m1, TEST)$class
  pred
}

AIrun = function(num.feat=1, resamp=FALSE){
  
  USETRAIN = TRAINDATA
  usesites = sites
  
  if(resamp == TRUE) USETRAIN = USERAND
  
  selFeat = numeric()
  predres = numeric()
  for(cvi in 1:length(unique(usesites))){
    print(cvi)
    
    test = unlist(runCV(usesites, min=TRUE)[cvi])
    #if(length(test) < 20){
    #  print(c("leave out: ", cvi))
    #  next}
    
    predLIST = vector(mode='list', length(num.feat))
    intPred = numeric()
    for(repi in 1:repet){
      
      CV_strat_inner = runCV(usesites, min=FALSE);
      train = unlist(CV_strat_inner[-cvi])
      
      ### calculate predictions for multiple PRS
      
      varsel = OR2(data=USETRAIN, out2=curflag, train=train)
      
      selFeat = rbind(selFeat, t(varsel) ) 	### record selected features
      
      listcounter = 0
      for(curn in num.feat){
        listcounter = listcounter + 1
        curvarsel = which(abs(varsel) >= sort(abs(varsel), decreasing=TRUE)[curn])[1:curn]
        #curvarsel = which(varsel >= sort(varsel, decreasing=TRUE)[curn])[1:curn]
        #curvarsel = which(varsel <= sort(varsel)[curn])[1:curn]
        
        TRAIN = data.frame(USETRAIN[train,][,curvarsel])
        TEST = data.frame(USETRAIN[test,][,curvarsel])
        
        colnames(TEST) = colnames(TRAIN)
        
        pred = AImodel(TRAIN, TEST, train)
        
        if(length(predLIST[[listcounter]]) == 0){
          predLIST[[listcounter]] = pred
        }else{
          predLIST[[listcounter]] = cbind(predLIST[[listcounter]], pred)
        }
      }	## end different variable numbers
    }	## end repeats
    
    cvORres = lapply(predLIST, function(x){					### calculate ORs
      
      testmedout = rep(NA, nrow(x))
      testmedfreq = apply(x,1,function(x){length(which(x==2))/length(x)})
      testmedout[testmedfreq<0.5] = 0
      testmedout[testmedfreq>=0.5] = 1
      testmed = testmedout
      
      print(length(intersect(
        which(is.na(testmed)==TRUE), which(curflag[test]==1)))/
          length( which(curflag[test]==1)))
      
      if(length(unique(testmed))==1){
        curcoeff = 1
      }else{
        if(length(unique(testmed[which(is.na(testmed)==FALSE)])) == 1 || length(unique(curflag[test])) == 1){
          curcoeff = 1
        }else{
          curcoeff = fisher.test(testmed, curflag[test])$estimate
        }
      }
    })
    predres = rbind(predres, cvORres)
  }
  return( list(predres, selFeat) )
}


############################
############################
############################ RUN
### Step 1 - run non-permuted data between 0 and 30 features

#num.feat = c(1,2,3,4,5,6,7,8,9,10,15,30,50,75,100,125,150,175,200)
num.feat = c(1,3,5,10,15,30,50,75,100,125,150)
#num.feat = c(15)
#predMed = numeric()
#predStepTot = list()

library(doParallel)
cl <- makeCluster(threads)
registerDoParallel(cl)

predSteps=foreach(permi=1:repet) %dopar% {
  require(MASS)
  AIrun(num.feat=num.feat, resamp=FALSE)
}

predStepTot = list()
predMed = numeric()
for(i in 1:length(predSteps)) {
  predStep = predSteps[[i]]
  predStepTot[[i]] = predStep[[2]]
  ORextract = matrix(unlist(predStep[[1]]), ncol=length(num.feat))
  out = apply(ORextract,2,function(x){ median(x) })
  predMed = rbind(predMed, out)
}

res = apply(predMed,2,median)
names(res) = num.feat

dir.create(file.path(path, output_name), showWarnings = FALSE)

write.table(t(format(round(res,3), nsmall = 3)),paste0(path,"/",output_name,"/ORs_",output_name,"_",ze,".tsv"), quote = F,col.names = T, row.names = F, sep = "\t")
}
end_time <- Sys.time()
print(end_time - start_time)
