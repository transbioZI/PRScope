cl <- makeCluster(30)
registerDoParallel(cl)

previous_features_discovery = c()

for(FEATURE_SIZE in features_sizes) {
  if (!file.exists(paste0(repository_path,"/",experiment_path,"/models"))){
    dir.create(file.path(paste0(repository_path,"/",experiment_path,"/models")))
  }
  if (!file.exists(paste0(repository_path,"/", experiment_path,"/feat_",FEATURE_SIZE))){
    dir.create(file.path(repository_path,"/", experiment_path, paste0("/feat_",FEATURE_SIZE)))
  } 
}

for(FEATURE_SIZE in features_sizes) {
  for(yy in c(1:num_it)) {
    CV_DISCOVERY = createFolds(c(1:nrow(DISCOVERY_DATA)), k=5)
    CV_RES_DISCOVERY = matrix(NA, nrow(DISCOVERY_DATA), 2)
    CV_RES_DISCOVERY_1 = matrix(NA, nrow(DISCOVERY_DATA), 2)
    CV_RES_DISCOVERY_2 = matrix(NA, nrow(DISCOVERY_DATA), 2)
    
    for(cvi in 1:length(CV_DISCOVERY)) {
      train = unlist(CV_DISCOVERY[-cvi])
      test  = unlist(CV_DISCOVERY[cvi])
      
      FEATURES_FOR_CV = sort_by_correlation(DISCOVERY_DATA[train,], as.numeric(DISCOVERY_GROUP[train]))
      
      SELECTED_FEATURES_FOR_CV = as.vector(FEATURES_FOR_CV[1:FEATURE_SIZE,])
      
      SELECTED_FEATURES_FOR_CV = sort(match(SELECTED_FEATURES_FOR_CV[!is.na(SELECTED_FEATURES_FOR_CV)],colnames(DISCOVERY_DATA)))
      
      ## build clustering on training data
      ## save cluster solution to align clustering across folds
      cluster_B_CV = clustering(DISCOVERY_DATA[train,], DISCOVERY_GROUP[train], SELECTED_FEATURES_FOR_CV, FEATURE_SIZE, 100)
      
      cluster_B_CVmatch = rep(NA, nrow(DISCOVERY_DATA))
      cluster_B_CVmatch[train][which(DISCOVERY_GROUP[train] == 2)] = cluster_B_CV
      if(cvi > 1){
        if(cor.test(cluster_B_CVmatch, cluster_B_CV_RES)$estimate < 0){
          cluster_B_CV = abs(cluster_B_CV-3)
        }
      }else{
        cluster_B_CV_RES = cluster_B_CVmatch
      }
      
      #### predict for both clusters and store output
      predCV_clu1 = train_and_predict(TRAIN=DISCOVERY_DATA[train,], 
                                      VALIDATION=DISCOVERY_DATA[test,], 
                                      GROUP=DISCOVERY_GROUP[train], 
                                      CLUSTER_GROUP=cluster_B_CV, 
                                      TARGET_CLUSTER=1, 
                                      FEATURES = SELECTED_FEATURES_FOR_CV, 
                                      REPEATS=repeat_)
      
      
      predCV_clu2 = train_and_predict(TRAIN=DISCOVERY_DATA[train,], 
                                      VALIDATION=DISCOVERY_DATA[test,], 
                                      GROUP=DISCOVERY_GROUP[train], 
                                      CLUSTER_GROUP=cluster_B_CV, 
                                      TARGET_CLUSTER=2, 
                                      FEATURES = SELECTED_FEATURES_FOR_CV, 
                                      REPEATS=repeat_)
      
      CV_RES_DISCOVERY[test,1] = predCV_clu1[[3]]    # Biotype 1
      CV_RES_DISCOVERY[test,2] = predCV_clu2[[3]]    # Biotype 2
      CV_RES_DISCOVERY_1[test,1] = predCV_clu1[[1]]  # Biotype 1 probability from biotype 1 model
      CV_RES_DISCOVERY_1[test,2] = predCV_clu1[[2]]  # Biotype 2 probability from biotype 1 model
      CV_RES_DISCOVERY_2[test,1] = predCV_clu2[[1]]  # Biotype 1 probability from biotype 2 model
      CV_RES_DISCOVERY_2[test,2] = predCV_clu2[[2]]  # Biotype 2 probability from biotype 2 model
    }
    
    c_x = CV_RES_DISCOVERY[,1]
    c_x_prob = CV_RES_DISCOVERY_1[,2]
    c_y = CV_RES_DISCOVERY[,2]
    c_y_prob = CV_RES_DISCOVERY_2[,2]
    if(length(previous_features_discovery) != 0) {
      if(abs(cor(c_x_prob,previous_features_discovery)) < abs(cor(c_y_prob,previous_features_discovery))) {
        c_x = CV_RES_DISCOVERY[,2]
        c_y = CV_RES_DISCOVERY[,1]
        c_x_prob = CV_RES_DISCOVERY_2[,2]
        c_y_prob = CV_RES_DISCOVERY_1[,2]
      } 
    }
    
    previous_features_discovery = c_x_prob
    
    result = "OR\tp_value\tconf\ttype\tBiotype\n"
    x = fisher.test(table(c_x, c_y))
    result = paste0(result,format(x$estimate,nsmall = 2),"\t",format(x$p.value,nsmall = 2),"\t",format(x$conf.int[1],nsmall = 2),"-",format(x$conf.int[2],nsmall = 2),'\t',"DISCOVERY",'\t',"Inverse","\n")
    x = fisher.test(c_x, DISCOVERY_GROUP)
    result = paste0(result,format(x$estimate,nsmall = 2),"\t",format(x$p.value,nsmall = 2),"\t",format(x$conf.int[1],nsmall = 2),"-",format(x$conf.int[2],nsmall = 2),'\t',"DISCOVERY",'\t',"Biotype A","\n")
    x = fisher.test(c_y, DISCOVERY_GROUP)
    result = paste0(result,format(x$estimate,nsmall = 2),"\t",format(x$p.value,nsmall = 2),"\t",format(x$conf.int[1],nsmall = 2),"-",format(x$conf.int[2],nsmall = 2),'\t',"DISCOVERY",'\t',"Biotype B","\n")
    
    location = paste0(repository_path,"/", experiment_path,"/",paste0("feat_",FEATURE_SIZE),"/",paste0("results_discovery_feat_",FEATURE_SIZE,"_iter_",yy,".txt"))
    
    writeLines(result,location)
    location = paste0(repository_path,"/", experiment_path,"/",paste0("feat_",FEATURE_SIZE),"/",paste0("discovery_cv_cluster_assignment_feat_",FEATURE_SIZE,"_iter_",yy,".tsv"))
    tab = data.frame(id = DISCOVERY_IDs, c1 = c_x,  c2 = c_y,c1_prob=c_x_prob,c2_prob=c_y_prob, group = DISCOVERY_GROUP) 
    write.table(tab,location,sep = "\t",row.names = F, quote = F)
  }
}

previous_features_validation = c()

FEATURES_FOR_VALIDATION = sort_by_correlation(DISCOVERY_DATA, as.numeric(DISCOVERY_GROUP))

for(FEATURE_SIZE in features_sizes) {
  
  result = "OR\tp_value\tconf\ttype\tBiotype\n"
  
  SELECTED_FEATURES = as.vector(FEATURES_FOR_VALIDATION[1:FEATURE_SIZE,])
  
  SELECTED_FEATURES = sort(match(SELECTED_FEATURES[!is.na(SELECTED_FEATURES)],colnames(DISCOVERY_DATA)))
  
  cluster_B = clustering(DISCOVERY_DATA, DISCOVERY_GROUP, SELECTED_FEATURES, FEATURE_SIZE, 100)
  
  for(yy in c(1:num_it)) {
    result = "OR\tp_value\tconf\ttype\tBiotype\n"
    
    ################################### create clustering on full discovery data and predict in validation
   
    cluster_Bmatch = rep(NA, nrow(DISCOVERY_DATA))
    cluster_Bmatch[which(DISCOVERY_GROUP == 2)] = cluster_B
    if(cor.test(cluster_Bmatch, cluster_B_CV_RES)$estimate < 0){
      cluster_B = abs(cluster_B-3)
    }
    
    pred_validation_clu1 = train_and_predict(TRAIN=DISCOVERY_DATA, 
                                    VALIDATION=VALIDATION_DATA, 
                                    GROUP=DISCOVERY_GROUP, 
                                    CLUSTER_GROUP=cluster_B, 
                                    TARGET_CLUSTER=1, 
                                    FEATURES = SELECTED_FEATURES, 
                                    REPEATS=repeat_)
    
    
    pred_validation_clu2 = train_and_predict(TRAIN=DISCOVERY_DATA, 
                                    VALIDATION=VALIDATION_DATA, 
                                    GROUP=DISCOVERY_GROUP, 
                                    CLUSTER_GROUP=cluster_B, 
                                    TARGET_CLUSTER=2, 
                                    FEATURES = SELECTED_FEATURES, 
                                    REPEATS=repeat_)
    
    c_x = pred_validation_clu1[[3]]
    c_x_prob = pred_validation_clu1[[2]]
    c_y = pred_validation_clu2[[3]]
    c_y_prob = pred_validation_clu2[[2]]
    
    switch_clusters = FALSE
    
    if(length(previous_features_validation) != 0) {
      if(abs(cor(c_x_prob,previous_features_validation)) < abs(cor(c_y_prob,previous_features_validation))) {
        c_x = pred_validation_clu2[[3]]
        c_y = pred_validation_clu1[[3]]
        c_x_prob = pred_validation_clu2[[2]]
        c_y_prob = pred_validation_clu1[[2]]

        switch_clusters = TRUE
      } 
    }
    
    previous_features_validation = c_x_prob
    
    # save models
    if(switch_clusters == TRUE) {
      saveRDS(pred_validation_clu1[[4]],paste0(repository_path,"/",experiment_path,"/models/model_feat_",FEATURE_SIZE,"_iter_",yy,"_2",".RDS"))
      saveRDS(pred_validation_clu2[[4]],paste0(repository_path,"/",experiment_path,"/models/model_feat_",FEATURE_SIZE,"_iter_",yy,"_1",".RDS"))
    } else {
      saveRDS(pred_validation_clu1[[4]],paste0(repository_path,"/",experiment_path,"/models/model_feat_",FEATURE_SIZE,"_iter_",yy,"_1",".RDS"))
      saveRDS(pred_validation_clu2[[4]],paste0(repository_path,"/",experiment_path,"/models/model_feat_",FEATURE_SIZE,"_iter_",yy,"_2",".RDS"))
    }
    x = fisher.test(table( c_x, c_y ))
    result = paste0(result,format(x$estimate,nsmall = 2),"\t",format(x$p.value,nsmall = 2),"\t",format(x$conf.int[1],nsmall = 2),"-",format(x$conf.int[2],nsmall = 2),'\t',"Validation",'\t',"Inverse","\n")
    x = fisher.test(table( c_x, VALIDATION_GROUP)) 
    result = paste0(result,format(x$estimate,nsmall = 2),"\t",format(x$p.value,nsmall = 2),"\t",format(x$conf.int[1],nsmall = 2),"-",format(x$conf.int[2],nsmall = 2),'\t',"Validation",'\t',"Biotype A","\n")
    x = fisher.test(table( c_y, VALIDATION_GROUP))    
    result = paste0(result,format(x$estimate,nsmall = 2),"\t",format(x$p.value,nsmall = 2),"\t",format(x$conf.int[1],nsmall = 2),"-",format(x$conf.int[2],nsmall = 2),'\t',"Validation",'\t',"Biotype B","\n")
    
    tab = data.frame(id = VALIDATION_IDs, c1 = c_x,  c2 = c_y, c1_prob=c_x_prob,c2_prob=c_y_prob, group = VALIDATION_GROUP)
    location = paste0(repository_path,"/", experiment_path,"/",paste0("feat_",FEATURE_SIZE),"/",paste0("validation_cluster_assignment_feat_",FEATURE_SIZE,"_iter_",yy,".tsv"))
    
    write.table(tab,location,sep = "\t",row.names = F, quote = F)
    
    studies = sapply(FEATURES_FOR_VALIDATION[1:FEATURE_SIZE,1], function(x) {
      str_split(x,"_",n=3)[[1]][3]
    })
    
    location = paste0(repository_path,"/", experiment_path,"/",paste0("feat_",FEATURE_SIZE),"/",paste0("selected_studies_feat_",FEATURE_SIZE,"_iter_",yy,".tsv"))
    
    write.table(studies,location,sep = "\t",row.names = F)
    
    location = paste0(repository_path,"/", experiment_path,"/",paste0("feat_",FEATURE_SIZE),"/",paste0("results_validation_feat_",FEATURE_SIZE,"_iter_",yy,".txt"))
    
    writeLines(result,location)
    cat(result)
  }
}

stopCluster(cl)
