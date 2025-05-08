collect_res = list()
df_discovery = list()
df_validation = list()

for(a in c(1:length(features_sizes))) {
  feat = features_sizes[a]
  fea = paste0(experiment_path,"/feat_",feat)
  res = calculate_feature_statistics(fea,feat,num_it)
  df_discovery_res = res[[1]]
  df_validation_res = res[[2]]

  collect_res[[a]] = read_results(fea,feat,num_it)

  df_discovery_res$c1_median = apply(df_discovery_res[,paste0("c1_",c(1:num_it)),drop =FALSE],1,median)
  df_discovery_res$c2_median = apply(df_discovery_res[,paste0("c2_",c(1:num_it)),drop =FALSE],1,median)

  df_validation_res$c1_median = apply(df_validation_res[,paste0("c1_",c(1:num_it)),drop =FALSE],1,median)
  df_validation_res$c2_median = apply(df_validation_res[,paste0("c2_",c(1:num_it)),drop =FALSE],1,median)


  df_discovery_res$c1_prop_median = apply(df_discovery_res[,paste0("c1_prob_",c(1:num_it)),drop =FALSE],1,median)
  df_discovery_res$c2_prop_median = apply(df_discovery_res[,paste0("c2_prob_",c(1:num_it)),drop =FALSE],1,median)

  df_validation_res$c1_prop_median = apply(df_validation_res[,paste0("c1_prob_",c(1:num_it)),drop =FALSE],1,median)
  df_validation_res$c2_prop_median = apply(df_validation_res[,paste0("c2_prob_",c(1:num_it)),drop =FALSE],1,median)

  df_discovery[[length(df_discovery)+1]] = df_discovery_res
  df_validation[[length(df_validation)+1]] = df_validation_res
}

switch_clusters = function(feat, path, iter) {
  for(f in feat) {
    for(i in c(1:iter)) {
      p = paste0(path,"/feat_",f)
      ite_ = paste0(p,"/results_feat_",f,"_iter_",i,".txt")
      dt = read.table(ite_,header = T,sep = "\t")
      if(i == 1) {
        dt = dt[c(1,3,2,4,6,5,7,9,8),]
      } else {
        dt = dt[c(1,3,2),]
      }
      write.table(dt,ite_,sep = "\t",row.names = F, quote = F)
      
      ite_discovery = paste0(p,"/discovery_cv_cluster_assignment_feat_",f,"_iter_",i,".tsv")
      ite_validation = paste0(p,"/validation_cluster_assignment_feat_",f,"_iter_",i,".tsv")
      df_discovery_res = read.table(ite_discovery,header = T,sep = "\t")
      df_validation_res = read.table(ite_validation,header = T,sep = "\t")
      colnames(df_discovery_res) = c("id",	"c2", "c1",	"c2_prob",	"c1_prob", "group")
      write.table(df_discovery_res,ite_discovery,sep = "\t",row.names = F, quote = F)
      if(i == 1) {
        colnames(df_validation_res) = c("id",	"c2", "c1",	"c2_prob",	"c1_prob", "group")
        write.table(df_validation_res,ite_validation,sep = "\t",row.names = F, quote = F)
      }
    }
  }
}
