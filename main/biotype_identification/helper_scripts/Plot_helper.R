calculate_feature_statistics = function(feat_path,feat,num_it) {
  ite_discovery = paste0(feat_path,"/discovery_cv_cluster_assignment_feat_",feat,"_iter_",1,".tsv")
  ite_validation = paste0(feat_path,"/validation_cluster_assignment_feat_",feat,"_iter_",1,".tsv")

  df_discovery = read.table(ite_discovery,header = T,sep = "\t")
  df_validation = read.table(ite_validation,header = T,sep = "\t")


  colnames(df_discovery)[2:ncol(df_discovery)] = paste0(colnames(df_discovery)[2:ncol(df_discovery)],"_1")
  colnames(df_validation)[2:ncol(df_validation)] = paste0(colnames(df_validation)[2:ncol(df_validation)],"_1")
  if(num_it > 1) {
    for(iter in c(2:num_it)) {
      ite_discovery = paste0(feat_path,"/discovery_cv_cluster_assignment_feat_",feat,"_iter_",iter,".tsv")
      ite_validation = paste0(feat_path,"/validation_cluster_assignment_feat_",feat,"_iter_",iter,".tsv")
      df_b = read.table(ite_discovery,header = T,sep = "\t")
      df_p = read.table(ite_validation,header = T,sep = "\t")
      colnames(df_b)[2:ncol(df_b)] = paste0(colnames(df_b)[2:ncol(df_b)],"_",iter)
      colnames(df_p)[2:ncol(df_p)] = paste0(colnames(df_p)[2:ncol(df_p)],"_",iter)
      df_discovery = merge(df_discovery,df_b,by = "id")
      df_validation = merge(df_validation,df_p,by = "id")
    }
  }
  ls = list()
  ls[[1]] = df_discovery
  ls[[2]] = df_validation
  return(ls)
}

read_results = function(feat_path,feat,num_it) {
  ite_ = paste0(feat_path,"/results_discovery_feat_",feat,"_iter_",1,".txt")
  df_ = read.table(ite_,header = T,sep = "\t")

  a = df_[1,1:2]
  ab = df_[2:3,1:2]

  ite_ = paste0(feat_path,"/results_validation_feat_",feat,"_iter_",1,".txt")
  df_ = read.table(ite_,header = T,sep = "\t")
  c = df_[1,1:2]
  cd = df_[2:3,1:2]

  ax = list(a,ab,c,cd,1)
  df_end = t(do.call(rbind,ax))
  colnames(df_end) = c("Discovery_Inverse","Discovery_C1","Discovery_C2","Validation_Inverse","Validation_C1","Validation_C2","iter")
  if(num_it > 1) {
    for(iter in c(2:num_it)) {
      ite_ = paste0(feat_path,"/results_discovery_feat_",feat,"_iter_",iter,".txt")
      df_ = read.table(ite_,header = T,sep = "\t")

      ite_1 = paste0(feat_path,"/results_validation_feat_",feat,"_iter_",iter,".txt")
      df_1 = read.table(ite_1,header = T,sep = "\t")

      a = df_[1,1:2]
      ab = df_[2:3,1:2]

      c = df_1[1,1:2]
      cd = df_1[2:3,1:2]

      ax = list(a,ab,c,cd,iter)
      df_ = t(do.call(rbind,ax))
      colnames(df_) = c("Discovery_Inverse","Discovery_C1","Discovery_C2","Validation_Inverse","Validation_C1","Validation_C2","iter")
      df_end = rbind(df_end,df_)
    }
  }
  return(df_end)
}

calculate_stat_for_res = function(dfs, iter) {
  dt = data.frame(label = character(), or = numeric(), feature = numeric(), p_value = numeric(), minus_log10_p_value = numeric(), iter = numeric(), color = character(),datatype = character(), significant = character())
  count = 1
  for(df in dfs) {
    df = as.data.frame(df)
    for(i in 1:iter) {
      dfx = df[df$iter == i,]
      dt = rbind(dt,list("C1",as.numeric(dfx["OR","Discovery_C1"]),features_sizes[count], as.numeric(dfx["p_value","Discovery_C1"]), -log10(p.adjust(as.numeric(dfx["p_value","Discovery_C1"]), method = "BH")),i, "red", "Discovery", check_sig(as.numeric(dfx["p_value","Discovery_C1"]))))
      dt = rbind(dt,list("C2",as.numeric(dfx["OR","Discovery_C2"]),features_sizes[count], as.numeric(dfx["p_value","Discovery_C2"]), -log10(p.adjust(as.numeric(dfx["p_value","Discovery_C2"]), method = "BH")),i, "green", "Discovery", check_sig(as.numeric(dfx["p_value","Discovery_C2"]))))
      dt = rbind(dt,list("Inverse",as.numeric(dfx["OR","Discovery_Inverse"]),features_sizes[count], as.numeric(dfx["p_value","Discovery_Inverse"]), -log10(p.adjust(as.numeric(dfx["p_value","Discovery_Inverse"]), method = "BH")),i,"blue", "Discovery", check_sig(as.numeric(dfx["p_value","Discovery_Inverse"]))))

      dt = rbind(dt,list("C1",as.numeric(dfx["OR","Validation_C1"]),features_sizes[count], as.numeric(dfx["p_value","Validation_C1"]), -log10(p.adjust(as.numeric(dfx["p_value","Validation_C1"]), method = "BH")),i, "red", "Validation", check_sig(as.numeric(dfx["p_value","Validation_C1"]))))
      dt = rbind(dt,list("C2",as.numeric(dfx["OR","Validation_C2"]),features_sizes[count], as.numeric(dfx["p_value","Validation_C2"]), -log10(p.adjust(as.numeric(dfx["p_value","Validation_C2"]), method = "BH")),i, "green", "Validation", check_sig(as.numeric(dfx["p_value","Validation_C2"]))))
      dt = rbind(dt,list("Inverse",as.numeric(dfx["OR","Validation_Inverse"]),features_sizes[count], as.numeric(dfx["p_value","Validation_Inverse"]), -log10(p.adjust(as.numeric(dfx["p_value","Validation_Inverse"]), method = "BH")),i,"blue", "Validation", check_sig(as.numeric(dfx["p_value","Validation_Inverse"]))))

    }
    count = count +1
  }
  colnames(dt) = c("label","or","feature","p_value","minus_log10_p_value","iter","color","datatype", "significant")
  return(dt)
}

define_groups = function(df,c1_column_name,c2_column_name,feat) {
  
  gr = rep(NA,nrow(df))
  
  gr[which(df[,c1_column_name] == 1 & df[,c2_column_name] == 1 & df$group_1 == 1)] = 1 #"both_right"
  gr[which(df[,c1_column_name] == 2 & df[,c2_column_name] == 2 & df$group_1 == 2)] = 1 #"both_right"
  
  gr[which(df[,c1_column_name] == 1 & df[,c2_column_name] == 1 & df$group_1 == 2)] = 2 #"both_wrong"
  gr[which(df[,c1_column_name] == 2 & df[,c2_column_name] == 2 & df$group_1 == 1)] = 2 #"both_wrong"
  
  gr[which(df[,c1_column_name] == 2 & df[,c2_column_name] == 1 & df$group_1 == 1)] = 3#"unique_c2"
  gr[which(df[,c1_column_name] == 1 & df[,c2_column_name] == 2 & df$group_1 == 1)] = 4#"unique_c1"
  
  gr[which(df[,c1_column_name] == 2 & df[,c2_column_name] == 1 & df$group_1 == 2)] = 4#"unique_c1"
  gr[which(df[,c1_column_name] == 1 & df[,c2_column_name] == 2 & df$group_1 == 2)] = 3#"unique_c2"
  
  end_df = data.frame(gr)
  end_df$feat = feat
  colnames(end_df) = c("group","feature")
  end_df$id = df$id
  end_df$outcome = df$group_1
  end_df$c1 = df[,c1_column_name]
  end_df$c2 = df[,c2_column_name]
  end_df$main_group = NA
  end_df[which(df$group_1 == 1),]$main_group = "HC"
  end_df[which(df$group_1 == 2),]$main_group = "SCZ"
  return(end_df)
}

calc_NagelkerkeR2_vector = function(vec,group) {
  m = glm((group - 1)~ vec, family = binomial())
  return(NagelkerkeR2(m)$R2)
}

plot_adj_heatmap = function(df,col1,col2, title,th) {

  res = list()
  for(i in c(1:length(features_sizes))) {
    df_b = df[[i]]
    cl = list()
    dx = define_groups(df_b,col1,col2, features_sizes[i])

    res[[i]] = dx
  }

  row_adj = list()
  for(i in c(1:(length(features_sizes)))) {
    adj_row = c()
    for(j in c(1:(length(features_sizes)))) {
      gr1 = as.numeric(res[[i]]$group)
      gr2 = as.numeric(res[[j]]$group)
      x = RRand(gr1,gr2)
      adj_row = c(adj_row,x$adjRand)
    }
    row_adj[[i]] = adj_row
  }
  data = do.call(rbind,row_adj)
  
  melted_cormat <- reshape2::melt(data)
  pa = wes_palettes %>% names()
  pal = wes_palette(name = pa[20], n = 40, type="continuous")
  colnames(melted_cormat) = c("X1","X2","value")
  melted_cormat[melted_cormat < th] = 0
  p = ggplot(data = melted_cormat, aes(X2, X1, fill = value))+
    geom_tile(color = "black")+coord_cartesian(expand = FALSE) +
    theme_bw() %+%
  theme(panel.background = element_rect(fill = "grey90"),
        panel.grid.major = element_line(linetype = 3, colour = "grey60"),
        axis.text = element_text(colour = 1, size = 10),
        axis.title = element_text(colour = 1, size = 12),
        legend.background = element_blank(),
        legend.key = element_blank(),
       legend.position = "right")+
       scale_fill_gradientn(colours = pal, breaks = c(0,th,1),labels = c("< threshold",paste0("threshold > ",th),"1"))+
       #scale_fill_gradient2(midpoint = th,high = scales::muted("red"),low = scales::muted("blue"))+
       guides(fill = guide_legend(title.position = "top",direction = "vertical",
                               title.theme = element_text(angle = 0, size = 12, colour = "black"),
                               barheight = .5, barwidth = .95,
                               title.hjust = 0.5, raster = FALSE,
                               title = "ARI"))+scale_color_hue(direction = -1, h.start=90)+
                               scale_x_continuous(labels = features_sizes, breaks = c(1:length(features_sizes)))+
                               scale_y_continuous(labels = features_sizes, breaks = c(1:length(features_sizes)))


    #geom_text(aes(label = round(value, digits = 2)), color = "black", size = 4) +
    #scale_fill_gradient2(low = "white", high = "red",mid ="white",
    #                     midpoint = th, limit = c(0,1), space = "Lab",
    #                     name="Adjusted\nRand Index") +
    #theme(axis.text=element_text(size=8))+
    #scale_x_continuous(labels = ls, breaks = breaks) +
    #scale_y_continuous(labels = ls, breaks = breaks)
  p = ggpar(p, xlab = "Number of selected studies", ylab = "Number of selected studies", x.text.angle = 90)

  return(p)
}

plot_adj_heatmap_for_k_means = function(path,th) {

  res = list()
  for(i in c(1:length(features_sizes))) {
    res[[i]] = readLines(paste0(path,"/feat_",features_sizes[i],"/kmeans_feat_",features_sizes[i]))
  }

  row_adj = list()
  for(i in c(1:(length(features_sizes)))) {
    adj_row = c()
    for(j in c(1:(length(features_sizes)))) {
      gr1 = as.numeric(res[[i]])
      gr2 = as.numeric(res[[j]])
      x = RRand(gr1,gr2)
      adj_row = c(adj_row,x$adjRand)
    }
    row_adj[[i]] = adj_row
  }
  data = do.call(rbind,row_adj)

  melted_cormat <- reshape2::melt(data)
  pa = wes_palettes %>% names()
  pal = wes_palette(name = pa[20], n = 40, type="continuous")
  colnames(melted_cormat) = c("X1","X2","value")
  melted_cormat[melted_cormat < th] = 0
  p = ggplot(data = melted_cormat, aes(X2, X1, fill = value))+
    geom_tile(color = "black")+coord_cartesian(expand = FALSE) +
    theme_bw() %+%
  theme(panel.background = element_rect(fill = "grey90"),
        panel.grid.major = element_line(linetype = 3, colour = "grey60"),
        axis.text = element_text(colour = 1, size = 10),
        axis.title = element_text(colour = 1, size = 12),
        legend.background = element_blank(),
        legend.key = element_blank(),
       legend.position = "right")+
       scale_fill_gradientn(colours = pal, breaks = c(0,th,1),labels = c("< threshold",paste0("threshold: ",th),"1"))+
       #scale_fill_gradient2(midpoint = th,high = scales::muted("red"),low = scales::muted("blue"))+
       guides(fill = guide_legend(title.position = "top",direction = "vertical",
                               title.theme = element_text(angle = 0, size = 12, colour = "black"),
                               barheight = .5, barwidth = .95,
                               title.hjust = 0.5, raster = FALSE,
                               title = "ARI"))+scale_color_hue(direction = -1, h.start=90)+
                               scale_x_continuous(labels = features_sizes, breaks = c(1:length(features_sizes)))+
                               scale_y_continuous(labels = features_sizes, breaks = c(1:length(features_sizes)))
  p = ggpar(p, xlab = "Number of selected studies", ylab = "Number of selected studies", x.text.angle = 90)

  return(p)
}

plot_heatmap_for_counts = function(data) {
  melted_cormat <- reshape2::melt(data)
  pa = wes_palettes %>% names()
  pal = wes_palette(name = pa[20], n = 40, type="continuous")
  colnames(melted_cormat) = c("X1","X2","value")
  melted_cormat[melted_cormat > 50] = max(melted_cormat$value)
  p = ggplot(data = melted_cormat, aes(X2, X1, fill = value))+
    geom_tile(color = "black")+coord_cartesian(expand = FALSE) +
    theme_bw() %+%
    theme(panel.background = element_rect(fill = "grey90"),
          panel.grid.major = element_line(linetype = 3, colour = "grey60"),
          axis.text = element_text(colour = 1, size = 10),
          axis.title = element_text(colour = 1, size = 12),
          legend.background = element_blank(),
          legend.key = element_blank(),
          legend.position = "right")+
    scale_fill_gradientn(colours = pal, breaks = c(0,50,max(melted_cormat$value)),labels = c("0",paste0("threshold: ",50),max(melted_cormat$value)))+
    #scale_fill_gradient2(midpoint = 0,high = scales::muted("red"),low = scales::muted("blue"))+
    guides(fill = guide_legend(title.position = "top",direction = "vertical",
                               title.theme = element_text(angle = 0, size = 12, colour = "black"),
                               barheight = .5, barwidth = .95,
                               title.hjust = 0.5, raster = FALSE,
                               title = "ARI"))+scale_color_hue(direction = -1, h.start=90)+
    scale_x_continuous(labels = features_sizes, breaks = c(1:length(features_sizes)))+
    scale_y_continuous(labels = features_sizes, breaks = c(1:length(features_sizes)))
  p = ggpar(p, xlab = "Number of selected studies", ylab = "Number of selected studies", x.text.angle = 90)
  
  return(p)
}

change_labels = function(x,rep,new) {
  lab = x
  for(i in c(1:length(rep))) {
    lab = str_replace(lab, rep[i], new[i])
  }
  return(lab)
}

or_and_pvalue = function(data) {
  ina = c()
  for(i in unique(data$feature)) {
    for(j in unique(data$label)) {
      ors = data[data$feature == i & data$label == j,]$or
      median_or = median(ors)
      which(data$feature == i & data$label == j & data$or == median_or)
      ina = c(ina,which(data$feature == i & data$label == j & data$or == median_or))
    }
  }
  return(ina)
}

or_plot = function(training_test, feat, plot_b) {
  plot_this = abc[abc$datatype ==training_test,]
  plot_this = plot_this[plot_this$label !="Inverse", ]
  plot_this$label = change_labels(plot_this$label,c("C1","C2"), c("Biotype 1","Biotype 2"))
  plot_this = plot_this[or_and_pvalue(plot_this),]
  options(repr.plot.width =9, repr.plot.height =9) 
  feat = c(0, feat)
  #breaks = seq(from = 0, to = max(feat), by = 5)
  #breaks[1] = 1
  or_plot = ggline(plot_this, x = "feature", y = "or", linetype = "solid", shape=NA, size = 6,
                   color = "label",
                   add.params = list(linetype = "dotted"))+ ylim(0.5,4.2)+
    scale_color_manual(values=group.colors)+
    #geom_hline(aes(yintercept =1), linetype = "dotted", color="chartreuse3", linewidth=0.5)+
    theme(axis.text=element_text(size=40,color="black", face = "bold"), axis.title = element_blank(),
          panel.border = element_blank(),panel.background = element_rect(fill = plot_b,
                                                                         colour = plot_b,
                                                                         size = 0.5, linetype = "solid"),
          panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                          colour = "gray55"), 
          panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                          colour = "gray55"),
          
          strip.background = element_rect(fill = "transparent", colour = "transparent"),
          plot.background = element_rect(fill = plot_b),legend.key.size = unit(5,"line"),
          legend.background = element_blank(), legend.position = "none", axis.line =element_blank(),
          legend.box.background = element_blank(), legend.text = element_text(size=45,color="black"),
          legend.key = element_blank())+guides(color = guide_legend(override.aes = list(size = 10))) 
    #scale_x_discrete(labels = as.character(ls), breaks = c(1:length(ls)))
  or_plot = ggpar(or_plot,legend.title = "", xlab =  "number of selected studies", ylab = "odds ratio", legend = c(0.80,0.95)) # 
  return(or_plot)
}

density_plot = function(training_test,sumstat, feat, path, title) {
  if (!file.exists(paste0(path,"/density"))){
    dir.create(file.path(path, "/density"))
  }
  path = paste0(path,"/density")
  for(l in c(1:length(feat))) {
    c1 = which(training_test[[l]]$c1_median == 2 & training_test[[l]]$group_1 ==2)
    c2 = which(training_test[[l]]$c2_median == 2 & training_test[[l]]$group_1 ==2)
    he = which(training_test[[l]]$group_1 ==1)
    
    plots = list()
    for(x in c(1:ncol(sumstat))) {
      rs = paste0("########################\n",
      paste0(title," number of GWASs ", feat[l]),"\n",
      colnames(sumstat)[x],"\n",
      "C1 vs C2","\n",
      paste0("sig :",format.pval(t.test(sumstat[c1,x],sumstat[c2,x])$p.value)),"\n",
      t.test(sumstat[c1,x],sumstat[c2,x])$estimate[1],"-",t.test(sumstat[c1,x],sumstat[c2,x])$estimate[2],"\n",
      "C1 vs HC","\n",
      paste0("sig :",format.pval(t.test(sumstat[c1,x],sumstat[he,x])$p.value)),"\n",
      t.test(sumstat[c1,x],sumstat[he,x])$estimate[1],"-",t.test(sumstat[c1,x],sumstat[he,x])$estimate[2],"\n",
      "C2 vs HC","\n",
      paste0("sig :",format.pval(t.test(sumstat[c2,x],sumstat[he,x])$p.value)),"\n",
      t.test(sumstat[c2,x],sumstat[he,x])$estimate[1],"-",t.test(sumstat[c2,x],sumstat[he,x])$estimate[2],"\n")
      cat(rs)
      c1_density = sumstat[c1,x]
      c2_density = sumstat[c2,x]
      h_density = sumstat[he,x]
      p_va = str_split(colnames(sumstat)[x],"_",3)[[1]][2]
      stu = str_split(colnames(sumstat)[x],"_",3)[[1]][3]
      df = data.frame(x = c(c1_density,c2_density,h_density), 
                      group = c(rep("C1",length(c1_density)),rep("C2",length(c2_density)),rep("H",length(h_density))),
                      p_val = p_va,
                      feat = feat[l])
      mu <- ddply(df, "group", summarise, grp.mean=mean(x))
      
      p<-ggplot(df, aes(x=x, fill=group))+scale_fill_manual(values=c("blue", "red", "green"))+
        geom_density(alpha=0.7)+ ggtitle(paste0("p-value=",p_va))+theme(plot.title = element_text(size = 10))+
        geom_vline(data=mu, aes(xintercept=grp.mean),
                   linetype="dashed")
      
      plots[[x]] = p
    }
    res = annotate_figure(ggarrange(plotlist = plots, common.legend = T,legend="bottom"), top = text_grob(paste0(title,", Number of GWASs used as feature = ",feat[l]), size = 12))
    res %>% export_plot(paste0(path,"/density_",stu,"_feat_",feat[l],".png"))
  }
}

mean_plot = function(training_test,sumstat, feat, title, best_prs, prs_name, df = NULL, plot_b = NULL, y_min = NULL,y_max = NULL) {
  if(is.null(df)) {
    c1_means_total = c()
    c2_means_total = c()
    healthy_total = c()
    
    col_name = colnames(sumstat)
    nn = c()
    if(length(prs_name) == 1) {
      nn = col_name[grepl(prs_name,col_name)]
    } else {
      for(p in prs_name) {
        nn = c(nn, col_name[grepl(p,col_name)])
      }
    }
    
    ss = sumstat[,nn,drop = F]
    x = prcomp(ss,center = T, scale = T)
    
    cor_with_dis = cor(rowSums(ss),x$x[,1])
    
    if(cor_with_dis < 0 ) {
      pcaa = x$x[,1]*(-1)
    } else {
      pcaa = x$x[,1]  
    }
    pcaa = scale(pcaa)
    for(l in c(1:length(feat))) {
      c1_m = c()
      c2_m = c()
      h_m = c()
      c1 = which(training_test[[l]]$c1_median == 2 & training_test[[l]]$group_1 ==2)
      c2 = which(training_test[[l]]$c2_median == 2 & training_test[[l]]$group_1 ==2)
      he = which(training_test[[l]]$group_1 ==1)
      for(xi in c(1:3)) {
        c1_means = c()
        c2_means = c()
        h_means = c()
        min_clust_size = min(length(c1), length(c2))
        
        #c1 = sample(c1, min_clust_size)
        #c2 = sample(c2, min_clust_size)
        c1 = c1
        c2 = c2
        if(dim(ss)[2] == 1) {
          c1_mean = mean(ss[c1,1]) # mean(ss[c1,x])
          c2_mean = mean(ss[c2,1])
          h_mean = mean(ss[he,1])
        } else {
          
          c1_mean = mean(pcaa[c1]) # mean(ss[c1,x])
          c2_mean = mean(pcaa[c2])
          h_mean = mean(pcaa[he])
        }
        c1_m = c(c1_m,c1_mean)
        c2_m = c(c2_m,c2_mean)
        h_m = c(h_m,h_mean)
      }
      c1_means = median(c1_m)
      c2_means = median(c2_m)
      h_means = median(h_m)
      c1_means_total = c(c1_means_total,c1_means)
      c2_means_total = c(c2_means_total,c2_means)
      healthy_total = c(healthy_total,h_means)
    }
    df = data.frame(x = feat, y = c(c1_means_total,c2_means_total,healthy_total), z = c(rep("Biotype 1",length(c1_means_total)),rep("Biotype 2",length(c2_means_total)),rep("Healthy",length(healthy_total))))
  }
  df$x = as.numeric(df$x)
  #df$y = scale(as.numeric(df$y))
  df$z = as.character(df$z)
  df$x = as.factor(df$x)
  or_plot = ggline(df, x = "x", y ="y", color ="z",linetype = "solid", shape=NA, size = 6,) +
    scale_color_manual(values=group.colors)+ ylim(y_min,y_max)+
    theme(axis.text=element_text(size=40,color="black", face = "bold"), axis.title = element_blank(),
          panel.border = element_blank(),panel.background = element_rect(fill = plot_b,
                                                                         colour = plot_b,
                                                                         size = 0.5, linetype = "solid"),
          panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                          colour = "gray55"), 
          panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                          colour = "gray55"),
          
          strip.background = element_rect(fill = "transparent", colour = "transparent"),
          plot.background = element_rect(fill = plot_b),legend.key.size = unit(5,"line"),
          legend.background = element_blank(), legend.position = "none", axis.line =element_blank(),
          legend.box.background = element_blank(), legend.text = element_text(size=45,color="black"),
          legend.key = element_blank())+guides(color = guide_legend(override.aes = list(size = 10))) 
    
  
  or_plot = ggpar(or_plot,legend.title = "", xlab =  "number of selected studies", legend = c(0.80,0.85),
                  ylab = "odds ratio")
  #or_plot = ggpar(or_plot,legend.title = "", xlab =  "number of selected studies",
   #               ylab = "odds ratio")
  return(or_plot)
}

pvalue_plot = function(training_test,feat, plot_b) {
  plot_this = abc[abc$label %in% c("C1","C2"),]
  plot_this = plot_this[plot_this$datatype ==training_test,]
  plot_this$label = change_labels(plot_this$label,c("C1","C2"), c("Biotype 1","Biotype 2"))
  plot_this = plot_this[or_and_pvalue(plot_this),]
  feat = c(0, feat)
  breaks = seq(from = 0, to = max(feat), by = 5)
  breaks[1] = 1
  options(repr.plot.width =9, repr.plot.height =9) 
  
  pvalue_plot = ggline(plot_this, x = "feature", y = "minus_log10_p_value", linetype = "solid", shape=NA, size = 6,
                       color = "label", add.params = list(linetype = "dotted"))+
    geom_hline(aes(yintercept =-log10(0.05)), linetype = "solid", color="red", linewidth=2)+
    scale_color_manual(values=group.colors)+
    theme(axis.text=element_text(size=40,color="black", face = "bold"), axis.title = element_blank(),
          panel.border = element_blank(),panel.background = element_rect(fill = plot_b,
                                                                         colour = plot_b,
                                                                         size = 0.5, linetype = "solid"),
          panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                          colour = "gray55"), 
          panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                          colour = "gray55"),
          
          strip.background = element_rect(fill = "transparent", colour = "transparent"),
          plot.background = element_rect(fill = plot_b),
          legend.background = element_blank(), legend.position = "none", axis.line =element_blank(),
          legend.box.background = element_blank(), legend.text = element_text(size=16,color="black"),
          legend.key = element_blank())#+
    #scale_x_discrete(labels = as.character(ls), breaks = c(1:length(ls)))+ ylim(-0.2,max(plot_this$minus_log10_p_value)+1)
  #pvalue_plot = ggpar(pvalue_plot,legend.title = "" ,xlab ="" ,ylab="", legend = "none", y.text.angle = 90)
  pvalue_plot =  ggpar(pvalue_plot,legend.title = "", xlab = "number of selected studies", 
                     ylab = "-log10 p-value (BH corrected)",  y.text.angle = 90)
  return(pvalue_plot)
}

explained_var_plot = function(training_test,feat) {
  plot_this = df_explained_var[df_explained_var$datatype ==training_test,]
  plot_this = plot_this[plot_this$label %in% c("C1","C2"),]
  plot_this$label = change_labels(plot_this$label,c("C1","C2"), c("Biotpye 1","Biotype 2"))
  feat = c(0, feat)
  breaks = seq(from = 0, to = max(feat), by = 5)
  breaks[1] = 1
  options(repr.plot.width =9, repr.plot.height =9) 
  expressed_variance = ggline(plot_this, x = "feature", y = "nagelkerkeR2", linetype = "solid", shape=19,
                              color = "label", add ="median", add.params = list(linetype = "dotted"))+
    scale_color_manual(values=group.colors)+
    theme(axis.text=element_text(size=10, color = "black"),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
          plot.background = element_rect(fill = "transparent", colour = "transparent"),
          panel.background = element_rect(fill = "transparent", colour = "gray"),
          panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
          strip.background = element_rect(fill = "transparent", colour = "transparent"),
          legend.background = element_blank(),
          legend.box.background = element_blank(),
          legend.key = element_blank())#+
    #scale_x_discrete(labels = as.character(ls), breaks = c(1:length(ls))) +ylim(-0.005,max(plot_this$nagelkerkeR2)+0.01)
  expressed_variance =ggpar(expressed_variance,legend.title = "", xlab = "number of selected studies", ylab = "Nagelkerke R2", y.text.angle = 90, legend = c(0.9,0.9))
  return(expressed_variance)
}

calculate_counts = function(data,cluster,feat) {
  res = list()
  for(i in feat) {
    df_b = data[[which(i==features_sizes)]]
    dx = define_groups(df_b,paste0(cluster,"_median"),paste0(cluster,"_median"), features_sizes[i])
    res[[which(feat == i)]] = dx
  }
  dtx = do.call(cbind,res)
  
  group_switch_count <- function(vec) {
    sum(sapply(1:(length(vec) - 1), function(i) vec[i] + vec[i + 1] == 3))
  }
  
  
  group_s_counts = apply(dtx[,grepl("^group",colnames(dtx))], 1, group_switch_count)
  return(group_s_counts)
}

check_sig = function(value) {
  if(value < 0.05) {
    return(19)
  } else {
    return(21)
  }
}

create_beautiful_radarchart <- function(data, color = "#00AFBB",
                                        vlabels = colnames(data), vlcex = 1.5,
                                        caxislabels = NULL, title = NULL,legend_label = NULL, ...){
  radarchartcirc(
    data, axistype = 1, maxmin = T,
    # Customize the polygon
    pcol = color, #pfcol = scales::alpha(color, 0.3), 
    pfcol = NA,
    plwd = 2, plty = c(1),
    # Customize the grid
    cglcol = "grey", cglty = 1, cglwd = 0.8, centerzero = T,
    # Customize the axis
    axislabcol = "black",
    # Variable labels
    vlcex = vlcex, vlabels = vlabels,
    caxislabels = caxislabels, title = title, ...
  )
  if(is.null(legend_label)) {
    legend_label = c("Cluster 1","Cluster 2")
  }
  op <- par(cex = 0.6)
  legend(0.9,1.4,
         legend=legend_label,pt.cex = 1.2, cex = 1.5,
         pch=c(16),
         col=color,
         lty=c(1), bty = "n", title = NA)
}

create_beautiful_radarchart_fill <- function(data, color = "#00AFBB",
                                        vlabels = colnames(data), vlcex = 1.2, pfcol = scales::alpha(color, 0.3),
                                        caxislabels = NULL, title = NULL,legend_label = NULL, vlabel_col = NULL, ...){
  radarc(
    data, axistype = 1, maxmin = T,
    # Customize the polygon
    pcol = color, pfcol = pfcol, plwd = 2, plty = c(1),
    # Customize the grid
    cglcol = c("gray35","gray35","gray4","gray35","gray35"), cglty = c(1,2,1,1), cglwd = c(0.8,3,0.8,0.8), centerzero = F,
    # Customize the axis
    axislabcol = "black",
    # Variable labels
    vlcex = vlcex, vlabels = vlabels,calcex = 1.5,
    caxislabels = caxislabels, title = title,vlabel_col =vlabel_col, ...
  )
  if(is.null(legend_label)) {
    legend_label = c("Biotype 1","Biotype 2")
  }
  op <- par(cex = 0.6)
  #legend(0.8,1.4, cex = 1.3,
  #       legend=legend_label,
  #       pch=c(16),
  #       col=color,
  #       lty=c(1), bty = "n", title = NA)
}

radarc = function (df, axistype = 0, seg = 2, pty = 32, pcol = 1:8, plty = 1:6, 
          plwd = 1, pdensity = NULL, pangle = 45, pfcol = NA, cglty = 1, 
          cglwd = 1, cglcol = "navy", axislabcol = "blue", title = "", 
          maxmin = TRUE, na.itp = TRUE, centerzero = FALSE, vlabels = NULL, 
          vlcex = NULL, caxislabels = NULL, calcex = NULL, paxislabels = NULL, 
          palcex = NULL,vlabel_col = NULL ,...) 
{
  if (!is.data.frame(df)) {
    cat("The data must be given as dataframe.\n")
    return()
  }
  if ((n <- length(df)) < 3) {
    cat("The number of variables must be 3 or more.\n")
    return()
  }
  if (maxmin == FALSE) {
    dfmax <- apply(df, 2, max)
    dfmin <- apply(df, 2, min)
    df <- rbind(dfmax, dfmin, df)
  }
  #par(bg = "aliceblue")
  
  plot(c(-1.02, 1.02), c(-1.02, 1.02), type = "n", frame.plot = FALSE, 
       axes = FALSE, xlab = "", ylab = "", asp = 1, cex.main = 1.5,
       ...)
  theta <- seq(90, 450, length = n + 1) * pi/180
  theta <- theta[1:n]
  xx <- cos(theta)
  yy <- sin(theta)
  theta2 <- (90 + 0:120 * 3) * pi/180
  xxx <- cos(theta2)
  yyy <- sin(theta2)
  CGap <- ifelse(centerzero, 0, 1)
  for (i in 0:seg) {
    polygon(xxx * (i + CGap)/(seg + CGap), yyy * (i + CGap)/(seg + 
                                                               CGap), lty = cglty[i+1], lwd = cglwd[i+1], border = cglcol[i+1] )
    if (axistype == 1 | axistype == 3) 
      CAXISLABELS <- paste(i/seg * 100, "(%)")
    if (axistype == 4 | axistype == 5) 
      CAXISLABELS <- sprintf("%3.2f", i/seg)
    if (!is.null(caxislabels) & (i < length(caxislabels))) 
      CAXISLABELS <- caxislabels[i + 1]
    if (axistype == 1 | axistype == 3 | axistype == 4 | axistype == 
        5) {
      #if (is.null(calcex)) 
      #  text(0, (i + CGap)/(seg + CGap), substitute(paste(bold(CAXISLABELS))), 
      #       col = axislabcol)
      #else text(0, (i + CGap)/(seg + CGap), substitute(paste(bold(CAXISLABELS))), 
      #          col = axislabcol, cex = calcex)
    }
  }
  if (centerzero) {
    arrows(0, 0, xx * 1, yy * 1, lwd = cglwd[1], lty = cglty[1], 
           length = 0, col = cglcol[1])
  }
  else {
    arrows(xx/(seg + CGap), yy/(seg + CGap), xx * 1, yy * 
             1, lwd = cglwd[1], lty = cglty[1], length = 0, col = cglcol[1])
  }
  PAXISLABELS <- df[1, 1:n]
  if (!is.null(paxislabels)) 
    PAXISLABELS <- paxislabels
  if (axistype == 2 | axistype == 3 | axistype == 5) {
    if (is.null(palcex)) 
      text(xx[1:n], yy[1:n], PAXISLABELS, col = axislabcol)
    else text(xx[1:n], yy[1:n], PAXISLABELS, col = axislabcol, 
              cex = palcex)
  }
  
  VLABELS <- vlabels
  
  if (is.null(vlcex)) {
      #text(xx * 1.15, yy * 1.15, VLABELS, col = vlabel_col)
      for(au in c(1:ncol(df))) {
        text.col =par("col") 
        if(!is.null(VLABELS)) {
          if(vlabel_col[au] == "darkblue") {
            text.col = "black"
          }
        }
        if(is.null(VLABELS)) {
          #text(xx[au] * 1.15, yy[au]  * 1.15, VLABELS[au])
          legend(xx[au] * 1.08 -0.08, yy[au] * 1.08 +0.07,"",text.font=2, adj = 0.9, text.col = text.col)
        } else {
          #text(xx[au] * 1.15, yy[au]  * 1.15, VLABELS[au])
          #legend(xx[au] * 1.08 -0.08, yy[au] * 1.08 +0.07, VLABELS[au],text.font=2, box.col = "transparent", bg = "transparent", adj = 1.1, text.col = text.col, cex = 2.5)
        }
        
      }
    }
    
  else {
      #text(xx * 1.15, yy * 1.15, VLABELS, cex = vlcex, col = vlabel_col)
      for(au in c(1:ncol(df))) {
        text.col =par("col") 
        if(!is.null(VLABELS)) {
          if(vlabel_col[au] == "darkblue") {
            text.col = "black"
          }
        }
        if(is.null(VLABELS)) {
          #text(xx[au] * 1.15, yy[au]  * 1.15, VLABELS[au])
          legend(xx[au] * 1.08-0.08, yy[au] * 1.08+0.07,"",text.font=2, adj = 1, text.col = text.col, cex = 1.5)
        } else {
          #text(xx[au] * 1.15, yy[au]  * 1.15, VLABELS[au])
          legend(xx[au] * 1.08-0.05, yy[au] * 1.08+0.05, VLABELS[au],text.font=2, box.col = vlabel_col[au], bg = vlabel_col[au], adj = 1.0, text.col = text.col, cex = 2.0)
        }
        
      }
    }
  series <- length(df[[1]])
  SX <- series - 2
  if (length(pty) < SX) {
    ptys <- rep(pty, SX)
  }
  else {
    ptys <- pty
  }
  if (length(pcol) < SX) {
    pcols <- rep(pcol, SX)
  }
  else {
    pcols <- pcol
  }
  if (length(plty) < SX) {
    pltys <- rep(plty, SX)
  }
  else {
    pltys <- plty
  }
  if (length(plwd) < SX) {
    plwds <- rep(plwd, SX)
  }
  else {
    plwds <- plwd
  }
  if (length(pdensity) < SX) {
    pdensities <- rep(pdensity, SX)
  }
  else {
    pdensities <- pdensity
  }
  if (length(pangle) < SX) {
    pangles <- rep(pangle, SX)
  }
  else {
    pangles <- pangle
  }
  if (length(pfcol) < SX) {
    pfcols <- rep(pfcol, SX)
  }
  else {
    pfcols <- pfcol
  }
  for (i in 3:series) {
    xxs <- xx
    yys <- yy
    scale <- CGap/(seg + CGap) + (df[i, ] - df[2, ])/(df[1, 
    ] - df[2, ]) * seg/(seg + CGap)
    if (sum(!is.na(df[i, ])) < 3) {
      cat(sprintf("[DATA NOT ENOUGH] at %d\n%g\n", i, df[i, 
      ]))
    }
    else {
      for (j in 1:n) {
        if (is.na(df[i, j])) {
          if (na.itp) {
            left <- ifelse(j > 1, j - 1, n)
            while (is.na(df[i, left])) {
              left <- ifelse(left > 1, left - 1, n)
            }
            right <- ifelse(j < n, j + 1, 1)
            while (is.na(df[i, right])) {
              right <- ifelse(right < n, right + 1, 1)
            }
            xxleft <- xx[left] * CGap/(seg + CGap) + 
              xx[left] * (df[i, left] - df[2, left])/(df[1, 
                                                         left] - df[2, left]) * seg/(seg + CGap)
            yyleft <- yy[left] * CGap/(seg + CGap) + 
              yy[left] * (df[i, left] - df[2, left])/(df[1, 
                                                         left] - df[2, left]) * seg/(seg + CGap)
            xxright <- xx[right] * CGap/(seg + CGap) + 
              xx[right] * (df[i, right] - df[2, right])/(df[1, 
                                                            right] - df[2, right]) * seg/(seg + CGap)
            yyright <- yy[right] * CGap/(seg + CGap) + 
              yy[right] * (df[i, right] - df[2, right])/(df[1, 
                                                            right] - df[2, right]) * seg/(seg + CGap)
            if (xxleft > xxright) {
              xxtmp <- xxleft
              yytmp <- yyleft
              xxleft <- xxright
              yyleft <- yyright
              xxright <- xxtmp
              yyright <- yytmp
            }
            xxs[j] <- xx[j] * (yyleft * xxright - yyright * 
                                 xxleft)/(yy[j] * (xxright - xxleft) - xx[j] * 
                                            (yyright - yyleft))
            yys[j] <- (yy[j]/xx[j]) * xxs[j]
          }
          else {
            xxs[j] <- 0
            yys[j] <- 0
          }
        }
        else {
          xxs[j] <- xx[j] * CGap/(seg + CGap) + xx[j] * 
            (df[i, j] - df[2, j])/(df[1, j] - df[2, j]) * 
            seg/(seg + CGap)
          yys[j] <- yy[j] * CGap/(seg + CGap) + yy[j] * 
            (df[i, j] - df[2, j])/(df[1, j] - df[2, j]) * 
            seg/(seg + CGap)
        }
      }
      if (is.null(pdensities)) {
        polygon(xxs, yys, lty = pltys[i - 2], lwd = plwds[i - 
                                                            2], border = pcols[i - 2], col = pfcols[i - 
                                                                                                      2])
      }
      else {
        polygon(xxs, yys, lty = pltys[i - 2], lwd = plwds[i - 
                                                            2], border = pcols[i - 2], density = pdensities[i - 
                                                                                                              2], angle = pangles[i - 2], col = pfcols[i - 
                                                                                                                                                         2])
      }
      points(xx * scale, yy * scale, pch = ptys[i - 2], 
             col = pcols[i - 2])
    }
  }
}
