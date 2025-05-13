print_study_list = function(columns) {

  studies = unique(sapply(columns, function(x) {
    str_split(x,"_",3)[[1]][3]
  }))

  complete_results = data.frame(matrix(ncol = 4, nrow = length(studies)))

  colnames(complete_results) <- c('study_id', 'reported_trait', 'efo_trait', 'publication_title')

  complete_results$study_id = studies

  complete_results$efo_trait = sapply(studies, function(x) {
    if(grepl("GCST", x, fixed = TRUE)) {
      get_studies(study_id = x)@studies$reported_trait
    }  else {
      x
    }

  })

  complete_results$reported_trait = sapply(studies, function(x) {
    if(grepl("GCST", x, fixed = TRUE)) {
      get_studies(study_id = x)@studies$reported_trait
    }  else  {
      x
    }
  })

  complete_results$publication_title = sapply(studies, function(x) {
    if(grepl("GCST", x, fixed = TRUE)) {
      get_studies(study_id = x)@publications$title
    } else{
      x
    }
  })

  return(complete_results)
}

find_corr = function(data,thr) {
  cor_heatmap = round(cor(data),3)

  melted_cormat <- reshape2::melt(cor_heatmap)

  melted_cormat$V1 = sapply(str_split(melted_cormat$Var1,"_",3),function(x) {
    x[3]
  })
  melted_cormat$V2 = sapply(str_split(melted_cormat$Var2,"_",3),function(x) {
    x[3]
  })
  melted_cormat <- melted_cormat[melted_cormat$value >= thr | melted_cormat$value <= -1*thr,]
  melted_cormat <- melted_cormat[melted_cormat$V1 != melted_cormat$V2,]

  return(melted_cormat)
}

calc_NagelkerkeR2 = function(data,group) {
  n2 = c()
  for(i in c(1:ncol(data))) {
    m = glm((group - 1)~ data[,i], family = binomial())
    n2 = c(n2,NagelkerkeR2(m)$R2)
  }

  return(n2)
}

calc_NagelkerkeR2_vector = function(vec,group) {
  m = glm((group - 1)~ vec, family = binomial())
  return(NagelkerkeR2(m)$R2)
}

sort_by_correlation = function(data,group, n = 10){
  studies = unique(sapply(colnames(data), function(x) {
    str_split(x,"_",n=3)[[1]][3]
  }))
  cor_val = c()
  ll = list()
  for(st in c(1:length(studies))) {
    sub_data = data[,grepl(studies[st],colnames(data)), drop=FALSE]
    r2 = abs(cor(sub_data, group))
    sub_data = sub_data[,order(r2, decreasing = T),drop = F]
    ll[[st]] = colnames(sub_data)[order(r2, decreasing = T)][1:n]
    cor_val = c(cor_val,max(r2))
  }

  lx = do.call(rbind,ll)
  return(lx[order(cor_val,decreasing = T),])
}

sort_by_correlation_pca = function(data,group, n = 10){
  studies = unique(sapply(colnames(data), function(x) {
    str_split(x,"_",n=3)[[1]][3]
  }))
  cor_val = c()
  ll = list()
  for(st in c(1:length(studies))) {
    sub_data = data[,grepl(studies[st],colnames(data)), drop=FALSE]
    pc_res = prcomp(sub_data)
    sd_res = sort(pc_res$sdev^2/sum(pc_res$sdev^2), decreasing = T)
    total = 0
    count = 0
    for(sdr in sd_res) {
      if(total < 0.9) {
        total = total + sdr
        count = count +1
      }
    }

    m = glm((group - 1)~ pc_res$x[,1], family = binomial())
    mm = max(-log10(coef(summary(m))[,4][-1]))
    #mm = NagelkerkeR2(m)$R2
    r2 = abs(cor(sub_data, group))
    sub_data = sub_data[,order(r2, decreasing = T),drop = F]
    ll[[st]] = colnames(sub_data)[order(r2, decreasing = T)][1:n]
    cor_val = c(cor_val,mm)
  }

  lx = do.call(rbind,ll)
  return(lx[order(cor_val,decreasing = T),])
}

get_first_cor = function(DATA,features, group, n) {
  pv_sel_orig = cor(DATA[,features], as.numeric(group) )
  return(names(pv_sel_orig[order(sort( abs(pv_sel_orig), decreasing=TRUE)),][1:n]))
}

filter_cor = function(prs_data,cor_thres,iter) {

  prs = colnames(prs_data)
  #First remove one from each pair with a correlation of higher than cor_thres
  #The number of variables remaining depends on the choice of which
  #variable of a pair of highly correlated variables gets removed.
  #Therefore, let's find the solution that maximizes the number of variables
  #over iter iterations
  set.seed(342781793)
  seeds <- sample(1:100000000, iter)

  cors <- cor(prs_data)

  find_max = foreach(j=1:length(seeds)) %dopar% {
    res_vec <- vector()
    set.seed(seeds[j])
    new_seeds <- sample(1:100000000, length(prs))
    for(i in 1:length(prs)) {
      problems <- which(abs(cors[,i]) > cor_thres)
      problems <- problems[!(problems %in% res_vec) & problems != i]
      #Note that the diagonal contains 1
        if(length(problems) > 1) {
          rid <- i
          res_vec[length(res_vec)+1] <- rid
        } else if (length(problems) == 1) {
          random_choice <- c(i, problems)
          set.seed(new_seeds[i])
          rid <- sample(random_choice, 1)
          res_vec[length(res_vec)+1] <- rid
        }
    }
    res_vec
  }

  #Find the vector with the shortest length
  remove <- NULL
  for(i in 1:length(find_max)) {
    if(i == 1) {
      remove <- find_max[[i]]
    } else if (length(find_max[[i]]) < length(remove)) {
      remove <- find_max[[i]]
    }
  }

  #Exclude the variables
  if(length(remove) != 0 ) {
    prs_data <- prs_data[, -remove]
  }
  print(paste0("removed prs: ", length(remove)))
  #Output
return(prs_data)
}

check_multi_prs_data_rows = function(data, quantiles, title) {
  plot(density(data[1,]), ylim = c(0,1), main = paste0("Density-",title," (Rows)"), xlim=c(-10,10))
  n = dim(data)[1]-1
    for(i in 1:n){
    lines(density(c(data[i,])))
  }

  quantile_n = quantiles
  res = lapply(c(1:nrow(data)), function(x) {
    quantile(data[x,],probs = seq(0,1,1/quantile_n))

  })
  ls = do.call(rbind, res)
  plot(abs(apply(ls, 2, mean)), pch = 19, type ="b", lty = 2, main = paste0("Mean of quantiles-",title," (Rows)") ,xaxt='n')
  axis(1, at=c(1:(quantile_n+1)),labels=paste0(seq(0,1,1/quantile_n)*100,"%"), las=2)
  ls = scale(ls, center = T, scale = F)
  plot(density(ls[,1]), ylim = c(0,2), main = paste0("Quantile density-",title," (Rows)"), xlim=c(-5,5))
  n = dim(ls)[2]-1
    for(i in 1:n){
    lines(density(c(ls[,i])))
  }
}

check_multi_prs_data_columns = function(data, quantiles, title) {
  plot(density(data[,1]), ylim = c(0,1), main = paste0("Density-",title," (Columns)"), xlim=c(-10,10))
  n = dim(data)[2]-1
    for(i in 1:n){
    lines(density(c(data[,i])))
  }

  quantile_n = quantiles
  res = lapply(c(1:ncol(data)), function(x) {
    quantile(data[,x],probs = seq(0,1,1/quantile_n))

  })
  ls = do.call(rbind, res)
  plot(abs(apply(ls, 2, mean)), pch = 19, type ="b", lty = 2, main = paste0("Mean of quantiles-",title," (Columns)"), xaxt='n')
  axis(1, at=c(1:(quantile_n+1)),labels=paste0(seq(0,1,1/quantile_n)*100,"%"), las=2)
  ls = scale(ls, center = T, scale = F)
  plot(density(ls[,1]), ylim = c(0,2), main = paste0("Quantile density-",title," (Columns)"), xlim=c(-5,5))
  n = dim(ls)[2]-1
    for(i in 1:n){
    lines(density(c(ls[,i])))
  }
}

export_plot = function(plot,path,width = 900, height = 900) {
  ggexport(plot, filename = path, width = width, height = height)
}

create_heatmap = function(data,th) {
  pa = wes_palettes %>% names()
  pal = wes_palette(name = pa[20], n = 3, type="continuous")
  cor_heatmap = round(cor(data),3)
  melted_cormat <- reshape2::melt(cor_heatmap)
  colnames(melted_cormat) = c("A1","A2","value")
  melted_cormat$A1 = sapply(str_split(melted_cormat[,1],"_",4),function(x) {
    paste0(x[3]," : ",x[2])
  })
  melted_cormat$A2 = sapply(str_split(melted_cormat[,2],"_",4),function(x) {
    paste0(x[3]," : ",x[2])
  })
  #melted_cormat[abs(melted_cormat$value) < th ,]$value = 0
  #melted_cormat <- melted_cormat[melted_cormat$A1 != melted_cormat$A2,]

  p = ggplot(data = melted_cormat, aes(A2, A1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                         midpoint = 0, limit = c(-1,1), space = "Lab",
                         name="Pearson\nCorrelation") +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 90, vjust = 1,
                                     size = 8, hjust = 1),
                                     axis.title.x = element_blank(),
                                     axis.title.y = element_blank())+
    coord_fixed()
 return(p)
}

