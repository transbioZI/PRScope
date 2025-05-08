train_and_predict = function(TRAIN, GROUP, CLUSTER_GROUP, TARGET_CLUSTER, FEATURES, VALIDATION, REPEATS=101){

  # MODEL CREATION : BIOTYPE<TARGET_CLUSTER>, CONTROL
  MODEL_FEATURES = train_models(TRAIN,GROUP, CLUSTER_GROUP, TARGET_CLUSTER, FEATURES, REPEATS)
  
  # PREDICTION
  RESULT = predict_models(MODEL_FEATURES, VALIDATION, FEATURES, REPEATS)
  
  return(RESULT)
}

train_models = function(DATA, GROUP, CLUSTER_GROUP, TARGET_CLUSTER, FEATURES, REPEATS = 101, N_TREES = 1000, THREADS = 20) {
  
  HC = which(GROUP == 1)
  CASE = which(GROUP == 2)
  BIOTYPE = CASE[which(CLUSTER_GROUP==TARGET_CLUSTER) ]
  
  stopifnot(length(HC) > length(BIOTYPE))
  
  LABELS = rep(0, length(BIOTYPE) * 2)
  LABELS[1:length(BIOTYPE)] = 1
  
  # MODEL CREATION : BIOTYPE<biotype_number>, CONTROL
  MODEL_FEATURES = foreach(repi=1:REPEATS) %dopar% {
    require(ranger)
    CONTROL = sample(HC, length(BIOTYPE) ) # downsampling
    TRAINING_DATA = data.frame( X = DATA[c(BIOTYPE, CONTROL),  FEATURES],Y = as.factor(LABELS) )
    ranger(dependent.variable.name = "Y", data=TRAINING_DATA, importance="impurity", num.trees = N_TREES, num.threads = THREADS, classification = TRUE)
  }
  
  return(MODEL_FEATURES)
}

predict_models = function(MODEL_FEATURES, VALIDATION, FEATURES, REPEATS = 101) {
  
  # PREDICTION
  PREDICTION_RESULTS_FEATURES_SIZES = numeric()
  for(repi in 1:REPEATS){
    MODEL = MODEL_FEATURES[[repi]]
    VALIDATION_DATA = data.frame( X = VALIDATION[,FEATURES] )
    PREDICTION_RESULTS = predict(MODEL, data=VALIDATION_DATA)$predictions
    PREDICTION_RESULTS_FEATURES_SIZES = cbind(PREDICTION_RESULTS_FEATURES_SIZES, as.numeric(PREDICTION_RESULTS))
  }
  
  # BIOTPYE MEMBERSHIP
  RESULT = list()
  RESULT[[1]] = apply(PREDICTION_RESULTS_FEATURES_SIZES, 1, function(x) {length(x[x==1])/length(x)}) # BIOYPE_1 probability
  RESULT[[2]] = apply(PREDICTION_RESULTS_FEATURES_SIZES, 1, function(x) {length(x[x==2])/length(x)}) # BIOYPE_2 probability
  RESULT[[3]] = ceiling(apply(PREDICTION_RESULTS_FEATURES_SIZES, 1, median)) # assignment
  RESULT[[4]] = MODEL_FEATURES # if you want to save models
  return(RESULT)
}

# z-score transformation using median, 
# perform repeated k means clustering
# identify consensus cluster solution
clustering = function(DATA, GROUP, FEATURES, FEATURE_SIZE, REPEATS = 100, KMEANS_RESULT_SAVE_PATH = NA){
  
  # 2 is cases
  # 1 is HC

  Z_TRANSFORMATION = apply(DATA,2,function(x){ ( x[which(GROUP == 2)] - median(x[which(GROUP == 1)]) ) /sd(x[which(GROUP == 1)])})
  
  #SELECTED_FEATURES = as.vector(FEATURES[1:FEATURE_SIZE,])

  #SELECTED_FEATURES = sort(match(SELECTED_FEATURES[!is.na(SELECTED_FEATURES)],colnames(DATA)))

  DATA_SELECTED_FEATURES = Z_TRANSFORMATION[,FEATURES]

  CLUSTER_SUBJECTS = foreach(repclu=1:REPEATS) %dopar% {
    ksol = flexclust::kcca(DATA_SELECTED_FEATURES, k=2)
    flexclust::clusters(ksol)
  }
  
  CLUSTER_SUBJECTS_MATRIX = matrix(0, length(which(GROUP == 2)), length(which(GROUP == 2)))
  
  for(repclu in c(1:REPEATS)){
    CLUSTER_RESULT = CLUSTER_SUBJECTS[[repclu]]
    CLUSTER_SUBJECTS_MATRIX[which(CLUSTER_RESULT == 1),which(CLUSTER_RESULT == 1)] = CLUSTER_SUBJECTS_MATRIX[which(CLUSTER_RESULT == 1),which(CLUSTER_RESULT == 1)] + 1
    CLUSTER_SUBJECTS_MATRIX[which(CLUSTER_RESULT == 2),which(CLUSTER_RESULT == 2)] = CLUSTER_SUBJECTS_MATRIX[which(CLUSTER_RESULT == 2),which(CLUSTER_RESULT == 2)] + 1
  }

  CLUSTER = kmeans(CLUSTER_SUBJECTS_MATRIX,2)
  CLUSTER_FINAL = rep(0, length(which(GROUP == 2)))
  CLUSTER_FINAL[which(CLUSTER$cluster == 1)] = 1
  CLUSTER_FINAL[which(CLUSTER$cluster == 2)] = 2

  if(!is.na(KMEANS_RESULT_SAVE_PATH)) {
    writeLines(as.character(CLUSTER_FINAL),paste0(KMEANS_RESULT_SAVE_PATH,"/kmeans_feat_",length(FEATURE_SIZE)))
  }
  
  return(CLUSTER_FINAL)
}
