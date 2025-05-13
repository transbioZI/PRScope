#Filter by correlation
#14.05.2024 Joonas Naamanka

library(data.table)
library(tidyverse)



###############################################
#Arguments


dir <- args[1]  #Directory to the prs data
file <- args[2] #Name of the prs data file
output_name <- args[3] #Name of the output file
cor_thres <- args[4] #The threshold for the correlation filtering, a value between 0 and 1
iter <- args[5] #Number of iterations, e.g. 10000


PRS_prsice <- fread(paste0(dir, file), data.table = F)


############################################
#PRS variables
prs <- names(PRS_prsice)


prs_data <- PRS_prsice[,prs]
#First remove one from each pair with a correlation of higher than cor_thres
#The number of variables remaining depends on the choice of which
#variable of a pair of highly correlated variables gets removed.
#Therefore, let's find the solution that maximizes the number of variables
#over iter iterations

find_max <- list()
set.seed(342781793)
seeds <- sample(1:100000000, iter)


cors <- cor(prs_data)


for(j in 1:length(seeds)) {
  find_max[[j]] <- vector()
  set.seed(seeds[j]) 
  new_seeds <- sample(1:100000000, length(prs))
for(i in 1:length(prs)) {
problems <- which(abs(cors[,i]) > cor_thres)
problems <- problems[!(problems %in% find_max[[j]]) & problems != i]

if(length(problems) > 1) {
rid <- i
find_max[[j]][length(find_max[[j]])+1] <- rid
} else if (length(problems) == 1) {
random_choice <- c(i, problems)
set.seed(new_seeds[i]) 
rid <- sample(random_choice, 1)
find_max[[j]][length(find_max[[j]])+1] <- rid
}
}
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

prs_data <- prs_data[, -remove]

#Output

fwrite(prs_data, paste0(dir, output_name), sep = "\t")


