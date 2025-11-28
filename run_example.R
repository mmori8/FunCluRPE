source("FunCluRPE.R")
source("data_generation.R")

# generating the data, two replications
data_scenario1 <- scenario1(replications = 2)
coefi <- data_scenario1$coefficients
labels <- data_scenario1$true_labels

# calling the function
output <- lapply(coefi, FunCluRPE, G=2)

# confusion matrices of final clusterings 
mapply(function(X,Y) {table(X$final_cluster_output, Y)}, X=output, Y=labels, SIMPLIFY = FALSE)
