# MetabolicGeneClusterCoexpressionAnalysis

#Enhancing gene regulatory network inference through data integration with markov random fields 
Home of the GRACE (Gene Regulatory network inference ACcuracy Enhancement) algorithm 


*Summary* 


General usage 


coexpression matrix

MGCC assumes 



Specific parameters and datasets 
n.cpus <- 2 
beta <- 1 # equal emphasis on precision (atrm) and recall (bp coreg pairs) - introducing (novel) combined f-measure (physical + functional evidence, singularity handling, minimum size) 
b.normalize_precision <- TRUE 
b.jaccard = TRUE 
n.sets <- 20 
lambda.gridSearch <- c(0.01,seq(0.5,2.5,0.5)) 
th.percentage.hyperparameters = 0.99 # grid search 
max.call = 200 # simulated annealing 
n.models <- 100 # 100 models 
n.sample_size <- 0.632 # as in traditional bootrapping approaches 



cite:
