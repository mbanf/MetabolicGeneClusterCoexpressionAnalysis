# PlantClusterFinder probabilistic high confidence gene cluster prediction algorithm

#This is the probabilistic metabolic gene cluster high confidence prediction algorithm based on co-expression analysis used in Schlapfer et al. to predict high confidence clusters based on co-expressed metabolic enzymes. 


*General usage*

th.pcc_distribution = 0.99 # set the co-expression matrix top threshold. In the paper we only use co-expression relationships between enzymes in the top 0.01 percentile 

# estimate the minimal for co-expression levels each dataset based on 
th.coexpressionLevels_statistical_significance = 0.2 

# threshold for high confidence prediction based on probabilistic ranking. Check different values to 
th.pvalue_clusters = 1e-2

# matrix should be prepared and saved as R matrix in a rds file first  

# build coexpressionmatrix from scratch - rows genes, columns treatments


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

Schlapfer P, Zhang P, Chuan W, Kim T, Chae L, Dreher K, Nilo-Poyanco R, Arvind Chavali, and Rhee SY. (2017) Genome-wide prediction of metabolic enzymes, pathways, and gene clusters in plants. Plant Physiology 173(4):2041-2059.
