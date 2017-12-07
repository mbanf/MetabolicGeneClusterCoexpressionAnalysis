# PlantClusterFinder probabilistic high confidence gene cluster prediction algorithm

#This is the probabilistic metabolic gene cluster high confidence prediction algorithm based on co-expression analysis used in Schlapfer et al. to predict high confidence clusters based on co-expressed metabolic enzymes. 


We performed coexpression analyses for four species (Arabidopsis, soybean, rice, and tomato) based on coexpression data sets provided by ATTED-II (Aoki et al., 2016; Supplemental Table S6). To identify gene clusters with coexpressed enzyme pairs, we defined a genome-wide coexpression threshold, based on Pearson’s correlation coefficient (Usadel et al., 2009), beyond which a gene pair would be considered coexpressed. We used the correlation coefficient of the upper 99th percentile of the coexpression distribution for each species’ coex- pression matrix.
The genome size and the number of experimental samples influence the selection of a coexpression threshold (Usadel et al., 2009). Therefore, we applied Student’s t tests with multiple hypothesis test correction (Usadel et al., 2009) to ensure that all 99th percentile-based Pearson’s correlation coefficient thresholds exceeded minimal coexpression values for mean- ingful statistical significance. Information about the number of genes and experimental samples was acquired from ATTED-II (Aoki et al., 2016; Supplemental Table S6).
For soybean, rice, and tomato, ATTED-II gene identifiers were mapped to Phytozome version 10 gene identifiers. To do this, we retrieved protein sequences from GenBank using RefSeq protein identifiers corresponding to each ENTREZ identifier (via GenBank’s internal mapping file ftp://ftp. ncbi.nih.gov/gene/DATA/gene2refseq.gz). Then, protein sequences were BLASTed against the Phytozome protein sequences, all splicing variants included, using BLAST+ 2.2.28+ (Camacho et al., 2009). All perfect matches (100% sequence identity and end-to-end length coverage) were kept. For a query sequence without a perfect match, the best hit was retained if the alignment identity was equal to or above 99% and the ratio of alignment over query length exceeded 95%. To eliminate uncertainty, all queries with more than one perfect hit were removed to generate the final identifier mapping (Supplemental Table S12). For Arabidopsis, we acquired a mapping table from ATTED-II. We limited our analyses to these four species because their enzyme coverage, with respect to identifier mapping, exceeded 50%. Subsequent analyses were performed only on the mapped genes. We applied the coexpression threshold derived from the original coexpression matrix described above to identify coexpressed gene pairs.
From all metabolic gene cluster predictions, we only retained clusters with at least one coexpressed enzyme pair (Supplemental Fig. S8). To further rank the retained clusters, we assigned a rank for each cluster:

![Alt text](/probabilisticranking.png?raw=true "Ranking threshold")

Here, k denotes the number of coexpressed metabolic gene pairs with a prob- ability pcoex=1 = 0.01 and n 2 k is the number of noncoexpressed metabolic gene pairs, each given a probability of pcoex=0 = 0.99. The rank of a cluster is defined as 1 minus the binomial probability of a metabolic gene cluster’s composition of coexpressed and noncoexpressed enzyme pairs, without considering the order of their occurrence in the cluster. The lower this binomial probability, the less likely a cluster’s composition has occurred by chance, and hence the higher its rank. Therefore, we reward the number of coexpressed metabolic gene pairs per cluster while penalizing the number of metabolic gene pairs without coex- pression. In addition to disregarding the order of coexpressed enzyme pairs, the binomial coefficient is used to normalize the rank with respect to the cluster size (the total number of metabolic gene pairs).



*General usage*

th.pcc_distribution = 0.99 # set the co-expression matrix top threshold. In the paper we only use co-expression relationships between enzymes in the top 0.01 percentile 

# estimate the minimal for co-expression levels each dataset based on 
th.coexpressionLevels_statistical_significance = 0.2 

As in Usadel 2009. 

What is signiﬁcant and what is relevant?The ‘signiﬁcance’ or P value of a co-expression score refersto the likelihood of obtaining the same or a betterco-expression score than the observed one by chance alone.In order to assess this likelihood, correlation scores can beconverted into test values.These can then be compared withstatistical tables or can be converted into P values. Thecorrelation coefﬁcients can be transformed by the following formula to values that approximately follow a t distributionwith n - 2 degrees of freedom:

ss = (r x sqrt( n - 1)) / sqrt(1 - r^2)

where the number of samples is indicated by n, and r is theobserved co-expression score (Fisher 1915). Therefore thelikelihood of observing a given correlation coefﬁcient r canbe calculated from the transformed value. In Excel thefollowing formula could be used to directly convertcorrelation coefﬁcients into P values.

P = TDIST(ABS(R/sqrt((1-r^2) / (n - 2))), (n-2),2)

Nevertheless, to analytically calculate P values severalassumptions have to be met. For Pearson’s correlationanalysis, the data have to be normally distributed for eachgene and bivariate normally distributed for gene pairs,which might not hold true for all gene pairs. Also, thehypothesis of independence of experimental conditionscannot be essentially satisﬁed because the experimentswere conducted for a particular biological purpose, with theresult that, for PCC P value calculations, there are devia-tions from the necessary behaviour of the data. Anotherimportant factor to be considered when using P values isthat when one queries with one gene about 22 000 statisticaltests are made (for the ~22 000 other probe sets on theATH1 array). Clearly using a P value cut-off of 5% wouldyield many genes by chance (22 000 x 0.05, or 1100 to beprecise), thus the P values need to be corrected for multipletesting. This can be achieved by the simple (and very con-servative) Bonferroni correction, which multiplies theresulting P value by the number of tests that are made (e.g.22 000). Still, in cases where many arrays are used tocompute the co-expression scores, a correlation coefﬁcientas little as 0.2 can become very signiﬁcant. Therefore, espe-cially when considering a large number of samples, a signiﬁ-cant correlation might not be of practical importance. Forthis reason, we do not recommend using less conservativeP value correction methods (such as false discoveryrate control), unless these are combined with thresholdsimposed on r2, which assesses the variance in commonbetween two genes in question. This value may be calcu-lated simply by squaring the obtained r value and is indeedimmediately returned by CressExpress. This r2value as ameasure of shared variance varies on a scale of 0 (no sharedvariance) to 1 (100% variance in common). It follows thatfor lower correlation values such as 0.3 the resulting r2would only be 0.09 or 9% shared variance, which might notbe biologically relevant.



# threshold for high confidence prediction based on probabilistic ranking. Check different values to 
th.pvalue_clusters = 1e-2

![Alt text](/ranking.jpeg?raw=true "Ranking threshold")

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
