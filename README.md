# QTLNetwork4BiCount

GLMGEX: A tool for epistasis and gene-environment interaction association analysis in discrete traits

# Package installation

The current GitHub version of QTLNetwork4BiCount can be installed via:
```
library(devtools)
install_github("jojolubi/QTLNetwork4BiCount")
```
If your Rstudio does not have the devtools package, you can download the zip file and install it locally：
```
install.packages("path/to/package_name", repos = NULL, type = "source")
```
snpStats can be downloaded as follows：
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("snpStats")
```
# Example
```
#Load packages ----------------------------------------------------------------

library(QTLNetwork4BiCount)
library(snpStats)
library(qqman)

# Load data --------------------------------------------------------------------

# Read the genotype file
# It is simulation data
# Data must have been charged
# The genotype data is consistent with PhenoGenoMap.Rdata
# Here we only show how to read the plink file
bedfile <- "./inst/test.bed"
bimfile <- "./inst/test.bim"
famfile <- "./inst/test.fam"
genotypeData <- read.plink(bedfile, bimfile, famfile)

# Read in pca data, which can be used as a covariate
# Here we only show how to read in pca files
pca <- read.table("./inst/test.eigenvec")
eigval <- read.table("./inst/test.eigenval")

# Read the file, which is simulation data
# It contains genotype data, phenotype data, and map files
load("./inst/PhenoGenoMap.RData")
G <- snp_Data$genotypes
phe <- snp_Data$phenotype
map <- snp_Data$map

# For proportional data
# Convert the y range from [0, 1] to the real number domain
# Note: the logit function cannot handle sample values of 0 and 1
y <- logit(phe$y_rate)

# Convert genotype matrix SnpMatrix format to standardized matrix form
Ga <- define_Ga(G)

# Convert genotype matrix SnpMatrix format to 0,1,2 matrix form
Ga <- define_Ga(Ga,standardization = F)

# 1D scan --------------------------------------------------------------------

# Take 0-1 binomial distribution as an example
# Gaussian or Poisson distribution, the operation is similar
# When the binomial distribution response variable is the number of successes, weights are required
# When the distribution is negative binomial distribution, theta are required

# In the case of multiple environments
# Calculate p-value
# First input form
y <- phe$y_binomial
cov <- phe[,2:3]
env <- phe$environment
rae <- single.snp.score(phenotype = y,covariates = cov, environment = env,genodata = Ga,map = map,family = "binomial")

# Check how many threads there are
library(parallel)
detectCores()

# Second input form
rae <- single.snp.score(phenotype = "y_binomial",covariates = "cov1+cov2", environment = "environment",data = phe,
                        genodata = Ga,map = map,family = "binomial",parallel = T,num_threads = 4,filename = "gwasp1.txt")

# QQ plot
p1 <- rae$p_values
qqman::qq(p1)

# Manhattan plot
gwas1d <- data.frame(rae$gwas1D[,-4])
gwas1d[,2] <- as.numeric(gwas1d[,2])
gwas1d[,3] <- as.numeric(gwas1d[,3])
gwas1d[,4] <- as.numeric(gwas1d[,4])
qqman::manhattan(gwas1d)

# Screen for significant markers
# The default Bonferroni correction is 0.5/m (m is the number of markers)
rae1 <- single.snp.estimate(phenotype = "y_binomial", covariates = "cov1+cov2", environment = "environment", data = phe,genodata = Ga,map = map,
                            Pvalues = p1, Bonferr = 0.5,family = "binomial")
# Screen for significant markers and estimate parameters
# You can choose forward selection method or lasso method for screening
rae1 <- single.snp.estimate(phenotype = y, covariates = cov, environment = env, genodata = Ga,map = map,Pvalues = p1, Bonferr = 0.5,
                            family = "binomial",selection = T,estimate = T)

# In the case of single environment
# Calculate p value
phe1 <- phe[1:200,]
ra <- single.snp.score(phenotype = "y_binomial",covariates = "cov1+cov2",data = phe1,genodata = Ga,map = map,family = "binomial")
ps1 <- ra$p_values

# Screen for significant markers and estimate parameters
ra1 <- single.snp.estimate(phenotype = "y_binomial",covariates = "cov1+cov2",data = phe1,genodata = Ga,map = map,Pvalues = ps1, Bonferr = 0.5,
                           family = "binomial",lasso = T,estimate = T)

# 2D scan --------------------------------------------------------------------

# In the case of multiple environments
# Significant markers obtained by 1D scanning screening
sig <- rae1[["gwas1D"]][,1]

# Calculate p-value
# SignificantASNP is not required if 1D scan gives no results
raae <- epi.snp.score(significantASNP = sig, phenotype = y, covariates = cov, environment = env,genodata = Ga,map = map,family = "binomial",
                      parallel = T, num_threads = 4,filename = "gwasp2.txt")
p2 <- raae$p_values

# Screen for significant markers and estimate parameters
raae1 <- epi.snp.estimate(significantASNP = sig,phenotype = y, covariates = cov, environment = env, genodata = Ga,map = map,Pvalues = p2, Bonferr = 0.5,
                          family = "binomial",selection = T,estimate = T)

# In the case of single environment
# Significant markers obtained by 1D scanning screening
sigs <- ra1[["gwas1D"]][,1]

#Calculate p-value
raa <- epi.snp.score(significantASNP = sigs, phenotype = "y_binomial",covariates = "cov1+cov2",data = phe1,genodata = Ga,map = map,family = "binomial",
                     parallel = T, num_threads = 4,filename = "gwasp2s.txt")
ps2 <- raa$p_values

#Screen for significant markers and estimate parameters
raa1 <- epi.snp.estimate(significantASNP = sigs,phenotype = "y_binomial",covariates = "cov1+cov2",data = phe1,genodata = Ga,map = map,Pvalues = ps2,
                         Bonferr = 0.5, family = "binomial",lasso = T,estimate = T)
```
