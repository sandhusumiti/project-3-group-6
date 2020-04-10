# This code was written to run limma on microarray results on
# group 6 chemicals: 3-methylcholanthrene, fluconazole, and pirinixic acid

# Load packages and set working directory
library(limma)
setwd("~/Documents/Applications 2/Project_3")

# Read in dataframe with array_id and chemical columns for group 6
samples <- read.csv('group_6_mic_info.csv',as.is=TRUE)

# Read in the full RMA normalized matrix of all experiments
rma <- read.table('liver-normalization-rma.txt',
  sep='\t',
  as.is=TRUE,
  header=TRUE,
  row.names=1,
)

# Subset the full expression matrix to just the chemicals in this comparison 
# (3-methylcholanthrene, fluconazole, and pirinixic acid)
rma.subset <- rma[paste0('X',samples$array_id)]
# Subset of only 3-methylcholanthrene and control
rma.3methyl <- rma.subset[paste0('X',samples$array_id[samples$chemical=='3-METHYLCHOLANTHRENE' | (samples$chemical=='Control' & samples$vehicle=='CMC_.5_%')])]
# Subset of only fluconazole and control
rma.fluco <- rma.subset[paste0('X',samples$array_id[samples$chemical=='FLUCONAZOLE' | (samples$chemical=='Control' & samples$vehicle=='CORN_OIL_100_%')])]
# Subset of only pirinixic acid and control
rma.piri <- rma.subset[paste0('X',samples$array_id[samples$chemical=='PIRINIXIC_ACID' | (samples$chemical=='Control' & samples$vehicle=='CMC_.5_%')])]

# Construct a design matrix modeling treatment vs control for use by limma
# First matrix, 3-methylcholanthrene
design_3methyl <- model.matrix(
  ~factor(
    samples$chemical[samples$vehicle=='CMC_.5_%'],
    levels=c('Control','3-METHYLCHOLANTHRENE')
  )
)
colnames(design_3methyl) <- c('Intercept','3-METHYLCHOLANTHRENE') # Add column names
# Second matrix, fluconazole
design_fluco <- model.matrix(
  ~factor(
    samples$chemical[samples$vehicle=='CORN_OIL_100_%'],
    levels=c('Control','FLUCONAZOLE')
  )
)
colnames(design_fluco) <- c('Intercept','FLUCONAZOLE') # Add column names
# Third matrix, pirinixic acid
design_piri <- model.matrix(
  ~factor(
    samples$chemical[samples$vehicle=='CMC_.5_%'],
    levels=c('Control','PIRINIXIC_ACID')
  )
)
colnames(design_piri) <- c('Intercept','PIRINIXIC_ACID') # Add column names

# Run limma in three steps:
# 1 - run lmFit to fit the data to the design matrix
# 2 - run eBayes to get statistics from microarray
# 3 - run topTable to get the top differentially expressed genes

# 3-methylcholanthrene
fit_3methyl <- lmFit(rma.3methyl, design_3methyl)
fit_3methyl <- eBayes(fit_3methyl)
t <- topTable(fit_3methyl, coef=2, n=nrow(rma.3methyl), adjust='BH')
sig <- subset(t, adj.P.Val < 0.05)
print(nrow(sig)) # 58 significant genes
write.csv(t,'updated_3methylcholanthrene_limma_results.csv') # write out the results to file
# fluconazole
fit_fluco <- lmFit(rma.fluco, design_fluco)
fit_fluco <- eBayes(fit_fluco)
t_fluco <- topTable(fit_fluco, coef=2, n=nrow(rma.fluco), adjust='BH')
sig <- subset(t_fluco, adj.P.Val < 0.05)
print(nrow(sig)) # 8761 significant genes
write.csv(t_fluco,'updated_fluconazole_limma_results.csv') # write out the results to file
# pirinixic acid
fit_piri <- lmFit(rma.piri, design_piri)
fit_piri <- eBayes(fit_piri)
t_piri <- topTable(fit_piri, coef=2, n=nrow(rma.piri), adjust='BH')
sig <- subset(t_piri, adj.P.Val < 0.05)
print(nrow(sig)) # 1997 significant genes
write.csv(t_piri,'updated_pirinixic_acid_limma_results.csv') # write out the results to file
