# This code uses the output from the run_limma.R code to create figures from the data

# Set the working directory and load packages
setwd("~/Documents/Applications 2/Project_3/limma_results")
library("ggplot2")
library(ggpubr)
library(dplyr)
library(ggrepel)

# Create histograms of fold change values from the significant DE genes from each of your analyses
# load limma data
# 3-methylcholanthrene
methyl <- read.csv('updated_3methylcholanthrene_limma_results.csv',as.is=TRUE)
# fluconazole
fluco <- read.csv('updated_fluconazole_limma_results.csv', as.is=TRUE)
# pirinixic acid
piri <- read.csv('updated_pirinixic_acid_limma_results.csv', as.is=TRUE)

# subset the data to get the signficiant DE genes (adjusted p-value <0.05) 
# 3-methylcholanthrene
sig_methyl <- subset(methyl, adj.P.Val < 0.05)
# fluconazole
sig_fluco <- subset(fluco, adj.P.Val < 0.05)
# pirinixic acid
sig_piri <- subset(piri, adj.P.Val < 0.05)

# Make the histogram of the log fold change
# 3-methylcholanthrene
m <- ggplot(sig_methyl, aes(x=logFC)) + 
  geom_histogram(color='black', fill='lightblue', binwidth=0.1) +
  labs(title="3-methylcholanthrene", x="Log Fold Change", y="Count")
# fluconazole
f <- ggplot(sig_fluco, aes(x=logFC)) + 
  geom_histogram(color='black', fill='pink', binwidth=0.17) +
  labs(title="Fluconazole", x="Log Fold Change", y="Count")
# pirinixic acid
p <- ggplot(sig_piri, aes(x=logFC)) + 
  geom_histogram(color='black', fill='lightgreen', binwidth=0.25) +
  labs(title="Pirinixic acid", x="Log Fold Change", y="Count")

# Display histograms from each chemical
ggarrange(m, f, p, 
          labels = c("A", "B", "C"),
          ncol = 3, nrow = 1)


# Create volcano plots of fold change vs nominal p-value
# Calculate the log10 p-value
sig_methyl$P.Value <- log10(sig_methyl$P.Value) * -1
sig_fluco$P.Value <- log10(sig_fluco$P.Value) * -1
sig_piri$P.Value <- log10(sig_piri$P.Value) * -1
# Read in affymap
map <- read.csv('../refseq_affy_map.csv', as.is=TRUE)
map <- map %>% distinct(PROBEID, .keep_all=TRUE)
map <- rbind(map, c('none', '1373778_at', '1373778_at', 'none'))

# Get the top ten most up- and down-regulated genes in each chemical
# 3-methylcholanthrene
meth_up <- top_n(sig_methyl, 10, logFC)
meth_down <- top_n(sig_methyl, -10, logFC)
# Get the gene names for each probe
labels <- rbind(meth_up, meth_down)
gene <- which(map$PROBEID %in% labels$X)
for (x in 1:nrow(labels)) {
  for (probe in gene) {
    if (labels$X[x] == map$PROBEID[probe]) {
      labels$gene[x] <- map$SYMBOL[probe]
    }
  }
}
# Make chart
pth <- log10(0.05)*-1
m <- ggplot(sig_methyl, aes(x=logFC, y=P.Value)) +
  geom_point() + xlim(-7,7) +
  geom_hline(yintercept = pth, linetype="dashed") +
  geom_text_repel(data=labels,aes(label=gene),direction="both")+
  geom_vline(xintercept = 1.5, linetype="dashed") +
  geom_vline(xintercept = -1.5, linetype="dashed") +
  geom_point(data=subset(sig_methyl, logFC < -1.5), col="mediumblue") + 
  geom_point(data=subset(sig_methyl, logFC > 1.5), col="green4") +
  labs(title="3-methylcholanthrene", x="Log Fold Change", y="-log10(P-Value)")

# Fluconazole
fluco_up <- top_n(sig_fluco, 10, logFC)
fluco_down <- top_n(sig_fluco, -10, logFC)
# Get the gene names for each probe
labels <- rbind(fluco_up, fluco_down)
gene <- which(map$PROBEID %in% labels$X)
for (x in 1:nrow(labels)) {
  for (probe in gene) {
    if (labels$X[x] == map$PROBEID[probe]) {
      labels$gene[x] <- map$SYMBOL[probe]
    }
  }
}
# Make chart
f <- ggplot(sig_fluco, aes(x=logFC, y=P.Value)) +
  geom_point() + xlim(-7,7) +
  geom_hline(yintercept = pth, linetype="dashed") +
  geom_text_repel(data=labels,aes(label=gene),direction="both")+
  geom_vline(xintercept = 1.5, linetype="dashed") +
  geom_vline(xintercept = -1.5, linetype="dashed") +
  geom_point(data=subset(sig_fluco, logFC < -1.5), col="mediumblue") + 
  geom_point(data=subset(sig_fluco, logFC > 1.5), col="green4") +
  labs(title="Fluconazole", x="Log Fold Change", y="-log10(P-Value)")
# pirinixic acid
piri_up <- top_n(sig_piri, 10, logFC)
piri_down <- top_n(sig_piri, -10, logFC)
# Get the gene names for each probe
labels <- rbind(piri_up, piri_down)
gene <- which(map$PROBEID %in% labels$X)
for (x in 1:nrow(labels)) {
  for (probe in gene) {
    if (labels$X[x] == map$PROBEID[probe]) {
      labels$gene[x] <- map$SYMBOL[probe]
    }
  }
}
# Make chart
p <- ggplot(sig_piri, aes(x=logFC, y=P.Value)) +
  geom_point() + xlim(-7,7) +
  geom_hline(yintercept = pth, linetype="dashed") +
  geom_text_repel(data=labels,aes(label=gene),direction="both")+
  geom_vline(xintercept = 1.5, linetype="dashed") +
  geom_vline(xintercept = -1.5, linetype="dashed") +
  geom_point(data=subset(sig_piri, logFC < -1.5), col="mediumblue") + 
  geom_point(data=subset(sig_piri, logFC > 1.5), col="green4") +
  labs(title="Pirinixic acid", x="Log Fold Change", y="-log10(P-Value)")

# Display the volcano plots from each chemical
ggarrange(m, f, p, 
          labels = c("A", "B", "C"),
          ncol = 3, nrow = 1)

''' Old scatter plots that were not volcano plots
# Create scatter plots of fold change vs nominal p value
# 3-methylcholanthrene
m <- ggplot(sig_methyl, aes(x=P.Value, y=logFC)) +
  geom_point(color="black", fill="lightblue", shape=21) +
  labs(title="3-methylcholanthrene", x="Nominal P-Value", y="Log Fold Change")
# fluconazole
f <- ggplot(sig_fluco, aes(x=P.Value, y=logFC)) +
  geom_point(color="black", fill="pink", shape=21) +
  labs(title="Fluconazole", x="Nominal P-Value", y="Log Fold Change")
# pirinixic acid
p <- ggplot(sig_piri, aes(x=P.Value, y=logFC)) +
  geom_point(color="black", fill="lightgreen", shape=21) +
  labs(title="Pirinixic acid", x="Nominal P-Value", y="Log Fold Change")'''