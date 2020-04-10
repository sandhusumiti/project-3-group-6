# This code compares the results from limma analysis and deseq analysis

# Set working directory and load packages
setwd("~/Documents/Applications 2/Project_3")
library(dplyr)
library(ggplot2)
library(ggpubr)

# Read in the data files from limma for each chemical
methyl <- read.csv('limma_results/3methylcholanthrene_limma_results.csv',as.is=TRUE)
fluco <- read.csv('limma_results/fluconazole_limma_results.csv', as.is = TRUE)
piri <- read.csv('limma_results/pirinixic_acid_limma_results.csv', as.is=TRUE)

# Create map data frame from the reference Affy map
map <- read.csv('refseq_affy_map.csv', as.is=TRUE)

# Get the probe id from the limma results for each chemical
methyl_probe <- methyl$X
fluco_probe <- fluco$X
piri_probe <- piri$X

# Find probe id in affy map, add this information to the row in limma results to column "REFSEQID"
# 3-methylcholanthrene
for (row in 1:length(methyl_probe)){
  refseq <- list(map$REFSEQ[which(map$PROBEID == methyl_probe[row])])
  methyl$REFSEQID[row] <- refseq
}
# Fluconazole
for (row in 1:length(fluco_probe)){
  refseq <- list(map$REFSEQ[which(map$PROBEID == fluco_probe[row])])
  fluco$REFSEQID[row] <- refseq
}
# Pirinixic acid
for (row in 1:length(piri_probe)){
  refseq <- list(map$REFSEQ[which(map$PROBEID == piri_probe[row])])
  piri$REFSEQID[row] <- refseq
}

# Expand rows with multiple refseqid and delete rows with no refseqids
# Make empty data frames for each chemical
# 3-methylcholanthrene
maptomethyl <- methyl[FALSE, names(methyl) %in% c("X", "REFSEQID", "logFC", "AveExpr")]   # Subset only probeid and refseqid
as.data.frame(maptomethyl)
# fluconazole
maptofluco <- fluco[FALSE, names(fluco) %in% c("X", "REFSEQID", "logFC", "AveExpr")]   # Subset only probeid and refseqid
as.data.frame(maptofluco)
# pirinixic acid
maptopiri <- piri[FALSE, names(piri) %in% c("X", "REFSEQID", "logFC", "AveExpr")]   # Subset only probeid and refseqid
as.data.frame(maptopiri)

# Fill in the empty data frames
# 3-methylcholanthrene
for (row in 1:nrow(methyl)) {
  if (methyl[row, "P.Value"] < 0.05) {
    if(abs(methyl[row, "logFC"]) > 1.5) {
      methyl_probe<- methyl[row,"X"] # Get probe id
      refseq <- methyl[row, "REFSEQID"] # Get refseqid
      logfc <- methyl[row, "logFC"]
      aveexpr <- methyl[row, "AveExpr"]
      if (length(refseq[[1]])>0) { # if the entry is NOT character(0)
        for (id in refseq[[1]]) { # For each refseqid in the list
          id_row <- data.frame(methyl_probe, id, logfc, aveexpr) # initialize a row with probe and id
          names(id_row)<-c("X","REFSEQID", "logFC", "AveExpr") # Label the row
          maptomethyl <- rbind(maptomethyl, id_row) # add the row to the dataframe
        }
      }
    }
  }
}
# fluconazole
for (row in 1:nrow(fluco)) {
  if (fluco[row, "P.Value"] < 0.05) {
    if(abs(fluco[row, "logFC"]) > 1.5) {
      fluco_probe<- fluco[row,"X"] # Get probe id
      refseq <- fluco[row, "REFSEQID"] # Get refseqid
      logfc <- fluco[row, "logFC"]
      aveexpr <- fluco[row, "AveExpr"]
      if (length(refseq[[1]])>0) { # if the entry is NOT character(0)
        for (id in refseq[[1]]) { # For each refseqid in the list
          id_row <- data.frame(fluco_probe, id, logfc, aveexpr) # initialize a row with probe and id
          names(id_row)<-c("X","REFSEQID", "logFC", "AveExpr") # Label the row
          maptofluco <- rbind(maptofluco, id_row) # add the row to the dataframe
        }
      }
    }
  }
}
# Pirinixic acid
for (row in 1:nrow(piri)) {
  if (piri[row, "P.Value"] < 0.05) {
    if(abs(methyl[row, "logFC"]) > 1.5) {
      piri_probe<- piri[row,"X"] # Get probe id
      refseq <- piri[row, "REFSEQID"] # Get refseqid
      logfc <- piri[row, "logFC"]
      aveexpr <- piri[row, "AveExpr"]
      if (length(refseq[[1]])>0) { # if the entry is NOT character(0)
        for (id in refseq[[1]]) { # For each refseqid in the list
          id_row <- data.frame(piri_probe, id, logfc, aveexpr) # initialize a row with probe and id
          names(id_row)<-c("X","REFSEQID", "logFC", "AveExpr") # Label the row
          maptopiri <- rbind(maptopiri, id_row) # add the row to the dataframe
        }
      }
    }
  }
}

# Collapse the new dataframes based on refseqid, taking the median of the logfoldchange
# 3-methylcholanthrene
methylmap <- maptomethyl %>% group_by(REFSEQID) %>% summarise(logFC=median(logFC), AveExpr=median(AveExpr))
# fluconazole
flucomap <- maptofluco %>% group_by(REFSEQID) %>% summarise(logFC=median(logFC), AveExpr=median(AveExpr))
# pirinixic acid
pirimap <- maptopiri %>% group_by(REFSEQID) %>% summarise(logFC=median(logFC), AveExpr=median(AveExpr))

# Read in the DESeq data for each chemical
d_methyl <- read.csv("deseq_mode1.csv", as.is = TRUE)
d_fluco <- read.csv("deseq_mode2.csv", as.is = TRUE)
d_piri <- read.csv("deseq_mode3.csv", as.is=TRUE)

# Subset only the significantly differentially expressed genes, nominal p-value < 0.05
d_methyl <- subset(d_methyl, pvalue < 0.05)
d_methyl <- subset(d_methyl, abs(log2FoldChange) >1.5)
d_fluco <- subset(d_fluco, pvalue < 0.05)
d_fluco <- subset(d_fluco, abs(log2FoldChange) >1.5)
d_piri <- subset(d_piri, pvalue < 0.05)
d_piri <- subset(d_piri, abs(log2FoldChange) >1.5)

# Calculate concordance for each chemical
# Obtain counts of the number of differentially expressed genes that match between limma and deseq
# 3-methylcholanthrene
methyl_count <- 0 # initialize count
for (row in 1:nrow(d_methyl)) {
  id_deseq <- d_methyl[row, "X"] # get refseqid
  logfc <- d_methyl[row, "log2FoldChange"] # get deseq log2foldchange
  limma_row <- methylmap[which(methylmap$REFSEQID == id_deseq), ] # extract the matching info from limma results
  if (nrow(limma_row) != 0) { # if there is a match between deseq and limma
    if(sign(logfc) == sign(limma_row[["logFC"]])) { # And if the logfoldchange values are both positive or negative
      methyl_count <- methyl_count+1 # Add one to the count
    }
  }
}
# Fluconazole
fluco_count <- 0 # initialize count
for (row in 1:nrow(d_fluco)) {
  id_deseq <- d_fluco[row, "X"] # get refseqid
  logfc <- d_fluco[row, "log2FoldChange"] # get deseq log2foldchange
  limma_row <- flucomap[which(flucomap$REFSEQID == id_deseq), ] # extract the matching info from limma results
  if (nrow(limma_row) != 0) { # if there is a match between deseq and limma
    if(sign(logfc) == sign(limma_row[["logFC"]])) { # And if the logfoldchange values are both positive or negative
      fluco_count <- fluco_count+1 # Add one to the count
    }
  }
}
# Pirinixic acid
piri_count <- 0 # initialize count
for (row in 1:nrow(d_piri)) {
  id_deseq <- d_piri[row, "X"] # get refseqid
  logfc <- d_piri[row, "log2FoldChange"] # get deseq log2foldchange
  limma_row <- pirimap[which(pirimap$REFSEQID == id_deseq), ] # extract the matching info from limma results
  if (nrow(limma_row) != 0) { # if there is a match between deseq and limma
    if(sign(logfc) == sign(limma_row[["logFC"]])) { # And if the logfoldchange values are both positive or negative
      piri_count <- piri_count+1 # Add one to the count
    }
  }
}

# Calculate concordance for each chemical
# 3-methylcholanthrene
methyl_concordance <- (2 * methyl_count)/(nrow(d_methyl) + nrow(methylmap)) # 0.2303143
# fluconazole
fluco_concordance <- (2 * fluco_count)/(nrow(d_fluco) + nrow(flucomap)) # 0.5534543
# pirinixic acid
piri_concordance <- (2 * piri_count)/(nrow(d_piri) + nrow(pirimap)) # 0.5214058

# Make scatter plot with treatment on x-axis and concordance of y-axis (figure 2a)
# make a table of the number of degs and concordance from deseq data
scatter <- data.frame("Concordance"=c(methyl_concordance, fluco_concordance, piri_concordance), 
                      "Treatment"=c(nrow(d_methyl), nrow(d_fluco), nrow(d_piri)))
chems <- c("3-methylcholanthrene", "Fluconazole", "Pirinixic acid")
plot1 <- ggplot(scatter, aes(x=Treatment, y=Concordance)) +
  geom_point(color="black", fill=c('lightblue', 'pink', 'lightgreen'), shape=21, size=3) +
  stat_smooth(method = "lm", linetype="dashed", se=FALSE, color="black") +
  geom_text(label=chems, nudge_x = c(75, 0, -60), nudge_y = c(0, -0.01, 0)) +
  stat_cor(aes(label=paste(..rr.label..))) +
  labs(title="RNA-Seq Analysis", x="Number of differentially expressed genes", y="Concordance")

# Make the same table for limma data
scatter2 <- data.frame("Concordance"=c(methyl_concordance, fluco_concordance, piri_concordance), 
                      "Treatment"=c(nrow(methylmap), nrow(flucomap), nrow(pirimap)))
plot2 <- ggplot(scatter2, aes(x=Treatment, y=Concordance)) +
  geom_point(color="black", fill=c('lightblue', 'pink', 'lightgreen'), shape=21, size=3)+
  geom_text(label=chems, nudge_x = c(12, -8, 5), nudge_y = c(0, 0, 0.007)) +
  stat_smooth(method = "lm", linetype="dashed", se=FALSE, color="black") +
  stat_cor(aes(label=paste(..rr.label..))) +
  labs(title="Microarray Analysis", x="Number of differentially expressed genes", y="Concordance")

# Display scatter plots side by side
ggarrange(plot1, plot2, 
          labels = c("A", "B"),
          ncol = 2, nrow = 1)

# Find median of the base mean from deseq results
# 3-methylcholanthrene
methyl_med <- median(d_methyl$baseMean)
m_deseq_above <- subset(d_methyl, baseMean > methyl_med) # genes above median
m_deseq_below <- subset(d_methyl, baseMean < methyl_med) # genes below median
# fluconazole
fluco_med <- median(d_fluco$baseMean)
f_deseq_above <- subset(d_fluco, baseMean > fluco_med) # genes above median
f_deseq_below <- subset(d_fluco, baseMean < fluco_med) # genes below median
# pirinixic acid
piri_med <- median(d_piri$baseMean)
p_deseq_above <- subset(d_piri, baseMean > piri_med) # genes above median
p_deseq_below <- subset(d_piri, baseMean < piri_med) # genes below median

# find median of the average expression from limma results
# 3-methylcholanthrene
methyl_limmamed <- median(methylmap$AveExpr)
m_limma_above <- subset(methylmap, AveExpr > methyl_limmamed)
m_limma_below <- subset(methylmap, AveExpr < methyl_limmamed)
# fluconazole
fluco_limmamed <- median(flucomap$AveExpr)
f_limma_above <- subset(flucomap, AveExpr > fluco_limmamed)
f_limma_below <- subset(flucomap, AveExpr < fluco_limmamed)
# pirinixic acid
piri_limmamed <- median(pirimap$AveExpr)
p_limma_above <- subset(pirimap, AveExpr > piri_limmamed)
p_limma_below <- subset(pirimap, AveExpr < piri_limmamed)

# Recompute concordance for each of the above and below groups
# concordance for below groups
# 3-methylcholanthrene
mcount <- 0
for (row in 1:nrow(m_deseq_below)) {
  id_deseq <- m_deseq_below[row, "X"]
  logfc <- m_deseq_below[row, "log2FoldChange"]
  below_row <- m_limma_below[which(m_limma_below$REFSEQID == id_deseq), ]
  if (nrow(below_row) != 0) {
    if(sign(logfc) == sign(below_row[["logFC"]])) {
      mcount <- mcount+1
    }
  }
}
m_below_concord <- (2*mcount)/(nrow(m_deseq_below) + nrow(m_limma_below)) # 0.1120797
# fluconazole
fcount <- 0
for (row in 1:nrow(f_deseq_below)) {
  id_deseq <- f_deseq_below[row, "X"]
  logfc <- f_deseq_below[row, "log2FoldChange"]
  below_row <- f_limma_below[which(f_limma_below$REFSEQID == id_deseq), ]
  if (nrow(below_row) != 0) {
    if(sign(logfc) == sign(below_row[["logFC"]])) {
      fcount <- fcount+1
    }
  }
}
f_below_concord <- (2*fcount)/(nrow(f_deseq_below) + nrow(f_limma_below)) # 0.3616364
# pirinixic acid
pcount <- 0
for (row in 1:nrow(p_deseq_below)) {
  id_deseq <- p_deseq_below[row, "X"]
  logfc <- p_deseq_below[row, "log2FoldChange"]
  below_row <- p_limma_below[which(p_limma_below$REFSEQID == id_deseq), ]
  if (nrow(below_row) != 0) {
    if(sign(logfc) == sign(below_row[["logFC"]])) {
      pcount <- pcount+1
    }
  }
}
p_below_concord <- (2*pcount)/(nrow(p_deseq_below) + nrow(p_limma_below)) # 0.3516524

# Concordance for above groups
# 3-methylcholanthrene
macount <- 0
for (row in 1:nrow(m_deseq_above)) {
  id_deseq <- m_deseq_above[row, "X"]
  logfc <- m_deseq_above[row, "log2FoldChange"]
  above_row <- m_limma_above[which(m_limma_above$REFSEQID == id_deseq), ]
  if (nrow(above_row) != 0) {
    if(sign(logfc) == sign(above_row[["logFC"]])) {
      macount <- macount+1
    }
  }
}
m_above_concord <- (2*macount)/(nrow(m_deseq_above) + nrow(m_limma_above)) # 0.2515567
# fluconazole
facount <- 0
for (row in 1:nrow(f_deseq_above)) {
  id_deseq <- f_deseq_above[row, "X"]
  logfc <- f_deseq_above[row, "log2FoldChange"]
  above_row <- f_limma_above[which(f_limma_above$REFSEQID == id_deseq), ]
  if (nrow(above_row) != 0) {
    if(sign(logfc) == sign(above_row[["logFC"]])) {
      facount <- facount+1
    }
  }
}
f_above_concord <- (2*facount)/(nrow(f_deseq_above) + nrow(f_limma_above)) # 0.5152451
# pirinixic acid
pacount <- 0
for (row in 1:nrow(p_deseq_above)) {
  id_deseq <- p_deseq_above[row, "X"]
  logfc <- p_deseq_above[row, "log2FoldChange"]
  above_row <- p_limma_above[which(p_limma_above$REFSEQID == id_deseq), ]
  if (nrow(above_row) != 0) {
    if(sign(logfc) == sign(above_row[["logFC"]])) {
      pacount <- pacount+1
    }
  }
}
p_above_concord <- (2*pacount)/(nrow(p_deseq_above) + nrow(p_limma_above)) # 0.482746

# Make a barplot of all of these values
# make into a table
bar <- data.frame("Concordance"=c(m_above_concord, m_below_concord, methyl_concordance,
                                  f_above_concord, f_below_concord, fluco_concordance,
                                  p_above_concord, p_below_concord, piri_concordance), 
                  "Analysis"=c("Above", "Below", "Overall", "Above", "Below", "Overall", "Above", "Below", "Overall"),
                  "Chemical"=c("3-methylcholanthrene","3-methylcholanthrene","3-methylcholanthrene",
                               "Fluconazole","Fluconazole","Fluconazole",
                               "Pirinixic Acid","Pirinixic Acid","Pirinixic Acid"))
# Plot data with a side-by-side bar chart
ggplot(bar, aes(fill=Analysis, x=Chemical, y=Concordance)) + 
  geom_bar(stat="identity", position="dodge") +
  labs(title="Above and Below Median Concordance")
