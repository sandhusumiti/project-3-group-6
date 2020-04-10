library(ggplot2)
library(ggrepel)
library(DESeq2)
library(ggpubr)# Create scatter plots of fold change vs nominal p value


c1<-read.table("SRR1177963.txt",header=TRUE)
head(c1)
#write.csv(c1,'SRR1177963.csv')

c2<-read.table("SRR1177964.txt",header=TRUE)
head(c2)
#write.csv(c2,'SRR1177964.csv')

c3<-read.table("SRR1177965.txt",header=TRUE)
head(c3)
#write.csv(c3,'SRR1177965.csv')

c4<-read.table("SRR1177997.txt",header=TRUE)
head(c4)
#write.csv(c4,'SRR1177997.csv')

c5<-read.table("SRR1177999.txt",header=TRUE)
head(c5)
#write.csv(c5,'SRR1177999.csv')

c6<-read.table("SRR1178002.txt",header=TRUE)
head(c6)
#write.csv(c6,'SRR1178002.csv')

c7<-read.table("SRR1178014.txt",header=TRUE)
head(c7)
#write.csv(c7,'SRR1178014.csv')

c8<-read.table("SRR1178021.txt",header=TRUE)
head(c8)
#write.csv(c8,'SRR1178021.csv')

c9<-read.table("SRR1178047.txt",header=TRUE)
head(c9)
#write.csv(c9,'SRR1178047.csv')


counts<-merge(c1,c2, by=c("Geneid"))
counts<-merge(counts,c3, by=c("Geneid"))
counts<-merge(counts,c4, by=c("Geneid"))
counts<-merge(counts,c5, by=c("Geneid"))
counts<-merge(counts,c6, by=c("Geneid"))
counts<-merge(counts,c7, by=c("Geneid"))
counts<-merge(counts,c8, by=c("Geneid"))
counts<-merge(counts,c9, by=c("Geneid"))

counts<- counts[c(1,7,13,19,25,31,37,43,49,55)]
box_data <- counts[-c(1)]

box_data<- subset(box_data,rowSums(final_data==0)==0)

box_data$SRR1177963 <- log(box_data$SRR1177963)
box_data$SRR1177964 <- log(box_data$SRR1177964)
box_data$SRR1177965 <- log(box_data$SRR1177965)
box_data$SRR1177997 <- log(box_data$SRR1177997)
box_data$SRR1177999 <- log(box_data$SRR1177999)
box_data$SRR1178002 <- log(box_data$SRR1178002)
box_data$SRR1178014 <- log(box_data$SRR1178014)
box_data$SRR1178021 <- log(box_data$SRR1178021)
box_data$SRR1178047 <- log(box_data$SRR1178047)


boxplot(box_data,las=2,ylab="log(counts)",
        ylim=c(-1,23),
        par(mar = c(9, 4, 0, 2)+ 0.1),
        col = c("lightgreen","lightgreen","lightgreen",
                "lightblue","lightblue","lightblue",
                "pink","pink","pink"),
        at = c(1,2,3,4,5,6,7,8,9))
        #names = c("Pirinixic_acid_1","Pirinixic_acid_2","Pirinixic_acid_3",
         #         "3-methylcholanthrene_1","3-methylcholanthrene_2","3-methylcholanthrene_3",
          #        "Fluconazole_1","Fluconazole_2","Fluconazole_3"))

mtext("Samples", side = 1, line = 7, cex = 1, font = 2)
legend(1, 22, legend=c("Pirinixic_acid_1", "3-methylcholanthrene_2","Fluconazole_1"),
       fill=c("lightgreen", "lightblue","pink"), cex=0.8)
#write.csv(counts,"counts.csv")


colnames(counts) <- c("Geneid","SRR1177963","SRR1177964",
                    "SRR1177965","SRR1177997",
                    "SRR1177999","SRR1178002",
                    "SRR1178014","SRR1178021","SRR1178047")

write.csv(counts,"counts.csv")



#data2<-read.csv("first_file.csv",header=TRUE)
#data2<-data2[c(2,8)]
#head(data2)

controls<-read.csv("control_counts.csv",header=TRUE)

controls<-controls[c(1,5:7,22:24)]
head(controls)

final_data <- merge(counts,controls,by=c("Geneid"))

#write.csv(final_data,'group_6_rna_counts.csv')

#installing required packages
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager",lib="/usr4/bf527/sandhu08")
#BiocManager::install("DESeq2",lib="/usr4/bf527/sandhu08")




#mode AhR
#load counts
#cnts<- read.csv('group_6_rna_counts.csv',row.names=1)

#filter out rows that have any zeros for funzies
mode1 <- final_data[c(1,5:7,11:13)]
write.csv(mode1,"cnts_mode1.csv")
cnts_mode1<-read.csv('cnts_mode1.csv',row.names=2)
cnts_mode1<-cnts_mode1[-c(1)]
#write.csv(cnts_mode1,'cnts.csv')#,row.names=1)

cnts_mode1<- subset(cnts_mode1,rowSums(final_data==0)==0)
#cnts <- cnts[c(-1)]

#sample information
info_mode1<- read.csv('group_6_rna_info.csv')
info_mode1 <- info_mode1[c(1:3,10:12),]
#info <- info[c(1:4,13:15),]
#create the DeSeq object

dds_mode1<-DESeqDataSetFromMatrix(
  countData=cnts_mode1,
colData=info_mode1,
  design=~mode_of_action
)

#relevel mode_of_action as factor
dds_mode1$mode_of_action <- relevel(dds_mode1$mode_of_action,ref='Control')
#info$mode_of_action <- relevel(info$mode_of_action, ref = "Control")

#run DESeq
dds_mode1 <- DESeq(dds_mode1)
res_mode1 <- results(dds_mode1, contrast = c('mode_of_action','AhR','Control'))
res_mode1 <- lfcShrink(dds_mode1, coef=2)

res_mode1 <- res_mode1[order(res_mode1$padj),]

#write out DE results
write.csv(res_mode1,'deseq_mode1.csv')

#write out matrix of normalized counts

dds_mode1 <- estimateSizeFactors(dds_mode1)
write.csv(counts(dds_mode1,normalized=TRUE),'deseq_norm_counts_mode1.csv')

sig_genes1 <- res_mode1[which(res_mode1$padj<0.05),]
write.csv(sig_genes1,'sig_genes1.csv')



#plot(sig_genes1$log2FoldChange, sig_genes1$pvalue)
plot_data1<-read.csv("sig_genes1.csv", head=TRUE)



plot_data1 <- plot_data1[!is.na(plot_data1$pvalue),]

plot_data1 <-plot_data1[order(plot_data1$log2FoldChange),]
plot_data1["log_pval"]<-(-log10(plot_data1$pvalue))
#labels1= rbind(head(plot_data1),tail(plot_data1))
#labels1["gene"]<-c("Mme","Lox","Slc13a3","Cyp3a9",
                   #"Hspb1","Adh7","Ugt2b1","Sult2a6",
                   #"Ddhd1","Oat","Ugt1a7c","Cyp1a2")

down_pval1<-plot_data1[which(plot_data1$pvalue<0.05 & plot_data1$log2FoldChange<(-1.5)),]
down_pval1["gene"]<-c("Mme")
up_pval1<-plot_data1[which(plot_data1$pvalue<0.05 & plot_data1$log2FoldChange>1.5),]
up_pval1["gene"]<-c("Ugt1a7c","Cyp1a2")
                   
vol_m<- ggplot(plot_data1,aes(x=log2FoldChange,y=log_pval))+
  geom_point()+xlim(-7,7)+
  geom_vline(xintercept = (-1.5), linetype="dashed") +
  geom_vline(xintercept = 1.5, linetype="dashed") +
  #scale_y_continuous(breaks = round(seq(min(-log(plot_data1$pvalue)), max(-log(plot_data1$pvalue)), by = 50),1))+
  geom_hline(yintercept = 0.05)+
  geom_text_repel(data=down_pval1,aes(label=gene))+
  geom_text_repel(data=up_pval1,aes(label=gene))+
  
  labs(title="3-methylcholanthrene", x="Log Fold Change", y="Nominal P-value",color = "Genes")+

  geom_point(data=down_pval1, mapping=aes(x = log2FoldChange, y = log_pval, color = "mediumblue"))+
  geom_point(data=up_pval1,mapping=aes(x = log2FoldChange, y = log_pval, color = "green4"))+
  scale_color_manual(values = c("green4","mediumblue"),
                     labels = c("Up-regulated","Down-regulated")) +
  scale_shape_identity()


#mode CAR/PXR
#load counts
#cnts<- read.csv('group_6_rna_counts.csv',row.names=1)

#filter out rows that have any zeros for funzies
mode2 <- final_data[c(1,8:10,14:16)]
write.csv(mode2,"cnts_mode2.csv")
cnts_mode2<-read.csv('cnts_mode2.csv',row.names=2)
cnts_mode2<-cnts_mode2[-c(1)]
cnts_mode2<- subset(cnts_mode2,rowSums(final_data==0)==0)
#cnts <- cnts[c(-1)]

#sample information
info_mode2<- read.csv('group_6_rna_info.csv')
info_mode2 <- info_mode2[c(4:6,13:15),]
#info <- info[c(1:4,13:15),]
#create the DeSeq object

dds_mode2<-DESeqDataSetFromMatrix(
  countData=cnts_mode2,
  colData=info_mode2,
  design=~mode_of_action
)

#relevel mode_of_action as factor
dds_mode2$mode_of_action <- relevel(dds_mode2$mode_of_action,ref='Control')
#info$mode_of_action <- relevel(info$mode_of_action, ref = "Control")

#run DESeq
dds_mode2 <- DESeq(dds_mode2)
res_mode2 <- results(dds_mode2, contrast = c('mode_of_action','CAR/PXR','Control'))
res_mode2 <- lfcShrink(dds_mode2, coef=2)
res_mode2 <- res_mode2[order(res_mode2$padj),]

#write out DE results
write.csv(res_mode2,'deseq_mode2.csv')

#write out matrix of normalized counts

dds_mode2 <- estimateSizeFactors(dds_mode2)
write.csv(counts(dds_mode2,normalized=TRUE),'deseq_norm_counts_mode2.csv')

sig_genes2 <- res_mode2[which(res_mode2$padj<0.05),]
write.csv(sig_genes2,'sig_genes2.csv')



#plot(sig_genes2$log2FoldChange, sig_genes2$pvalue)
plot_data2<-read.csv("sig_genes2.csv", head=TRUE)
plot_data2 <- plot_data2[!is.na(plot_data2$pvalue),]

plot_data2 <-plot_data2[order(plot_data2$log2FoldChange),]
plot_data2["log_pval"]<-(-log10(plot_data2$pvalue))
labels2= rbind(head(plot_data2),tail(plot_data2))
labels2["gene"]<-c("Stac3","Nrep","St8sia1","Aldh1b1",
                  "Slc13a4","Atp2b2","Lifr","Orm1",
                  "Abcc3","Cyp3a23/3a1","Slc5a1","Cited4")

down_pval2<-plot_data2[which(plot_data2$pvalue<0.05 & plot_data2$log2FoldChange<(-1.5)),]
up_pval2<-plot_data2[which(plot_data2$pvalue<0.05 & plot_data2$log2FoldChange>1.5),]

vol_f<- ggplot(plot_data2,aes(x=log2FoldChange,y=log_pval))+
  geom_point()+xlim(-7,7)+
  geom_vline(xintercept = (-1.5), linetype="dashed") +
  geom_vline(xintercept = 1.5, linetype="dashed") +
  #scale_y_continuous(breaks = round(seq(min(-log(plot_data2$pvalue)), max(-log(plot_data2$pvalue)), by = 50),1))+
  geom_hline(yintercept = 0.05,linetype="dashed")+
  geom_text_repel(data=labels2,aes(label=gene))+
  geom_point(data=down_pval2, mapping=aes(x = log2FoldChange, y = log_pval, color = "mediumblue"))+
  geom_point(data=up_pval2,mapping=aes(x = log2FoldChange, y = log_pval, color = "green4"))+
  labs(title="Fluconazole", x="Log Fold Change", y="Nominal P-value",color = "Genes")+
  scale_color_manual(values = c("green4","mediumblue"),
                     labels = c("Up-regulated","Down-regulated")) +
  scale_shape_identity()

#mode PPARA
#load counts
#cnts<- read.csv('group_6_rna_counts.csv',row.names=1)

#filter out rows that have any zeros for funzies
mode3 <- final_data[c(1:4,11:13)]
write.csv(mode3,"cnts_mode3.csv")
cnts_mode3<-read.csv('cnts_mode3.csv',row.names=2)
cnts_mode3<-cnts_mode3[-c(1)]
cnts_mode3<- subset(cnts_mode3,rowSums(final_data==0)==0)
#cnts <- cnts[c(-1)]

#sample information
info_mode3<- read.csv('group_6_rna_info.csv')
info_mode3 <- info_mode3[c(7:12),]
#info <- info[c(1:4,13:15),]
#create the DeSeq object

dds_mode3<-DESeqDataSetFromMatrix(
  countData=cnts_mode3,
  colData=info_mode3,
  design=~mode_of_action
)

#relevel mode_of_action as factor
dds_mode3$mode_of_action <- relevel(dds_mode3$mode_of_action,ref='Control')
#info$mode_of_action <- relevel(info$mode_of_action, ref = "Control")

#run DESeq
dds_mode3 <- DESeq(dds_mode3)
res_mode3 <- results(dds_mode3, contrast = c('mode_of_action','PPARA','Control'))
res_mode3 <- lfcShrink(dds_mode3, coef=2)
res_mode3 <- res_mode3[order(res_mode3$padj),]

#write out DE results
write.csv(res_mode3,'deseq_mode3.csv')

#write out matrix of normalized counts

dds_mode3 <- estimateSizeFactors(dds_mode3)
write.csv(counts(dds_mode3,normalized=TRUE),'deseq_norm_counts_mode3.csv')

sig_genes3 <- res_mode3[which(res_mode3$padj<0.05),]
write.csv(sig_genes3,'sig_genes3.csv')





plot_data3<-read.csv("sig_genes3.csv", head=TRUE)
#plot(plot_data3$log2FoldChange, plot_data3$pvalue)
plot_data3 <- plot_data3[!is.na(plot_data3$pvalue),]

plot_data3 <-plot_data3[order(plot_data3$log2FoldChange),]
plot_data3["log10_p_val"]<-(-log10(plot_data3$pvalue))
labels= rbind(head(plot_data3),tail(plot_data3))
labels["gene"]<-c("Apoa4","Cyp2c7","Sult2a1","Dhrs7l1",
                  "Cyp2c11","Cyp3a9","Hdc","LOC292543",
                  "Ugt1a2","Ptprn2","Aqp7","Fabp3")

down_pval<-plot_data3[which(plot_data3$pvalue<0.05 & plot_data3$log2FoldChange<(-1.5)),]
up_pval<-plot_data3[which(plot_data3$pvalue<0.05 & plot_data3$log2FoldChange>1.5),]

vol_p <- ggplot(plot_data3,aes(x=log2FoldChange,y=log10_p_val))+
  geom_point()+xlim(-7,7)+
  #scale_y_continuous(breaks = round(seq(min(plot_data3$log10_p_val), max(plot_data3$log10_p_val)), by = 50),1)+
  geom_vline(xintercept = (-1.5), linetype="dashed") +
  geom_vline(xintercept = 1.5, linetype="dashed") +
  geom_hline(yintercept = 0.05, linetype="dashed")+
  geom_text_repel(data=labels,aes(label=gene),direction="both")+
  geom_point(data=down_pval, mapping=aes(x = log2FoldChange, y = log10_p_val, color = "mediumblue"))+
    

  geom_point(data=up_pval,mapping=aes(x = log2FoldChange, y = log10_p_val, color = "green4"))+
  labs(title="Pirinixic acid", x="Log Fold Change", y="Nominal P-value",color = "Genes")+
  scale_color_manual(values = c("green4","mediumblue"),
                     labels = c("Up-regulated","Down-regulated")) +
  scale_shape_identity()


#deseq_data <- read.csv("example_deseq_results.csv")

#data sorted by adj p-val
#deseq_data <- deseq_data[order(deseq_data$padj),]
#dim(deseq_data)

#genes with padj <0.05
#sig_genes <- deseq_data[which(deseq_data$padj<0.05),]
#dim(sig_genes)
#hist(sig_genes$log2FoldChange)

#norm_data <- read.csv("example_deseq_norm_counts.csv")

#plot(sig_genes$log2FoldChange, sig_genes$pvalue)

#multiqc data: /projectnb/bf528/users/group6/project_2/analysis/counts/multiqc_data

#ggplot(sig_methyl, aes(x=logFC, y=P.Value)) +
  #geom_point() + xlim(-7,7) +
  #geom_hline(yintercept = 0.05, linetype="dashed") +
  #geom_vline(xintercept = min(meth_up$logFC), linetype="dashed") +
  #geom_vline(xintercept = max(meth_down$logFC), linetype="dashed") +
  #geom_point(data=meth_down, col="mediumblue") + 
  #geom_point(data=meth_up, col="green4")

# 3-methylcholanthrene
m <- qplot(plot_data1$log2FoldChange,
      geom="histogram",
      fill=I("lightblue"), 
      col=I("black"), 
      main = "3-methylcholanthrene", 
      xlab = "Log Fold Change",
      ylab="Count")

# fluconazole

f <- qplot(plot_data2$log2FoldChange,
                geom="histogram",
                fill=I("pink"), 
                col=I("black"), 
                main = "Fluconazole",xlim = c(-4,4), 
                xlab = "Log Fold Change",
                ylab="Count")

# pirinixic acid
p <- qplot(plot_data3$log2FoldChange,
           geom="histogram",
           fill=I("lightgreen"), 
           col=I("black"), 
           main = "Pirinixic acid", 
           xlab = "Log Fold Change",
           ylab="Count")

# Display scatter plots
ggarrange(m, f, p, 
          labels = c("A", "B", "C"),
          ncol = 3, nrow = 1)

ggarrange(vol_m, vol_f, vol_p, 
          labels = c("A", "B", "C"),
          ncol = 3, nrow = 1,
          common.legend = TRUE)
