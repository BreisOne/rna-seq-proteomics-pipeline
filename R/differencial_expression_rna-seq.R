# Load libraries for DESeq2

libraries <- c('ggpubr', 'dendsort','EnsDb.Hsapiens.v75','marray','apeglm','VennDiagram','clusterProfiler','GOplot','org.Hs.eg.db','apeglm','DESeq2', 'RColorBrewer', 'pheatmap', 'tidyverse','scales','ggrepel')
lapply(libraries,library, character.only = TRUE)

#Create counmatrix and metadata dataframe
BJ_RNA_seq <- df[3:60625,2:7]
rownames(BJ_RNA_seq) <- df[3:60625,1]
label <- c("C1","C2", "C3","KO1","KO2","KO3")
replicates <- c("1","2","3","1","2","3")
condition <- factor(c("Control","Control", "Control","Knockout", "Knockout","Knockout"))
colnames(BJ_RNA_seq)<- label
BJ_Metadata <- data.frame(condition, replicates)
rownames(BJ_Metadata) <- label

#Match metada and countmatrix and reorder the columns of the countmatrix
reorder_idx <- match(rownames(BJ_Metadata), colnames(BJ_RNA_seq))
BJ_RNA_seq_ro <- BJ_RNA_seq[,reorder_idx]

#Create the DESeq2 object
BJ_DESeq <- DESeqDataSetFromMatrix(countData = BJ_RNA_seq_ro,
                                    colData = BJ_Metadata,
                                    design = ~ condition)

####Quality control of samples####

#Determine the size factors to use for normalization
BJ_DESeq <- estimateSizeFactors(BJ_DESeq)

#Extract the normalized counts
BJ_normalized <- counts(BJ_DESeq, normalized = TRUE)

#Transform the normalized counts 
BJ_DESeq_VS <- vst(BJ_DESeq, blind=TRUE)
BJ_DESeq_rlog <- rlog(BJ_DESeq, blind=FALSE)

#Extract the matrix of transformed counts
vsd_mat_BJ <- assay(BJ_DESeq_VS)
rlog_mat_BJ <- assay(BJ_DESeq_rlog)
#Compute the correlation values between samples
vsd_cor_BJ <- cor(vsd_mat_BJ) 

# Plot the heatmap correlation
pheatmap(vsd_cor_BJ, annotation = select(BJ_Metadata, condition))

# Plot the PCA of PC1 and PC2
BJ_PCA <- plotPCA(BJ_DESeq_VS, intgroup=c("condition", "replicates"), returnData=TRUE)
percentVar <- round(100 * attr(BJ_PCA, "percentVar"))
ggplot(BJ_PCA, aes(PC1, PC2, color=condition, shape=replicates)) +
  geom_point(size=5) +
  ggtitle("PCA plot - RNA-seq DEG genes")+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+
  theme_bw()+
  theme(
    axis.title=element_text(size=12,face="bold")
  )

####Differential expression analysis####
BJ_Analysis <- DESeq(BJ_DESeq)

#Plot dispersions
plotDispEsts(BJ_Analysis)

# Extract the results of the differential expression analysis
BJ_res <- results(BJ_Analysis, 
                     contrast = c("condition","Knockout","Control"),
                     alpha = 0.05,
                    lfcThreshold = 0)
# Shrink the log2 fold change estimates to be more accurate
BJ_res <- lfcShrink(BJ_Analysis, 
                    contrast =  c("condition","Knockout", "Control"),
                    res = BJ_res,
                    type = "ashr")

# Save results as a data frame and Short by padj
BJ_res_all <- data.frame(BJ_res)
BJ_res_all <- BJ_res_all[order(BJ_res_all$padj),]

# Subset the results to only return the significant genes with p-adjusted values less than 0.05
BJ_res_sig <- subset(BJ_res_all, padj < 0.05)
BJ_res_sig <- subset(BJ_res_sig, abs(log2FoldChange)>=1.5)

####Anotacion genica####
BJ_res_sig$symbol <- mapIds(org.Hs.eg.db,
                            keys=rownames(BJ_res_sig),
                            column="SYMBOL",
                            keytype="ENSEMBL")

BJ_res_sig$TXBIOTYPE <- mapIds(EnsDb.Hsapiens.v75,
                            keys=rownames(BJ_res_sig),
                            column="TXBIOTYPE",
                            keytype="GENEID")

BJ_res_sig_biotype <- BJ_res_sig[BJ_res_sig$TXBIOTYPE=="protein_coding",]
BJ_res_sig_biotype <- BJ_res_sig_biotype[!is.na(BJ_res_sig_biotype$TXBIOTYPE),]

BJ_res_all$symbol <- mapIds(org.Hs.eg.db,
                            keys=rownames(BJ_res_all),
                            column="SYMBOL",
                            keytype="ENSEMBL")

BJ_res_all$TXBIOTYPE <- mapIds(EnsDb.Hsapiens.v75,
                               keys=rownames(BJ_res_all),
                               column="TXBIOTYPE",
                               keytype="GENEID")

BJ_res_all <- BJ_res_all[complete.cases(BJ_res_all),]
BJ_res_all_biotype <- BJ_res_all[!is.na(BJ_res_all$TXBIOTYPE)&BJ_res_all$TXBIOTYPE=="protein_coding",]

####Visualisation of results####
# Create MA plot
plotMA(BJ_res)

# Generate logical column 
BJ_res_all_biotype <- BJ_res_all_biotype %>% mutate(threshold_padj = padj < 0.05,Overexpressed = log2FoldChange >= 1.5,Underexpressed = log2FoldChange <= -1.5)
BJ_res_all_biotype <- BJ_res_all_biotype[order(BJ_res_all_biotype$padj),]

# Create the volcano plot
ggplot(BJ_res_all_biotype, aes(x = log2FoldChange, y = -log10(padj), color = interaction(threshold_padj,Overexpressed,Underexpressed), label = symbol)) + 
  geom_point() +
  scale_colour_manual(values = c("gray","gray","gray","#DC3220","gray","#005AB5"))+
  geom_label_repel(data = head(BJ_res_all_biotype,40), aes(label=symbol), max.overlaps = Inf) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "#AA4499")+
  geom_vline(xintercept = 1.5, linetype="dashed", color = "#AA4499")+
  geom_vline(xintercept = -1.5, linetype="dashed", color = "#AA4499")+
  scale_x_continuous(limits = c(-7, 7), breaks = seq(from = -7.5, to = 7.5, by = 1.5))+
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") + 
  theme_bw()+
  theme(legend.position = "none", 
        plot.title = element_text(size = rel(1.5), hjust = 0.5), 
        axis.title = element_text(size = rel(1.25)))
