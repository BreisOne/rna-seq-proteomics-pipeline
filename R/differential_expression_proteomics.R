
# Load the libraries
libraries <- c('RColorBrewer','readr', 'pheatmap', 'tidyverse','scales','ggrepel','ggplot2','readxl','DEP')
lapply(libraries,library, character.only = TRUE)

# Helper function
save_plot <- function(plot,name, w, h){
  pdf(file=paste0("./results/figures/",name), width=w, height=h)
  plot
}

####GENE NAMES ANOTATION####
# Import the data from the Thermo software Excel file and select the columns of interest.

LFQ_proteomics <- read_excel("./data/proteomics.xlsx")

data.matrix <- LFQ_proteomics[,c(4:5,7,9,24:29)]
colnames(data.matrix)[5:10] <- c("Control_1","Control_2","Control_3", "ALMS1KO_1","ALMS1KO_2","ALMS1KO_3")
view(data.matrix)

# Export the accession numbers to annotate the gene names.The annotation was made with the batch search on the UNIPROT website.

write.csv(data.matrix$Accession,file="UNIPROT_IDS",row.names = FALSE)

#Import list of mapped genes with UNIPROT proteinID

MappedProteins <- read_delim("./data/MappedProteins.txt", 
                             "\t", escape_double = FALSE, trim_ws = TRUE)

# Perform a join to merge the protein.ID with the gene.name in the original dataframe

data.joined <- left_join(data.matrix, MappedProteins, by = c("Accession"= "From"))
data.joined <- data.joined[,c(11, 1:10)]
colnames(data.joined)[1:2] <-c("Gene_names", "protein_IDs")

####DEP: Differential Expression analysis of Proteomics data.####
#Generate the experimental design dataframe

#####Data preprocessing#####
# Check for duplicates in the gene names.
data.joined$Gene_names %>% duplicated() %>% any()

# See exactly which gene names are duplicated.
data.joined %>% group_by(Gene_names) %>% summarize(frequency = n()) %>% 
  arrange(desc(frequency)) %>% filter(frequency > 1)

#Generate unique identifiers for each Gene/Protein
data_unique <- make_unique(data.joined, "Gene_names", "protein_IDs", delim = ";")

#Validate that all names are unique.
data_unique$name %>% duplicated() %>% any()

# Generate the dataframe of the experimental design.
l <- colnames(data_unique)[6:11]
c <- c("Control","Control","Control", "ALMS1KO","ALMS1KO","ALMS1KO")
r <- c("1","2","3","1","2","3")
colnames <- c("label", "condition","replicate")
exp_design <- data.frame(l,c,r, stringsAsFactors = FALSE)
colnames(exp_design) <- colnames

# Generate a SummarisedExperiment object from the experimental design
exp_Cols <- c(6:11)

#Generate the SummarizedExperiment object.
data_se <- make_se(data_unique,exp_Cols, exp_design)

#####Explotarory data analysis and normalization#####
# Plot a barplot of the protein identification overlap between samples
plot_frequency(data_se)

# Filter for proteins that are identified in all replicates of at least one condition
data_filt <- filter_missval(data_se, thr = 1)

# Plot a barplot of the number of identified proteins per samples
plot_numbers(data_filt)

# Normalize the data
data_norm <- normalize_vsn(data_filt)

# Visualize normalization by boxplots for all samples before and after normalization
plot_normalization(data_filt, data_norm)

# Plot a heatmap of proteins with missing values
plot_missval(data_norm)

#####Imputation and DE#####
# Impute missing data using random draws from a Gaussian distribution centered around a minimal value (for MNAR)
data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)

# Plot intensity distributions before and after imputation
plot_imputation(data_norm, data_imp)

# Differential enrichment analysis  based on linear models and empherical Bayes statistics

# Test every sample versus control
data_diff <- test_diff(data_imp, type = "control", control = "Control")

# Denote significant proteins based on user defined cutoffs
dep <- add_rejections(data_diff, alpha = 0.05, lfc=0)

# Plot the first and second principal components
plot_pca(dep, x = 1, y = 2, n = 500, point_size = 4)

# Plot the Pearson correlation matrix
plot_cor(dep, significant = TRUE, lower = 0, upper = 1, pal = "Reds")

#####Vizualization of DE analysis results#####

# Plot a heatmap of all significant proteins with the data centered per protein
heatmap_proteome <- plot_heatmap(dep, type = "centered", col_limit = 3, show_row_names = FALSE,
             indicate = c("condition", "replicate"), kmeans = FALSE,
             k = 6)

# Plot a volcano plot for the contrast "ALMS1KO_vs_Control""
plot_volcano(dep, contrast = "ALMS1KO_vs_Control", label_size = 1, add_names = FALSE)

# Plot a barplot for USP15 and IKBKG
plot_single(dep, proteins = c("TGFBI","S100A6"))

plot_single(dep, proteins = "TGFBI", type = "centered")

# Plot a frequency plot of significant proteins for the different conditions
plot_cond(dep)

#####Results Exportation#####
# Generate a results table
data_results <- get_results(dep)
# Generate a wide data.frame
Proteome_BJ <- get_df_wide(dep)
write.csv2(Proteome_BJ, file="./results/dataframes/TableS2_Proteome_Results.csv",row.names = FALSE,)

####PAPER FIGURES####

#####Volcanoplot#####
colnames(data_results)[c(1,3,4,7)]<- c("gene.names", "pvalue", "padj","log2FC")
# Generate logical column 
data_results <- data_results %>% mutate(threshold_padj = padj < 0.05,Overexpressed = log2FC >= 0,Underexpressed = log2FC <= 0)

#Order by padj
data_results <- data_results[order(data_results$padj),]
data_results <- data_results[!(data_results$gene.names==data_results$ID),]

# Create volcano plot
Proteome_volcano_plot <- ggplot(data_results, aes(x = log2FC, y = -log10(padj), color = interaction(threshold_padj, Overexpressed,Underexpressed), label = gene.names)) + 
                                geom_point() +
                                scale_colour_manual(values = c("gray","#DC3220","gray", "#005AB5"))+
                                geom_label_repel(data = head(data_results,40), aes(label=gene.names), max.overlaps = Inf) +
                                geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "#AA4499")+
                                annotate(geom="text", x=1.3, y=-6, label="-Log(0.05)", color="red")+
                                scale_x_continuous(limits = c(-7, 7), breaks = seq(from = -7.5, to = 7.5, by = 1.5))+
                                scale_y_continuous(limits = c(0,14), breaks = seq(from = 0, to = 14, by= 2))+
                                xlab("log2 fold change") + 
                                ylab("-log10 adjusted p-value") + 
                                theme_bw()+
                                theme(legend.position = "none", 
                                      plot.title = element_text(size = rel(1.5), hjust = 0.5), 
                                      axis.title = element_text(size = rel(1.25)))

Proteome_volcano_plot

save_plot(Proteome_volcano_plot,"Fig1D_proteome_volcano_plot.pdf",9,5)
dev.off()

#####Heatmap Plot#####
colnames(Proteome_BJ)[c(2:7)] <- c("C1","C2","C3","KO1","KO2","KO3")

heat_colors <- colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100)

colours_heatmap_annotation = list(Condition = c(Control = "#6395ED", Knockout = "#FF4500"),
                                  Samples = c("1" = "#FACBB3", "2" = "#B3CDE3", "3" = "#CCEBC5"))

#Create new metadata
# Create metadata dataframe
label <- c("C1","C2", "C3","KO1","KO2","KO3")
Samples <- c("1","2","3","1","2","3")
Condition <- factor(c("Control","Control", "Control","Knockout", "Knockout","Knockout"))
colnames(BJ_RNA_seq)<- label
BJ_Metadata <- data.frame(Condition, Samples)
rownames(BJ_Metadata) <- label

#Create dendrogram for samples
mat_cluster_cols <- hclust(dist(t(Proteome_BJ[!is.na(Proteome_BJ$gene.names)& Proteome_BJ$ALMS1KO_vs_Control_p.adj<0.05,c(2:7)])))
plot(mat_cluster_cols, main = "Unsorted Dendrogram", xlab = "", sub = "")

proteome_heatmap <- pheatmap(Proteome_BJ[!is.na(Proteome_BJ$gene.names)& Proteome_BJ$ALMS1KO_vs_Control_p.adj<0.05,c(2:7)],
                             color = heat_colors,
                             cluster_rows = TRUE,
                             cluster_cols = mat_cluster_cols,
                             show_rownames = FALSE,
                             legend_breaks = c(-1.4,-0.75,0,0.75,1.4),
                             annotation_col = select(BJ_Metadata, Samples, Condition),
                             annotation_colors = colours_heatmap_annotation,
                             scale = "row",
                             main = "Proteome BJ-5TA TGF-B 24h")
proteome_heatmap

save_plot(proteome_heatmap,"Fig4A_proteome_heatmap.pdf",9,5)
dev.off()