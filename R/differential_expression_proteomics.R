library(readr)
library(readxl)
library(ggplot2)
library(tidyverse)

#Load the librarys
libraries <- c('RColorBrewer', 'pheatmap', 'tidyverse','scales','ggrepel','readxl','DEP')
lapply(libraries,library, character.only = TRUE)

#Importamos los datos del excel del software de thermo y seleccionamos las columnas de interes.
DPGDVP340140_3525_proteins <- read_excel("C:/Users/Brise/OneDrive/Tesis/KO BJ-5TA/Proteoma TGF-B/Resultados/DPGDVP340140_3525 proteins.xlsx")

data.matrix <- DPGDVP340140_3525_proteins[,c(4:5,7,9,24:29)]
colnames(data.matrix)[5:10] <- c("Control_1","Control_2","Control_3", "ALMS1KO_1","ALMS1KO_2","ALMS1KO_3")
view(data.matrix)
#Exportamos los numeros de acceso para anotar los gene names

write.csv(data.matrix$Accession,file="UNIPROT_IDS",row.names = FALSE)

#Impostarmos el la lista de genes mapeados con el Protein ID de Uniprot.

MappedProteins <- read_delim("C:/Users/Brise/OneDrive/Tesis/KO BJ-5TA/Proteoma TGF-B/Resultados/MappedProteins.txt", 
                             "\t", escape_double = FALSE, trim_ws = TRUE)
View(MappedProteins)

# Hacemos un join para mezclar el protein.ID con el gene.name  en el dataframe original

data.joined <- left_join(data.matrix, MappedProteins, by = c("Accession"= "From"))
data.joined <- data.joined[,c(11, 1:10)]
colnames(data.joined)[1:2] <-c("Gene_names", "protein_IDs")

#DEP: Differential Enrichment analysis of Proteomics data.
# Generamos el data frame del diseño experimental
install.packages("BiocManager")
BiocManager::install("DEP")
library(DEP)

# Check for duplicates in the gene names.
data.joined$Gene.names %>% duplicated() %>% any()

# See exactly which gene names are duplicated.
data.joined %>% group_by(gene.names) %>% summarize(frequency = n()) %>% 
  arrange(desc(frequency)) %>% filter(frequency > 1)

#Generate unique identifiers for each Gene/Protein
data_unique <- make_unique(data.joined, "Gene.names", "protein.IDs", delim = ";")

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

# Plot a heatmap of all significant proteins with the data centered per protein
heatmap_proteome <- plot_heatmap(dep, type = "centered", col_limit = 3, show_row_names = FALSE,
             indicate = c("condition", "replicate"), kmeans = FALSE,
             k = 6)

#Heatmap Plot
heat_colors <- colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100)

colours_heatmap_annotation = list(
  condition = c(Control = "#6395ED", Knockout = "#FF4500"),
  replicates = c("1" = "#FACBB3", "2" = "#B3CDE3", "3" = "#CCEBC5"))

colnames(Proteome_BJ)[c(2:7)] <- c("C1","C2","C3","KO1","KO2","KO3")

#Create dendrogram for samples
mat_cluster_cols <- hclust(dist(t(Proteome_BJ[!is.na(Proteome_BJ$gene_names)& Proteome_BJ$`ALMS1KO_vs_Control_p,adj`<0.05,c(2:7)])))
plot(mat_cluster_cols, main = "Unsorted Dendrogram", xlab = "", sub = "")

#Reverse de order in the dendrogram for samples
sort_hclust <- function(...) as.hclust(rev(as.dendrogram(...)))
mat_cluster_cols <- sort_hclust(mat_cluster_cols)
plot(mat_cluster_cols, main = "Sorted Dendrogram", xlab = "", sub = "")


heat <- pheatmap(Proteome_BJ[!is.na(Proteome_BJ$gene_names)& Proteome_BJ$`ALMS1KO_vs_Control_p,adj`<0.05,c(2:7)],
                 color = heat_colors,
                 cluster_rows = TRUE,
                 cluster_cols = mat_cluster_cols,
                 show_rownames = FALSE,
                 legend_breaks = c(-1.4,-0.75,0,0.75,1.4),
                 annotation_col = select(BJ_Metadata, replicates, condition),
                 annotation_colors = colours_heatmap_annotation,
                 scale = "row",
                 main = "Proteome BJ-5TA TGF-B 24h"
)
heat

# Plot a volcano plot for the contrast "ALMS1KO_vs_Control""
plot_volcano(dep, contrast = "ALMS1KO_vs_Control", label_size = 1, add_names = FALSE)

# Plot a barplot for USP15 and IKBKG
plot_single(dep, proteins = c("TGFBI","S100A6"))

plot_single(dep, proteins = "TGFBI", type = "centered")

# Plot a frequency plot of significant proteins for the different conditions
plot_cond(dep)

# Generate a results table
data_results <- get_results(dep)
view(data_results)
write.csv(data_results, file="data_total_BJ-5TAKO.csv",row.names = FALSE,)
# Number of significant proteins
data_significant <- data_results %>% filter(significant)
view(data_significant)
write.csv(data_significant, file="data_significant2_BJ-5TAKO.csv",row.names = FALSE,)
# Column names of the results table
colnames(data_results)
# Generate a wide data.frame
df_wide <- get_df_wide(dep)
write.csv(df_wide, file="df_wide_BJ-5TAKO.csv",row.names = FALSE,)
# Generate a long data.frame
df_long <- get_df_long(dep)
write.csv(df_long, file="df_long_BJ-5TAKO.csv",row.names = FALSE,)
# Save analyzed data
save(data_se, data_norm, data_imp, data_diff, dep, file = "data.RData")

####Volcanoplot####
colnames(data_results)[c(1,3,4,7)]<- c("gene.names", "pvalue", "padj","log2FC")
# Generate logical column 
data_results <- data_results %>% mutate(threshold_padj = padj < 0.05,Overexpressed = log2FC >= 1.5,Underexpressed = log2FC <= -1.5)

#Order by padj
data_results <- data_results[order(data_results$padj),]
data_results <- data_results[!(data_results$gene.names==data_results$ID),]

# Create volcano plot
ggplot(data_results, aes(x = log2FC, y = -log10(padj), color = interaction(threshold_padj, Overexpressed,Underexpressed), label = gene.names)) + 
  geom_point() +
  scale_colour_manual(values = c("gray","gray","gray","#DC3220","gray", "#005AB5"))+
  geom_label_repel(data = head(data_results,40), aes(label=gene.names), max.overlaps = Inf) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "#AA4499")+
  geom_vline(xintercept = 1.5, linetype="dashed", color = "#AA4499")+
  geom_vline(xintercept = -1.5, linetype="dashed", color = "#AA4499")+
  annotate(geom="text", x=1.3, y=-6, label="-Log(0.05)", color="red")+
  scale_x_continuous(limits = c(-7, 7), breaks = seq(from = -7.5, to = 7.5, by = 1.5))+
  scale_y_continuous(limits = c(0,14), breaks = seq(from = 0, to = 14, by= 2))+
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") + 
  theme_bw()+
  theme(legend.position = "none", 
        plot.title = element_text(size = rel(1.5), hjust = 0.5), 
        axis.title = element_text(size = rel(1.25)))
