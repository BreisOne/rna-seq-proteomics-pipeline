#Load the libraries
libraries <- c('rstatix','org.Hs.eg.db','AnnotationDbi','RColorBrewer', 'pheatmap','ggVennDiagram', 'tidyverse','scales','ggrepel','readxl', 'VennDiagram', 'ggvenn','ggpubr')
lapply(libraries,library, character.only = TRUE)

#Helper function###
save_plot <- function(plot,name, w, h){
  png(file=paste0("./results/figures/",name), width=w, height=h, units = "in", res = 300)
  plot
}

##### Load proteome results ####

Proteome_BJ <- read_excel("./results/dataframes/TableS2_Proteome_Results.xlsx")

Whole_proteome_gene_list <- unlist(Proteome_BJ[, "gene_names"])

Transcriptome_BJ <- read_excel("./results/dataframes/TableS1_Transcriptome_Results.xlsx")

Whole_transcriptome_gene_list <- unlist(Transcriptome_BJ[!is.na(Transcriptome_BJ$symbol), "symbol"])

#### Overlapping between genes (transcriptome) and protein-associated genes (proteome)

x_2 <- list(Transcriptome= Whole_transcriptome_gene_list,
            Proteome= Whole_proteome_gene_list)

overlap <- calculate.overlap(x_2)

whole_overlap_plot <-ggVennDiagram(x_2,label_alpha = 0,color = "black")+
                                   scale_fill_gradient(low="#F4FAFE",high = "#4981BF")

save_plot(whole_overlap_plot,"RNA_seq_proteome_overlap_plot.png",7,7)
dev.off()
#### Overlapping between DEGs in transcriptome and proteome

DEG_overlap_transcriptome <- unlist(unique(Transcriptome_BJ[!is.na(Transcriptome_BJ$symbol)&Transcriptome_BJ$padj<0.05, "symbol"]))
DEG_overlap_proteome <- unlist(unique(Proteome_BJ[!is.na(Proteome_BJ$gene_names)&Proteome_BJ$ALMS1KO_vs_Control_p.adj<0.05, "gene_names"]))

Overlap_sig <- list(Transcriptome = DEG_overlap_transcriptome,
                    Proteome = DEG_overlap_proteome)

Overlap_sig_cal <- calculate.overlap(Overlap_sig)
 
transcriptome.genes.commons <- Transcriptome_BJ[Transcriptome_BJ$symbol %in% Overlap_sig_cal$a3 & Transcriptome_BJ$padj<0.05,]
proteome.genes.commons <- Proteome_BJ[Proteome_BJ$gene_names %in% Overlap_sig_cal$a3 & Proteome_BJ$ALMS1KO_vs_Control_p.adj<0.05,]

Mix.genes.commons <- transcriptome.genes.commons[,c(8,3)]
colnames(Mix.genes.commons) <- c("Gene_names", "Transcriptome_Log2FC")

Mix2.genes.commons <- proteome.genes.commons[,c(8,18)]
colnames(Mix2.genes.commons) <- c("Gene_names","Proteome_Log2FC")

Mix3.genes.commons <- merge(Mix.genes.commons, Mix2.genes.commons, by = "Gene_names" , all = FALSE)
mix.final.genes.commons <- Mix3.genes.commons[,c(2,3)]

row.names(mix.final.genes.commons) <- Mix3.genes.commons$Gene_names
colnames(mix.final.genes.commons) <- c("Trasncriptome", "Proteome")

### Scatter plot of Log2 FC of common genes in transcriptome and proteome

expression_plot <- ggplot(mix.final.genes.commons, aes(x=`Proteome`, y=`Trasncriptome`)) +
                              geom_point(color="#4657AD")+
                              geom_smooth(method=lm, color="#4657AD", fill="#4657AD")+
                              geom_hline(yintercept = 0, linetype="dashed", color = "darkred")+
                              geom_vline(xintercept = 0, linetype="dashed", color = "darkred")+
                              xlim(-6,6)+
                              ylim(-6.5,6.5)+
                              xlab("Proteome Log2 Fold Change")+
                              ylab("Trasnciptome Log2 Fold Change")+
                              ggtitle("Correlation Transcriptome vs proteome")+
                              stat_cor(method = "pearson", label.x = -5, label.y = 4)+
                              theme_bw()+
                              theme(
                                plot.title = element_text(hjust = 0.5)
                              )

save_plot(expression_plot,"Cor_plot_genes.png",7,7)
dev.off()

#### Overlapping between DEG in transcriptome and proteome

overlapping_genes_venn <- ggVennDiagram(Overlap_sig,label_alpha = 0,color = "black")+
                                        scale_fill_gradient(low="#F4FAFE",high = "#4981BF")

save_plot(overlapping_genes_venn,"overlap_DEG_plot.png",7,7)
dev.off()

#### Fisher's exact test Overlapping genes adj p-value < 0.05 ####

x <- seq(0, 100, by = 1) #X axis sequence

m <- 158 #DE genes in proteome
n <- 1712 #DE genes in rna-seq
k <- 75 #DE genes overlapping

density <- dhyper(x, m, n, k)
density_frame <- data.frame(x, density)
colnames(density_frame)[1] <- "Genes"

Common_DEGs_fisher_plot<- ggplot(density_frame, aes(x=Genes, y=density))+
                                  geom_line( size = 1, color = "#4657AD")+
                                  geom_vline(xintercept = 75,size = 1, linetype="dashed", color = "darkred")+
                                  xlab("Genes")+
                                  ylab("Density")+
                                  ggtitle("Fisher's exact test Overlapping genes adj p-value < 0.05")+
                                  theme_classic()+
                                  theme(
                                    plot.title = element_text(hjust = 0.5))

save_plot(Common_DEGs_fisher_plot,"Fishers_exact_test_overlapping_genes.png",7,7)
dev.off()