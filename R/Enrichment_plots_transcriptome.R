libraries <- c('readxl', 'RColorBrewer', 'tidyverse')
lapply(libraries,library, character.only = TRUE)

#Helper function
save_plot <- function(plot,name, w, h){
  png(file=paste0("./results/figures/",name), width=w, height=h, units = "in", res = 300)
  plot
}

####GO terms Overlapping genes
GO_terms_transcriptome <- read_excel("./results/dataframes/TableS3_GO_terms_and_pathways_transcriptome_1712_genes.xlsx", 
                                   sheet = "GO_terms")

GO_terms_transcriptome <- GO_terms_transcriptome %>% 
                                mutate(ordering = as.numeric(factor(Category, levels = c("CC","MF","BP")))+ `Adjusted P-value`)  
  
Enrich_plot_go_terms_transcriptome <- ggplot(GO_terms_transcriptome[GO_terms_transcriptome$`Adjusted P-value`<0.003,],
                                              aes(x=reorder(Term, desc(ordering)),
                                                y=`-Log(adj p-value)`, 
                                                fill= Category)) +
                                            geom_bar(
                                              stat="identity",
                                              color="black",
                                              position="dodge",
                                              alpha=.8, 
                                              width=.9) +
                                            ggtitle("TOP 20 GO terms for trancriptome gene") +
                                            theme_classic()+
                                            theme()+
                                            coord_flip()+
                                            xlab(" ")+
                                            geom_hline(yintercept = 1.3, color = "red", linetype = "dashed", size = 1)

Enrich_plot_go_terms_transcriptome

save_plot(Enrich_plot_go_terms_transcriptome,"Fig2B_transcriptome_go_terms_enrich_plot.png",9,5)
dev.off()

####KEGG Pathway
Pathways_transcriptome <- read_excel("./results/dataframes/TableS3_GO_terms_and_pathways_transcriptome_1712_genes.xlsx", 
                                sheet = "Pathways")

Pathways_transcriptome <- Pathways_transcriptome %>% 
                            mutate(ordering = as.numeric(factor(Database)) + `Adjusted P-value`)

Enrich_plot_pathways_transcriptome <-ggplot(Pathways_transcriptome[Pathways_transcriptome$`Adjusted P-value`<0.01,],
                                        aes(
                                          x=reorder(Term, desc(ordering)),
                                          y=`-Log(adj p-value)`, 
                                          fill=Database)
                                          )+
                                        geom_bar(
                                          stat="identity",
                                          color="black",
                                          position="dodge",
                                          alpha=.8, 
                                          width=.9
                                        ) +
                                        ggtitle("Top significant pathways for transcriptome genes") +
                                        theme_classic()+
                                        theme()+
                                        coord_flip()+
                                        xlab(" ")+
                                        geom_hline(yintercept = 2, color = "red", linetype = "dashed", size = 1)

Enrich_plot_pathways_transcriptome

save_plot(Enrich_plot_pathways_transcriptome,"Fig2C_transcriptome_pathways_enrich_plot.png",9,5)
dev.off()
