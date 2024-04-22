libraries <- c('readxl', 'RColorBrewer', 'tidyverse')
lapply(libraries,library, character.only = TRUE)

#Helper function
save_plot <- function(plot,name, w, h){
  png(file=paste0("./results/figures/",name), width=w, height=h, units = "in", res = 300)
  plot
}

####ENRICHMENT ANALYSIS TRANSCRIPTOME####
#####GO terms results#####
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

#####Pathways results#####

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

####ENRICHMENT ANALYSIS PROTEOME####

#####GO terms results#####
GO_terms_proteome <- read_excel("./results/dataframes/TableS4_GO_terms_and_pathways_proteome_158_genes.xlsx", 
                                sheet = "GO_terms")

GO_terms_proteome <- GO_terms_proteome %>% 
                        mutate(ordering = as.numeric(factor(Category, levels = c("CC","MF","BP")))+ `Adjusted P-value`)  

Enrich_plot_go_terms_proteome <- ggplot(GO_terms_proteome[GO_terms_proteome$`Adjusted P-value`<0.001611,],
                                            aes(
                                               x=reorder(Term, desc(ordering)),
                                               y=`-Log(adj p-value)`, 
                                               fill= Category, 
                                               #group=Term,
                                               #label=Term
                                            )
                                            ) +
                                              geom_bar(
                                                stat="identity",
                                                color="black",
                                                position="dodge",
                                                alpha=.8, 
                                                width=.9
                                              ) +
                                              ggtitle("TOP 20 GO terms for Proteome gene") +
                                              theme_classic()+
                                              theme()+
                                              coord_flip()+
                                              xlab(" ")+
                                              geom_hline(yintercept = 2, color = "red", linetype = "dashed", size = 1)

Enrich_plot_go_terms_proteome

save_plot(Enrich_plot_go_terms_proteome,"Fig4B_proteome_go_terms_enrich_plot.png",9,5)
dev.off()

#####Pathways Results#####
Pathways_proteome <- read_excel("./results/dataframes/TableS4_GO_terms_and_pathways_proteome_158_genes.xlsx", 
                                sheet = "Pathways")

Pathways_proteome <- Pathways_proteome %>% 
                                  mutate(ordering = as.numeric(factor(Database)) + `Adjusted P-value`)

Enrich_plot_pathways_proteome <- ggplot(Pathways_proteome[Pathways_proteome$`Adjusted P-value`<0.01,],
                                          aes(
                                            x=reorder(Term, desc(ordering)),
                                            y=`-Log(adj p-value)`, 
                                            fill=Database
                                          )
                                          ) +
                                          geom_bar(
                                            stat="identity",
                                            color="black",
                                            position="dodge",
                                            alpha=.8, 
                                            width=.9
                                          ) +
                                          ggtitle("Top significant pathways for proteome genes") +
                                          theme_classic()+
                                          theme()+
                                          coord_flip()+
                                          xlab(" ")+
                                          geom_hline(yintercept = 2, color = "red", linetype = "dashed", size = 1)

Enrich_plot_pathways_proteome

save_plot(Enrich_plot_pathways_proteome,"Fig4C_proteome_pathways_enrich_plot.png",9,5)
dev.off()

####ENRICHMENT ANALYSIS OVERLAP GENES/PROTEINS####
#####Go terms results#####
GO_terms_overlapping <- read_excel("./results/dataframes/TableS6_GO_terms_and_pathways_overlapping_75_genes.xlsx", 
                                   sheet = "GO_terms")

GO_terms_overlapping <- GO_terms_overlapping %>% 
  mutate(ordering = as.numeric(factor(Category, levels = c("CC","MF","BP")))+ `Adjusted P-value`) 

Enrich_plot_go_terms_overlap_genes <- ggplot(GO_terms_overlapping,
                                                  aes(
                                                    x=reorder(Term, dplyr::desc(ordering)),
                                                    y=`-Log(adj p-value)`, 
                                                    fill= Category
                                                  )
                                                  ) +
                                                  geom_bar(
                                                    stat="identity",
                                                    color="black",
                                                    position="dodge",
                                                    alpha=.8, 
                                                    width=.9
                                                  ) +
                                                  ggtitle("GO terms for overlapping genes") +
                                                  theme_classic()+
                                                  theme(axis.text = element_text(size = 11))+
                                                  coord_flip()+
                                                  xlab(" ")+
                                                  geom_hline(yintercept = 1.3, color = "red", linetype = "dashed", size = 1)

save_plot(Enrich_plot_go_terms_overlap_genes,"Fig6D_overlap_genes_go_terms_enrich_plot.png",11,5.25)
dev.off()

#####Pathways results#####
Pathways_overlapping <- read_excel("./results/dataframes/TableS6_GO_terms_and_pathways_overlapping_75_genes.xlsx", 
                                   sheet = "Pathways")

Pathways_overlapping <- Pathways_overlapping %>% 
                            mutate(ordering = as.numeric(factor(Database)) + `Adjusted P-value`)

Enrich_plot_pathways_overlap_genes <- ggplot(Pathways_overlapping[Pathways_overlapping$`Adjusted P-value`<0.01,],
                                            aes(
                                              x=reorder(Term, dplyr::desc(ordering)),
                                              y=`-Log(adj p-value)`, 
                                              fill=Database, 
                                              #group=Term,
                                              #label=Term
                                            )
                                            ) +
                                            geom_bar(
                                              stat="identity",
                                              color="black",
                                              position="dodge",
                                              alpha=.8, 
                                              width=.9
                                            ) +
                                            ggtitle("Pathways for overlapping genes") +
                                            theme_classic()+
                                            theme(text = element_text(color = "black"))+
                                            coord_flip()+
                                            xlab(" ")+
                                            geom_hline(yintercept = 2, color = "red", linetype = "dashed", size = 1)

Enrich_plot_pathways_overlap_genes

save_plot(Enrich_plot_pathways_overlap_genes,"Fig6E_overlap_genes_pathways_enrich_plot.png",10,7.3)
dev.off()

