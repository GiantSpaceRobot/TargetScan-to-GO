## ---------------------------
##
## Script name: TargetScan-to-GO-analysis.R
##
## Purpose of script: Given raw TargetScan output, process data and analyse using EnrichR
##
## Author: Dr. Paul Donovan
##
## Date Created: 04-10-19
##
## Email: pauldonovan@rcsi.com
##
## ---------------------------
##
## Notes: 
##   
##
## ---------------------------

#install.packages("enrichR")
#install.packages("ggplot2")
library("enrichR")
library("ggplot2")
library("DOSE")
library("enrichplot")
library("clusterProfiler")
#library("gridExtra")

### Use command line arguments
args = commandArgs(trailingOnly=TRUE)
my.input <- args[1]
my.output <- args[2]

my.results <- read.csv(my.input, 
                       sep = "\t", 
                       header = F)

my.genes <- as.character(my.results$V1)
dbs <- c("Enrichr_Libraries_Most_Popular_Genes",
         "GO_Molecular_Function_2018", 
         "GO_Cellular_Component_2018", 
         "GO_Biological_Process_2018",
         "Chromosome_Location_hg19",
         "OMIM_Disease",
         "OMIM_Expanded",
         "ENCODE_Histone_Modifications_2015",
         "CORUM",
         "InterPro_Domains_2019",
         "KEA_2015", # Kinase Enrichment Analysis
         "miRTarBase_2017",
         "NCI-Nature_2016", # Pathway Interaction Database
         "NURSA_Human_Endogenous_Complexome", # Protein - protein interactions and their networks
         "Panther_2016", # protein analysis through evolutionary relationships
         "Pfam_Domains_2019",
         "Pfam_InterPro_Domains",
         "Phosphatase_Substrates_from_DEPOD",
         "PPI_Hub_Proteins",
         "Rare_Diseases_GeneRIF_Gene_Lists",
         "Reactome_2016",
         "TargetScan_microRNA_2017",
         "WikiPathways_2019_Human",
         "KEGG_2019_Human")
my.enrichR <- enrichr(my.genes, databases = dbs)

my.summary <- printEnrich(data = my.enrichR, 
            file = paste0(my.output, "_Summary.txt"), 
            sep = "\t", 
            columns = c(1:9))

write.table(x = my.enrichR[["Enrichr_Libraries_Most_Popular_Genes"]], file = paste0(my.output, "_Enrichr_Libraries_Most_Popular_Genes.txt"))
write.table(x = my.enrichR[["GO_Molecular_Function_2018"]], file = paste0(my.output, "_GO_Molecular_Function_2018.txt"))
write.table(x = my.enrichR[["GO_Cellular_Component_2018"]], file = paste0(my.output, "_GO_Cellular_Component_2018.txt"))
write.table(x = my.enrichR[["GO_Biological_Process_2018"]], file = paste0(my.output, "_GO_Biological_Process_2018.txt"))
write.table(x = my.enrichR[["Chromosome_Location_hg19"]], file = paste0(my.output, "_Chromosome_Location_hg19.txt"))
write.table(x = my.enrichR[["OMIM_Disease"]], file = paste0(my.output, "_OMIM_Disease.txt"))
write.table(x = my.enrichR[["OMIM_Expanded"]], file = paste0(my.output, "_OMIM_Expanded.txt"))
write.table(x = my.enrichR[["ENCODE_Histone_Modifications_2015"]], file = paste0(my.output, "_ENCODE_Histone_Modifications_2015.txt"))
write.table(x = my.enrichR[["CORUM"]], file = paste0(my.output, "_CORUM.txt"))
write.table(x = my.enrichR[["InterPro_Domains_2019"]], file = paste0(my.output, "_InterPro_Domains_2019.txt"))
write.table(x = my.enrichR[["KEA_2015"]], file = paste0(my.output, "_KEA_2015.txt"))
write.table(x = my.enrichR[["miRTarBase_2017"]], file = paste0(my.output, "_miRTarBase_2017.txt"))
write.table(x = my.enrichR[["NCI-Nature_2016"]], file = paste0(my.output, "_NCI-Nature_2016.txt"))
write.table(x = my.enrichR[["NURSA_Human_Endogenous_Complexome"]], file = paste0(my.output, "_NURSA_Human_Endogenous_Complexome.txt"))
write.table(x = my.enrichR[["Panther_2016"]], file = paste0(my.output, "_Panther_2016.txt"))
write.table(x = my.enrichR[["Pfam_Domains_2019"]], file = paste0(my.output, "_Pfam_Domains_2019.txt"))
write.table(x = my.enrichR[["Pfam_InterPro_Domains"]], file = paste0(my.output, "_Pfam_InterPro_Domains.txt"))
write.table(x = my.enrichR[["Phosphatase_Substrates_from_DEPOD"]], file = paste0(my.output, "_Phosphatase_Substrates_from_DEPOD.txt"))
write.table(x = my.enrichR[["PPI_Hub_Proteins"]], file = paste0(my.output, "_PPI_Hub_Proteins.txt"))
write.table(x = my.enrichR[["Rare_Diseases_GeneRIF_Gene_Lists"]], file = paste0(my.output, "_Rare_Diseases_GeneRIF_Gene_Lists.txt"))
write.table(x = my.enrichR[["Reactome_2016"]], file = paste0(my.output, "_Reactome_2016.txt"))
write.table(x = my.enrichR[["TargetScan_microRNA_2017"]], file = paste0(my.output, "_TargetScan_microRNA_2017.txt"))
write.table(x = my.enrichR[["WikiPathways_2019_Human"]], file = paste0(my.output, "_WikiPathways_2019_Human.txt"))
write.table(x = my.enrichR[["KEGG_2019_Human"]], file = paste0(my.output, "_KEGG_2019_Human.txt"))

df.mpg <- my.enrichR[["Enrichr_Libraries_Most_Popular_Genes"]]
df.mf <- my.enrichR[["GO_Molecular_Function_2018"]]
df.cc <- my.enrichR[["GO_Cellular_Component_2018"]]
df.bp <- my.enrichR[["GO_Biological_Process_2018"]]
df.chr <- my.enrichR[["Chromosome_Location_hg19"]]
df.disease <- my.enrichR[["OMIM_Disease"]]
df7.disease.exp <- my.enrichR[["OMIM_Expanded"]]
df8.mods <- my.enrichR[["ENCODE_Histone_Modifications_2015"]]
df9.CORUM <- my.enrichR[["CORUM"]] 
df10.InterPro <- my.enrichR[["InterPro_Domains_2019"]]
df11.kea <- my.enrichR[["KEA_2015"]] 
df12.miRTarBase <- my.enrichR[["miRTarBase_2017"]]
df13.NCI <- my.enrichR[["NCI-Nature_2016"]]
df14.NURSA <- my.enrichR[["NURSA_Human_Endogenous_Complexome"]]
df15.Panther <- my.enrichR[["Panther_2016"]]
df16.Pfam <- my.enrichR[["Pfam_Domains_2019"]]
df17.Pfam.InterPro <- my.enrichR[["Pfam_InterPro_Domains"]]
df18.Phos <- my.enrichR[["Phosphatase_Substrates_from_DEPOD"]]
df19.PPI <- my.enrichR[["PPI_Hub_Proteins"]]
df20.RareDisease <- my.enrichR[["Rare_Diseases_GeneRIF_Gene_Lists"]]
df21.reactome <- my.enrichR[["Reactome_2016"]]
df22.TargetScan <- my.enrichR[["TargetScan_microRNA_2017"]]
df23.WikiPathways <- my.enrichR[["WikiPathways_2019_Human"]]
df24.KEGG <- my.enrichR[["KEGG_2019_Human"]]

df.to.plot <- function(my.dataframe, analysis = "") {
  sub.df <- head(my.dataframe, n = 10)
  sub.df <- sub.df[order(-sub.df$Adjusted.P.value),]
  sub.df$log10.of.padj <- abs(log10(sub.df$Adjusted.P.value))
  sub.df <- subset(sub.df, log10.of.padj>1.3)
  number.of.terms <- nrow(sub.df)
  text.size <- ((2/number.of.terms)+1.18)^3 #### Create bespoke text size for each plot
  level_order <- as.list(sub.df$Term) # Create factor of terms to reorder ggplot2
  my.plot <- 
    ggplot(data=sub.df, aes(x=factor(Term, level = level_order), y=log10.of.padj)) +
      geom_bar(stat="identity", fill="darkred") +  #steelblue
      geom_text(
        aes(label=Term), 
        hjust = 1,
        color="black", 
        size=text.size) +
      theme_minimal() +
      theme(axis.text.y=element_blank(), 
            axis.ticks.y=element_blank(), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) +
      scale_x_discrete(name = "") +
      scale_y_continuous(name = "Log10 of Adjusted P value") +
      labs(title=analysis) +
      coord_flip()
  return(my.plot)
}

pdf(file = paste0(my.output, "_EnrichR_Barplots.pdf"))
df.to.plot(df.mpg, "Enrichr Libraries Most Popular Genes")
df.to.plot(df.mf, "GO Molecular Function 2018")
df.to.plot(df.cc, "GO Cellular Component 2018")
df.to.plot(df.bp, "GO Biological Process 2018")
df.to.plot(df.chr, "Chromosome Location hg19")
df.to.plot(df.disease, "OMIM Disease")
df.to.plot(df7.disease.exp, "OMIM Expanded")
df.to.plot(df8.mods, "ENCODE Histone Modifications 2015")
df.to.plot(df9.CORUM, "CORUM") 
df.to.plot(df10.InterPro, "InterPro Domains 2019")
df.to.plot(df11.kea, "KEA 2015") 
df.to.plot(df12.miRTarBase, "miRTarBase 2017")
df.to.plot(df13.NCI, "NCI-Nature 2016")
df.to.plot(df14.NURSA, "NURSA Human Endogenous Complexome")
df.to.plot(df15.Panther, "Panther 2016")
df.to.plot(df16.Pfam, "Pfam Domains 2019")
df.to.plot(df17.Pfam.InterPro, "Pfam InterPro Domains")
df.to.plot(df18.Phos, "Phosphatase Substrates from DEPOD")
df.to.plot(df19.PPI, "PPI Hub Proteins")
df.to.plot(df20.RareDisease, "Rare Diseases GeneRIF Gene Lists")
df.to.plot(df21.reactome, "Reactome 2016")
df.to.plot(df22.TargetScan, "TargetScan microRNA 2017")
df.to.plot(df23.WikiPathways, "WikiPathways 2019 Human")
df.to.plot(df24.KEGG, "KEGG 2019 Human")
dev.off()

### Network plot generation
target.genes <- as.character(my.input$Target.gene)
symbols <- mapIds(org.Hs.eg.db, keys = target.genes, keytype = "SYMBOL", column="ENTREZID")
symbols.na.omit <- na.omit(symbols)
edo <- enrichKEGG(symbols.na.omit, pvalueCutoff = 0.05, qvalueCutoff = 0.05)
#barplot(edo)
edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
pdf(file = paste0(my.output, "_Network.pdf"))
cnetplot(edox, 
         circular = F, 
         foldChange=NULL, 
         showCategory = 5, 
         colorEdge = T)
dev.off()


