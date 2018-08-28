#################################################################################
####################### e d g e R - ZIKV Wash 31/07/2017 ########################
#################################################################################

source("https://bioconductor.org/biocLite.R")
library(edgeR)
library(DESeq2)
library(ShortRead)
library(Rsubread)
library(readr)
library(RColorBrewer)
library(gplots)
library(org.Hs.eg.db)
library(statmod)
library(GO.db)
library(clusterProfiler)
library(pheatmap)

library(rgl)
library(amap)
library(e1071)
library(scales)
library(cluster)

library(bios2mds)

# install.packages("bios2mds", dependencies = TRUE, repos = c(CRAN="https://cran.r-project.org/"))


#calling the file without setting the first collunm as row's name
countdata_DGE_ZIKV <- read.table("gene_count_ZIKV.txt",          # arquivo
                                  sep=",",                       # usar espaço como separador
                                  header = TRUE,                 # usar a primeira linha como head
                                  #row.names = 1,                # usar a primeira coluna como nome para a linha inteira
                                  check.names = FALSE)           # tira a PORRA do 'X' do inicio no nome das colunas


#removing the '.bam' at the end of each colunm name
colnames(countdata_DGE_ZIKV) <- gsub("\\.[sb]am$", "", colnames(countdata_DGE_ZIKV))

# Ajustando o sexo dos pacientes
#setting the sex of the patients 
sex <- c(rep("U", 56))

# Usando repetições para gerar fatores representando condições
condition <- c(rep("CTL", 8), rep("ZIKA", 8), rep("ZIKA", 8), rep("CTL", 8), rep("ZIKA", 8), rep("CTL", 8), rep("ZIKA", 8))


#trnsform the condition and the replicates in a matrix
group_list <- list('sex'=sex,'condition'=condition)

#transform the condition and the replicates into a factor
group <- factor(paste0(group_list$sex, ".", 
                       group_list$condition))

(group)

#ceating the edgeR dataframe
dgelist <- DGEList(counts = countdata_DGE_ZIKV[, 2:ncol(countdata_DGE_ZIKV)], 
                   genes = countdata_DGE_ZIKV$Geneid, 
                   group = group)

#creating a matrix that stores the gene_id based on the gene symbol
ann <- toTable(org.Hs.egSYMBOL)

ann2 <- toTable(org.Hs.egUNIPROT)

ann3 <- toTable(org.Hs.egPATH)
ann

#createing a matrix that takes the symbols used by the dataframe and introducing the gene_id to each one
m <- match(dgelist$genes$genes, 
           ann$symbol)

dgelist$genes$entrez <- ann$gene_id[m]

#createing a matrix that takes the gene_id used by the dataframe and introducing the uniprot_id to each one
m2 <- match(dgelist$genes$entrez,
            ann2$gene_id)

dgelist$genes$uniprot <- ann2$uniprot_id[m2]

#createing a matrix that takes the gene_id used by the dataframe and introducing the KEGG_id to each one
m3 <- match(dgelist$genes$entrez,
            ann3$gene_id)

dgelist$genes$kegg <- ann3$path_id[m3]


#this step slects genes with cpm that are greater tham 0.5 in at least 6 libraryes (inside repetitions)
keep <- rowSums(cpm(dgelist) >= 0.5) >= 6

#setting the design of the experiment according to the technical repetitions set above
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

# Modificar aqui, e colocar o número igual aos fatores sex (88) e condition (ver CTL e ZIKA)
col <- c(rep(2, 8), rep(4, 8), rep(4, 8), rep(2, 8), rep(4, 8), rep(2, 8), rep(4, 8))

# 10() 11(h) (12, Zika) 19(h) 20(h) (24, Zika) 33(m) 34(h) (36, Zika) 43(m) (48, Zika)
# 2 ct
# 4 zika

forma <- c(rep(23, 8), rep(22, 8), rep(22, 8), rep(23, 8), rep(22, 8), rep(23, 8), rep(22, 8))

# 22 sintomatico
# 23 controle

# Abaixo, reconfigurar a margem no próprio Plots, abaixo
list_MDS <- plotMDS(dgelist, col=col, pch=forma, top = 200)  # Error in plot.new() : figure margins too large

#write.csv(list_MDS$distance.matrix, file="list_MDS_modified.csv")


################# MDS WITH PERCENTAGE VARIATIONS OF  EACH PC #################
# Modificaçao feita em08/08/2017 

dist <- read.table("list_MDS_modified.csv", #arquivo
                   sep=",", #usar espaço como separador
                   header = TRUE, #usar a primeira linha como head
                   row.names = 1, #usar a primeira coluna como nome para a linha inteira
                   check.names = FALSE)  #tira a PORRA do 'X' do inicio no nome das colunas

#mds_pca <- mmds(as.matrix(dist), pc = 5)
#mds_pca


#plot(mds_pca$coord$PC1, mds_pca$coord$PC2, pch=forma, bg = col, cex = 1.3)

############################################################################

#creating the contrast between CHIKV and CTL
con_ZIKA <- makeContrasts(U.ZIKA - U.CTL, levels=design)


#### TUTORIAL ZHOU AND ROBINSON, 2014 ####


dgelist <- calcNormFactors(dgelist)

#dgelist <- estimateGLMRobustDisp(dgelist, design)


dgelist <- estimateGLMCommonDisp(dgelist, design = design)
dgelist <- estimateGLMTrendedDisp(dgelist,design=design)
dgelist <- estimateGLMTagwiseDisp(dgelist, design = design, prior.df = 10)  

fit <- glmFit(dgelist, design = design)
lr <- glmLRT(fit, coef=2, contrast = con_ZIKA)  # Substituir con_CHIKV

##########################################

## Comentei para não precisar dar o comando:
# summary(decideTests(lr))  # ERRO dims [product 169616] do not match the length of object [17]

## Comentei para não precisar dar o comando:
# plotMD(lr)

#########-------------------------------------------------------------------------

#salvando algumas informações umportantes do dataframe em uma tabela separada. 
DE <- topTags(lr, n=nrow(dgelist))$table
nrow(DE)
DE$AveLogCPM <- lr$AveLogCPM
DE$padj <- p.adjust(DE$PValue, method = "fdr", n = length(DE$PValue))

# 
write.table(as.data.frame(DE), sep = ",", "zika_geneCountDE_7samples.txt", row.names = F)


#aqui eu filtro o dado pelo log do foldchange acima de 1.5e abaixo de -1.5
DE_subset <- subset(DE, logFC >= 1.0 | logFC <= -1.0)
nrow(DE_subset)

#aqui eu estou filtrando pelo 
DE_subset <- subset(DE_subset, padj <= 0.01)
nrow(DE_subset)

#separando o dataset em genes up e down regulados
DE_down <- subset(DE_subset, logFC <= -1.0)         # Downregulated genes
DE_up <- subset(DE_subset, logFC >= 1.0)            # Upregulated genes


## Transformar cada um numa tabela para 
## aproveitar melhor a lista de dados.

# Tabela de contagem de genes 'upregulados'
write.table(as.data.frame(DE_up), sep = ",", "zika_UPREGULATED_7samples.txt", row.names = F)

# Tabela de contagem de genes 'downregulados'
write.table(as.data.frame(DE_down), sep = ",", "zika_DOWNREGULATED_7samples.txt", row.names = F)

  # MA PLOT
with(DE, plot(AveLogCPM, logFC, pch=21, xlim=c(-2,18), axes=FALSE, main="ZIKV vs CTL"))
axis(1, at=c(-2,0,2,4,6,8,10,12,14,16,18))
axis(2, at=c(-10,-8,-6,-4,-2,0,2,4,6,8))
box()
with(subset(DE_subset, logFC <= -1.0), points(AveLogCPM, logFC, pch=21, col="black", bg="steelblue4"))
with(subset(DE_subset, logFC >= 1.0), points(AveLogCPM, logFC,pch=21, col="black", bg="tomato3"))
legend("topright", legend = c("Non-DE", "Up", "Down"), pch=c(21,21,21), col="black", pt.bg = c("white","tomato3","steelblue4"))

head(DE_up$genes)

write.csv(DE_up$genes, "genes_UPrregulados.txt")
write.csv(DE_down$genes, "genes_DOWNrregulados.txt")


