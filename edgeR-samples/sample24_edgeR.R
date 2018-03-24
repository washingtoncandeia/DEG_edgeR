##---------------------------------------------------------
# Análise de contagem de genes RNA-seq, edgeR
# Versão Modificada para análises de replicatas separadas
# Data: 16 AGO 2017
# Washington C. Araujo
##-------------------------------------------------------
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
library(rgl)
library(amap)
library(e1071)
library(scales)
library(cluster)
library(bios2mds)


# install.packages("bios2mds", dependencies = TRUE, repos = c(CRAN="https://cran.r-project.org/"))


#calling the file without setting the first collunm as row's name
countdata_DGE_ZIKV <- read.table("gene_count_ZIKV.txt", #arquivo
                                 sep=",", #usar espaço como separador
                                 header = TRUE, #usar a primeira linha como head
                                 row.names = 1, #usar a primeira coluna como nome para a linha inteira
                                 check.names = FALSE)  #tira a PORRA do 'X' do inicio no nome das colunas


# Limpeza fornecida por Diego
cpm_log <- cpm(countdata_DGE_ZIKV, log = TRUE)
median_log2_cpm <- apply(cpm_log, 1, median)
hist(median_log2_cpm)
expr_cutoff <- -1
abline(v = expr_cutoff, col = "red", lwd = 3)
sum(median_log2_cpm > expr_cutoff)
data_clean_zika <- countdata_DGE_ZIKV[median_log2_cpm > expr_cutoff, ]

# Modificação feita: colocar data_clean_zika)
#removing the '.bam' at the end of each colunm name
colnames(data_clean_zika) <- gsub("\\.[sb]am$", "", colnames(data_clean_zika))

# Ajustando o sexo dos pacientes
#setting the sex of the patients 
sex <- c(rep("CTL", 8),
         rep("F_12",8),
         rep("F_24",8),
         rep("CTL",8),
         rep("F_36",8),
         rep("F_48",8))

# Usando repetições para gerar fatores representando condições
condition <- c(rep("CTL", 8), 
               rep("ZIKA", 8),
               rep("ZIKA", 8),
               rep("CTL", 8), 
               rep("ZIKA", 8), 
               rep("ZIKA", 8))



#trnsform the condition and the replicates in a matrix
group_list <- list('sex'=sex,'condition'=condition)

#transform the condition and the replicates into a factor
group <- factor(paste0(group_list$sex, ".", 
                       group_list$condition))

head(group)

#ceating the edgeR dataframe
dgelist <- DGEList(counts = data_clean_zika, 
                   genes = row.names(data_clean_zika), 
                   group = group)

# Verificando a coluna genes
head(dgelist$genes)

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
################################# Design do experimento #################################
#setting the design of the experiment according to the technical repetitions set above
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

# Modificar aqui, e colocar o número igual aos fatores sex (88) e condition (ver CTL e ZIKA)
col <- c(rep(2, 8), rep(4, 8), rep(4, 8), rep(2, 8), rep(4, 8), rep(4, 8))

# 10(CTL) (12, Zika) (24, Zika) 33(CTL) 3(36, Zika) 43(CTL) (48, Zika)
# 2 vermelho - controles
# 4 azul - ZIKA virus

forma <- c(rep(23, 8), rep(22, 8), rep(22, 8), rep(23, 8), rep(22, 8), rep(22, 8))

# 22 sintomatico
# 23 controle

# Abaixo, reconfigurar a margem no próprio Plots, abaixo
list_MDS <- plotMDS(dgelist, col=col, pch=forma, top = nrow(dgelist))  # Error in plot.new() : figure margins too large

#write.csv(list_MDS$distance.matrix, file="list_MDS_modified.csv")

############################################################################
#creating the contrast between CHIKV and CTL

#### TUTORIAL ZHOU AND ROBINSON, 2014 ####
dgelist <- estimateGLMRobustDisp(dgelist, design)

#dgelist <- estimateGLMCommonDisp(dgelist, design = design)
#dgelist <- estimateGLMTrendedDisp(dgelist,design=design)
#dgelist <- estimateGLMTagwiseDisp(dgelist, design = design, prior.df = 10)  

#setting the design of the experiment according to the technical repetitions set above

##------------------- Análise de Zika Virus - Replicata 24 -------------------##

con_zika_24 <- makeContrasts(F_24.ZIKA - CTL.CTL, levels=design)
fit <- glmFit(dgelist, design = design)
lr_24 <- glmLRT(fit, coef=2, contrast = con_zika_24)
DE_24 <- topTags(lr_24, n=nrow(dgelist))$table
nrow(DE_24)
DE_24$AveLogCPM <- lr_24$AveLogCPM
DE_24$padj <- p.adjust(DE_24$PValue, method = "fdr", n = length(DE_24$PValue))

# Padj filter 
DE_24_padj <- subset(DE_24, padj <= 0.05)
nrow(DE_24_padj)

#logFC filter
DE_24_FC <- subset(DE_24_padj, logFC >= 1.0 | logFC <= -1.0)   # Observar se deixa 1.5
nrow(DE_24_FC)


# Análise de Expressão Diferencial

############################################################################

#salvando algumas informações umportantes do dataframe em uma tabela separada. 
DE_24 <- topTags(lr_24, n=nrow(dgelist))$table
nrow(DE_24)
DE_24$AveLogCPM <- lr_24$AveLogCPM
DE_24$padj <- p.adjust(DE_24$PValue, method = "fdr", n = length(DE_24$PValue))

# Escrever para ZIKV replicata 12
write.table(as.data.frame(DE_24), sep = ",", "DE_24_ZIKV.txt", row.names = F)

#aqui eu filtro o dado pelo log do foldchange acima de 2.0 e abaixo de -2.0
DE_subset_24 <- subset(DE_24, logFC >= 1.0 | logFC <= -1.0)  
nrow(DE_subset_24)

#aqui eu estou filtrando pelo 
DE_subset_24 <- subset(DE_subset_24, padj <= 0.05)
nrow(DE_subset_24)

#separando o dataset em genes up e down regulados
DE_down_24 <- subset(DE_subset_24, logFC <= -1.0)         # Downregulated genes
DE_up_24 <- subset(DE_subset_24, logFC >= 1.0)            # Upregulated genes

## Transformar cada um numa tabela para 
## aproveitar melhor a lista de dados.

# Tabela de contagem de genes 'upregulados'
write.table(as.data.frame(DE_up_24), sep = ",", "DE_up_24_count_ZIKV_19-11-2017.txt", row.names = F)

# Tabela de contagem de genes 'downregulados'
write.table(as.data.frame(DE_down_24), sep = ",", "DE_down_24_count_ZIKV_19-11-09-2017.txt", row.names = F)

# MA PLOT
summary(DE_24$AveLogCPM)

pdf("24_MA_DE-19-11-2017.pdf")
with(DE_24, plot(AveLogCPM, logFC, pch=21, xlim=c(-0.5,18), axes=FALSE, main="ZIKV Sample 24 vs CTL"))
axis(1, at=c(-0.5,0.5,1.55,2.55,3.5,4.55,5.5,6.5,7.55,8.5,9.55,10.5,12,14,16,18))
axis(2, at=c(-10,-8,-6,-4,-2,0,2,4,6,8,10,12))
box()
with(subset(DE_subset_24, logFC <= -1.0), points(AveLogCPM, logFC, pch=21, col="black", bg="steelblue4"))
with(subset(DE_subset_24, logFC >= 1.0), points(AveLogCPM, logFC,pch=21, col="black", bg="tomato3"))
legend("topright", legend = c("Non-DE", "Up", "Down"), pch=c(21,21,21), col="black", pt.bg = c("white","tomato3","steelblue4"))
dev.off()

# Tabela Up
#write.table(DE_up_24$genes, "24_genes_DE_UP_26-09-2017.txt", row.names = T, col.names = T)

# Tabela Down
#write.table(DE_down_24$genes, "24_genes_DE_DOWN_26-09-2017.txt", row.names = T, col.names = T)

write(DE_up_24$genes, "24_Lista_DE_UPREGULATED_novembro.txt")

write(DE_down_24$genes, "24_Lista_DE_DOWNREGULATED_novembro.txt")

