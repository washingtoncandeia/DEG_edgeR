##-----------------------------------
# Código Diego - 23 de março de 2018
# Modificar 
##-----------------------------------

library(edgeR)
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
library(RedeR)
library(STRINGdb)
library(bios2mds)
library(scatterplot3d)
library(gclus)
library(graphics)
library(plyr)
library(pheatmap)
library(Glimma)
library(dplyr)

featureCounts_table <- read.table("~/Dropbox/DualSeq_IMT/ALL/gene_count.txt", #arquivo
                                  sep="\t", #usar espaço como separador
                                  header = TRUE, #usar a primeira linha como head
                                  row.names = 1, #usar a primeira coluna como nome para a linha inteira
                                  check.names = FALSE)  #tira o 'X' do inicio no nome das colunas

colnames(featureCounts_table) <- gsub("\\.[sb]am$", "", colnames(featureCounts_table))

#selecting only the columns of the experiment + the controls

featureCounts_table<- featureCounts_table %>% dplyr:: select(grep("Chr", names(featureCounts_table)),
                                                             grep("Start", names(featureCounts_table)),
                                                             grep("End", names(featureCounts_table)),
                                                             grep("Strand", names(featureCounts_table)),
                                                             grep("Length", names(featureCounts_table)),
                                                             grep("10_", names(featureCounts_table)), #control
                                                             grep("11_", names(featureCounts_table)), #control
                                                             grep("19_", names(featureCounts_table)), #control
                                                             grep("20_", names(featureCounts_table)), #control
                                                             grep("33_", names(featureCounts_table)), #control
                                                             grep("34_", names(featureCounts_table)), #control
                                                             #grep("43_", names(featureCounts_table)), #control
                                                             grep("12_", names(featureCounts_table)), #ZIKV
                                                             grep("24_", names(featureCounts_table)), #ZIKV
                                                             grep("36_", names(featureCounts_table)), #ZIKV
                                                             grep("48_", names(featureCounts_table))  #ZIKV
                                                             )


#importing the condition table
peridotConditions <- read.table("~/Dropbox/DualSeq_IMT/ALL/peridotConditions.txt",
                                sep = "\t",
                                header = T,
                                check.names = F)


peridotConditions$condition <- as.factor(peridotConditions$condition)

peridotConditions <- subset(peridotConditions, peridotConditions$sample %in% colnames(featureCounts_table))

#countdata_DGE_CHIKV <- featureCounts_table[,c(6:ncol(featureCounts_table))]

### REMOVE GENES LOWLY EXPRESSED 
# 
# rpkm <- rpkm.default(featureCounts_table[6:ncol(featureCounts_table)], featureCounts_table$Length)
# keep.exprs <- rowSums(rpkm>0.05) >=3
# data_clean <- featureCounts_table[keep.exprs,, ]


CountTableLess5 = as.data.frame(apply(featureCounts_table[,c(6:ncol(featureCounts_table))], c(1, 2), function(x){
  if(x < 5){
    x = 0
  }else{
    x = x
  }
}))

#Sum all columns in each line
peridotRow = rowSums(CountTableLess5)

#Create dataframe only with the rows q have had sum of columns other than 0
CountTableWithout0 = CountTableLess5[!peridotRow == 0,]

#Create empty numeric list
baseMeanConditions = replicate(length(levels(peridotConditions$condition)), numeric())

#Calculate the baseMeans
for(i in 1:length(levels(peridotConditions$condition))){
  for(j in 1:length(rownames(CountTableWithout0))){
    baseMeanConditions[[i]] = c(baseMeanConditions[[i]], 
                                sum(CountTableWithout0[j, peridotConditions$condition == levels(peridotConditions$condition)[i]])/length(CountTableWithout0[j, peridotConditions$condition == levels(peridotConditions$condition)[i]]))
  }
}


#Sum baseMeans of each condition
baseMean = numeric()
for(i in 1:length(levels(peridotConditions$condition))){
  if(length(baseMean) == 0){
    baseMean = c(baseMean, baseMeanConditions[[i]])
  }else{
    baseMean = baseMean+baseMeanConditions[[i]]
  }
}

#Calculate the baseMean general
baseMean = baseMean/length(levels(peridotConditions$condition))

#Dataframe with baseMeans
dfBaseMeans = data.frame(baseMeanConditions, row.names = rownames(CountTableWithout0))

#Colnames of basemans with "baseMeanX" where X mean the name of condition
colnames(dfBaseMeans) = paste("baseMean", levels(peridotConditions$condition), sep = "")

#Create empty numeric list
onethirdCounts = replicate(length(levels(peridotConditions$condition)), numeric())

#Create empty dataframe
dfOneThird = data.frame(row.names = rownames(CountTableWithout0))

#Verify baseMeans with counts > 0 and which are have more of one third of samples counts
calcOneThird = function(condition){
  for(i in 1:length(rownames(CountTableWithout0))){
    increment = 0
    for(j in colnames(CountTableWithout0[,peridotConditions$condition == levels(peridotConditions$condition)[condition]])){
      if(CountTableWithout0[i, j] != 0){
        increment = increment + 1
      }
    }
    if(increment > length(CountTableWithout0[j, peridotConditions$condition == levels(peridotConditions$condition)[condition]])/3){
      onethirdCounts[[condition]] = c(onethirdCounts[[condition]], 1)
    }else{
      onethirdCounts[[condition]] = c(onethirdCounts[[condition]], 0)
    }
  }
  name = paste("onethirdCounts", condition, sep = "")
  dfOneThird[,name] <<- onethirdCounts[[condition]]
}

#Run function of calcOneThird
for(i in 1:length(levels(peridotConditions$condition))){
  calcOneThird(i)
}

#Sum the rows in dfOneThird
sumDfOneThird = rowSums(dfOneThird)

#Get names filtred for sums != 0
namesFiltered = rownames(dfOneThird[!sumDfOneThird == 0,])

#Create dfBaseMeans2 with filtred results
dfBaseMeans2 = dfBaseMeans[namesFiltered,]

#Create final dataframe with all filters
CountTableFinal = CountTableWithout0[namesFiltered,]



cpm <- cpm(CountTableFinal)
keep.exprs <- rowSums(cpm>2) >= 3
data_clean <- CountTableFinal[keep.exprs,, ]

sex <- c(rep("A",8),    #10 CTL
         rep("A",8),    #11 CTL
         rep("A",8),    #19 CTL
         rep("A",8),    #20 CTL
         rep("A",8),    #33 CTL
         rep("A",8),    #34 CTL
         #rep("A",8),    #43 CTL
         rep("A",8),    #12 ZIKV
         rep("A",8),    #24 ZIKV
         rep("A",8),    #36 ZIKV
         rep("A",8)     #48 ZIKV
         )


condition <- c(rep("CTL",8),     #10 CTL
               rep("CTL",8),     #11 CTL
               rep("CTL",8),     #19 CTL
               rep("CTL",8),     #20 CTL
               rep("CTL",8),     #33 CTL
               rep("CTL",8),     #34 CTL
               #rep("CTL",8),     #43 CTL
               rep("ZIKV",8),    #12 ZIKV
               rep("ZIKV",8),    #24 ZIKV
               rep("ZIKV",8),    #36 ZIKV
               rep("ZIKV",8)     #48 ZIKV
               )

group_list <- list('sex'=sex,'condition'=condition)

group <- factor(paste0(group_list$sex, ".", 
                       group_list$condition))


dgelist <- DGEList(counts = data_clean, 
                   genes = row.names(data_clean), 
                   group = group)


dgelist <- calcNormFactors(dgelist)

design <- model.matrix(~0+group)

colnames(design) <- gsub("group", "", colnames(design))

contr.matrix <- makeContrasts(CHIKV = A.ZIKV-A.CTL,
                              levels=design
                              )

v <- voom(dgelist, design, plot=TRUE)

vfit <- lmFit(v, design) 
vfit <- contrasts.fit(vfit, contrasts=contr.matrix) 
efit <- eBayes(vfit) 
plotSA(efit)

tfit <- treat(vfit, lfc=0)
dt <- decideTests(tfit)

summary(dt)


##### selecionando os DEGs por biblioteca e colocando-os em variaveis
# para identificar apenas os diferencialmente expressos estatisticamente precisa colocar um filtro no adj.P.val para <= 0.05

CHIKV_ALL <- topTreat(tfit, coef = "CHIKV", n=Inf)
CHIKV_DE <- subset(CHIKV_ALL, CHIKV_ALL$adj.P.Val <= 0.05)



######PLOTS

#MDS

col <- c(rep("#303030",48), #control
         rep("#D9B561",32)  #CHIKV
         )


mds <- plotMDS(dgelist, 
               method = "logFC",
               gene.selection = "pairwise",
               col="black",
               bg = col,
               pch = 21,
               top = nrow(dgelist), 
               cex = 1.5,
               axes = T)

box()
abline(h = 0, v = 0, lty = 3)
legend("bottomleft", 
       legend = c("Coltrol", "ZIKV"), 
       pch = 16, 
       col = c("#303030", "#D9B561"))

write.csv(mds$distance.matrix, "~/Dropbox/DualSeq_IMT/ZIKV/mds.csv")

dist <- read.table("~/Dropbox/DualSeq_IMT/ZIKV/mds.csv",
                   header = T,
                   row.names = 1,
                   check.names = FALSE,
                   sep = ",")



mds_pca <- mmds(as.matrix(dist), pc=5)

plot(mds_pca$coord$PC1,
     mds_pca$coord$PC2,
     #main = "MDS CHIKV"
     pch=21,
     bg = col,
     cex = 1.5,
     xlab = ("PC1 (14.6%)"),
     ylab = ("PC2 (10%)"),
     axes = T)

box()
abline(h = 0, v = 0, lty = 3)
legend("top", 
       legend = c("Coltrol", "ZIKV"), 
       pch = 16, 
       col = c("#303030", "#D9B561"))




#HEATMAP
Lib48 <- as.matrix(rowMeans(data_clean[,73:80]))
Lib36 <- as.matrix(rowMeans(data_clean[,65:72]))
Lib24 <- as.matrix(rowMeans(data_clean[,57:64]))
Lib12 <- as.matrix(rowMeans(data_clean[,49:56]))
Lib34 <- as.matrix(rowMeans(data_clean[,41:48]))
Lib33 <- as.matrix(rowMeans(data_clean[,33:40]))
Lib20 <- as.matrix(rowMeans(data_clean[,25:32]))
Lib19 <- as.matrix(rowMeans(data_clean[,17:24]))
Lib11 <- as.matrix(rowMeans(data_clean[,9:16]))
Lib10 <- as.matrix(rowMeans(data_clean[,1:8]))


mean_count_group <- cbind(Lib10, Lib11, Lib19, Lib20, Lib33, Lib34, 
                          Lib12, Lib24, Lib36,Lib48)

mean_count_group <- as.data.frame(mean_count_group)

colnames(mean_count_group) <- c("10", "11", "19", "20", "33", "34",
                                "12", "24", "36", "48")



rld <- log2(mean_count_group+0.1)
COR <- cor(as.matrix(rld))



mat_col = data.frame(Condition=factor(c("Control", #10
                                        "Control", #11
                                        "Control", #18
                                        "Control", #19
                                        "Control", #20
                                        "Control", #22
                                        "ZIKA", #23
                                        "ZIKA", #31
                                        "ZIKA", #32
                                        "ZIKA" #33
                                        )
                                      )
                     )

row.names(mat_col) = colnames(mean_count_group)

mat_colors <- list(Condition = c('Control' = "#303030", 
                                 'ZIKA' = "#D9B561" 
                                  )
                   )





pheatmap(COR,
         color = colorRampPalette(c("white","gray30"))(n = 300),   #colors from RColorBrewer,
         show_rownames = T,
         cluster_rows = T,
         cluster_cols = T,
         # breaks =  c(seq(-5.2,2.33,length=100),                            # for blue
         #             seq(2.34,9.86,length=100),                         # for black
         #             seq(9.87,17.5,length=100)
         # ),
         annotation_col = mat_col,
         annotation_colors = mat_colors,
         treeheight_row = 0,
         legend = T
)




rld2 <- log2(data_clean+0.1)
COR2 <- cor(as.matrix(rld2))


pheatmap(COR2,
         color = colorRampPalette(c("white","gray30"))(n = 300),   #colors from RColorBrewer,
         show_rownames = F,
         cluster_rows = T,
         cluster_cols = T
         )






########### enriquecimentio

GO_ZIKA <- enrichGO(CHIKV_DE$genes,
                    keyType = "SYMBOL",
                    OrgDb = "org.Hs.eg.db",
                    ont = "BP",
                    pAdjustMethod = 'fdr'
                    )
GO_ZIKA_symp <- simplify(GO_ZIKA, cutoff = 0.7)


SYMBOL <- toTable(org.Hs.egSYMBOL)

kegg_zika <- enrichKEGG(gene = subset(SYMBOL$gene_id, SYMBOL$symbol %in% CHIKV_DE$genes),
                        keyType = "ncbi-geneid",
                        organism = "hsa",
                        qvalueCutoff = 0.05,
                        pAdjustMethod = "fdr"
                        )





splited_col <- strsplit(GO_ZIKA_symp@result$geneID, '/')

GO_teste <- enrichGO(unlist(splited_col[183]),
                    keyType = "SYMBOL",
                    OrgDb = "org.Hs.eg.db",
                    ont = "BP",
                    pAdjustMethod = 'fdr'
)


##### Volcan plot

with(subset(CHIKV_ALL, logFC >= -6), plot(logFC,-log10(adj.P.Val), pch=16, axes=T, xlab = "Log2(FC)", ylab = "-log10(P-value adjusted)", main = " "))
with(subset(subset(CHIKV_ALL, adj.P.Val <= 0.01), logFC < 0), points(logFC,-log10(adj.P.Val), pch=16, col="steelblue"))
with(subset(subset(CHIKV_ALL, adj.P.Val <= 0.01), logFC > 0), points(logFC,-log10(adj.P.Val),pch=16, col="tomato3"))
axis(2, las=1, at=c(0,30,60,90,120))
axis(1, las=1, at=c(-5,0,5))
box()

#####reder

net <- cnetplot(GO_ZIKA_symp, categorySize="geneNum", showCategory = 95)
net <- cnetplot(GO_ZIKA_symp[95], categorySize="geneNum", showCategory = 1)

rdp <- RedPort()
calld(rdp)
addGraph(rdp, net)
resetd(rdp)



