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

featureCounts_table <- read.table("gene_count.txt", #arquivo
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
                                                             grep("33_", names(featureCounts_table)), #control
                                                             grep("12_", names(featureCounts_table)), #ZIKV
                                                             grep("24_", names(featureCounts_table)), #ZIKV
                                                             grep("36_", names(featureCounts_table)), #ZIKV
                                                             grep("48_", names(featureCounts_table))  #ZIKV
                                                             )


#importing the condition table
peridotConditions <- read.table("peridotConditions.txt",
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

sex <- c(rep("F",8),    #10 CTL
         rep("F",8),    #33 CTL
         rep("F_12",8), #12 GBS
         rep("F_24",8), #24 ZIKV
         rep("F_36",8), #36 ZIKV
         rep("F_48",8)  #48 ZIKV
         )


condition <- c(rep("CTL",8),     #10 CTL
               rep("CTL",8),     #33 CTL
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

contr.matrix <- makeContrasts(F_12 = F_12.ZIKV-F.CTL,
                              F_24 = F_24.ZIKV-F.CTL,
                              F_36 = F_36.ZIKV-F.CTL,
                              F_48 = F_48.ZIKV-F.CTL,
                              levels=design
                              )

v <- voom(dgelist, design, plot=TRUE)

vfit <- lmFit(v, design) 
vfit <- contrasts.fit(vfit, contrasts=contr.matrix) 
efit <- eBayes(vfit) 
plotSA(efit)

tfit <- treat(vfit, lfc=1.5)
dt <- decideTests(tfit)

summary(dt)


##### selectin DEGs por biblioteca

DE_12 <- topTreat(tfit, coef = "F_12", n=Inf)
DE_24 <- topTreat(tfit, coef = "F_24", n=Inf)
DE_36 <- topTreat(tfit, coef = "F_36", n=Inf)
DE_48 <- topTreat(tfit, coef = "F_48", n=Inf)  


