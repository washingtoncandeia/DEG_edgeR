# Airway smooth muscle cells
# Tutoriais:
# http://bioconductor.org/packages/release/data/experiment/vignettes/airway/inst/doc/airway.html
# https://github.com/Bioconductor/CSAMA2015/blob/master/materials/labs/2_Tuesday/rnaseqCSAMA.md
# https://www.bioconductor.org/help/course-materials/2015/BioC2015/bioc2015rnaseq.html#count
# https://f1000research.com/articles/4-1070/v2
# Tutorial RNA0-seq 
library(airway)
library(GEOquery)

# Obendo informações de amostras com GEOquery:
dir <- system.file('extdata', package = 'airway')
geofile <- file.path(dir, 'GSE52778_series_matrix.txt')
gse <- getGEO(filename = geofile)

# Agora, organizando os dados.
pdata <- pData(gse)[ , grepl('ch1', names(pData(gse)))]
names(pdata) <- c('treatment', 'tissue', 'ercc_mix', 'cell', 'celltype')
pdataclean <- data.frame(treatment = sub("treatment: (.*)", "\\1", pdata$treatment),
                         cell = sub("cell line: (.*)", "\\1", pdata$cell),
                         row.names = rownames(pdata))

pdataclean$dex <- ifelse(grepl('Dex', pdataclean$treatment), 'trt', 'untrt')
pdataclean$albut <- ifelse(grepl('Albut', pdataclean$treatment), 'trt', 'untrt')
pdataclean$SampleName <- rownames(pdataclean)
pdataclean$treatment <- NULL

# A informação que conecta as informações de amostras advindas de GEO com
# aquelas de SRA run id pode ser baixada usando Send to: File no NCBI.
srafile <- file.path(dir, 'SraRunInfo_SRP033351.csv')
srp <- read.csv(srafile)
srp

srpsmall <- srp[ , c("Run","avgLength","Experiment","Sample","BioSample","SampleName")]
srpsmall

# Estes dois data.frames são juntos e então nós fazemos um subset para apenas as amostras
# não tratadas com o albuterol (estas amostras não foram incluidas na publicação).
coldata <- merge(pdataclean, srpsmall, by = 'SampleName')
rownames(coldata) <- coldata$Run
coldata <- coldata[coldata$albut == 'untrt', ]
coldata$albut <- NULL
coldata

# Por fim, a tabela de amostras foi salva para um arquivo CSV para checagem futura.
# Este arquivo está incluso no diretório inst/extdata do pacote.
write.csv(coldata, file = 'sample_table.csv')
