library("limma")
source("metaAnalyze.R")

library("biomaRt")
ensembl = useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host = "uswest.ensembl.org")
View(listDatasets(ensembl))
View(listFilters(ensembl))
q = listAttributes(ensembl)
View(q[grep("affy|agilent|illumina", q$name),])



#### AGILENT
#### 
parentFolderPath = "./Data_in"

studyName = "GSE100191"
x = read.maimages(files = dir(paste0(parentFolderPath, "/", studyName)), 
                  path = paste0(parentFolderPath, "/", studyName), 
                  source = "agilent", green.only = TRUE)
y = neqc(x, status = x$genes$ControlType, negctrl = -1, regular = 0)

colnames(y) = unname(fixGeoSampleNames(colnames(y)))

rownames(y$E) = y$genes$ProbeName
agilent_ensembl = c("agilent_wholegenome_4x44k_v2", "entrezgene")
mapping = getBM(attributes = agilent_ensembl, 
                mart = ensembl,
                # filters = "agilent_wholegenome_4x44k_v2",
                # values = rownames(y$E), 
                verbose = F,
                uniqueRows = T)
mapping = getGeneProbeMappingDirect(featureDf = mapping, 
                                    geneColname = "entrezgene", 
                                    probeColname = "agilent_wholegenome_4x44k_v2")

exprsByGene = calcExprsByGeneEmat(y$E, mapping)

if (any(is.na(exprsByGene))) {
  warning(sprintf('Imputing missing expression values for study %s.', studyName))
  resultImputed = impute.knn(exprsByGene)
  exprsByGene = resultImputed$data
}
# colnames(exprsByGene) = fixGeoSampleNames(colnames(exprsByGene))

eset <- new("ExpressionSet", exprs = exprsByGene)

saveRDS(eset, "./Data_in/GEO/GSE100191.rds")


#==================

studyName = "GSE68215"
# x = SCAN_TwoColor(inFilePattern = paste0(parentFolderPath, "/", studyName, "/*.txt.gz"))

x = read.maimages(files = dir(paste0(parentFolderPath, "/", studyName)),
                  path = paste0(parentFolderPath, "/", studyName),
                  source = "agilent", green.only = FALSE)

# y = x[,seq(1,72,2)]

x = new("EListRaw", list(E = x$R, 
                         Eb = x$Rb, 
                         targets = x$targets, 
                         genes = x$genes, 
                         source = x$source))

y = neqc(x, status = x$genes$ControlType, negctrl = -1, regular = 0)
colnames(y) = unname(fixGeoSampleNames(colnames(y)))

rownames(y$E) = y$genes$ProbeName
agilent_ensembl = c("agilent_wholegenome_4x44k_v1", "entrezgene")
mapping = getBM(attributes = agilent_ensembl, 
                mart = ensembl,
                verbose = F,
                uniqueRows = T)
mapping = getGeneProbeMappingDirect(featureDf = mapping, 
                                    geneColname = "entrezgene", 
                                    probeColname = "agilent_wholegenome_4x44k_v1")

exprsByGene = calcExprsByGeneEmat(y$E, mapping)
any(is.na(exprsByGene))

if (any(is.na(exprsByGene))) {
  warning(sprintf('Imputing missing expression values for study %s.', studyName))
  resultImputed = impute.knn(exprsByGene)
  exprsByGene = resultImputed$data
}
# colnames(exprsByGene) = fixGeoSampleNames(colnames(exprsByGene))

eset <- new("ExpressionSet", exprs = exprsByGene)

saveRDS(eset, "./Data_in/GEO/GSE68215.rds")


#=================

studyName = "GSE78068"
x = read.maimages(files = dir(paste0(parentFolderPath, "/", studyName)), 
                  path = paste0(parentFolderPath, "/", studyName), 
                  source = "agilent", green.only = TRUE)

y = neqc(x, status = x$genes$ControlType, negctrl = -1, regular = 0)
colnames(y) = unname(fixGeoSampleNames(colnames(y)))

rownames(y$E) = y$genes$ProbeName
agilent_ensembl = c("agilent_wholegenome_4x44k_v1", "entrezgene")
mapping = getBM(attributes = agilent_ensembl, 
                mart = ensembl,
                verbose = F,
                uniqueRows = T)
mapping = getGeneProbeMappingDirect(featureDf = mapping, 
                                    geneColname = "entrezgene", 
                                    probeColname = "agilent_wholegenome_4x44k_v1")

exprsByGene = calcExprsByGeneEmat(y$E, mapping)
any(is.na(exprsByGene))

if (any(is.na(exprsByGene))) {
  warning(sprintf('Imputing missing expression values for study %s.', studyName))
  resultImputed = impute.knn(exprsByGene)
  exprsByGene = resultImputed$data
}
# colnames(exprsByGene) = fixGeoSampleNames(colnames(exprsByGene))

eset <- new("ExpressionSet", exprs = exprsByGene)
saveRDS(eset, "./Data_in/GEO/GSE78068.rds")

### ILLUMINA ==========

studyName = "GSE37107"
q = read.ilmn(files = "GSE37107_non-normalized.txt.gz", path = paste0(parentFolderPath),
              probeid = "PROBE_ID", expr = "RIT",  other.columns = "Detection")
esetGEO = getGEO(filename = "./Data_in/GEO/GSE37107_series_matrix.txt.gz")
# pl = getGEO("GPL6947") 

y = neqc(q)
colnames(y) = esetGEO$geo_accession

mapping = getBM(attributes = c("illumina_humanht_12_v3", "entrezgene"), 
                mart = ensembl,
                verbose = F,
                uniqueRows = T)
mapping = getGeneProbeMappingDirect(featureDf = mapping, 
                                    geneColname = "entrezgene", 
                                    probeColname = "illumina_humanht_12_v3")

exprsByGene = calcExprsByGeneEmat(y$E, mapping)
any(is.na(exprsByGene))

if (any(is.na(exprsByGene))) {
  warning(sprintf('Imputing missing expression values for study %s.', studyName))
  resultImputed = impute.knn(exprsByGene)
  exprsByGene = resultImputed$data
}
# colnames(exprsByGene) = fixGeoSampleNames(colnames(exprsByGene))

eset <- new("ExpressionSet", exprs = exprsByGene)
saveRDS(eset, "./Data_in/GEO/GSE37107.rds")

##================

studyName = "GSE48348"
w = read.csv( "./Data_in/GEO/GSE48348_non-normalized.txt.gz", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
for (i in 1:(1468/2)) {
  cat("Substituting sample", i, "\n")
  w[1, 2*i] = paste0("Sample_S", i)
}


write_delim(w, path =  "./Data_in/GEO/GSE48348_non-normalized_corr.txt", delim = "\t", col_names = FALSE)


q = read.ilmn(files = "GSE48348_non-normalized_corr.txt", path = paste0(parentFolderPath),
              probeid = "ID_REF", expr = "Sample",  other.columns = "Detection")
esetGEO = getGEO(filename = "./Data_in/GEO/GSE48348_series_matrix.txt.gz")

y = neqc(q)
colnames(y) = esetGEO$geo_accession

mapping = getBM(attributes = c("illumina_humanht_12_v3", "entrezgene"), 
                mart = ensembl,
                verbose = F,
                uniqueRows = T)
mapping = getGeneProbeMappingDirect(featureDf = mapping, 
                                    geneColname = "entrezgene", 
                                    probeColname = "illumina_humanht_12_v3")

exprsByGene = calcExprsByGeneEmat(y$E, mapping)
any(is.na(exprsByGene))

if (any(is.na(exprsByGene))) {
  warning(sprintf('Imputing missing expression values for study %s.', studyName))
  resultImputed = impute.knn(exprsByGene)
  exprsByGene = resultImputed$data
}
# colnames(exprsByGene) = fixGeoSampleNames(colnames(exprsByGene))

eset <- new("ExpressionSet", exprs = exprsByGene)
saveRDS(eset, "./Data_in/GEO/GSE48348.rds")


#===========
studyName = "GSE12051"
# x = read.maimages(files=dir(paste0(parentFolderPath, "/", studyName)), path=paste0(parentFolderPath, "/", studyName),
#                   source="smd", columns = list(E = "Beadstudio_intensity"))

q = read.ilmn(files = dir(paste0(parentFolderPath, "/", studyName)), 
              path = paste0(parentFolderPath, "/", studyName),
              probeid = "ProbeSet_name", 
              expr = "Beadstudio_intensity")

mode(q$E) <- "numeric"
y = backgroundCorrect.matrix(q$E, method = "normexp", normexp.method = "mle", offset = 16)
y = log2(y)
y = normalizeQuantiles(y)


library("illuminaHumanv1.db")

x <- illuminaHumanv1ENTREZID
mapped_probes <- mappedkeys(x)
xx = as.list(x[mapped_probes])
z = xx[rownames(q)]
mapping = data.frame(unlist(z))
mapping$ID = rownames(mapping)

mapping = getGeneProbeMappingDirect(featureDf = mapping, 
                                    geneColname = "unlist.z.", 
                                    probeColname = "ID")


exprsByGene = calcExprsByGeneEmat(y, mapping)
any(is.na(exprsByGene))
exprsByGene = exprsByGene[complete.cases(exprsByGene),]

esetGEO = getGEO(filename = "./Data_in/GEO/GSE12051_series_matrix.txt.gz")

eset <- new("ExpressionSet", exprs = exprsByGene)

names = esetGEO$geo_accession#[which(esetGEO$description %in% colnames(eset))]
colnames(eset) = names



saveRDS(eset, "./Data_in/GEO/GSE12051.rds")

#=========
studyName = "GSE58795"
platformInfo = "gpl10379hursta2a520709custommmpm.cdf.gzcdf"
require(platformInfo, character.only=TRUE)
cwd = setwd(file.path(parentFolderPath, studyName))
# cdfFile = dir(pattern="cdf.gz")
# env = make.cdf.env(filename=cdfFile, compress=TRUE) #hursta2a520709cdf # gpl10379hursta2a520709custommmpm.cdf.gzcdf

require("gpl10379hursta2a520709custommmpm.cdf.gzcdf", character.only = TRUE)
celfiles = affy::list.celfiles()
eset = justRMA(filenames = celfiles,
               cdfname = platformInfo,
               normalize = TRUE,
               background = TRUE,
               verbose = TRUE,
               destructive = FALSE)

setwd(cwd)

sampleNames(eset) = fixGeoSampleNames(sampleNames(eset))

pl = getGEO("GPL10379")
  
geneColname = "EntrezGeneID"
probeColname = "ID" 

featureDf = pl@dataTable@table
idx = sapply(featureDf, is.factor)
featureDf[idx] = lapply(featureDf[idx], as.character) # convert factors into characters




mapping = getGeneProbeMappingDirect(featureDf, geneColname = geneColname, probeColname = probeColname)
exprsByGene = calcExprsByGene(eset, mapping)
sum(is.na(exprsByGene))

if (any(is.na(exprsByGene))) {
  warning(sprintf('Imputing missing expression values for study %s.', studyName))
  resultImputed = impute.knn(exprsByGene)
  exprsByGene = resultImputed$data
}
# colnames(exprsByGene) = fixGeoSampleNames(colnames(exprsByGene))

eset = ExpressionSet(assayData = exprsByGene, 
                     phenoData = phenoData(eset), 
                     experimentData = experimentData(eset))

saveRDS(eset, "./Data_in/GEO/GSE58795.rds")


#========
# SYNOVIUM

studyName = "GSE39340"
q = read.ilmn(files = "GSE39340_non_normalized.txt.gz", path = paste0(parentFolderPath),
              probeid = "ID_REF", expr = "AVG_Signal",  other.columns = "Detection")
esetGEO = getGEO(filename = "./Data_in/GEO/GSE39340_series_matrix.txt.gz")
# pl = getGEO("GPL6947") 

y = neqc(q)
colnames(y) = esetGEO$geo_accession

mapping = getBM(attributes = c("illumina_humanht_12_v4", "entrezgene"), 
                mart = ensembl,
                verbose = F,
                uniqueRows = T)
mapping = getGeneProbeMappingDirect(featureDf = mapping, 
                                    geneColname = "entrezgene", 
                                    probeColname = "illumina_humanht_12_v4")

exprsByGene = calcExprsByGeneEmat(y$E, mapping)
sum(is.na(exprsByGene))
exprsByGene = exprsByGene[complete.cases(exprsByGene),]

# if (any(is.na(exprsByGene))) {
#   warning(sprintf('Imputing missing expression values for study %s.', studyName))
#   resultImputed = impute.knn(exprsByGene)
#   exprsByGene = resultImputed$data
# }
# colnames(exprsByGene) = fixGeoSampleNames(colnames(exprsByGene))

eset <- new("ExpressionSet", exprs = exprsByGene)

saveRDS(eset, "./Data_in/GEO/GSE39340.rds")







############################

RG <- read.maimages(files = dir(paste0(parentFolderPath, "/", studyName)),
                    path = paste0(parentFolderPath, "/", studyName),
                    source = "genepix", 
                    other.columns = c("F635 SD","B635 SD","F532 SD",
                                    "B532 SD","B532 Mean","B635 Mean","F Pixels","B Pixels"))
RGmodel <- kooperberg(RG)


#=========
studyName = "GSE21537"
x = read.maimages(files = dir(paste0(parentFolderPath, "/", studyName)), 
                  path = paste0(parentFolderPath, "/", studyName), 
                  source = "genepix.mean", green.only = FALSE)

x = new("EListRaw", list(E = x$R, Eb = x$Rb, targets = x$targets, genes = x$genes, source = x$source))
colnames(x) = unname(fixGeoSampleNames(colnames(x)))
x$genes$ID = gsub(" ", "", x$genes$ID)

pl = getGEO("GPL7768")
featureDf = pl@dataTable@table
idx = sapply(featureDf, is.factor)
featureDf[idx] = lapply(featureDf[idx], as.character) # convert factors into characters
mapping = getGeneProbeMappingDirect(featureDf, 
                                    geneColname = "Entrez Gene", 
                                    probeColname = "ID")

y =  limma::backgroundCorrect(x,  method = "normexp", normexp.method = "mle", offset = 16)
y$E = log2(y$E)
y$E = normalizeBetweenArrays(y$E, method = "quantile")

rownames(y$E)  = y$genes$ID
exprsByGene = calcExprsByGeneEmat(y$E, mapping)
sum(is.na(exprsByGene))

if (any(is.na(exprsByGene))) {
  warning(sprintf('Imputing missing expression values for study %s.', studyName))
  resultImputed = impute.knn(exprsByGene)
  exprsByGene = resultImputed$data
}

eset <- new("ExpressionSet", exprs = exprsByGene)

saveRDS(eset, "./Data_in/GEO/GSE21537.rds")


## ====
studyName = "GSE3698"

x = read.maimages(files = dir(paste0(parentFolderPath, "/", studyName)), 
                  path = paste0(parentFolderPath, "/", studyName), 
                  source = "genepix.mean", green.only = TRUE)
colnames(x) = unname(fixGeoSampleNames(colnames(x)))
pl = getGEO("GPL3050")
featureDf = pl@dataTable@table
idx = sapply(featureDf, is.factor)
featureDf[idx] = lapply(featureDf[idx], as.character) # convert factors into characters
mapping = getGeneProbeMappingAnno(featureDf, 
                                  dbName = 'org.Hs.egREFSEQ2EG', 
                                  interName = "GB_LIST")


y =  limma::backgroundCorrect(x,  method = "normexp", normexp.method = "mle", offset = 16)
y$E = log2(y$E)
y$E = normalizeBetweenArrays(y$E, method = "quantile")

rownames(y$E) = y$genes$Name

exprsByGene = calcExprsByGeneEmat(y$E, mapping)
sum(is.na(exprsByGene))
exprsByGene = exprsByGene[complete.cases(exprsByGene),]

if (any(is.na(exprsByGene))) {
  warning(sprintf('Imputing missing expression values for study %s.', studyName))
  resultImputed = impute.knn(exprsByGene)
  exprsByGene = resultImputed$data
}

eset <- new("ExpressionSet", exprs = exprsByGene)

saveRDS(eset, "./Data_in/GEO/GSE3698.rds")


## ====

studyName = "GSE47726"
q = read.ilmn(files = "GSE47726_non-normalized.txt.gz", path = paste0(parentFolderPath),
              probeid = "ID_REF", expr = "SAMPLE",  other.columns = "Detection")
esetGEO = getGEO("GSE47726")
# pl = getGEO("GPL6947") 

y = neqc(q)
colnames(y) = esetGEO$GSE47726_series_matrix.txt.gz$geo_accession
# HumanWG6_V2
# options(error = rlang::entrace)
mapping = getBM(attributes = c("illumina_humanwg_6_v3", "hgnc_symbol"), 
                mart = ensembl,
                verbose = F,
                uniqueRows = T)
mapping = getGEO("GPL8234")
mapping = mapping@dataTable@table[,1:2]

mapping = getGeneProbeMappingDirect(featureDf = mapping, 
                                    geneColname = "SYMBOL", 
                                    probeColname = "ID")

exprsByGene = calcExprsByGeneEmat(y$E, mapping)
sum(is.na(exprsByGene))
exprsByGene = exprsByGene[complete.cases(exprsByGene),]


# if (any(is.na(exprsByGene))) {
#   warning(sprintf('Imputing missing expression values for study %s.', studyName))
#   resultImputed = impute.knn(exprsByGene)
#   exprsByGene = resultImputed$data
# }
# colnames(exprsByGene) = fixGeoSampleNames(colnames(exprsByGene))

eset <- new("ExpressionSet", exprs = exprsByGene)

saveRDS(eset, "./Data_in/GEO/GSE47726.rds")


## ====
# GSE47727
studyName = "GSE47727"
q = read.ilmn(files = "GSE47727_non_normalized.txt.gz", path = paste0(parentFolderPath),
              probeid = "ID_REF", expr = "Sample")
esetGEO = getGEO("GSE47727")
# pl = getGEO("GPL6947") 
q$E = q$E[complete.cases(q$E),]

y =  limma::backgroundCorrect(q,  method = "normexp", normexp.method = "mle", offset = 0)
y$E = log2(y$E)
y$E = normalizeBetweenArrays(y$E, method = "quantile")

colnames(y$E) = esetGEO$GSE47727_series_matrix.txt.gz$geo_accession

# plotDensity(y$E)
# plot(prcomp(t(y$E))$x)

# HumanHT-12 V3.0
mapping = getBM(attributes = c("illumina_humanht_12_v3", "entrezgene_id"), 
                mart = ensembl,
                verbose = F,
                uniqueRows = T)
mapping = getGeneProbeMappingDirect(featureDf = mapping, 
                                    geneColname = "entrezgene_id", 
                                    probeColname = "illumina_humanht_12_v3")

exprsByGene = calcExprsByGeneEmat(y$E, mapping)
sum(is.na(exprsByGene))
exprsByGene = exprsByGene[complete.cases(exprsByGene),]

meta = pData(esetGEO$GSE47727_series_matrix.txt.gz)
meta = meta[, c("geo_accession", "title", "source_name_ch1", "characteristics_ch1", "characteristics_ch1.1",
                "platform_id", "data_row_count", "age (yrs):ch1", "gender:ch1")]
meta$study = "GSE47727"

eset <- new("ExpressionSet", exprs = exprsByGene)
pData(eset) = meta
saveRDS(eset, "./Data_in/GEO/GSE47727.rds")

## ====
# GSE47728
studyName = "GSE47728"
q = read.ilmn(files = "GSE47728_non_normalized.txt.gz", path = paste0(parentFolderPath),
              probeid = "ID_REF", expr = "Sample")
esetGEO = getGEO("GSE47728")
# pl = getGEO("GPL6947") 
unique(which(is.na(q$E), arr.ind = T)[,1])
q$E = q$E[complete.cases(q$E),]

y =  limma::backgroundCorrect(q,  method = "normexp", normexp.method = "mle", offset = 0)
y$E = log2(y$E)
y$E = normalizeBetweenArrays(y$E, method = "quantile")

colnames(y$E) = esetGEO$GSE47728_series_matrix.txt.gz$geo_accession

# plotDensity(y$E)
# plot(prcomp(t(y$E))$x)

# HumanHT-12 V3.0
mapping = getBM(attributes = c("illumina_humanht_12_v4", "entrezgene_id"), 
                mart = ensembl,
                verbose = F,
                uniqueRows = T)
mapping = getGeneProbeMappingDirect(featureDf = mapping, 
                                    geneColname = "entrezgene_id", 
                                    probeColname = "illumina_humanht_12_v4")

exprsByGene = calcExprsByGeneEmat(y$E, mapping)
sum(is.na(exprsByGene))
exprsByGene = exprsByGene[complete.cases(exprsByGene),]

meta = pData(esetGEO$GSE47728_series_matrix.txt.gz)
# meta = meta[, c("geo_accession", "title", "source_name_ch1", "characteristics_ch1", "characteristics_ch1.1", "platform_id", "data_row_count", "age (yrs):ch1", "gender:ch1")]
meta$study = "GSE47728"


eset <- new("ExpressionSet", exprs = exprsByGene)
pData(eset) = meta
saveRDS(eset, "./Data_in/GEO/GSE47728.rds")

###############
####### Others
###############

##### GSE29746
eset = getGEO("GSE29746")
emat = exprs(eset[[1]])
meta = pData(eset[[1]])
meta$class = case_when(
  meta$`disease status:ch1` == "RASF (rheumatoid arthritis synovial fibroblasts)" ~ "RA",
  meta$`disease status:ch1` == "HSF (healthy synovial fibroblasts)" ~ "Healthy",
  meta$`disease status:ch1` == "OASF (osteoarthritis arthritis synovial fibroblasts)" ~ "OA"
)

colnames(meta)[!colnames(meta) %in% colnames(anno)][1] = "sample"
meta[,colnames(anno)  %in% colnames(meta)]
meta$study = "GSE29746"

featureDf = eset[[1]]@featureData@data
mapping = getGeneProbeMappingDirect(featureDf, 
                                    geneColname = "GENE", 
                                    probeColname = "ID")
y = neqc(emat, status = featureDf$CONTROL_TYPE, negctrl = "neg", regular = "FALSE", offset = 0)
exprsByGene = calcExprsByGeneEmat(y, mapping)
sum(is.na(exprsByGene))
exprsByGene = exprsByGene[complete.cases(exprsByGene),]
emat = exprsByGene

save(emat, meta, file = "./Data_in/GEO/GSE29746.Rdata")

eset <- new("ExpressionSet", exprs = exprsByGene)
pData(eset) = meta

saveRDS(eset, "./Data_in/GEO/GSE29746.rds")

# GSE19821

studyName = "GSE19821"
q = getGEO("GSE19821")

GSE19821_A = exprs(q[[1]])


# # GSE2053
# 
# studyName = "GSE2053"
# q = getGEO(filename = "Data_in/GSE2053_family.soft.gz")
# 
# meta = foreach(sample = names(q@gsms)) %do% {
#   q@gsms[[sample]]@header$characteristics_ch1
# }
# 
# int = foreach(sample = names(q@gsms)) %do% {
#   m = q@gsms[[sample]]@dataTable@table[,"SIGNAL_RAW"]
#   names(m) = q@gsms[[sample]]@dataTable@table[,"ID_REF"]
#   m
# }
# intB = foreach(sample = names(q@gsms)) %do% {
#   m = q@gsms[[sample]]@dataTable@table[,"BKD_RAW"]
#   names(m) = q@gsms[[sample]]@dataTable@table[,"ID_REF"]
#   m
# }
# names(int) = names(intB) = names(meta) = gsm
# 
# meta = do.call(rbind, meta)
# int = do.call(cbind, int)
# intB = do.call(cbind, intB)
# 
# gpl = q@gpls$GPL1740
# mapping = gpl@dataTable@table[, c("ID", "GB_ACC")]
# colnames(mapping) = c('probeSet', 'geneId')
# class(mapping[,1]) = class(mapping[,2]) = "character"
# mapping[mapping == ""] = NA
# mapping = mapping[complete.cases(mapping),]
# rownames(int) = rownames(intB) = NULL
# 
# mapping1 = getGeneProbeMappingAnno(mapping, dbName = 'org.Hs.egACCNUM2EG', interName = 'GB_ACC')
# 
# 
# targets = data.frame(colnames(int))
# rownames(targets) = targets[,1]
# x = new("EListRaw", list(E = int, 
#                          #Eb = intB, 
#                          targets = targets, 
#                          genes = gpl@dataTable@table, 
#                          source = "SMD Homo sapiens 37.6K Print_843"))
# 
# y =  limma::backgroundCorrect(x, method = "normexp", normexp.method = "mle", offset = 16)
# y$E = log2(y$E)
# y = normalizeBetweenArrays(y, method = "quantile")
# 
# rownames(y$E) = y$genes$ID
# cat("First 5 probe ids: ", rownames(y$E)[1:5], "\n")
# exprsByGene = calcExprsByGeneEmat(y$E, mapping)
# cat("Number of missing values =", sum(is.na(exprsByGene)), "\n")
# eset <- new("ExpressionSet", exprs = exprsByGene)
# message("Correction is done!")
# eset = fixMissingData(eset)

