


### GSE8361 
eset = getGEO("GSE8361")[[1]]
emat = eset@assayData$exprs
meta = pData(eset)
meta$class = case_when(
  meta$source_name_ch1 == "Peripheral blood cells in healthy child" ~ "Healthy",
  meta$source_name_ch1 == "Peripheral blood cells in patient with systemic onset juvenile idiopathic arthritis" ~ "sJIA",
  meta$source_name_ch1 == "Peripheral blood cells in patient with poly-articular type juvenile idiopathic arthritis" ~ "polyJIA"
)

featureDf = eset@featureData@data
mapping = getGeneProbeMappingDirect(featureDf, 
                                    geneColname = "Entrez Gene ID", 
                                    probeColname = "ID")
exprsByGene = calcExprsByGeneEmat(emat, mapping)
sum(is.na(exprsByGene))
emat = exprsByGene

save(emat, meta, file = "./Data_in/GEO/GSE8361.Rdata")

eset <- new("ExpressionSet", exprs = exprsByGene)
pData(eset) = meta
saveRDS(eset, "./Data_in/GEO/GSE8361.rds")

load("./Data_in/GEO/GSE8361.Rdata")
ematGSE8361 = emat
metaGSE8361 = meta

plotDensity(exprsByGene)

ematGSE8361[is.na(ematGSE8361)] = 0
plot(prcomp(t(ematGSE8361))$x, col = as.numeric(as.factor(metaGSE8361$class)), pch = 16, cex = 2)

boxplot(exprsByGene)

### GSE11083 --- 
eset = getGEO("GSE11083")[[1]]
emat = eset@assayData$exprs
meta = pData(eset)
meta$tissue = case_when(
  meta$source_name_ch1 == "neutrophil" | meta$source_name_ch1 == "Neutrophil"  ~ "Neutrophil",
  meta$source_name_ch1 == "PBMC" ~ "PBMC"
)

meta$class = sapply(meta$title, function(x) strsplit(x[1], " ")[[1]][2])
meta$class[meta$class %in% c("control", "Control")] = "Healthy"
meta$class[meta$class == "JIA"] = "polyJIA"


y = backgroundCorrect.matrix(emat, method = "normexp", normexp.method = "mle", offset = 16)
y = log2(y)
y = normalizeQuantiles(y)

plotDensity(y)


featureDf = eset@featureData@data
# mapping = getGeneProbeMappingAnno(featureDf, 
#                                   dbName = 'org.Hs.egACCNUM2EG',
#                                   interName = "GB_ACC")
mapping = getGeneProbeMappingDirect(featureDf, 
                                    geneColname = "ENTREZ_GENE_ID", 
                                    probeColname = "ID")
exprsByGene = calcExprsByGeneEmat(y, mapping)
sum(is.na(exprsByGene))
emat = exprsByGene

save(emat, meta, file = "./Data_in/GEO/GSE11083.Rdata")

eset <- new("ExpressionSet", exprs = exprsByGene)
pData(eset) = meta
saveRDS(eset, "./Data_in/GEO/GSE11083.rds")

load("./Data_in/GEO/GSE11083.Rdata")
ematGSE11083 = emat
metaGSE11083 = meta

plot(prcomp(t(ematGSE11083))$x, col = as.numeric(as.factor(meta$class)), pch = 16)


### GSE13849 --- 
eset = getGEO("GSE13849")[[1]]
emat = eset@assayData$exprs
meta = pData(eset)

meta$class = sapply(meta$title, function(x) strsplit(x, "_")[[1]][2])
meta$subclass = case_when(
  meta$class == "CTRL" ~ "Healthy",
  meta$class == "PolyRF-" ~ "polyJIA_RF-",
  meta$class == "PolyRF+" ~ "polyJIA_RF+"
)
meta$class = ifelse(meta$subclass == "Healthy", "Healthy", "polyJIA")

featureDf = eset@featureData@data
mapping = getGeneProbeMappingDirect(featureDf, 
                                    geneColname = "ENTREZ_GENE_ID", 
                                    probeColname = "ID")
exprsByGene = calcExprsByGeneEmat(emat, mapping)
sum(is.na(exprsByGene))
emat = exprsByGene

save(emat, meta, file = "./Data_in/GEO/GSE13849.Rdata")

eset <- new("ExpressionSet", exprs = exprsByGene)
pData(eset) = meta
saveRDS(eset, "./Data_in/GEO/GSE13849.rds")

load("./Data_in/GEO/GSE13849.Rdata")
ematGSE13849 = emat
metaGSE13849 = meta

plot(prcomp(t(ematGSE13849))$x, col = as.numeric(as.factor(metaGSE13849$subclass)), pch = 16)


### GSE13501 --- 
eset = getGEO("GSE13501")[[1]]
emat = eset@assayData$exprs
meta = pData(eset)

meta$class = case_when(
  meta$source_name_ch1 == "CTRL" ~ "Healthy",
  meta$source_name_ch1 == "oligoarthritis" ~ "oligoJIA",
  meta$source_name_ch1 == "ERA" ~ "ERA",
  meta$source_name_ch1 == "RF- polyarthritis" ~ "polyJIA",
  meta$source_name_ch1 == "Systemic" ~ "sJIA"
)
meta$subclass = ifelse(meta$class == "polyJIA", "polyJIA_RF-", meta$class)

featureDf = eset@featureData@data
mapping = getGeneProbeMappingDirect(featureDf, 
                                    geneColname = "ENTREZ_GENE_ID", 
                                    probeColname = "ID")
exprsByGene = calcExprsByGeneEmat(emat, mapping)
sum(is.na(exprsByGene))
emat = exprsByGene

save(emat, meta, file = "./Data_in/GEO/GSE13501.Rdata")

eset <- new("ExpressionSet", exprs = exprsByGene)
pData(eset) = meta
saveRDS(eset, "./Data_in/GEO/GSE13501.rds")

load("./Data_in/GEO/GSE13501.Rdata")
ematGSE13501 = emat
metaGSE13501 = meta

plot(prcomp(t(emat))$x, col = as.numeric(as.factor(meta$subclass)), pch = 16)


### GSE15645 --- 
eset = getGEO("GSE15645")[[1]]
emat = eset@assayData$exprs
meta = pData(eset)


meta$class = case_when(
  meta$`disease state:ch1` == "polyarticular JIA" ~ "polyJIA",
  meta$`disease state:ch1` == "control" ~ "Healthy"
)
meta$subclass = case_when(
  meta$`symptom:ch1` == "active" ~ "polyJIA_active",
  meta$`symptom:ch1` == "CRM" ~ "polyJIA_CRM",
  meta$`symptom:ch1` == "CR" ~ "polyJIA_CR",
  meta$`symptom:ch1` == "control" ~ "Healthy"
)

featureDf = eset@featureData@data
mapping = getGeneProbeMappingDirect(featureDf, 
                                    geneColname = "ENTREZ_GENE_ID", 
                                    probeColname = "ID")
exprsByGene = calcExprsByGeneEmat(emat, mapping)
sum(is.na(exprsByGene))
emat = exprsByGene

save(emat, meta, file = "./Data_in/GEO/GSE15645.Rdata")

eset <- new("ExpressionSet", exprs = exprsByGene)
pData(eset) = meta
saveRDS(eset, "./Data_in/GEO/GSE15645.rds")

load("./Data_in/GEO/GSE15645.Rdata")
ematGSE15645 = emat
metaGSE15645 = meta

plot(prcomp(t(emat))$x, col = as.numeric(as.factor(meta$subclass)), pch = 16)



### GSE17755 --- 

eset = readRDS("./Data_in/GEO/GSE17755.rds")
meta = read.csv("./Data_in/GSE17755_meta.csv", stringsAsFactors = F, row.names = 1)
emat = eset@assayData$exprs
pData(eset) = meta


save(emat, meta, file = "./Data_in/GEO/GSE17755.Rdata")
saveRDS(eset, "./Data_in/GEO/GSE17755.rds")

load("./Data_in/GEO/GSE17755.Rdata")
ematGSE17755 = emat
metaGSE17755 = meta

plot(prcomp(t(emat))$x, col = as.numeric(as.factor(meta$class)), pch = 16)

### GSE26112 --- 
eset = getGEO("GSE26112")[[1]]
emat = eset@assayData$exprs
meta = pData(eset)

meta$class = "polyJIA"
meta$subclass = case_when(
  meta$`flare status:ch1` == "no flare" ~ "polyJIA_inactive",
  meta$`flare status:ch1` == "flare" ~ "polyJIA_active"
)
meta$time = rep(c("t0", "t1"), 17)

featureDf = eset@featureData@data
mapping = getGeneProbeMappingDirect(featureDf, 
                                    geneColname = "Entrez_Gene_ID", 
                                    probeColname = "ID")

exprsByGene = calcExprsByGeneEmat(emat, mapping)
sum(is.na(exprsByGene))
emat = exprsByGene

save(emat, meta, file = "./Data_in/GEO/GSE26112.Rdata")

eset <- new("ExpressionSet", exprs = exprsByGene)
pData(eset) = meta
saveRDS(eset, "./Data_in/GEO/GSE26112.rds")

load("./Data_in/GEO/GSE26112.Rdata")
ematGSE26112 = emat
metaGSE26112 = meta

plot(prcomp(t(emat))$x, col = as.numeric(as.factor(meta$time)), pch = 16, cex = 2)


### GSE20307 --- 
eset = getGEO("GSE20307")[[1]]
emat = eset@assayData$exprs
meta = pData(eset)

meta$subclass = case_when(
  meta$`original diagnosis:ch1` == "RF- polyarticular JIA" ~ "polyJIA_RF-",
  meta$`original diagnosis:ch1` == "healthy control" ~ "Healthy",
  meta$`original diagnosis:ch1` == "oligoarticular JIA" ~ "oligoJIA",
  meta$`original diagnosis:ch1` == "systemic JIA" ~ "sJIA"
)
meta$class = ifelse(meta$subclass == "polyJIA_RF-", "polyJIA", meta$subclass)


featureDf = eset@featureData@data
mapping = getGeneProbeMappingDirect(featureDf, 
                                    geneColname = "ENTREZ_GENE_ID", 
                                    probeColname = "ID")

exprsByGene = calcExprsByGeneEmat(emat, mapping)
sum(is.na(exprsByGene))
emat = exprsByGene

save(emat, meta, file = "./Data_in/GEO/GSE20307.Rdata")

eset <- new("ExpressionSet", exprs = exprsByGene)
pData(eset) = meta
saveRDS(eset, "./Data_in/GEO/GSE20307.rds")

load("./Data_in/GEO/GSE20307.Rdata")
ematGSE20307 = emat
metaGSE20307 = meta

plot(prcomp(t(emat))$x, col = as.numeric(as.factor(meta$class)), pch = 16, cex = 2)


### GSE55319 --- 
eset = getGEO(GEO = "GSE55319", AnnotGPL = TRUE)

#
emat = eset[[1]]@assayData$exprs
meta = pData(eset[[1]])

meta$class = "polyJIA"
meta$time = ifelse(meta$`month of visit:ch1` == 0, "t0", "t1")

featureDf = eset[[1]]@featureData@data
mapping = getGeneProbeMappingDirect(featureDf, 
                                    geneColname = "Gene ID", 
                                    probeColname = "ID")

exprsByGene = calcExprsByGeneEmat(emat, mapping)
sum(is.na(exprsByGene))
emat = exprsByGene

save(emat, meta, file = "./Data_in/GEO/GSE55319_A.Rdata")

eset1 <- new("ExpressionSet", exprs = exprsByGene)
pData(eset1) = meta
saveRDS(eset1, "./Data_in/GEO/GSE55319_A.rds")

load("./Data_in/GEO/GSE55319_A.Rdata")
ematGSE55319_A = emat
metaGSE55319_A = meta

plot(prcomp(t(emat))$x, col = as.numeric(as.factor(meta$class)), pch = 16, cex = 2)

#
emat = eset[[2]]@assayData$exprs
meta = pData(eset[[2]])

meta$class = case_when(
  meta$`group:ch1` == "healthy control" ~ "Healthy",
  meta$`group:ch1` == "juvenile idiopathic arthritis patient" ~ "polyJIA"
)
meta$time = "t0"

featureDf = eset[[2]]@featureData@data
mapping = getGeneProbeMappingDirect(featureDf, 
                                    geneColname = "Gene ID", 
                                    probeColname = "ID")

exprsByGene = calcExprsByGeneEmat(emat, mapping)
sum(is.na(exprsByGene))
emat = exprsByGene

save(emat, meta, file = "./Data_in/GEO/GSE55319_B.Rdata")

eset1 <- new("ExpressionSet", exprs = exprsByGene)
pData(eset1) = meta
saveRDS(eset1, "./Data_in/GEO/GSE55319_B.rds")

load("./Data_in/GEO/GSE55319_B.Rdata")
ematGSE55319_B = emat
metaGSE55319_B = meta

plot(prcomp(t(emat))$x, col = as.numeric(as.factor(meta$class)), pch = 16, cex = 2)


### GSE67596 --- 
eset = getGEO("GSE67596")[[1]]
emat = eset@assayData$exprs
meta = pData(eset)

meta$tissue = case_when(
  meta$`cell type:ch1` == "neutrophil" ~ "Neutrophil",
  meta$`cell type:ch1` == "PBMC" ~ "PBMC"
)
meta$class = case_when(
  meta$`disease state:ch1` == "Polyrticular JIA" ~ "polyJIA",
  meta$`disease state:ch1` == "Pauciarticular JIA" ~ "oligoJIA",
  meta$`disease state:ch1` == "Healthy control" ~ "Healthy"
)


y = backgroundCorrect.matrix(emat, method = "normexp", normexp.method = "mle", offset = 16)
y = log2(y)
y = normalizeQuantiles(y)

plotDensity(y)

featureDf = eset@featureData@data
mapping = getGeneProbeMappingDirect(featureDf, 
                                    geneColname = "ENTREZ_GENE_ID", 
                                    probeColname = "ID")

exprsByGene = calcExprsByGeneEmat(y, mapping)
sum(is.na(exprsByGene))
emat = exprsByGene

save(emat, meta, file = "./Data_in/GEO/GSE67596.Rdata")

eset <- new("ExpressionSet", exprs = exprsByGene)
pData(eset) = meta
saveRDS(eset, "./Data_in/GEO/GSE67596.rds")

load("./Data_in/GEO/GSE67596.Rdata")
ematGSE67596 = emat
metaGSE67596 = meta

plot(prcomp(t(emat))$x, col = as.numeric(as.factor(meta$class)), pch = 16, cex = 2)


### GSE112057 --- 

emat = read.delim(file = "./Data_in/GSE112057_Normalized_dataset.txt.gz", header = TRUE, stringsAsFactors = F)
genes = emat$gene
emat = emat[,-1] %>% as.matrix()
rownames(emat) = genes

eset = getGEO("GSE112057")[[1]]
meta = pData(eset)

# all(sapply(meta$title, function(x) strsplit(x, "_")[[1]][1]) == colnames(emat))
colnames(emat) = meta$geo_accession

colnames(meta)[47] = "tissue"

meta$class = case_when(
  meta$`disease state (diagnosis):ch1` == "Crohn's Disease" ~ "IBD",
  meta$`disease state (diagnosis):ch1` == "Ulcerative Colitis" ~ "UC",
  meta$`disease state (diagnosis):ch1` == "Polyarticular JIA" ~ "polyJIA",
  meta$`disease state (diagnosis):ch1` == "Oligoarticular JIA" ~ "oligoJIA",
  meta$`disease state (diagnosis):ch1` == "Systemic JIA" ~ "sJIA",
  meta$`disease state (diagnosis):ch1` == "Control" ~ "Healthy"
)

save(emat, meta, file = "./Data_in/GEO/GSE112057.Rdata")

eset <- new("ExpressionSet", exprs = exprsByGene)
pData(eset) = meta
saveRDS(eset, "./Data_in/GEO/GSE112057.rds")

load("./Data_in/GEO/GSE112057.Rdata")
ematGSE112057 = emat
metaGSE112057 = meta

plot(prcomp(t(emat))$x, col = as.numeric(as.factor(meta$class)), pch = 16, cex = 2)


######
polyJIA_exprs = list(GSE8361 = ematGSE8361, 
                     GSE11083 = ematGSE11083, 
                     GSE13849 = ematGSE13849, 
                     GSE13501 = ematGSE13501, 
                     GSE15645 = ematGSE15645, 
                     GSE17755 = ematGSE17755, 
                     GSE26112 = ematGSE26112, 
                     GSE20307 = ematGSE20307, 
                     GSE55319_A = ematGSE55319_A,
                     GSE55319_B = ematGSE55319_B,
                     GSE67596 = ematGSE67596, 
                     GSE112057 = ematGSE112057)

polyJIA_meta = list(GSE8361 = metaGSE8361, 
                    GSE11083 = metaGSE11083, 
                    GSE13849 = metaGSE13849, 
                    GSE13501 = metaGSE13501, 
                    GSE15645 = metaGSE15645, 
                    GSE17755 = metaGSE17755, 
                    GSE26112 = metaGSE26112, 
                    GSE20307 = metaGSE20307, 
                    GSE55319_A = metaGSE55319_A,
                    GSE55319_B = metaGSE55319_B,
                    GSE67596 = metaGSE67596, 
                    GSE112057 = metaGSE112057)

save(polyJIA_exprs, polyJIA_meta, file = "./Data_out/polyJIA_data.Rdata")

rm(ematGSE8361, ematGSE11083, ematGSE13849, ematGSE13501, ematGSE15645, ematGSE17755, ematGSE26112, ematGSE20307, ematGSE55319_A, ematGSE55319_B, ematGSE67596, ematGSE112057, 
   metaGSE8361, metaGSE11083, metaGSE13849, metaGSE13501, metaGSE15645, metaGSE17755, metaGSE26112, metaGSE20307, metaGSE55319_A, metaGSE55319_B, metaGSE67596, metaGSE112057)


