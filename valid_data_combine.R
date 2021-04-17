source("metaAnalyze.R")
source("analysis_addon.R")
library("readxl")

bestGenes = read.csv("./Data_out/bestGEnes_DAS_corr.csv", stringsAsFactors = F)[,1]

sampleMeta = read.csv("./Metadata/sampleMeta.csv", stringsAsFactors = F, row.names = 1)
anno = read_excel("./Metadata/GSM_anno.xlsx", sheet = 1) %>% as.data.frame()

# gpl = getGEO("GPL91")
# 
# bestGenes[!bestGenes %in% gpl@dataTable@table$`Gene Symbol`]
# humanSyno(genes = c("DRAM1" , "SAMSN1", "SQOR" )) %>% unlist() %in%  gpl@dataTable@table$`Gene Symbol`

## Synovium GSE1919 -------------------------------------------------------------------------------
# Microarray
emat = readRDS("./Data_in/GEO/GSE1919.rds")
emat = exprs(emat)
colnames(emat) = fixGeoSampleNames(colnames(emat))
rownames(emat) = fixCustomCdfGeneIds(rownames(emat))
emat = geneId2Symbol(emat)

bestGenes[!bestGenes %in% rownames(emat)] # "DRAM1"  "SAMSN1" "SQOR" missing
emat.GSE1919 = emat

meta.GSE1919 = anno[anno$study == 'GSE1919',]
rownames(meta.GSE1919) = meta.GSE1919$sample

all((rownames(meta.GSE1919) == colnames(emat.GSE1919)))

# gpl = getGEO("GPL91")
# 
# bestGenes[!bestGenes %in% gpl@dataTable@table$`Gene Symbol`]
# humanSyno(genes = c("DRAM1" , "SAMSN1", "SQOR" )) %>% unlist() %in%  gpl@dataTable@table$`Gene Symbol`


# Synovium GSE3698???? ---------------------------------------------------------------------------
# Microarray
# gpl = getGEO("GPL3050")
# bestGenes[!bestGenes %in% gpl@dataTable@table$GENE_SYMBOL]
# humanSyno(genes = bestGenes[!bestGenes %in% gpl@dataTable@table$GENE_SYMBOL]) %>% unlist() %in%  gpl@dataTable@table$GENE_SYMBOL

## Synovium GSE89408 ----------------------------------------------------------------------------
# RNA-seq
load("./Data_out/validSyno.RData")
emat.GSE89408 = validSyno
meta.GSE89408 = sampleSynoValid

all(rownames(meta.GSE89408) == colnames(emat.GSE89408))

# ## Synovial fibroblasts. GSE29746 ----------
# # Microarray
# load(file = "./Data_in/GSE29746.Rdata")
# 
# emat.GSE29746 = exprsByGene
# meta.GSE29746 = meta

# Whole Blood. GSE17755 ------------------------------------------------------------------------
# Microarray
eset = readRDS("./Data_in/GEO/GSE17755.rds")
emat.GSE17755 = exprs(eset)
emat.GSE17755 = geneId2Symbol(emat.GSE17755)
meta.GSE17755 = read.csv("./Data_in/GSE17755_meta.csv", stringsAsFactors = F, row.names = 1)
all(rownames(meta.GSE17755) == colnames(emat.GSE17755))

s = intersect(rownames(meta.GSE17755), colnames(emat.GSE17755))
meta.GSE17755 = meta.GSE17755[s,]
emat.GSE17755 = emat.GSE17755[,s]

# PBMC. GSE15573 --------------------------------------------------------------------------------
# Microarray
load(file = "./Data_in/GSE15573.Rdata")
emat.GSE15573 = geneId2Symbol(exprsByGene)
meta.GSE15573 = meta

all((rownames(meta.GSE15573) == colnames(emat.GSE15573)))

# PBMC. GSE90081 --------------------------------------------------------------------------------
# RNA-seq

load("Data_out/validPBMC.RData")
emat.GSE90081 = validPBMC
meta.GSE90081 = samplePBMCValid

all(rownames(meta.GSE90081) == colnames(emat.GSE90081))

# # White Blood Cells. GSE117769 ------------------------------------------------------------------
# # RNAseq 
# load("./Data_out/GSE117769.RData")
# emat.GSE117769 = emat
# meta.GSE117769 = meta117769



emat.validList = list(GSE89408 = emat.GSE89408,
                      GSE1919 = emat.GSE1919,
                      GSE90081 = emat.GSE90081,
                      GSE17755 = emat.GSE17755, 
                      GSE15573 = emat.GSE15573)
meta.validList = list(GSE89408 = meta.GSE89408,
                      GSE1919 = meta.GSE1919,
                      GSE90081 = meta.GSE90081,
                      GSE17755 = meta.GSE17755, 
                      GSE15573 = meta.GSE15573)

# Check correct gene symbols
lapply(names(emat.validList), function(x) findGenes(iterGenes, emat.validList[[x]]))

q = findGenes(iterGenes, emat.validList[[1]])
rownames(emat.validList[[1]])[match(q[names(q) != ""], rownames(emat.validList[[1]]))] = names(q[names(q) != ""])

q = findGenes(iterGenes, emat.validList[[3]])
rownames(emat.validList[[3]])[match(q[names(q) != ""], rownames(emat.validList[[3]]))] = names(q[names(q) != ""])

# save(emat.validList, meta.validList, file = "./Data_out/validation_datasets_list.RData")

# load(file = "./Data_out/validation_datasets_list.RData")


plotHeatmap(emat.validList[[4]][rownames(emat.validList[[4]]) %in% bestGenes,], meta.validList[[4]][meta.validList[[4]]$class %in% c("RA", "Healthy"), ], rowNames = T)

##############  xCell Analysis ----------------------------------------------------------
library("xCell")
cellList = read.csv("./Data_in/cells_families.csv", stringsAsFactors = F)
cells = cellList[cellList$Type != "Epithelial", "Cells"]
cells = cells[!cells %in% c("Skeletal muscle", "Smooth muscle", "Myocytes", "Mesangial cells")]


rnaseq = c(TRUE, FALSE, TRUE, FALSE, FALSE)#, TRUE, FALSE)
names(rnaseq) = names(emat.validList)

# cellList = read.csv("./Data_in/cells_families.csv", stringsAsFactors = F)

xcell.validList = list()
for (name in names(emat.validList)) {
  message("Evaluating ", name, " for cell type scores...\n")
  q = xCell::xCellAnalysis(emat.validList[[name]], rnaseq = FALSE, parallel.sz = 3, parallel.type = "FORK")
  message("Computing ", name, " p-values...\n")
  p = xCellSignifcanceBetaDist(q)
  colnames(p) = colnames(q)
  xcell.validList[[name]] = list(xcell = q, pval = p)
}

# save(xcell.validList, file = "./Data_out/xcell.validList.RData")

