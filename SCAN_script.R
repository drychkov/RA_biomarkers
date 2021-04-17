require('foreach')
require('doParallel')
require("SCAN.UPC")


cl <- makePSOCKcluster(3)
registerDoParallel(cl, cores = 3)

# gseList="GSE34748"
# gseList = c('GSE11166', 'GSE14328', 'GSE34437', 'GSE50058', 'GSE72925')

parentFolderPath = "./Data_in/GEO"

studyMeta = read.csv("./Metadata/study_metadata.csv", stringsAsFactors = FALSE)
rownames(studyMeta) = studyMeta$study
studies = c("GSE68689", "GSE93272_A", "GSE93272_B")
studyMeta = studyMeta[studies,]
# studyName =studies[2]

# gseList = studyMeta$study
pkgName = NA
foreach(studyName = studyMeta$study, .packages = c("SCAN.UPC", "affy"), .verbose = TRUE) %dopar% {
    celFilePath = file.path(parentFolderPath, studyName, "*.CEL.gz")
    # celFilePath = file.path(parentFolderPath, studyName, "GSM635711.cel.gz")
    cat("Working on", studyName, "...", "\n")
    # pkgName = studyMeta[studyName, "platformInfo"]
    # cwd = setwd(file.path(parentFolderPath, studyName))
    files = affy::list.celfiles(file.path(parentFolderPath, studyName))
    if (is.na(pkgName)) {
      pkgName = InstallBrainArrayPackage(file.path(parentFolderPath, studyName, files[1]), "22.0.0", "hs", "entrezg")
      require(pkgName, character.only=TRUE)
    }
    
    outFilePath = file.path(parentFolderPath, "SCAN", paste0(studyName, "_SCAN.txt"))
    # normSet = SCANfast(celFilePath, probeSummaryPackage=pkgName, outFilePath)
    # normSet = SCAN(celFilePath, probeSummaryPackage = platformInfo, outFilePath, verbose = TRUE)
    normSet = SCAN(celFilePath, probeSummaryPackage = pkgName, outFilePath, verbose = TRUE)
    # normSet = SCAN(celFilePath, outFilePath)
    # normSet = SCAN(studyName, outFilePath)
    saveRDS(normSet, file = file.path(parentFolderPath, "SCAN", paste0(studyName, "_SCAN.rds")))
    
    pkgName = NA
    # setwd(cwd)
}
stopCluster(cl)

