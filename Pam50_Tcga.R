# Created Date: 18.Sept.2017
# Modified Date: 18.Sept.2017
# Owner: Samhita
# Modified By: Samhita
# Description: Dynamic script to perform PAM50 and prepare the data for further validation

# LIBRARIES REQUIRED
library(TCGAbiolinks)
library(genefu)
library(biomaRt)

# FUNCTIONS

## 1. Get the project ID
getProjectID = function(projectsList, tissue){
  projectID = projectsList[which(projectsList[[4]] == tissue),2]
  return(projectID)
}

## 2. Create tissue specific directory
createDirectory = function(subDir){
  mainDir = getwd()
  subDir = subDir
  dir.create(file.path(mainDir,subDir), showWarnings = FALSE)
  setwd(file.path(mainDir,subDir))
}

## 3. Data download and preparation
dataRetrieveAndPreparation = function(query){
  GDCdownload(query = query, method = "client")
  preparedData = GDCprepare(query)
  return(preparedData)
}

## 4. Get gene expression data
getGeneExp = function(data){
  geneExpData = SummarizedExperiment::assay(data)
  geneExpData = as.data.frame(t(geneExpData))
  return(geneExpData)
}

## 5. Get the annotation data from bioMart
getAnnotation = function(ensemblID){
  ensembl = useMart("ensembl")
  ensembleID = as.list(ensembleID)
  annotation = getBM(attributes = c("ensembl_gene_id", "external_gene_name", "entrezgene"),values = ensembleID, mart = ensembl)
  colnames(annotation)<-c("Probe", "Gene Symbol", "EntrezGene.ID")  
  annotation = subset(annotation, ensemblID %in% annotation$Probe)
  return(annotation)
}

## 6. Get the subtype classification with the given data
getPam50Subtypes = function(geneExp, annot){
  pam50Prediction = intrinsic.cluster.predict(sbt.model = pam50, data  = geneExp, annot = annot, do.mapping=TRUE, verbose=TRUE) 
  pam50Subtypes = as.data.frame(pam50Prediction$subtype)
  setDT(pam50Subtypes, keep.rownames = TRUE)[]
  colnames(pam50Subtypes) = c("Barcode","Classification_PAM50")
  return(pam50Subtypes)
}

## 7. Get the subtypes classified manually or from TCGA data
getTcgaSubtypes = function(tcgaData, tissue){
  # For breast cancer get the pam50 tcga classification directly from TCGA records
  if(tissue == "Breast"){
    classificationTcga = data.frame(preparedData$barcode,preparedData$subtype_PAM50.mRNA)
    colnames(classificationTcga) = c("Barcode","Classification_TCGA")
    levels(classificationTcga$Classification_TCGA) = c("Basal","Her2","LumA","LumB","Normal")
  }
  return(classificationtcga)
}

## 8. Merge both the classifications for further analysis and save to a file
saveSubtypeValidation = function(pam50, tcga, tissue){
  mergeData = merge(pam50, tcga, "Barcode") 
  setwd("..")
  filename = paste(tissue,"validationData","_")
  write.csv(mergeData, file = filename)
  return(mergeData)
}

## Main
tissue = readline(prompt = "Enter the tissue/ primary site: ")
tcgaProjects = getGDCprojects() # Get list of projects and their details
projectID = getProjectID(tcgaProjects,tissue) # Get the exact project ID
query = GDCquery(project = projectID$project_id, data.category = "Transcriptome Profiling", experimental.strategy = "RNA-Seq", workflow.type = "HTSeq - FPKM", legacy = FALSE) # Query the files based on the user requirements
createDirectory(projectID) # Create and set the path to a new directory for the particular tissue
preparedData = dataRetrieveAndPreparation(query) # Download the required data
geneExpData = getGeneExp(preparedData) # Get the gene expression data
ensemblGeneID = colnames(geneExpData)
annotationData = getAnnotation(ensemblGeneID) # Get the appropriate annotation as Entrez Gene ID
classificationPam50 = getPam50Subtypes(geneExpData, annotationData) # Get the subtype classification of the cancer
classificationTcga = getTcgaSubtypes(preparedData, tissue) # Get the subtype classification from tcga directly or manually curated
subtypeValidation = saveSubtypeValidation(classificationPam50, classificationTcga, tissue) # Merge both the classifications for further analysis and save to a file