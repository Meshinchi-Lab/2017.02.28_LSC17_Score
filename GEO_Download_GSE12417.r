#Jenny Smith

#March 13, 2017 


#Purpose: Download and analyze publically available mircroarray from NCBI GEO 

# References: 
#   Metzeler, K. H., Hummel, M., Bloomfield, C. D., Spiekermann, K., Braess, J., Sauerland, M.-C., … Buske, C. (2008). An 86-probe-set gene-expression signature predicts survival in cytogenetically normal acute myeloid leukemia. Blood, 112(10), 4193–4201. http://doi.org/10.1182/blood-2008-02-134411

# source("http://www.bioconductor.org/biocLite.R")
# biocLite("GEOquery")

library(parallel)
library(glmnet)
library(genefilter)
library(Biobase)
library(GEOquery)
library(stringr)


setwd("H:/RNA_seq_Analysis/2017.02.28_LSC17_Score/")

############Access GSE Class ##################

#gene series  gse12417 direct download using ftp  
gse12417 <- getGEO("gse12417", destdir = ".")


#save the object in order to avoid downloading again (takes ~2-3 min)
save(gse12417, file = "gse12417_CN-AML_series.RData")

#ExpressionSet class  
str(gse12417)

#print the datasets contained in the series 
names(gse12417)

#the number of 
length(gse12417)

#Lookat the information in each of the expression set matrices. 
#Note: these are already an "ExpressionSet" class. 
for (i in 1:length(gse12417)){
  print(head(gse12417[[i]]))
}


#Determine how the data was processed - access the data_processing column from each gene set. 
#Variance stabilizing normalization alogorithm. 
#probeset exprn vals calculated by median polish method. 
for (i in 1:length(gse12417)){
  print(unique(pData(gse12417[[i]])$data_processing))
}



#Print the column names
for (i in 1:length(gse12417)){
  print(i)
  print(colnames(pData(phenoData(gse12417[[i]]))))
}


###### Access Phenotype Data #############

#Access the phenotype data of each ExpressionSet and create a dataframe for the clinical data (from the "characteristics_ch1" column of each eset)
GetClinData <- function(ExpressionSet) {
  pheno_data <- pData(phenoData(ExpressionSet))
  
  clinData <- data.frame(pheno_data$title,
                         pheno_data$geo_accession,
                         str_split(as.character(pheno_data$characteristics_ch1), ";", simplify = TRUE))
  
  #add column names. 
  colnames(clinData) <- c("title", "geo_accession", "pheno_data", "Age_years", "OS_days", "vital_status")
  
  #split the column "pheno_data" into four additional columns. 
  clinData$pheno_data <- str_split_fixed(clinData$pheno_data, "[,()]", 4)
  colnames(clinData$pheno_data) <- c("cohort", "karyotype", "training_set", "FAB_group")
  
  #keep only the numeric information for the age, OS, and vital status.
  clinData$Age_years <- as.numeric(gsub(".+=([0-9]+) years", "\\1", clinData$Age_years))
  clinData$OS_days <- as.numeric(gsub(".+= ([0-9]+) days", "\\1", clinData$OS_days))
  clinData$vital_status <- as.numeric(gsub(".+: ([0,1])", "\\1", clinData$vital_status))

  return(clinData)
}


#call the function on the series matrix
gse12417_clinData_list <- lapply(gse12417, GetClinData)

head(gse12417_clinData_list)


#save the clinical data. 
write.csv(gse12417_clinData_list$`GSE12417-GPL570_series_matrix.txt.gz`, file = "gse12417_CN-AML_GPL570_testSet_ClinData.csv", row.names = FALSE)
write.csv(gse12417_clinData_list$`GSE12417-GPL96_series_matrix.txt.gz`, file = "gse12417_CN-AML_GPL96_trainingSet_clinData.csv",  row.names = FALSE)
write.csv(gse12417_clinData_list$`GSE12417-GPL97_series_matrix.txt.gz`, file = "gse12417_CN-AML_GPL97_trainingSet_clinData.csv", row.names = FALSE)


###### Access gene Expression Data #############

#Access gene expression data using exprs() and merge with official gene symbol for feature data using fData()
GetExpnData <- function(ExpressionSet){
  featureData <- fData(ExpressionSet) #feature data contains the probe IDs and the offical gene symbols 
  expnSet <- exprs(ExpressionSet) #Access the expression data. 
  expnSet <- merge(expnSet, featureData[,c(1,11)], by.x = 0, by.y = "ID") #merge df
  expnSet <- expnSet[,c(1, ncol(expnSet), 2:(ncol(expnSet)-1))] #rearrange columns
  
  return(expnSet)
}


gse12417_ExpnData_list <- lapply(gse12417, GetExpnData)


head(gse12417_ExpnData_list)


#write the expression data to a file
write.csv(gse12417_ExpnData_list$`GSE12417-GPL570_series_matrix.txt.gz`, file = "gse12417_CN-AML_GPL570_testSet_ExpnData.csv", row.names = FALSE)
write.csv(gse12417_ExpnData_list$`GSE12417-GPL96_series_matrix.txt.gz`, file = "gse12417_CN-AML_GPL96_trainingSet_ExpnData.csv", row.names = FALSE)
write.csv(gse12417_ExpnData_list$`GSE12417-GPL97_series_matrix.txt.gz`, file = "gse12417_CN-AML_GPL67_trainingSet_ExpnData.csv", row.names = FALSE)


#Need to Merge GLP96 and GPL97 (U133A and U133B chips). Will need to convert GSM# to ID (Patient 1-163)
patientIDs_A <- gse12417_clinData_list$`GSE12417-GPL96_series_matrix.txt.gz`[, 1:2]
patientIDs_B <- gse12417_clinData_list$`GSE12417-GPL97_series_matrix.txt.gz`[,1:2]


patientIDs_A_names <- gsub(".+ (P.+[0-9]{1,3}) .+", "\\1", patientIDs_A$title)
patientIDs_B_names <- gsub(".+ (P.+[0-9]{1,3}) .+", "\\1", patientIDs_B$title)



GSM_A <- colnames(gse12417_ExpnData_list$`GSE12417-GPL96_series_matrix.txt.gz`) 
GSM_B <- colnames(gse12417_ExpnData_list$`GSE12417-GPL97_series_matrix.txt.gz`)

gse12417_GPL96 <- gse12417_ExpnData_list$`GSE12417-GPL96_series_matrix.txt.gz`
gse12417_GPL97 <- gse12417_ExpnData_list$`GSE12417-GPL97_series_matrix.txt.gz`


colnames(gse12417_GPL96) <- c("probe_ID", "Gene_Symbol", patientIDs_A_names)
colnames(gse12417_GPL97) <- c("probe_ID", "Gene_Symbol", patientIDs_B_names)


gse12417_GPL96 <- rbind(GSM_A, gse12417_GPL96)
gse12417_GPL97 <- rbind(GSM_B, gse12417_GPL97)


gse12417_U133_geneChip <- rbind(gse12417_GPL96, gse12417_GPL97)
gse12417_U133_geneChip <- rbind(gse12417_U133_geneChip[22285,], gse12417_U133_geneChip)
gse12417_U133_geneChip <- gse12417_U133_geneChip[-22286, ]


dim(gse12417_U133_geneChip)
head(gse12417_U133_geneChip)



#save the file 
write.csv(gse12417_U133_geneChip, file = "gse12417_CN-AML_U133Combined_trainingSet_ExpnData.csv")
save(gse12417_U133_geneChip, file = "gse12417_CN-AML_U133Combined_trainingSet_ExpnData.RData")












