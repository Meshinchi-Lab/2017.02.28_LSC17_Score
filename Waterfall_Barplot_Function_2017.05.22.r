#Jenny Smith

#March 31,2017 

#Purpose: To create a waterfall (barplot) of patient expression data, or LSC17 scores. 
setwd(file.path(PROJHOME,"2017.02.28_LSC17_Score"))

GroupIDs <- function(clinData, col){
  #clindata has patient IDs as rownames. 
  #col is a chracter string of the column name in the clinical data with the factor/variable information. 
  list <- list()
  grps <- unique(clinData[,col])
  N <- length(grps)
  
  for (i in 1:length(grps)){
    if (grepl("[^a-zA-Z0-9 :]", grps[i])){
      grps[i] <- gsub("[^a-zA-Z0-9 :]", "\\.", grps[i]) #remove special characters and replace with a "."
    }
    IDs <- rownames(clinData[grepl(grps[i], clinData[,col]), ])
    list[[i]] <- IDs
  }
  
  names(list) <- grps
  return(list)
}

phenoVectors_MultipleGroups <- function(listofgoupsIDs){
  library(magrittr)
  #listofgoupsIDs contains a list with each item containing the IDs for each factor level.  
  #See GroupIDs function - this produced the input called  "listofgroupIDs"
  group <- names(listofgoupsIDs)
  
  vector <- NULL
  names <- NULL
  for (i in 1:length(listofgoupsIDs)){
    g <- group[i]
    vector <- c(vector, rep(g, length(listofgoupsIDs[[i]])))
    names <- c(names, listofgoupsIDs[[i]])
  }
  
  names(vector) <- names
  return(vector)
}


phenoVectors <- function(groupA, groupB){
  library(magrittr)
  #groupA and GroupB are character vectors with the patients IDs in each group
  g1 <- as.character(substitute(groupA))
  g2 <- as.character(substitute(groupB)) 
  
  vector <- c(rep(g1, length(groupA)), rep(g2, length(groupB)))
  names(vector) <- c(groupA, groupB)
  
  return(vector)
}

#Updated 4/6/18
waterfallPlot <- function(expnMatrix, geneName,phenoVector=NULL,unit,returnData=FALSE){
  #df is the gene expression data frame with patient IDs as a column names and genes as rownames. 
  #geneName is the vector of gene name
  #phenoVector is a character vector with patient IDs as names and the patietns status eg pos,neg
  #unit is a character vector for the 
  

  require(ggplot2)
  library(dplyr)
  
  #Select the gene of interest
  regex <- paste0("^", geneName, "$")
  expn <- expnMatrix %>%
    rownames_to_column("Gene") %>%
    filter(grepl(regex, Gene)) %>%
    gather(var,val, -Gene)
  
  #Create column for groups/phenos
  if(is.null(phenoVector)){
    
    phenoVector <- factor(ifelse(grepl("^RO|^BM", var), "NBM", "AML"), levels=c("NBM","AML"))
    expn <- expn %>% 
      mutate(Status=phenoVector)
  }else{
    
    #reorder if NBM is present in the phenovector
    if(any(grepl("NBM", phenoVector))){
      lvl <- c("NBM",unique(phenoVector)[-grep("NBM",unique(phenoVector))])
    }else{
      lvl <- unique(phenoVector)
    }
   
    phenoVector <- factor(phenoVector, levels=lvl)
    expn <- expn %>% 
      merge(.,data.frame(USI=names(phenoVector), Status=phenoVector) , by.x="var", by.y="USI") 
      # select(everything(), Status=y)
  }
  
  
  #arrange by group and value within group 
  expn <- expn %>% 
    arrange(Status,val)  %>%
    mutate(var=factor(var, levels = var))

  
  plot <- ggplot(data = expn, aes(x=var, y=val,fill=Status, color=Status)) +
    geom_bar(stat ="identity") +
    scale_fill_brewer(type="qual", palette = "Paired", aesthetics = c("colour", "fill")) +
    theme(plot.title = element_text(hjust = 0.5, size = 16),
          panel.background = element_rect(fill="white"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 16, color = "black"),
          axis.title = element_text(size = 14), 
          plot.margin = margin(r = 1, unit="cm"), 
          legend.text = element_text(size=11), 
          legend.title = element_blank()) +
    theme(legend.position="bottom") +
    guides(fill = guide_legend(nrow = 4,ncol=3, byrow = TRUE)) +
    guides(color = guide_legend(nrow = 4,ncol=3, byrow = TRUE)) +
    labs(x="Patient", y=unit, title=paste(geneName,"Expression"))
  
  
  if(returnData){
    res <- list("Data"=expn, "Plot"=plot)
  }else{
    res <- plot
  }
  
  return(res)
}









# waterfallPlot <- function(expnMatrix, geneName,phenoVector,unit,BM=FALSE){
#   #df is the gene expression data frame with patient IDs as a column names and genes as rownames. 
#   #geneName is the vector of gene name
#   #phenoVector is a character vector with patient IDs as names and the patietns status eg pos,neg
#   #unit is a character vector for the 
#   
#   library(colorspace)
#   library(RColorBrewer)
#   require(ggplot2)
#   library(magrittr)
#   
#   if(BM == TRUE){
#     BM <- rep("NBM", length(grep("^BM|^RO", colnames(expnMatrix)))) %>% setNames(grep("^BM|^RO", colnames(expnMatrix), value = TRUE))
#     phenoVector <- c(phenoVector, BM)
#   }else if (BM==FALSE){
#     phenoVector = phenoVector
#   }
#   
#   if (length(geneName) == 1){
#     expn <- expnMatrix[geneName, match(names(phenoVector), colnames(expnMatrix))] #correct order
#     
#     if (any(is.na(expn))){print("NAs Introduced. Check Column names match phenovector")}
#     
#     # expn <- expn[,which(! grepl("[[:alpha:]]", expn[1,]))] #special clause for GSE12417 epxn matrix
#     expn <- data.frame(t(expn), 
#                   Patient=names(expn),
#                   Status=phenoVector)
#     
#     colnames(expn)[1] <- unit
#   }else{
#     print("Please input only one gene at a time")
#   }
#   
#   N <- nrow(expn)
#   if (N < 50){
#     size=7
#   }else if (N >= 50 & N < 100){
#     size=5
#   }else if (N >= 100){
#     size=1
#   }
# 
#   plot <- ggplot(data = expn, aes(x=reorder(expn$Patient,expn[,unit]), y=expn[,unit], fill=Status)) +
#     geom_bar(stat ="identity") +
#     theme(plot.title = element_text(hjust = 0.5, size = 18),
#           panel.background = element_rect(fill="white"),
#           panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank(),
#           axis.text.x = element_text(angle = 45,hjust = 1, vjust = 0.5, size = size),
#           axis.text.y = element_text(size = 18, color = "black"),
#           axis.title = element_text(size = 16)) +
#     labs(x="Patient", y=unit, title="")
#   
#   # list <- list(expn, plot)
#   # names(list) <- c("expnData", "plo")
#   return(plot)
# }


