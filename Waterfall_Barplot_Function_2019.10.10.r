#Jenny Smith

#March 31,2017 

#Purpose: To create a waterfall (barplot) of patient expression data, or LSC17 scores. 

setwd(file.path(PROJHOME,"2017.02.28_LSC17_Score"))

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


# https://drsimonj.svbtle.com/ordering-categories-within-ggplot2-facets
#Example Data - when and how to order 
# ordered.data <- dat %>% 
#                 arrange(desc(Percent.Expressors_GT.10TPM_AML),
#                         gene_name,
#                         Cytogenetic.Category.1, desc(TPM)) %>% 
#                  mutate(Order=row_number(), 
#                         Labels=factor(Labels,levels=unique(Labels)))

#This function is prelim. It takes in a order data frame so that the facet variable
#will have the waterfall plot ordered be increasing expn, NOT USI (which looks a mess).
#However, I need to use eval() or similar with facet_wrap(~COLUMN) which I havent figured out yet.

waterfall_facets <- function(ordered.data,fill.col){
  ggplot(data = ordered.data,
         aes_string(x="Order", y="TPM",
                    fill=fill.col, 
                    color=fill.col)) +
    geom_bar(stat ="identity") +
    facet_wrap(~Labels, scales = 'free') +
    geom_hline(aes(yintercept = 10), 
               color="black",linetype="dashed") +
    scale_fill_brewer(type="qual",
                      palette = "Paired",
                      aesthetics = c("colour", "fill")) +
    theme(plot.title = element_text(hjust = 0.5, size = 16),
          strip.background = element_rect(fill="white"),
          strip.text = element_text(size=12),
          panel.background = element_rect(fill="white"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 14, color = "black"),
          axis.title = element_text(size = 14),
          plot.margin = margin(r = 1, t=1, l=1,unit="cm"),
          legend.text = element_text(size=11),
          legend.title = element_blank()) +
    theme(legend.position="bottom") +
    guides(fill = guide_legend(nrow = 4,ncol=3, byrow = TRUE)) +
    guides(color = guide_legend(nrow = 4,ncol=3, byrow = TRUE)) +
    labs(x="Patient", y="TPM", title="Expression of Plasma Membrane Associated Genes\n
               Highly Expressed in Pediatric AML")
}

#Example Ordered Data
#https://drsimonj.svbtle.com/ordering-categories-within-ggplot2-facets
dat <- full.TPM.matrix %>% 
  rownames_to_column("gene_id") %>%
  gather(USI,TPM,-gene_id) %>%
  mutate(Category=ifelse(TPM >= 10, 1,0)) %>%
  left_join(.,select(AMLGenes.stats,gene_id, gene_name, Percent.Expressors_GT.10TPM_AML,Median.TPM_AML),
            by="gene_id") %>% 
  left_join(., select(Groups, USI,Cytogenetic.Category.1:Categories),
            by=c("USI")) %>% 
  mutate_at(vars(Cytogenetic.Category.1,Cytogenetic.Category.2,Rare.Fusions),
            ~ifelse(. == "Unknown", "OtherAML", .)) %>%
  mutate(Cytogenetic.Category.1=factor(Cytogenetic.Category.1, levels=c("NBM","NBM.CD34",
                                                                        "inv.16.","MLL",
                                                                        "Normal","Other",
                                                                        "t.8.21.","OtherAML")),
         Cytogenetic.Category.2=factor(Cytogenetic.Category.2, levels=c("NBM","NBM.CD34",
                                                                        "del5q", "del7q" ,"del9q",
                                                                        "M6", "M7", "monosomy.7",
                                                                        "OtherAML")),
         Rare.Fusions=factor(Rare.Fusions, levels=c("NBM", "NBM.CD34","CBFA2T3.GLIS2",
                                                    "DEK.NUP214","NPM1.MLF1","NUP98.KDM5A",
                                                    "NUP98.NSD1","RBM15.MKL1",
                                                    "RUNX1.CBFA2T3","OtherAML")), 
         Labels=paste0(gene_name,"  (",Percent.Expressors_GT.10TPM_AML,"% Expn)")) %>%
  group_by(gene_name) %>%
  mutate(N=sum(Category == 1)) %>%
  ungroup()






