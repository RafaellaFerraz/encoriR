#' Get RBP Target Data
#'
#' @description Retrieve data of RBP-RNA interactions supported by the binding sites of RBPs derived from CLIP-seq data
#'
#'
#' @param assembly =[unique genome version]. hg38, mm10, dm6, ce10, sacCer3
#' @param geneType =[main gene type]: mRNA, lncRNA, pseudogene, circRNA, sncRNA
#' @param RBP =[protein name]: CBX7 ("all" for downloading all regulatory data)
#' @param clipExpNum =[integer]: minimum number of supporting CLIP-seq experiments
#' @param pancancerNum =[integer]: minimum number of Cancer types (Pan-Cancer, <=32)
#' @param target =[gene name]. e.g., TP53 ("all" for downloading all regulatory data)
#' @param cellType =[cell type]. e.g., HeLa ("all" for downloading all regulatory data)
#'
#' @return RBP	Name of RNA binding protein
#' @return geneID	ENSEMBL gene ID, e.g., ENSG00000148584
#' @return geneName	Name of gene, e.g., A1CF
#' @return geneType	Type of gene, e.g., protein_coding
#' @return clusterNum	The number of stacked peak regions
#' @return totalClipExpNum	The number of total CLIP-seq experiments
#' @return totalClipSiteNum	The number of total supporting binding sites from CLIP-seq experiments
#' @return clusterID	The ID of stacked peak regions
#' @return chromosome	The name of the chromosome, e.g., chr3, chrY
#' @return narrowStart	The starting position of the overlapped part of stacked target region (across target-predicting programs)
#' @return narrowEnd	The ending position of the overlapped part of stacked target region (across target-predicting programs)
#' @return broadStart	The starting position of the full range of stacked target region (across target-predicting programs)
#' @return broadEnd	The ending position of the full range of stacked target region (across target-predicting programs)
#' @return strand	Defines the strand. Either "+" or "-"
#' @return clipExpNum	Number of supporting CLIP-seq experiments
#' @return HepG2	The change (log2FC) of gene while knocking down RBP in HepG2 cell lines
#' @return K562	The change (log2FC) of gene while knocking down RBP in K562 cell lines
#' @return pancancerNum	Number of Cancer types (r<0, p-value<0.05)
#'
#' @export
#'
#' @examples
#' rbp_tp53 <- rbp_target(geneType = "mRNA", target = "TP53")
#' rbps <- rbp_target(geneType = "mRNA", target = c("TP53", "BRCA1"))

rbp_target <- function(assembly="hg38",
                       geneType,
                       RBP="all",
                       pancancerNum=0,
                       clipExpNum = 5,
                       target,
                       cellType="all"){
  rbps_out <- NULL
  for(gen in geneType){
    for(rbp in RBP){
      for (tar in target){
        for(cell in cellType){
          link <- paste("https://rna.sysu.edu.cn/encori/api/RBPTarget/?",
                         "assembly=",assembly,
                         "&geneType=",gen,
                         "&RBP=",rbp,
                         "&clipExpNum=",clipExpNum,
                         "&pancancerNum=",pancancerNum,
                         "&target=",tar,
                         "&cellType=",cell, sep = "")
          rbp_url <- utils::read.csv(url(link), comment.char = "#", sep = "\t", row.names = NULL)
          if(rbp_url[1,1] != "The target parameter haven't been set correctly! Or the input of target parameter is not available!"
             & !is.na(rbp_url[1,1])){
            rbps_out <- rbind(rbps_out, rbp_url)
          }

        }
      }
    }
  }
  if (is.null(rbps_out)){
    print("Data not available or incorrect parameters")
  } else{
    return(rbps_out)
  }
}
