#' Get Data for MiRNAs Cleavage Events
#'
#' @description Retrieve data for cleavage events of miRNAs on genes supported by degradome-seq data
#'
#' @param assembly =[genome version]. e.g., hg19
#' @param geneType =[main gene type]: mRNA, ncRNA
#' @param miRNA =[microRNA name]. e.g., hsa-miR-196a-5p ("all" for downloading all regulatory data)
#' @param ampdegraExpNum =[integer]: minimum number of supporting degradome-seq experiments
#' @param target =[gene name]. e.g., TP53 ("all" for downloading all regulatory data)
#' @param cellType =[cell type]. e.g., HeLa ("all" for downloading all regulatory data)
#'
#' @return A data frame object
#' @return miRNAid	MiRBase microRNA ID, e.g., MIMAT0002175
#' @return miRNAname	Name of microRNA, e.g., hsa-miR-485-5p
#' @return geneName	Name of gene, e.g., TP53
#' @return geneType	Type of gene, e.g., protein_coding
#' @return cleaveEventNum	Number of cleavage events
#' @return degraExpNum	Number of supporting degradome-seq experiments
#' @return degraSiteNum	Number of cleavage sites from supporting degradome-seq experiments
#' @return totalReads	Number of total reads from all samples
#' @return category	Defined by starscan program
#'
#' @export
#'
#' @examples
#'
#'dd_tp_brca <- degradome_RNA(geneType = "mRNA", target = c("TP53", "BRCA1"))
#'dd_h19 <- degradome_RNA(target = "H19", geneType = "ncRNA")


degradome_RNA <- function(assembly="hg19",
                          geneType,
                          miRNA="all",
                          ampdegraExpNum=1,
                          target,
                          cellType="all"){
  degradomes <- NULL
  for (gen in geneType){
    for (mir in miRNA){
      for(tar in target){
        for (cell in cellType){
          link <- paste("https://starbase.sysu.edu.cn/api/degradomeRNA/?",
                         "assembly=",assembly,
                         "&geneType=",gen,
                         "&miRNA=",mir,
                         "&ampdegraExpNum=",ampdegraExpNum,
                         "&target=",tar,
                         "&cellType=",cell, sep = "")
          degradome <- utils::read.csv(url(link), comment.char = "#", sep = "\t", row.names = NULL)
          if(degradome[1,1] != "The target parameter haven't been set correctly! Or the input of target parameter is not available!"
             & !is.na(degradome[1,1])){
            degradomes <- rbind(degradomes, degradome)
          }
        }
      }
    }
  }
  if (is.null(degradomes)){
    print("Data not available for selected targets")
  } else{
    BiocGenerics::colnames(degradomes)[1] <- "miRbase_ID"
    return(degradomes)
  }
}
