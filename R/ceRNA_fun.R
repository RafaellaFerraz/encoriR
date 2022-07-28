#' Get Data of the CeRNA Networks
#'
#' @description Retrieve the ceRNA networks from thousands of interactions of miRNA-targets supported by CLIP-seq data.
#'
#' @param assembly =[genome version]: hg19, mm10
#' @param geneType =[main gene type]: mRNA, lncRNA, pseudogene
#' @param ceRNA =[gene name]. e.g., MYC ("all" for downloading all regulatory data)
#' @param miRNAnum =[integer]: miminum number of miRNAs
#' @param family =[miRNA family name]. e.g., miR-200bc-3p/429,miR-513c-5p/514b-5p("all" for downloading all regulatory data)
#' @param pval =[decimal]: maximum p-value (<=0.01).
#' @param fdr =[decimal]: maximum false dicovery rate (<=0.01)
#' @param pancancerNum =[integer]: minimum number of Cancer types (Pan-Cancer, <=32)
#'
#' @return geneID	ENSEMBL gene ID, e.g., ENSG00000148584
#' @return geneName	Name of gene, e.g., A1CF
#' @return geneType	Type of gene, e.g., protein_coding
#' @return ceRNAid	ENSEMBL gene ID, e.g., ENSG00000018510
#' @return ceRNAname	Name of gene, e.g., AGPS
#' @return ceRNAgeneType	Type of gene, e.g., protein_coding
#' @return hitMiRNAFamilyNum	Number of hit miRNA families
#' @return hitMiRNAFamily	Name of hit miRNA family
#' @return pval	P value
#' @return fdr	False dicovery rate
#' @return pancancerNum	Number of Cancer types (r<0, p-value<0.05)
#'
#' @export
#'
#' @examples
#' ceRNA <- ceRNA_fun(geneType = "mRNA", family = "miR-200bc-3p/429,miR-513c-5p/514b-5p")
#'
ceRNA_fun <- function(assembly="hg19",
                  geneType,
                  ceRNA="all",
                  miRNAnum=2,
                  family,
                  pval=0.01,
                  fdr=0.01,
                  pancancerNum=0){
  options(useFancyQuotes = F)
  link <- paste("https://starbase.sysu.edu.cn/api/ceRNA/?",
                "assembly=",assembly,
                "&geneType=", geneType,
                "&ceRNA=",ceRNA,
                "&miRNAnum=",miRNAnum,
                "&family=", dQuote(family),
                "&pval=",pval,
                "&fdr=",fdr,
                "&pancancerNum=", pancancerNum, sep = "")

  ceRNAs <- utils::read.csv(url(link), comment.char = "#", sep = "\t", row.names = NULL)
  colnames(ceRNAs)[1] <- "miRbase_ID"
  return(ceRNAs)
}
