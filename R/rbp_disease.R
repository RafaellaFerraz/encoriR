#' Get Data for RBP and Somatic Gene Mutations
#'
#'
#'  @description Retrieve data for RBP-gene interactions and somatic mutations in human diseases
#'
#' @param assembly =[genome version]: hg19
#' @param RBP =[protein name]. e.g., CBX7 ("all" for downloading all regulatory data)
#' @param tissue =[tissue name or organ name]. e.g., adrenal gland, breast
#' @param disease =[cancers and rare diseases in human]. e.g., Rosai-Dorfman disease, Crohn disease, Proteus syndrome
#' @param target =[gene name]. e.g., MYC ("all" for downloading all regulatory data)
#'
#' @return RBP	Name of RNA binding protein, e.g., ACIN1
#' @return geneID	ENSEMBL gene ID, e.g., ENSG00000136997
#' @return geneName	Name of gene, e.g., MYC
#' @return tissue	Tissue name or organ name, e.g., adrenal gland, breast
#' @return diseaseNum	The number of diseases
#' @return diseases	Name of cancers or rare diseases in human
#' @return diseaseCosmicID	The comma-separated mutations (COSMIC ID) corresponded to diseases
#' @return cosmicNum	The number of mutations (COSMIC ID)
#' @return sampleNum	Number of samples
#' @return mutTypeNum	Number of mutation types
#' @return clipExpNum	Number of supporting CLIP-seq experiments
#' @return clipNum	Number of CLIP-seq sites
#'
#' @export
#'
#' @examples
#'
#' disease_myc <- rbp_disease(tissue = c("prostate","breast"),disease = "carcinoma", target = "MYC")
#'
rbp_disease <- function(assembly="hg19",
                        RBP="all",
                        tissue,
                        disease,
                        target){
  links <- paste("https://starbase.sysu.edu.cn/api/RBPDisease/?",
                 "assembly=",assembly,
                 "&RBP=",RBP,
                 "&tissue=",tissue,
                 "&disease=",disease,
                 "&target=",target, sep = "")
  rbps <- NULL
  for (link in links){
    rbp <- utils::read.csv(url(link), comment.char = "#", sep = "\t", row.names = NULL)
    rbps <- rbind(rbps, rbp)
  }
  return(rbps)
}
