#' Get Data for RBP and Somatic Gene Mutations
#'
#'
#'  @description Retrieve data for RBP-gene interactions and somatic mutations in human diseases
#'
#' @param assembly =[unique genome version]: hg19
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
  rbps <- NULL
  for (tis in tissue){
    for (dis in disease){
      for (tar in target){
        for (rbp in RBP){
          link <- paste("https://starbase.sysu.edu.cn/api/RBPDisease/?",
                        "assembly=",assembly,
                        "&RBP=",rbp,
                        "&tissue=",tis,
                        "&disease=",dis,
                        "&target=",tar, sep = "")
          rbp <- utils::read.csv(url(link), comment.char = "#", sep = "\t", row.names = NULL)
          if(rbp[1,1] != "The target parameter haven't been set correctly! Or the input of target parameter is not available!"
             & !is.na(rbp[1,1])){
            rbps <- rbind(rbps, rbp)
          }

        }
      }
    }
  }
  if (is.null(rbps)){
    print("Data not available or incorrect parameters")
  } else{
    return(rbps)
  }
}
