#'Get the Data for NcRNA-RNA Interaction Network.
#'
#' @description Retrieve the interaction network of ncRNA-RNA identified from high-throughput sequencing data of RNA-RNA interactome
#'
#' @param assembly =[genome version]: hg38, mm10
#' @param geneType =[main gene type]: mRNA, lncRNA, pseudogene, sncRNA, miRNA
#' @param RNA =[RNA name]. e.g., TP53 ("all" for downloading all regulatory data)
#' @param interNum =[integer]: miminum number of RNA-RNA interactions
#' @param expNum =[integer]: miminum number of experiments
#' @param cellType =[cell type]. e.g., HeLa ("all" for downloading all regulatory data)
#'
#' @return geneID	ENSEMBL gene ID, e.g., ENSG00000141510
#' @return geneName	Name of gene, e.g., TP53
#' @return geneType	Type of gene, e.g., protein_coding
#' @return pairGeneID	ENSEMBL gene ID, e.g., ENSG00000002834
#' @return pairGeneName	Name of gene, e.g., LASP1
#' @return pairGeneType	Type of gene, e.g., protein_coding
#' @return interactionNum	Number of RNA-RNA interactions
#' @return totalExpNum	Number of total experiments
#' @return totalSeqTypeNum	Number of total types of high-throughput sequencing
#' @return totalReadsNum	Number of total reads from all samples
#' @return interactionLocus	Pairing locus of RNA-RNA interactions
#' @return alignment	The consequential pairing of two genes
#' @return expNum	Number of supporting experiments
#' @return seqTypeNum	Number of types of high-throughput sequencing
#' @return readsNum	Number of reads
#' @return FreeEnergy	The minimum free energy of RNA-RNA pairs were calculated by RNAfold software
#' @return AlignScore(Smith-Waterman)	Calculated by Smith-Waterman algorithm. The analysis results with alignment score greater than 10 were kept and stored in ENCORI
#' @return pancancerNum	Number of Cancer types (r<0, p-value<0.05)
#'
#' @export
#'
#' @examples
#' int_tp_brca<- rna_rna(geneType = "mRNA", RNA = c("TP53", "BRCA1"))
#' int_h19 <- rna_rna(geneType = "lncRNA", RNA = "H19")
#'
rna_rna <- function(assembly="hg38",
                    geneType,
                    RNA,
                    interNum=1,
                    expNum=1,
                    cellType="all"){
  rna_interactions <- NULL
  for (gen in geneType){
    for (rna in RNA){
      for (cell in cellType){
        link <- paste("https://rnasysu.com/encori/api/RNARNA/?",
                      "assembly=",assembly,
                      "&geneType=",gen,
                      "&RNA=",rna,
                      "&interNum=",interNum,
                      "&expNum=",expNum,
                      "&cellType=",cell, sep = "")

        rna_int <- utils::read.csv(url(link), comment.char = "#", sep = "\t", row.names = NULL)

        if(rna_int[1,1] != "The RNA parameter haven't been set correctly! Or the input of RNA parameter is not available!"
           & !is.na(rna_int[1,1])){
          rna_interactions <- rbind(rna_interactions, rna_int)
        }

      }
    }
  }
  if (is.null(rna_interactions)){
    print("Data not available or incorrect parameters")
  } else{
    BiocGenerics::colnames(rna_interactions)[1] <- "miRbase_ID"
    return(rna_interactions)
  }
}
