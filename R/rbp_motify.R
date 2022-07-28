#' Filter the binding motifs of RBPs contained specific sequence pattern
#'
#' @description Retrieve binding motifs of RBPs contained specific sequence pattern, which were de novo identified from their CLIP-seq peak data by Homer.
#'
#' @param assembly =[genome version]. hg19, mm10, dm6, ce10, sacCer3
#' @param motif =[motif pattern]. e.g., RRAC
#' @param rankLimit =[integer]. Rank limit of RNA binding motifs (default is 10)
#'
#' @export
#'
#' @return RBP Name of RNA binding protein, e.g., ACIN1
#' @return DatasetID
#' @return MotifRank
#' @return IdentifiedMotif
#' @return QueryMotif
#' @return targetPeakNum
#' @return TargetPercentage
#' @return p.value
#' @return p.value.ln
#' @return MotifMatrix
#' @return Region
#' @return CellLine.Tissue
#' @return Properties PAR-CLIP, Differentiated Glioma stem cells, rep1
#' @return MainAccession
#' @return SubAccession
#' @return Citation
#'
#'
#' @examples
#'
#' mot_rrac <- rbp_motify(motif = "RRAC")
#'
rbp_motify <- function(assembly="hg19",
                       motif,
                       rankLimit=10){
  links <- paste("https://starbase.sysu.edu.cn/api/RBPMotifScan/?",
  "assembly=",assembly,
  "&motif=",motif,
  "&rankLimit=",rankLimit, sep = "")


  motifys <- NULL
  for (link in links){
    motify <- utils::read.csv(url(link), comment.char = "#", sep = "\t", row.names = NULL)
    motifys <- rbind(motifys, motify)
  }
  return(motifys)
}
