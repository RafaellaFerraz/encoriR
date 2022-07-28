#' Get the binding sites of CLIP-seq.
#'
#' @description Download the binding sites of CLIP-seq in bed format
#'
#' @param assembly =[genome version]. hg19, mm10, dm6, ce10, sacCer3
#' @param datasetID =[uniq dataset ID]. e.g., SBDH27
#'
#' @return chrom - The name of the chromosome (e.g. chr3, chrY, chr2_random) or scaffold (e.g. scaffold10671).
#' @return chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
#' @return chromEnd - The ending position of the feature in the chromosome or scaffold.
#' @return name - Defines the name of the BED line. This label is displayed to the left of the BED line in the Genome Browser window when the track is open to full display mode or directly to the left of the item in pack mode.
#' @return score - A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray).
#' @return strand - Defines the strand. Either "." (=no strand) or "+" or "-".
#'
#'
#' @export
#'
#' @examples
#'
#' bind <- binding_sites(datasetID = "SBDH27")
#'
binding_sites <- function(assembly="hg19",
                          datasetID){
  links <- paste("https://starbase.sysu.edu.cn/api/bindingSite/?",
                 "assembly=",assembly,
                 "&datasetID=",datasetID,sep = "")

  bindings <- NULL
  for (link in links){
    binding <- utils::read.table(url(link),header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="", comment.char = "#", row.names = NULL)
    bindings <- rbind(bindings, binding)
  }
  colnames(bindings) <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand")
  return(bindings)

}
