#' Get MiRNA Target Data
#'
#' @description Retrieve data of the miRNA-target interactions by intersecting the predicting target sites of miRNAs with binding sites of Ago protein
#'
#'
#' @param assembly genome version: hg19, mm10
#' @param geneType main gene type: mRNA, lncRNA, pseudogene, circRNA, sncRNA
#' @param miRNA microRNA name. e.g., hsa-miR-21-5p; "all" for downloading all regulatory data
#' @param clipExpNum integer: minimum number of supporting CLIP-seq experiments
#' @param ampdegraExpNum integer: minimum number of supporting degradome-seq experiments
#' @param pancancerNum integer: minimum number of Cancer types. Pan-Cancer, <=32
#' @param programNum integer: minimum number of target-predicting programs. <=7
#' @param program string: target-predicting programs. e.g., PITA, RNA22, miRmap, DIANA-microT, miRanda, PicTar and TargetScan
#' @param target gene name. e.g., TP53 "all" for downloading all regulatory data.
#' @param cellType cell type. e.g., HeLa. "all" for downloading all regulatory data.
#'
#' @return A data frame object
#' @return miRNAid	MiRBase microRNA ID, e.g., MIMAT0002175
#' @return miRNAname	Name of microRNA, e.g., hsa-miR-485-5p
#' @return geneID	ENSEMBL gene ID, e.g., ENSG00000141510
#' @return geneName	Name of gene, e.g., TP53
#' @return geneType	Type of gene, e.g., protein_coding
#' @return chromosome	The name of the chromosome, e.g., chr3, chrY
#' @return narrowStart	The starting position of the overlapped part of stacked target region (across target-predicting programs)
#' @return narrowEnd	The ending position of the overlapped part of stacked target region (across target-predicting programs)
#' @return broadStart	The starting position of the full range of stacked target region (across target-predicting programs)
#' @return broadEnd	The ending position of the full range of stacked target region (across target-predicting programs)
#' @return strand	Defines the strand. Either "+" or "-"
#' @return clipExpNum	Number of supporting CLIP-seq experiments
#' @return degraExpNum	Number of supporting degradome-seq experiments
#' @return RBP	Name of RNA binding protein
#' @return PITA	Number of target sites predicted by PITA
#' @return RNA22	Number of target sites predicted by RNA22
#' @return miRmap	Number of target sites predicted by miRmap
#' @return microT	Number of target sites predicted by microT
#' @return miRanda	Number of target sites predicted by miRanda
#' @return PicTar	Number of target sites predicted by PicTar
#' @return TargetScan	Number of target sites predicted by TargetScan
#' @return pancancerNum	Number of Cancer types (r<0, p-value<0.05)
#'
#' @export
#'
#' @examples
#'
#' tp53_vdr <-  mirna_target(target = c("TP53","VDR"), geneType = "mRNA")
#' h19 <- mirna_target(target = "H19", geneType = "lncRNA")

mirna_target <- function(assembly="hg19",
                         geneType,
                         miRNA="all",
                         clipExpNum=NULL,
                         ampdegraExpNum=NULL,
                         pancancerNum=NULL,
                         programNum=NULL,
                         program=NULL,
                         target,
                         cellType="all"){
  miRNA_out <- NULL
  for (genT in geneType){
    for (mir in miRNA){
      for (tar in target){
        for (cell in cellType){
          link <- paste("https://starbase.sysu.edu.cn/api/miRNATarget/?",
                         "assembly=", assembly,
                         "&geneType=", genT,
                         "&miRNA=", mir,
                         "&clipExpNum=", clipExpNum,
                         "&ampdegraExpNum=",ampdegraExpNum,
                         "&pancancerNum=",pancancerNum,
                         "&programNum=", programNum,
                         "&program=",program,
                         "&target=", tar,
                         "&cellType=",cell, sep = "")

          miRNA_url <- utils::read.csv(url(link), comment.char = "#", sep = "\t", row.names = NULL)
          if(miRNA_url[1,1] != "The target parameter haven't been set correctly! Or the input of target parameter is not available!"
             & !is.na(miRNA_url[1,1])){
            miRNA_out <- rbind(miRNA_out, miRNA_url)
          }
        }
      }
    }
  }
  if (is.null(miRNA_out)){
    print("Data not available for selected targets or incorrect parameters")
  } else {
    BiocGenerics::colnames(miRNA_out)[1] <- "miRbase_ID"
    return(miRNA_out)
  }
}
