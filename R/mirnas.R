#' Get elemnts which interact with miRNA
#'
#' @description Retrieve data of the miRNA-target interactions by intersecting the predicting target sites of miRNAs with binding sites of Ago protein
#'
#' @param assembly unique genome version: hg19, mm10
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
#' mRNA <-  targets_of_miRNAs(miRNA = c("hsa-miR-485-5p","hsa-miR-21-5p"), geneType = "mRNA")
#' lncRNA <- targets_of_miRNAs(miRNA = "hsa-miR-21-5p", geneType = "lncRNA")

targets_of_miRNAs <- function(assembly="hg19",
                              geneType,
                              miRNA,
                              clipExpNum=NULL,
                              ampdegraExpNum=NULL,
                              pancancerNum=NULL,
                              programNum=NULL,
                              program="None",
                              target="all",
                              cellType="all"){
  targets <- NULL
  for (mir in miRNA){
    for (tar in target){
      for (cell in cellType){
        for (genT in geneType){
          link <- paste("https://starbase.sysu.edu.cn/api/miRNATarget/?",
                        "assembly=",assembly,
                        "&geneType=",genT,
                        "&miRNA=",mir,
                        "&clipExpNum=",clipExpNum,
                        "&ampdegraExpNum=", ampdegraExpNum,
                        "&pancancerNum=",pancancerNum,
                        "&programNum=",programNum,
                        "&program=",program,
                        "&target=",tar,
                        "&cellType=",cell, sep = "")
          target_url <- utils::read.csv(base::url(description = link), comment.char = "#", sep = "\t", row.names = NULL)
          BiocGenerics::colnames(target_url)[1] <- "miRbase_ID"
          if(target_url[1,1] != "The miRNA parameter haven't been set correctly! Or the input of miRNA parameter is not available!" & !is.na(target_url[1,1])){
            targets <- rbind(targets, target_url)
          }
        }
      }
    }
  }
  if (is.null(targets)){
    print("Data not available for selected miRNAs or incorrect parameters")
  } else {
    return(targets)
  }
}

