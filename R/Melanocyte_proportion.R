#' Title
#'
#' @param bk.dat RNA-seq sequencing data with count values, where column names are gene symbols and row names are sample IDs.
#' @param out_dir Output results folder.
#' @param key The character string that correponds to malignant cells in cell.type.labels. Set to NULL if there are no malignant cells or the malignant cells between reference and mixture are from matched sample, in which case all cell types will be treated equally, this parameter was derived from the 'new.prism' function in BayesPrism package.
#'
#' @return data.frame
#' @export
#'
#' @examples

Melanocyte_proportion <- function(bk.dat=NULL, out_dir=NULL, key="Melanocyte_cluster4"){
    library(RColorBrewer)
	library(ggplot2)
	library(ggthemes)
	library(Seurat)
	library(BayesPrism)
	library(reshape2)
	library(factoextra)
	library(cluster)
	library(NbClust)
	library(ggpubr)
	print("----------make sure the colnames are gene symbol-----------")
	myPrism <- new.prism(
	  reference=sc.dat.filtered.pc, 
	  mixture=bk.dat,
	  input.type="count.matrix", 
	  cell.type.labels = cell.type.labels, 
	  cell.state.labels = cell.state.labels,
	  key=key,
	  outlier.cut=0.01,
	  outlier.fraction=0.1,
	)	
    bp.res <- run.prism(prism = myPrism, n.cores=50)	
    theta <- get.fraction(bp=bp.res,
                       which.theta="final",
                       state.or.type="type")
	theta<-data.frame(theta,check.names=F)
    theta$progression_index <- (1+(theta$Melanocyte_cluster4/(theta$Melanocyte_cluster3+theta$Melanocyte_cluster4)))*((theta$Melanocyte_cluster4+theta$Melanocyte_cluster3)/(theta$Melanocyte_cluster1+theta$Melanocyte_cluster2+theta$Melanocyte_cluster4+theta$Melanocyte_cluster3))
	write.csv(theta,file=paste0(out_dir,"theta.csv"))
	return(theta)
}