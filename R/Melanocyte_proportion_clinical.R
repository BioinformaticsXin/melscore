#' Title
#'
#' @param out_dir Output results folder.
#' @param theta The cell proportion matrix calculated from the 'Melanocyte_proportion' function.
#' @param Clinical Clinical data, consisting of one column with row names corresponding to sample names.
#' @param plot_var The results to be plotted，for example: plot_var=c("progression_index",'Melanocyte_cluster4', 'Melanocyte_cluster3', 'Melanocyte_cluster2', 'Melanocyte_cluster1')
#' @param tag The objects to be compared, for example: tag=list(c("Stage 0","Stage II"), c("Stage 0","Stage IV"), c("Stage I","Stage II"))
#' @param cols Set color, for example: cols=c("Metastatic"="#7876b1", "nevus"="#ffdc91", "Primary"="#6f99ad", "Normal"="#ee4c97")
#'
#' @return data.frame
#' @export 
#'
#' @examples

Melanocyte_proportion_clinical <- function(out_dir=NULL, theta=NULL, Clinical=NULL, plot_var=c("progression_index"), tag=NULL, cols=NULL){
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
	library(survival)
	library(survminer)
	library(tidyverse)
	
	print("----------make sure the colnames are gene symbol-----------")   
    theta<-theta[,c('Melanocyte_cluster4', 'Melanocyte_cluster3', 'Melanocyte_cluster2', 'Melanocyte_cluster1')]
	if (ncol(Clinical) != 1) {
		stop("The size of Clinical data is wrong. Please check inputs and selected regression type.")
	}
	else {
		select_sample <- intersect(rownames(Clinical),rownames(theta))
		name <- colnames(Clinical)[1]
		Clinical <- cbind(Clinical[select_sample,], theta[select_sample,])
		Clinical$progression_index <- (1+(Clinical$Melanocyte_cluster4/(Clinical$Melanocyte_cluster3+Clinical$Melanocyte_cluster4)))*((Clinical$Melanocyte_cluster4+Clinical$Melanocyte_cluster3)/(Clinical$Melanocyte_cluster1+Clinical$Melanocyte_cluster2+Clinical$Melanocyte_cluster4+Clinical$Melanocyte_cluster3))
		colnames(Clinical) <- c('Type',colnames(Clinical)[-1])
		if(is.numeric(Clinical$Type)){
		    colnames(Clinical) <- c('Value',colnames(Clinical)[-1])
		    Clinical <- na.omit(data.frame(Clinical, Type = rep("High_group", nrow(Clinical)), stringsAsFactors = FALSE))
		    Clinical[which(Clinical$Value < median(Clinical$Value)),]$Type <- "Low_group"
		}
		Clinical$Type <- factor(Clinical$Type, levels = unique(Clinical$Type))
		p_res <- list()
		for (y_label in plot_var){
			Clinical$plots <-Clinical[,y_label]
			p_res[[y_label]]<-tagplot(Clinical, tag=tag, y_label, cols=cols)
		}
		pdf(paste0(out_dir, name, "_boxplot.pdf"),width=3*length(plot_var), height=5)
		  p<-do.call("ggarrange", c(p_res, ncol=length(plot_var), nrow=1))
		  print(p)
		dev.off()
	}
	return (Clinical)
}





tagplot <- function(Clinical=NULL, tag=NULL, y_label=NULL, cols=NULL){
    if (is.null(tag)){
	    tag <- unique(Clinical$Type)
	}	
    p <- ggplot(Clinical[Clinical$Type %in% unique(unlist(tag)),], aes(x = Type, y = plots, fill = Type, color = Type)) +
		geom_boxplot(width = .8,show.legend = F,
					position = position_dodge(0.9),
					alpha = 0.5,
					outlier.color = 'grey50') +
		geom_point(position=position_jitterdodge(jitter.width = 0.5)) +
		theme_classic(base_size = 14) +
		theme(axis.text.x = element_text(angle = 45,hjust = 1,color = 'black'), legend.position = '') +
	   stat_compare_means(label = "p.signif", comparisons=tag) +
	   labs(x = "", y = y_label)
    if (!is.null(cols)){
        p <- p + scale_fill_manual(values = cols) + scale_color_manual(values = cols)
    }
	return (p)
}