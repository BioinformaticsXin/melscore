#' Title
#'
#' @param out_dir Output results folder.
#' @param theta The cell proportion matrix calculated from the 'Melanocyte_proportion' function.
#' @param Survival Survival data, consisting of two columns, with the first column representing time and the second column representing status.
#'
#' @return list
#' @export
#'
#' @examples

Melanocyte_proportion_survival <- function(out_dir=NULL, theta=NULL, Survival=NULL){
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
	if (ncol(Survival) != 2) {
		stop("The size of Survival data is wrong. Please check inputs and selected regression type.")
	}
	else {
	    name <- colnames(Survival)[1]
	    colnames(Survival) <- c('OS.time','OS')
		select_sample <- intersect(rownames(Survival),rownames(theta))
		Survival <- cbind(Survival[select_sample,], theta[select_sample,])
		Survival$progression_index <- (1+(Survival$Melanocyte_cluster4/(Survival$Melanocyte_cluster3+Survival$Melanocyte_cluster4)))*((Survival$Melanocyte_cluster4+Survival$Melanocyte_cluster3)/(Survival$Melanocyte_cluster1+Survival$Melanocyte_cluster2+Survival$Melanocyte_cluster4+Survival$Melanocyte_cluster3))
		#################################OS HR##################################
		HR_result <- c()
		for(label in c('Melanocyte_cluster4','Melanocyte_cluster3','Melanocyte_cluster2','Melanocyte_cluster1','progression_index')){
			this.result <- hazard_ratios(inputArr=Survival, out_prefix=NULL, variable_col=label, time_col="OS.time", status_col="OS", forest_plot=FALSE, width=7, height=7)
			HR_result <- rbind(HR_result, this.result)
		}	
		HR_result$HR <- log10(HR_result$HR)
		HR_result$HR_CI_025 <- log10(HR_result$HR_CI_025)
		HR_result$HR_CI_975 <- log10(HR_result$HR_CI_975)
		HR_result$p_col[HR_result$logrank_pvalue < 0.05 & HR_result$HR > 0] <- "Postive effect(P<0.05)"
		HR_result$p_col[HR_result$logrank_pvalue > 0.05 & HR_result$HR > 0] <- "Postive effect(P>=0.05)"
		HR_result$p_col[HR_result$logrank_pvalue >= 0.05 & HR_result$HR < 0] <- "Negtive effect(P>=0.05)"
		HR_result$p_col[HR_result$logrank_pvalue < 0.05 & HR_result$HR < 0] <- "Negtive effect(P<0.05)"
		HR_result$p[HR_result$logrank_pvalue >= 0.05] <- ""
		HR_result$p[HR_result$logrank_pvalue < 0.05] <- "*"
		p <- ggplot(HR_result)+
		  geom_hline(yintercept = 0, linewidth = 0.3)+
		  geom_linerange(aes(Object, ymin = HR_CI_025, ymax = HR_CI_975, color = p_col), show.legend = F)+
		  geom_point(aes(Object, HR, color = p_col)) +
		  geom_text(aes(Object, y = HR_CI_975 + 0.05, label = p, color = p_col), show.legend = F)+
		  scale_color_manual(name = "", values = c("Postive effect(P<0.05)" = "#d55e00", "Postive effect(P>=0.05)" = "#ffbd88", 
							 "Negtive effect(P<0.05)" = "#0072b2", "Negtive effect(P>=0.05)" = "#7acfff"))+
		  xlab("")+ ylab("log10(HR)")+ theme_bw()+ coord_flip()
		ggsave(p, file=paste0(out_dir,name,"_HR.pdf"), height = 2.5, width = 6)
		##############progression index与survival###################
		progression_group <- separatePatient(input=Survival, var_cols="progression_index", time="OS.time", event="OS", labels=c("progression_index-Low", "progression_index-High"), minprop=0.1, maxstat=FALSE, quantile_cutoff=0.5)
		Survival <- data.frame(Survival, group=progression_group$Categorize_Arr$progression_index, stringsAsFactors=FALSE)
		fit<-survfit(Surv(OS.time,OS)~group,data=Survival)
		pdf(paste0(out_dir,name,"_progression_group.pdf"), width=9, height=7, onefile=FALSE)
		p<-ggsurvplot(fit, data = Survival,
		   conf.int = FALSE,
		   pval = TRUE,
		   fun = "pct",
		   risk.table = TRUE,
		   size = 1,
		   linetype = "strata",
		   palette = rev(c(brewer.pal(7,"Set1"))[1:2]),
		   legend = "top")
		print(p)
		dev.off()
    }
	return (list(Survival,HR_result))
}



hazard_ratios <- function(inputArr=NULL, out_prefix=NULL, variable_col=NULL, time_col="OS_Time", status_col="OS_Status", forest_plot=FALSE, width=7, height=7){
  library(survminer)
  library(survival)
  formula_used <- NA
  if(length(variable_col)>1){
    formula_used <- as.formula(paste0("Surv(", time_col, ", ", status_col, ") ~ ", paste(variable_col, collapse="+")))
  }else{
    formula_used <- as.formula(paste0("Surv(", time_col, ", ", status_col, ") ~ ", variable_col))
  }
  print(formula_used)
  coxph_fit_res <- coxph(formula=formula_used, data=inputArr)
  if(forest_plot && !is.null(out_prefix)){
    outfile <- paste0(out_prefix, ".forest_plot.pdf")
    checkDir(dirname(outfile))
    print("Output file is")
    print(outfile)
    ggforest_obj <- ggforest(coxph_fit_res, data = inputArr)
    ggsave(filename=outfile, plot=ggforest_obj, width=width, height=height)
  }
  HR <- exp(coxph_fit_res$coefficients)
  HR_CI <- exp(confint(coxph_fit_res))
  res <- data.frame(Object=rownames(HR_CI), HR=HR, HR_CI, logrank_pvalue=summary(coxph_fit_res)$sctest["pvalue"], wald_pvalue=summary(coxph_fit_res)$waldtest["pvalue"], Likelihood_pvalue=summary(coxph_fit_res)$logtest["pvalue"])
  colnames(res)[3:4] <- c("HR_CI_025", "HR_CI_975")
  return(res)
}



separatePatient <- function(input=NULL, var_cols=NULL, time="time", event="event", labels=c("Low", "High"), minprop=0.1, maxstat=FALSE, quantile_cutoff=0.5){
  if(maxstat){
    library(survminer)
    res.cut <- surv_cutpoint(data=input, time=time, event=event, variables=var_cols, minprop=minprop, progressbar = TRUE)
    res.cat <- surv_categorize(res.cut, labels=labels)
    res.cat[, var_cols] <- factor(as.character(res.cat[, var_cols]), levels=labels)
    return(list(Categorize_Arr=res.cat, Cutpoint=res.cut))
  }
  if(!maxstat && !is.null(quantile_cutoff)){
    group_label <- rep(labels[1], dim(input)[1])
    cutpoint <- quantile(input[, var_cols], probs=quantile_cutoff)
    if(length(labels)==2){
      group_label[which(input[, var_cols]>cutpoint)] <- labels[2]
    }
    if(length(labels)==3){
      group_label[which(input[, var_cols]>cutpoint[1])] <- labels[2]
      group_label[which(input[, var_cols]>cutpoint[2])] <- labels[3]
    }
    res.cat <- data.frame(input[, c(time, event)], group_label=group_label, check.names=FALSE)
    colnames(res.cat)[3] <- var_cols
    res.cat[, var_cols] <- factor(as.character(res.cat[, var_cols]), levels=labels)
    return(list(Categorize_Arr=res.cat, Cutpoint=cutpoint))
  }
}