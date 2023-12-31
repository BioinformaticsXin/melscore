% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Melanocyte_proportion_clinical.R
\name{Melanocyte_proportion_clinical}
\alias{Melanocyte_proportion_clinical}
\title{Title}
\usage{
Melanocyte_proportion_clinical(
  out_dir = NULL,
  theta = NULL,
  Clinical = NULL,
  plot_var = c("progression_index"),
  tag = NULL,
  cols = NULL
)
}
\arguments{
\item{out_dir}{Output results folder.}

\item{theta}{The cell proportion matrix calculated from the 'Melanocyte_proportion' function.}

\item{Clinical}{Clinical data, consisting of one column with row names corresponding to sample names.}

\item{plot_var}{The results to be plotted，for example: plot_var=c("progression_index",'Melanocyte_cluster4', 'Melanocyte_cluster3', 'Melanocyte_cluster2', 'Melanocyte_cluster1')}

\item{tag}{The objects to be compared, for example: tag=list(c("Stage 0","Stage II"), c("Stage 0","Stage IV"), c("Stage I","Stage II"))}

\item{cols}{Set color, for example: cols=c("Metastatic"="#7876b1", "nevus"="#ffdc91", "Primary"="#6f99ad", "Normal"="#ee4c97")}
}
\value{
data.frame
}
\description{
Title
}
