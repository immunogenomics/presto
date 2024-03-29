% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pseudobulk.R
\name{pseudobulk_pairwise}
\alias{pseudobulk_pairwise}
\title{Pseudobulk pairwise}
\usage{
pseudobulk_pairwise(
  dge_formula,
  counts_df,
  meta_data,
  contrast_var,
  vals_test,
  verbose,
  min_counts_per_sample,
  present_in_min_samples
)
}
\arguments{
\item{dge_formula}{differential gene expression formula for DESeq2}

\item{counts_df}{counts matrix}

\item{meta_data}{data.frame of cell metadata}

\item{contrast_var}{cell metadata column to use for differential
gene expression}

\item{vals_test}{cell metadata columns}

\item{verbose}{verbose}

\item{min_counts_per_sample}{minimum counts per sample to include in
differential gene expression}

\item{present_in_min_samples}{minimum samples with gene counts to
include in differential gene expression}
}
\description{
Pseudobulk pairwise
}
