% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reactome.R
\name{fetch_reactome_api_genes}
\alias{fetch_reactome_api_genes}
\title{Download term - gene symbol mapping from Reactome}
\usage{
fetch_reactome_api_genes(pathways)
}
\arguments{
\item{pathways}{A character vector with Reactome patway IDs to get
corresponding genes from.}
}
\value{
A tibble with columns\code{term_id}, \code{accession_number} and
  \code{gene_symbol}.
}
\description{
Download term - gene symbol mapping from Reactome
}
\details{
This function interrogates Reactome API to get term-gene mapping for
  all pathways. This is considerable slower than
  \code{fetch_reactome_ensembl_genes}. Warning, occasionally, for some
  pathways, Reactome does not return gene symbol - only UniProt accession
  number is available.
}
