% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interactive_example.R
\name{fetch_terms_for_example}
\alias{fetch_terms_for_example}
\title{Create term data for interactive example}
\usage{
fetch_terms_for_example(de)
}
\arguments{
\item{de}{Differential expression results, use \code{yeast_de} data attached
to this package.}
}
\value{
A list of objects containing functional terms for GO and Reactome.
}
\description{
Create term data for interactive example
}
\examples{
\dontrun{
data(yeast_de)
term_data <- fetch_terms_for_example(yeast_de)
}
}
