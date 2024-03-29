% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reactome.R
\name{fetch_reactome_species}
\alias{fetch_reactome_species}
\title{List of available Reactome species}
\usage{
fetch_reactome_species(on_error = c("stop", "warn"))
}
\arguments{
\item{on_error}{A character vector specifying the error handling method. It
can take values `"stop"` or `"warn"`. The default is `"stop"`. `"stop"`
will halt the function execution and throw an error, while `"warn"` will
issue a warning and return `NULL`.}
}
\value{
A tibble with species names used by Reactome.
}
\description{
List of available Reactome species
}
\examples{
re <- fetch_reactome_species(on_error = "warn")
}
