#' Parse OBO file and return a tibble with key and value
#'
#' @param obo Obo file content as a character vector
#'
#' @return A tibble with term_id, key and value
parse_obo_file <- function(obo) {
  # Find index of start and end line of each term

  # Start lines
  starts <- stringr::str_which(obo, "\\[Term\\]")

  # Empty lines at end of each term
  blanks <- stringr::str_which(obo, "^$")
  blanks <- blanks[blanks > starts[1]]

  # No space at the end
  if(length(blanks) < length(starts))
    blanks <- c(blanks, length(obo) + 1)

  # End lines: ignore empty lines beyond terms
  ends <- blanks[1:length(starts)]

  # Parse each term
  purrr::map2_dfr(starts, ends, function(i1, i2) {
    obo_term <- obo[(i1 + 1):(i2 - 1)]
    trm <- obo_term |>
      stringr::str_split(":\\s", 2, simplify = TRUE)
    colnames(trm) <- c("key", "value")
    # assuming term_id is in the first line, if not, we are screwed
    tid <- trm[1, 2]
    cbind(trm, term_id = tid) |>
      as.data.frame()
  }) |>
    tibble::as_tibble()
}


#' Download GO term descriptions
#'
#' @param obo_file A URL or local file containing GO ontology, in OBO format.
#'
#' @return A tibble with term_id and term_name.
fetch_go_terms <- function(obo_file = "http://purl.obolibrary.org/obo/go.obo") {
  # Binding variables from non-standard evaluation locally
  key <- term_id <- value <- term_name <- NULL

  parsed <- readr::read_lines(obo_file) |>
    parse_obo_file()

  terms <- parsed |>
    dplyr::filter(key == "name") |>
    dplyr::select(term_id, term_name = value)

  alt_terms <- parsed |>
    dplyr::filter(key == "alt_id") |>
    dplyr::left_join(terms, by = "term_id") |>
    dplyr::select(term_id = value, term_name)

  dplyr::bind_rows(
    terms,
    alt_terms
  )
}



#' Find all species available from geneontology.org
#'
#' This function attempts to scrape HTML web page containing a table of
#' available species and corresponding file names. If the structure of the page
#' changes one day and the function stops working, go to
#' \url{http://current.geneontology.org/products/pages/downloads.html} and check
#' file names. The species designation used in this package is the GAF file name
#' without extension (e.g. for a file \file{goa_chicken.gaf} the designation is
#' \file{goa_chicken}).
#'
#' @param url URL of the Gene Ontology web page with downloads.
#'
#' @return A tibble with columns \code{species} and \code{designation}.
#' @import XML
#' @export
#'
#' @examples
#' go_species <- fetch_go_species()
fetch_go_species <- function(url = "http://current.geneontology.org/products/pages/downloads.html") {
  # Binding variables from non-standard evaluation locally
  species <- designation <- `Species/Database` <- File <- NULL

  assert_url_path(url)
  u <- httr::GET(url) |>
    httr::content("text", encoding = "UTF-8") |>
    XML::readHTMLTable(as.data.frame = TRUE)
  u[[1]] |>
    tibble::as_tibble() |>
    dplyr::mutate(
      species = `Species/Database` |>
        stringr::str_replace_all("\\n", "-") |>
        stringr::str_replace_all("\\s\\s+", " ") |>
        stringr::str_replace_all("(\\S)-", "\\1"),
      designation = File |>
        stringr::str_remove("\\..*$")
    ) |>
    dplyr::select(species, designation)
}


#' Download GO term gene mapping from geneontology.org
#'
#' @param species Species designation. Base file name for species file under
#'   \url{http://current.geneontology.org/annotations}. Examples are
#'   \file{goa_human} for human, \file{mgi} for mouse or \file{sgd} for yeast.
#'
#' @import assertthat
#' @return A tibble with columns \code{gene_symbol}, \code{uniprot_id} and \code{term_id}.
fetch_go_genes_go <- function(species) {
  # Binding variables from non-standard evaluation locally
  gene_synonym <- db_object_synonym <- gene_symbol <- symbol <- NULL
  uniprot_id <- db_id <- term_id <- go_term <- NULL

  gaf_file <- stringr::str_glue("http://current.geneontology.org/annotations/{species}.gaf.gz")
  assert_url_path(gaf_file)

  readr::read_tsv(gaf_file, comment = "!", quote = "", col_names = GAF_COLUMNS, col_types = GAF_TYPES) |>
    dplyr::mutate(gene_synonym = stringr::str_remove(db_object_synonym, "\\|.*$")) |>
    dplyr::select(gene_symbol = symbol, gene_synonym, db_id, term_id = go_term) |>
    dplyr::distinct()
}


#' Get functional term data from gene ontology
#'
#' Download term information (GO term ID and name) and gene-term mapping (gene
#' symbol and GO term ID) from gene ontology.
#'
#' @details This function relies on Gene Ontology's GAF files containing more
#'   generic information than gene symbols. Here, the third column of the GAF
#'   file (DB Object Symbol) is returned as \code{gene_symbol}, but, depending
#'   on the \code{species} argument it can contain other entities, e.g. RNA or
#'   protein complex names. Similarly, the eleventh column of the GAF file (DB
#'   Object Synonym) is returned as \code{gene_synonym}. It is up to the user to
#'   select the appropriate database.
#'
#' @param species Species designation. Examples are \file{goa_human} for human,
#'   \file{mgi} for mouse or \file{sgd} for yeast. Full list of available
#'   species can be obtained using \code{fetch_go_species} - column
#'   \code{designation}.
#'
#' @return A list with \code{terms} and \code{mapping} tibbles.
#' @export
#' @import assertthat
#'
#' @examples
#' go_data <- fetch_go_from_go("sgd")
fetch_go_from_go <- function(species) {
  assert_that(!missing(species), msg = "Argument 'species' is missing.")
  assert_species(species, "fetch_go_species")

  mapping <- fetch_go_genes_go(species)
  terms <- fetch_go_terms()

  list(
    terms = terms,
    mapping = mapping
  )
}



#' Download GO term gene mapping from Ensembl
#'
#' @param mart Object class \code{Mart} representing connection to BioMart
#'   database, created with, e.g., \code{useEnsembl}.
#'
#' @return A tibble with columns \code{ensembl_gene_id}, \code{gene_symbol} and
#'   \code{term_id}.
fetch_go_genes_bm <- function(mart) {
  # Binding variables from non-standard evaluation locally
  gene_symbol <- external_gene_name <- term_id <- go_id <- NULL

  biomaRt::getBM(
    attributes = c("ensembl_gene_id", "external_gene_name", "go_id"),
    mart = mart
  ) |>
    dplyr::rename(
      gene_symbol = external_gene_name,
      term_id = go_id
    ) |>
    dplyr::filter(term_id != "") |>
    tibble::as_tibble()
}



#' Get functional term data from Ensembl
#'
#' Download term information (GO term ID and name) and gene-term mapping
#' (gene ID, symbol and GO term ID) from Ensembl.
#'
#' @param mart Object class \code{Mart} representing connection to BioMart
#'   database, created with, e.g., \code{useEnsembl}.
#'
#' @return A list with \code{terms} and \code{mapping} tibbles.
#' @export
#' @import assertthat
#' @importFrom methods is
#'
#' @examples
#' \dontrun{
#' mart <- biomaRt::useEnsembl(biomart = "ensembl", dataset = "scerevisiae_gene_ensembl")
#' go_terms <- fetch_go_from_bm(mart)
#' }
fetch_go_from_bm <- function(mart) {
  assert_that(!missing(mart), msg = "Argument 'mart' is missing.")
  assert_that(is(mart, "Mart"))

  terms <- fetch_go_terms()
  mapping <- fetch_go_genes_bm(mart)

  list(
    terms = terms,
    mapping = mapping
  )
}

