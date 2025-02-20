#' URL of GO obo top dir
#'
#' @return A string with URL.
#' @noRd
get_go_obo_dir <- function() {
  getOption("GO_OBO_DIR", "http://purl.obolibrary.org/obo")
}

#' URL of GO gene association file
#'
#' @return A string with URL.
#' @noRd
get_go_obo_file <- function() {
  getOption("GO_OBO_FILE", "http://purl.obolibrary.org/obo/go.obo")
}

#' URL of GO species webpage
#'
#' @return A string with URL.
#' @noRd
get_go_species_url <- function() {
  getOption("GO_SPECIES_URL", "http://current.geneontology.org/products/pages/downloads.html")
}

#' URL of GO annotation server
#'
#' @return A string with URL.
#' @noRd
get_go_annotation_url <- function() {
  getOption("GO_ANNOTATION_URL", "http://current.geneontology.org/annotations")
}

#' URL of Ensembl's biomart
#'
#' @return A string with URL.
#' @noRd
get_biomart_url <- function() {
  getOption("ENSEMBL_BIOMART", "http://www.ensembl.org")
}

# See https://www.ensembl.org/info/data/biomart/biomart_restful.html for more details.
get_biomart_xml <- function(dataset) {
stringr::str_glue("
<?xml version='1.0' encoding='UTF-8'?>
<!DOCTYPE Query>
<Query virtualSchemaName = 'default' uniqueRows = '1' count='0' datasetConfigVersion='0.6' header='1' formatter='TSV' requestid='biomaRt'>
  <Dataset name = '{dataset}'>
    <Attribute name = 'ensembl_gene_id'/>
    <Attribute name = 'external_gene_name'/>
    <Attribute name = 'go_id'/>
    <Attribute name = 'go_linkage_type'/>
  </Dataset>
</Query>") |>
    stringr::str_replace_all("\n", "") |>
    stringr::str_replace_all("\\>\\s+\\<", "><")
}

#' Parse OBO file and return a tibble with key and value
#'
#' @param obo Obo file content as a character vector
#'
#' @return A tibble with term_id, key and value
#' @noRd
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
  ends <- blanks[seq_len(length(starts))]

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

#' Extract GO term IDs and names from parsed OBO data
#'
#' @param parsed OBO data parsed by \code{parse_obo_file}
#'
#' @return A tibble with term_id and term_name
#' @noRd
extract_obo_terms <- function(parsed) {
  # Binding variables from non-standard evaluation locally
  key <- term_id <- value <- name <- term_name <- namespace <- NULL

  terms <- parsed |>
    dplyr::filter(key %in% c("name", "namespace")) |>
    tidyr::pivot_wider(id_cols = term_id, names_from = key, values_from = value) |>
    dplyr::rename(term_name = name)

  alt_terms <- parsed |>
    dplyr::filter(key == "alt_id") |>
    dplyr::left_join(terms, by = "term_id") |>
    dplyr::select(term_id = value, term_name, namespace)

  dplyr::bind_rows(
    terms,
    alt_terms
  )
}


#' Download GO term descriptions
#'
#' @param use_cache Logical, if TRUE, the remote file will be cached locally.
#' @param on_error A character string indicating the error handling strategy:
#'   either "stop" to halt execution, "warn" to issue a warning and return
#'   `NULL` or "ignore" to return `NULL` without warnings. Defaults to "stop".
#'
#' @return A tibble with term_id and term_name.
#' @noRd
fetch_go_terms <- function(use_cache, on_error) {
  if(!assert_url_path(get_go_obo_dir(), on_error))
    return(NULL)

  obo_file <- get_go_obo_file()
  lpath <- cached_url_path("obo", obo_file, use_cache)
  readr::read_lines(lpath) |>
    parse_obo_file() |>
    extract_obo_terms()
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
#' @param on_error A character string indicating the error handling strategy:
#'   either "stop" to halt execution, "warn" to issue a warning and return
#'   `NULL` or "ignore" to return `NULL` without warnings. Defaults to "stop".
#'
#' @return A tibble with columns \code{species} and \code{designation}.
#' @export
#' @examples
#' go_species <- fetch_go_species(on_error = "warn")
fetch_go_species <- function(on_error = c("stop", "warn", "ignore")) {
  on_error <- match.arg(on_error)
  # Binding variables from non-standard evaluation locally
  species <- designation <- `Species/Database` <- File <- NULL

  url <- get_go_species_url()
  resp <- http_request(url, "")
  if(resp$is_error)
    return(catch_error("GO species website", resp, on_error))

  u <- resp$response |>
    httr2::resp_body_html() |>
    rvest::html_table()

  u[[1]] |>
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
#' @param use_cache Logical, if TRUE, the remote file will be cached locally.
#' @param on_error A character string indicating the error handling strategy:
#'   either "stop" to halt execution, "warn" to issue a warning and return
#'   `NULL` or "ignore" to return `NULL` without warnings. Defaults to "stop".
#'
#' @return A tibble with columns \code{gene_symbol}, \code{uniprot_id}, 
#'   \code{term_id} and \code{evidence}.
#' @noRd
fetch_go_genes_go <- function(species, use_cache, on_error) {
  # Binding variables from non-standard evaluation locally
  gene_id <- db_object_synonym <- symbol <- NULL
  db_id <- go_term <- evidence <- NULL

  url <- get_go_annotation_url()
  if(!assert_url_path(url, on_error))
    return(NULL)
  gaf_file <- stringr::str_glue("{url}/{species}.gaf.gz")

  lpath <- cached_url_path(stringr::str_glue("go_gaf_{species}"), gaf_file, use_cache)
  readr::read_tsv(lpath, comment = "!", quote = "", col_names = GAF_COLUMNS,
                  col_types = GAF_TYPES) |>
    dplyr::mutate(gene_id = stringr::str_remove(db_object_synonym, "\\|.*$")) |>
    dplyr::select(gene_symbol = symbol, gene_id, db_id, term_id = go_term, evidence) |>
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
#'   Object Synonym) is returned as \code{gene_id}. It is up to the user to
#'   select the appropriate database.
#'
#' @param species Species designation. Examples are \file{goa_human} for human,
#'   \file{mgi} for mouse or \file{sgd} for yeast. Full list of available
#'   species can be obtained using \code{fetch_go_species} - column
#'   \code{designation}.
#' @param use_cache Logical, if TRUE, the remote file will be cached locally.
#' @param on_error A character string indicating the error handling strategy:
#'   either "stop" to halt execution, "warn" to issue a warning and return
#'   `NULL` or "ignore" to return `NULL` without warnings. Defaults to "stop".
#'
#' @return A list with \code{terms} and \code{mapping} tibbles.
#' @importFrom assertthat assert_that
#' @noRd
fetch_go_from_go <- function(species, use_cache, on_error) {
  assert_that(!missing(species), msg = "Argument 'species' is missing.")
  assert_species(species, "fetch_go_species", on_error)

  mapping <- fetch_go_genes_go(species = species, use_cache = use_cache, on_error = on_error)
  if(is.null(mapping))
    return(NULL)

  terms <- fetch_go_terms(use_cache = use_cache, on_error = on_error)
  if(is.null(terms))
    return(NULL)

  list(
    terms = terms,
    mapping = mapping
  )
}



#' Download GO term gene mapping from Ensembl
#'
#' @param dataset Dataset you want to use. To see the different datasets
#'   available within a biomaRt you can e.g. do: mart <-
#'   biomaRt::useEnsembl(biomart = "ensembl"), followed by
#'   biomaRt::listDatasets(mart).
#' @param use_cache Logical, if TRUE, the remote data will be cached locally.
#' @param on_error A character string indicating the error handling strategy:
#'   either "stop" to halt execution, "warn" to issue a warning and return
#'   `NULL` or "ignore" to return `NULL` without warnings. Defaults to "stop".
#'
#' @return A tibble with columns \code{gene_id}, \code{gene_symbol},
#'   \code{term_id} and \code{evidence}.
#' @noRd
fetch_go_genes_bm <- function(dataset, use_cache, on_error) {
  xml <- get_biomart_xml(dataset) |>
    stringr::str_replace_all("\\s", "%20")
  biomart_path <- paste0(get_biomart_url(), "/biomart/martservice")
  if(!assert_url_path(biomart_path, on_error))
    return(NULL)

  req <- paste0(biomart_path, "?query=", xml)

  # Problems with cache, bfcneedsupdate returns error for this query
  # lpath <- cached_url_path(stringr::str_glue("biomart_{dataset}"), resp, use_cache)
  res <- readr::read_tsv(req, show_col_types = FALSE)
  if(ncol(res) == 4) {
    res |> rlang::set_names(c("gene_id", "gene_symbol", "term_id", "evidence"))
  } else {
    error_response("Problem with Biomart", on_error)
  }
}



#' Get functional term data from Ensembl
#'
#' Download term information (GO term ID and name) and gene-term mapping
#' (gene ID, symbol and GO term ID) from Ensembl.
#'
#' @param dataset Dataset you want to use. To see the different datasets
#'   available within a biomaRt you can e.g. do: mart <-
#'   biomaRt::useEnsembl(biomart = "ensembl"), followed by
#'   biomaRt::listDatasets(mart).
#' @param use_cache Logical, if TRUE, the remote file will be cached locally.
#' @param on_error A character string indicating the error handling strategy:
#'   either "stop" to halt execution, "warn" to issue a warning and return
#'   `NULL` or "ignore" to return `NULL` without warnings. Defaults to "stop".
#'
#' @importFrom assertthat assert_that is.string
#' @return A list with \code{terms} and \code{mapping} tibbles.
#' @noRd
fetch_go_from_bm <- function(dataset, use_cache, on_error) {
  assert_that(!missing(dataset), msg = "Argument 'dataset' is missing.")
  assert_that(is.string(dataset))

  mapping <- fetch_go_genes_bm(dataset, use_cache = use_cache, on_error = on_error)
  if(is.null(mapping))
    return(error_response("Could not retrieve mapping from Ensembl", on_error))
  terms <- fetch_go_terms(use_cache = use_cache, on_error)

  list(
    terms = terms,
    mapping = mapping
  )
}


#' Get Gene Ontology (GO) data
#'
#'
#' This function downloads term information (GO term ID and name) and gene-term
#' mapping (gene ID, symbol, and GO term ID) from either the Ensembl database
#' (using BioMart) or the Gene Ontology database (using GAF files), depending on
#' the provided argument.
#'
#' @details If \code{species} is provided, mapping from a Gene Ontology GAF file
#'   will be downloaded. GAF files contain more generic information than gene
#'   symbols. In this function, the third column of the GAF file (DB Object
#'   Symbol) is returned as \code{gene_symbol}, but, depending on the
#'   \code{species} argument it can contain other entities, e.g. RNA or protein
#'   complex names. Similarly, the eleventh column of the GAF file (DB Object
#'   Synonym) is returned as \code{gene_id}. It is up to the user to select
#'   the appropriate database.
#'
#'   Alternatively, if \code{dataset} is provided, mapping will be downloaded
#'   from Ensembl database. It will gene symbol and Ensembl gene ID.
#'
#' @param species (Optional) Species designation. Examples are \code{goa_human}
#'   for human, \code{mgi} for mouse, or \code{sgd} for yeast. Full list of
#'   available species can be obtained using \code{fetch_go_species} - column
#'   \code{designation}. This argument is used when fetching data from the Gene
#'   Ontology database.
#' @param dataset (Optional) A string representing the dataset passed to
#'   Ensebml's Biomart, e.g. 'scerevisiae_gene_ensembl'. To see the different
#'   datasets available within a biomaRt you can e.g. do: mart <-
#'   biomaRt::useEnsembl(biomart = "ensembl"), followed by
#'   biomaRt::listDatasets(mart).
#' @param use_cache Logical, if TRUE, the remote data will be cached locally.
#' @param on_error A character string indicating the error handling strategy:
#'   either "stop" to halt execution, "warn" to issue a warning and return
#'   `NULL` or "ignore" to return `NULL` without warnings. Defaults to "stop".
#'
#' @return A list with \code{terms} and \code{mapping} tibbles.
#' @export
#' @importFrom assertthat assert_that
#' @examples
#' # Fetch GO data from Ensembl
#' go_data_ensembl <- fetch_go(dataset = "scerevisiae_gene_ensembl", on_error = "warn")
#' # Fetch GO data from Gene Ontology
#' go_data_go <- fetch_go(species = "sgd", on_error = "warn")
fetch_go <- function(species = NULL, dataset = NULL, use_cache = TRUE,
                     on_error = c("stop", "warn", "ignore")) {
  on_error <- match.arg(on_error)

  assert_that(!(is.null(species) & is.null(dataset)),
              msg = "One of the arguments 'species' or 'dataset' must be supplied.")
  assert_that(is.null(species) | is.null(dataset),
              msg = "Only one of the arguments 'species' or 'dataset' must be supplied.")

  if (!is.null(species)) {
    fetch_go_from_go(species, use_cache = use_cache, on_error = on_error)
  } else {
    fetch_go_from_bm(dataset, use_cache = use_cache, on_error = on_error)
  }
}
