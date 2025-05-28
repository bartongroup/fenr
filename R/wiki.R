#' WikiPathways server base URL
#'
#' @return A string with URL. This can be changed by options(WIKI_BASE_URL = "different/url").
#' @noRd
get_wiki_url <- function() {
  getOption("WIKI_BASE_URL", "https://webservice.wikipathways.org")
}

#' List of available WikiPathways species
#'
#' @param on_error A character string indicating the error handling strategy:
#'   either "stop" to halt execution, "warn" to issue a warning and return
#'   `NULL` or "ignore" to return `NULL` without warnings. Defaults to "stop".
#'
#' @return A character vector with species names used by WikiPathways.
#' @export
#' @examples
#' spec <- fetch_wiki_species(on_error = "warn")
fetch_wiki_species <- function(on_error = c("stop", "warn", "ignore")) {
  on_error <- match.arg(on_error)

  resp <- http_request(get_wiki_url(), "listOrganisms", parameters = list(format = "json"))
  if(resp$is_error)
    return(catch_error("WikiPathways", resp, on_error))

  js <- httr2::resp_body_json(resp$response)
  tibble::tibble(designation = unlist(js))
}

#' Download pathway data from WikiPathways
#'
#' @param species Species name recognised by WikiPathways (see \code{fetch_wiki_species()})
#' @param on_error A character string indicating the error handling strategy:
#'   either "stop" to halt execution, "warn" to issue a warning and return
#'   `NULL` or "ignore" to return `NULL` without warnings. Defaults to "stop".
#'
#' @return A tibble with columns \code{gene_id} and \code{term_id}
#' @noRd
fetch_wiki_pathways <- function(species, on_error) {
  id <- name <- NULL

  resp <- http_request(get_wiki_url(), "listPathways",
                       parameters = list(format = "json", organism = species))
  if(resp$is_error)
    return(catch_error("WikiPathways", resp, on_error))

  js <- httr2::resp_body_json(resp$response)
  if(is(js[[1]], "character") && js[[1]] == "error")
    stop(stringr::str_glue("Cannot retrieve pathways from WikiPathways for species {species}."))

  js$pathways |>
    purrr::map(tibble::as_tibble) |>
    purrr::list_rbind() |>
    dplyr::select(term_id = id, term_name = name)
}


#' Helper function for `parse_wiki_gpml`
#'
#' @param s Input string
#' @param key Key to extract
#'
#' @return Value of the key
#' @noRd
extract_key <- function(s, key) {
  stringr::str_extract(s, paste0("(?<=", key, '=\")(.+?)(?=")'))
}

#' Parse GPML string provided by WikiPathways
#'
#' @param gpml GMPL string returned by WikiPathways query.
#' @param keys Keys to be extracted.
#'
#' @return A tibble with selected keys and values
#' @noRd
parse_wiki_gpml <- function(gpml, keys = c("TextLabel", "Type", "Database", "ID")) {
  gpml |>
    stringr::str_replace_all("\\n", "") |>
    stringr::str_extract_all("<DataNode.+?</DataNode>") |>
    purrr::pluck(1) |>
    purrr::map_dfr(function(s) {
      purrr::map_chr(keys, ~extract_key(s, .x)) |>
        purrr::set_names(keys)
    })
}


#' Download WikiPathways gene-pathway mapping using API
#'
#' @param pathways A character vector with pathway names.
#' @param databases Names of databases to use
#' @param types Names of types to use
#' @param on_error A character string indicating the error handling strategy:
#'   either "stop" to halt execution, "warn" to issue a warning and return
#'   `NULL` or "ignore" to return `NULL` without warnings. Defaults to "stop".
#'
#' @return A tibble with \code{term_id} and \code{gene_symbol}
#' @noRd
fetch_wiki_pathway_genes_api <- function(pathways, databases = NULL, types = NULL, on_error) {
  term_id <- TextLabel <-ID <- Type <- Database <- database <- type <- NULL

  raise_error <- FALSE
  pb <- progress::progress_bar$new(total = length(pathways))
  res <- purrr::map(pathways, function(pathway) {
    pb$tick()
    resp <- http_request(get_wiki_url(), "getPathway",
                         parameters = list(format = "json", pwId = pathway))
    if(resp$is_error) {
      catch_error("WikiPathways", resp, on_error)
      raise_error <<- TRUE
      return()
    }
    js <- httr2::resp_body_json(resp$response)
    parse_wiki_gpml(js$pathway$gpml) |>
      tibble::add_column(term_id = pathway)
  }) |>
    purrr::list_rbind()
  if(raise_error)
    return(NULL)

  res <- res |>
    dplyr::select(term_id, text_label = TextLabel, id = ID, type = Type, database = Database) |>
    dplyr::distinct()
  if(!is.null(databases)) {
    res <- res |>
      dplyr::filter(database %in% databases)
  }
  if(!is.null(types)) {
    res <- res |>
      dplyr::filter(type %in% types)
  }

  res
}


#' Get functional term data from WikiPathways
#'
#' Download term information (pathway ID and name) and gene-pathway mapping
#' (gene symbol and pathway ID) from WikiPathways.
#'
#' @details WikiPathways contain mapping between pathways and a variety of
#'   entities from various databases. Typically a gene symbol is returned in
#'   column \code{text_label} and some sort of ID in column \code{id}, but this
#'   depends on the species and databases used. For gene/protein enrichment,
#'   these should be filtered to contain gene symbols only. This can be done by
#'   selecting a desired databases and types. The default values for parameters
#'   \code{databases} and \code{types} attempt to select information from generic
#'   databases, but there are organism-specific databases not included in the
#'   selection. We suggest to run this function with \code{databases = NULL,
#'   types = NULL} to see what types and databases are available before making
#'   selection.
#'
#' @param species WikiPathways species designation, for example "Homo sapiens"
#'   for human. Full list of available species can be found using
#'   \code{fetch_wiki_species()}.
#' @param databases A character vector with database names to pre-filter mapping
#'   data. See details. Full result will be returned if NULL.
#' @param types A character vector with types of entities to pre-filter mapping
#'   data. See details. Full result will be returned if NULL.
#' @param on_error A character string indicating the error handling strategy:
#'   either "stop" to halt execution, "warn" to issue a warning and return
#'   `NULL` or "ignore" to return `NULL` without warnings. Defaults to "stop".
#'
#' @return A list with \code{terms} and \code{mapping} tibbles.
#' @export
#' @importFrom assertthat assert_that
#' @examples
#' wiki_data <- fetch_wiki("Bacillus subtilis", on_error = "warn")
fetch_wiki <- function(
    species,
    databases = c("Ensembl", "Entrez Gene", "HGNC", "HGNC Accession number", "Uniprot-TrEMBL"),
    types = c("GeneProduct", "Protein", "Rna", "RNA"),
    on_error = c("stop", "warn", "ignore")
  ) {
  on_error <- match.arg(on_error)
  assert_that(!missing(species), msg = "Argument 'species' is missing.")
  assert_species(species, "fetch_wiki_species", on_error)

  terms <- fetch_wiki_pathways(species, on_error)
  if(is.null(terms))
    return(NULL)
  mapping <- fetch_wiki_pathway_genes_api(terms$term_id, databases, types, on_error)
  list(
    terms = terms,
    mapping = mapping
  )
}
