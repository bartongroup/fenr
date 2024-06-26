# Gene Annotation File column names. GAF is used by GO and Reactome.
GAF_COLUMNS <- c(
  "db",
  "db_id",
  "symbol",
  "qualifier",
  "go_term",
  "db_ref",
  "evidence",
  "from",
  "aspect",
  "db_object_name",
  "db_object_synonym",
  "db_object_type",
  "taxon",
  "date",
  "assigned_by",
  "annotation_extension",
  "form_id"
)

GAF_TYPES <- rep("c", length(GAF_COLUMNS)) |>
  stringr::str_c(collapse = "")


#' Check if URL with a path is responding OK
#'
#' Stops with error message if the path is not accessible or server not
#' responding.
#'
#' @param url_path Full URL with a path, e.g. `https://reactome.org/download/current/ReactomePathways.txt`.
#' @param on_error A character string indicating the error handling strategy:
#'   either "stop" to halt execution, "warn" to issue a warning and return
#'   `NULL` or "ignore" to return `NULL` without warnings. Defaults to "stop".
#' @param timeout Timeout limit in seconds.
#'
#' @importFrom assertthat assert_that is.string
#' @return TRUE if assertion passed
#' @noRd
assert_url_path <- function(url_path, on_error = c("stop", "warn", "ignore"), timeout = 30) {
  on_error <- match.arg(on_error)
  assert_that(is.string(url_path))

  resp <- http_request(url_path, "", timeout = timeout)

  if(resp$is_error) {
    catch_error(stringr::str_glue("{url_path}"), resp, on_error)
    return(FALSE)
  }
  return(TRUE)
}



#' Check if species are valid
#'
#' @param species A string, species designation for a given database
#' @param fetch_fun A string, name of the function to retrieve available
#'   species, must return a tibble with a column \code{designation}.
#' @param on_error A character string indicating the error handling strategy:
#'   either "stop" to halt execution, "warn" to issue a warning and return
#'   `NULL` or "ignore" to return `NULL` without warnings. Defaults to "stop".
#'
#' @importFrom assertthat assert_that is.string
#' @return A tibble with valid species - a response from \code{fetch_fun}
#' @noRd
assert_species <- function(species, fetch_fun, on_error) {
  assert_that(is.string(species))
  assert_that(is.string(fetch_fun))
  f <- match.fun(fetch_fun)
  valid_species <- f(on_error = on_error)
  if(is.null(valid_species))
    return(NULL)
  assert_that(species %in% valid_species$designation,
    msg = stringr::str_glue("Invalid species {species}. Use {fetch_fun}() to find all available species.")
  )
  valid_species
}


#' Return a field for a given species
#'
#' Used to match species designation with a species ID.
#'
#' @param species Species designation
#' @param fetch_fun  A string with name of a function to retrieve available
#'   species, must return a tibble with a column \code{designation}.
#' @param col_name Column name in the tibble returned by \code{fetch_fun} to
#'   extract, e.g. \code{tax_id}
#' @param on_error A character string indicating the error handling strategy:
#'   either "stop" to halt execution, "warn" to issue a warning and return
#'   `NULL` or "ignore" to return `NULL` without warnings. Defaults to "stop".
#'
#' @return A value extracted from column \code{col_name} at row where
#'   \code{designation} = \code{species}.
#' @noRd
match_species <- function(species, fetch_fun, col_name, on_error) {
  designation <- NULL

  sp <- assert_species(species, fetch_fun, on_error)
  if(is.null(sp))
    return(NULL)
  sp |>
    dplyr::filter(designation == species) |>
    dplyr::pull(get(col_name)) |>
    unique()
}


#' Get features annotated with a given term.
#'
#' @param term_data An object class \code{fenr_terms}, created by
#'   \code{prepare_for_enrichment}.
#' @param term_id A string with a functional term ID.
#'
#' @return A character vector containing feature IDs annotated with the term ID.
#'
#' @importFrom methods is
#' @importFrom assertthat assert_that is.string
#' @export
#' @examples
#' \dontrun{
#' go_data <- fetch_go(species = "sgd")
#' go_terms <- prepare_for_enrichment(go_data$terms, go_data$mapping, feature = "gene_symbol")
#' feats <- get_term_features(go_terms, "GO:0000001")
#' }
get_term_features <- function(term_data, term_id) {
  # Check term_data class
  assert_that(is(term_data, "fenr_terms"))
  assert_that(is.string(term_id))

  term_data$term2feature[[term_id]]
}


#' Get terms annotating a given feature.
#'
#' @param term_data An object class \code{fenr_terms}, created by
#'   \code{prepare_for_enrichment}.
#' @param feature_id A string with a feature ID
#'
#' @return A character vector containing functional term IDs annotating given feature.
#' @importFrom assertthat assert_that is.string
#' @importFrom methods is
#' @export
#' @examples
#' \dontrun{
#' go_data <- fetch_go(species = "sgd")
#' go_terms <- prepare_for_enrichment(go_data$terms, go_data$mapping, feature = "gene_symbol")
#' trms <- get_feature_terms(go_terms, "GEM1")
#' }
get_feature_terms <- function(term_data, feature_id) {
  # Check term_data class
  assert_that(is(term_data, "fenr_terms"))
  assert_that(is.string(feature_id))

  term_data$feature2term[[feature_id]]
}



#' Test correctness of data structure returned by fetch_*
#'
#' @param re The data structure returned by fetch_*
#' @noRd
test_fetched_structure <- function(re) {
  testthat::expect_is(re, "list")
  testthat::expect_length(re, 2)
  testthat::expect_named(re)
  testthat::expect_equal(names(re), c("terms", "mapping"))
}


#' Test if returned terms tibble is consistent with expected selection
#'
#' @param returned A tibble, the terms data returned by  fetch_*.
#' @param expected A tibble, a selection of expected terms data.
#' @noRd
test_terms <- function(returned, expected) {
  # Binding variables from non-standard evaluation locally
  term_id <- term_name <- NULL

  testthat::expect_is(returned, "tbl")
  testthat::expect_contains(names(returned), c("term_id", "term_name"))

  returned <- returned |> dplyr::select(term_id, term_name)

  merged <- expected |>
    dplyr::left_join(returned, by = c("term_id", "term_name")) |>
    tidyr::drop_na()
  testthat::expect_equal(nrow(expected), nrow(merged))
}

#' Test if returned mapping is consistent with expected selection
#'
#' @param returned A tibble, the mapping returned by  fetch_*.
#' @param expected A tibble, a selection of expected mapping.
#' @param feature_id A character string, the name of the column containing the feature
#'   ID.
#' @noRd
test_mapping <- function(returned, expected, feature_id) {
  # Binding variables from non-standard evaluation locally
  term_id <- NULL

  testthat::expect_is(returned, "tbl")
  testthat::expect_contains(names(returned), c("term_id", feature_id))

  returned <- returned |> dplyr::select(term_id, tidyselect::all_of(feature_id))

  merged <- expected |>
    dplyr::left_join(returned, by = c("term_id", feature_id)) |>
    dplyr::distinct() |>
    tidyr::drop_na()
  testthat::expect_equal(nrow(expected), nrow(merged))
}

#' Catch and Handle Errors Based on Specified Action
#'
#' This function handles errors by either stopping execution or issuing a
#' warning, depending on the specified action. It takes a message and an action
#' choice (`"stop"` or `"warn"`) as inputs.
#'
#' @param server Name of the server.
#' @param resp Response from \code{http_request()}. A list with "status" and
#'   "desc" strings.
#' @param on_error A character string indicating the error handling strategy:
#'   either "stop" to halt execution, "warn" to issue a warning and return
#'   `NULL` or "ignore" to return `NULL` without warnings. Defaults to "stop".
#'
#' @return No return value if `on_error` is "stop". If `on_error` is "warn" or
#'   "ignore", the function returns `NULL` after issuing a warning.
#' @noRd
catch_error <- function(server, resp, on_error = c("stop", "warn", "ignore")) {
  on_error <- match.arg(on_error)

  st <- stringr::str_glue("Cannot access {server}. {resp$status}: {resp$description}.")

  error_response(st, on_error)
}

#' Error Response Handler
#'
#' This function handles errors by either stopping execution or issuing a
#' warning based on the specified action. It is designed to provide a
#' standardized way to handle error messages within functions.
#'
#' @param msg A character string specifying the error message to be displayed.
#' @param on_error A character string indicating the error handling strategy:
#'   either "stop" to halt execution, "warn" to issue a warning and return
#'   `NULL` or "ignore" to return `NULL` without warnings. Defaults to "stop".
#'
#' @return No return value if `on_error` is "stop". If `on_error` is "warn" or
#'   "ignore", the function returns `NULL` after issuing a warning.
#' @noRd
error_response <- function(msg, on_error = c("stop", "warn", "ignore")) {
  on_error <- match.arg(on_error)

  if(on_error == "stop") {
    stop(msg)
  } else if (on_error == "warn") {
    warning(msg, "\nNULL returned.", call. = FALSE)
    return(NULL)
  } else {
    return(NULL)
  }
}

#' Perform a HTTP request using httr2
#'
#' This function sends an HTTP request to a specified API endpoint and returns
#' various details about the response. It is designed to work with the httr2
#' package.
#'
#' @param base_url A string specifying the base URL of the API.
#' @param path A string specifying the path of the specific API endpoint.
#' @param parameters An optional named list of query parameters to be included
#'   in the request.
#' @param timeout Timeout limit in seconds.
#'
#' @return A list containing the following elements:
#' \itemize{
#'   \item{response}{The complete response object returned by the API.}
#'   \item{is_error}{Logical value indicating if the response was an error.}
#'   \item{status}{The HTTP status code of the response.}
#'   \item{description}{The description of the HTTP status.}
#' }
#'
#' @details The function constructs a request using the base URL and the path.
#'   If provided, query parameters are appended to the request. The function
#'   then performs the request and checks for errors. It returns a list
#'   containing the response, error status, HTTP status code, and the
#'   description of the status. In case of a time-out, when the server does not
#'   return any code, a list with NULL response, and fixed status and
#'   description are returned.
#' @noRd
http_request <- function(base_url, path, parameters = NULL, timeout = 30) {
  req <- httr2::request(base_url) |>
    httr2::req_url_path_append(path) |>
    httr2::req_timeout(timeout)

  if(!is.null(parameters)) {
    req <- req |>
      httr2::req_url_query(!!!parameters)
  }

  resp <- tryCatch(
    {req |>
        httr2::req_error(is_error = ~FALSE) |>
        httr2::req_perform()
    },
    error = function(err) NULL
  )

  if(is.null(resp)) {
    ans <- list(
      response = NULL,
      is_error = TRUE,
      status = "unknown",
      description = "Unknown error (probably timeout)"
    )
  } else {
    ans <- list(
      response = resp,
      is_error = httr2::resp_is_error(resp),
      status = httr2::resp_status(resp),
      description = httr2::resp_status_desc(resp)
    )
  }

  ans
}


#' Test a function for correct response when server is non-responsive.
#'
#' @details When the server is not responsive (using
#'   httr2::with_mocked_responses), the function should throw an error, if
#'   `on_error = "stop"` and throw a warning and return NULL, it `on_error =
#'   "warn"`.
#'
#'
#' @param fun A function to test
#' @param ... Parameters to pass to the function
#'
#' @return Test result
#' @noRd
test_unresponsive_server <- function(fun, ...) {
  testthat::expect_error(fun(..., on_error = "stop"))
  testthat::expect_warning({res <- fun(..., on_error = "warn")})
  testthat::expect_null(res)
}


#' A mock function used with httr2::with_mocked_responses() for testing an
#' unresponsive server behaviour.
#'
#' @param req A httr2 request
#'
#' @return Mocked response with error 500.
#' @noRd
mocked_500 <- function(req) {
  httr2::response(status_code = 500)
}
