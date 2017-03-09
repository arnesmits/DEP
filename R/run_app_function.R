#' proteomeR shiny apps
#'
#' \code{run_app} Launches an interactive shiny app for interactive analysis of proteomics data.
#'
#' @param app The name of the app. Currently supporting "LFQ" or "TMT".
#' @examples
#' run_app("LFQ")
#' run_app("TMT")
#' @export
run_app <- function(app) {
  # locate all the shiny apps that exist
  valid_apps <- list.files(system.file("shiny_apps", package = "proteomeR"))

  valid_apps_msg <-
    paste0(
      "Valid apps are: '",
      paste(valid_apps, collapse = "', '"),
      "'")

  # if an invalid app is given, throw an error
  if (missing(app) || !nzchar(app) ||
      !app %in% valid_apps) {
    stop(
      'Please run `run_app()` with a valid app as an argument.n',
      valid_apps_msg,
      call. = FALSE)
  }

  # find and launch the app
  appDir <- system.file("shiny_apps", app, package = "proteomeR")
  suppressWarnings(shiny::runApp(appDir, display.mode = "normal"))
}
