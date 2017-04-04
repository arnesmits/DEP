#' proteomeR shiny apps
#'
#' \code{run_app} Launches an interactive shiny app for interactive analysis of proteomics data.
#'
#' @param app The name of the app. Currently supporting "LFQ" or "TMT".
#' @examples
#' \dontrun{
#' # Run the app
#' run_app("LFQ")
#'
#' run_app("TMT")
#'
#' }
#' @export
run_app <- function(app) {
  # Locate all the shiny apps that exist
  valid_apps <- list.files(system.file("shiny_apps", package = "proteomeR"))

  valid_apps_msg <-
    paste0(
      "Valid apps are: '",
      paste(valid_apps, collapse = "', '"),
      "'")

  # Show error if an unvalid app-name is given
  if (missing(app) || !nzchar(app) ||
      !app %in% valid_apps) {
    stop(
      'Please run `run_app()` with a valid app as argument\n',
      valid_apps_msg,
      call. = FALSE)
  }

  # Launch the app
  appDir <- system.file("shiny_apps", app, package = "proteomeR")
  suppressWarnings(shiny::runApp(appDir, display.mode = "normal"))
}
