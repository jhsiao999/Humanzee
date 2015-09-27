#' @export

run_igraphScatter <- function() {
  appDir <- system.file("shiny-apps", "igraphScatter", package = "Humanzee")

  if (appDir == "") {
    stop("Could not find App directory. Try re-installing `Humanzee`.", call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal")
}