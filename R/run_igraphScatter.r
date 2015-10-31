#' Run Scatter shiny app
#' 
#' Call igraphScatter. The app was developed for Tung et al. (2016)
#' to faciliate exploratory data analysis.
#' 
#' @export
#' @examples
#' run_igraphScatter()
#' 
run_igraphScatter <- function() {
  appDir <- system.file("shiny-apps", "igraphScatter", package = "Humanzee")

  if (appDir == "") {
    stop("Could not find App directory. Try re-installing `Humanzee`.", call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal")
}