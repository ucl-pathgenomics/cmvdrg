#' Runs the CMV Resistance Phenotyping Shiny Application
#' 
#' @import shiny
#' 
#' @export

runShinyCMV = function() {
  appDir <- system.file("shiny-app", package = "cmvdrg")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `cmvdrg`.", call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal")
}