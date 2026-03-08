#' @title Run Shiny app
#' @description Runs the EventProjection Shiny app.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
runShinyApp <- function() {
  app_pkgs <- c(
    "shinyMatrix",
    "shinyFeedback",
    "shinyjs",
    "shinybusy",
    "readxl",
    "writexl",
    "DT",
    "prompter"
  )

  missing_pkgs <- app_pkgs[!vapply(
    app_pkgs, requireNamespace, logical(1), quietly = TRUE
  )]

  if (length(missing_pkgs) > 0) {
    stop(
      "To run the Shiny app, please install these optional packages first: ",
      paste(missing_pkgs, collapse = ", "),
      "\nExample:\ninstall.packages(c(",
      paste(sprintf('"%s"', missing_pkgs), collapse = ", "),
      "))",
      call. = FALSE
    )
  }

  shiny::shinyAppDir(system.file("shinyApp", package = "EventProjection"))
}

#' @title Run Shiny app
#' @description Backward-compatible wrapper for \code{runShinyApp()}.
#'
#' @export
runShinyApp_eventPred <- function() {
  runShinyApp()
}
