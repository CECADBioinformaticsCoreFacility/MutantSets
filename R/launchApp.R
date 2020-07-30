#' \code{MutantSets} package
#' 
#' See README on \href{https://github.com/RichardJActon/MutantSets/}{GitHub}
#' 
#' @docType package
#' @name MutantSets
NULL

utils::globalVariables(c(
	".", "AC", "AF", "ALT", "Amino_Acid_Length", "CHROM", "DP", "EFF",
	"ERRORS", "Genotype_Number", "ID", "Indiv", "NUMALT", "POS", "QA", "QR",
	"QUAL", "REF", "TYPE", "WARNINGS", "end", "feature", "gt_GT", "keep",
	"pc", "phase", "pos", "quantile", "start", "strand"
))


#' launches the shinyAppDemo app
#'
#' @export launchApp
#'
#' @return shiny application object
#'
#'
#' @import shiny
#' @import dplyr
#' @import shinydashboard
#'
launchApp <- function() {
	shinyApp(ui = ui, server = server)
}
