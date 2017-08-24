#' Parse OTU's ConsensusLineage label to hierarchical taxonomy
#'
#' @param lineage A vector of ConsensusLineage labels.
#'
#' @return A list, where each element of the list is a vector of classifications at that taxonomic level, and each element of the vector corresponds to an element in the input vector.
#'
#' @export
#'
#' @examples
lineage_2_taxonomy <- function(lineage) {
  library(data.table)
  llist <- tstrsplit(lineage, ";")
  lapply(llist, function(v) {
    sapply(v, function(w) {
      out <- sub("[[:alpha:]]__", "", w)
      if (out == "") {
        out <- NA_character_
      }
      out
    })
  })
}