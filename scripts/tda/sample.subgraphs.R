#' Mapper subgraph per data point
#'
#' @param df \code{data.frame} with 2 columns: \code{samples, vertices}
#' @param graf the top level graph with all the vertices
#'
#' @return list of subgraphs each corresponding to a sample, named by sample
#' @export
#'
#' @examples
sample.subgraphs <- function(df, graf) {
  samples <- df$samples
  vertices <- df$vertices
  uniq.samps <- unique(samples)
  sverts <- lapply(uniq.samps, function(s, vertices, samples) {
    vertices[samples == s]
  }, vertices = vertices, samples = samples)
  names(sverts) <- uniq.samps
  sgrafs <- lapply(sverts, function(vs, graf) {
    igraph::induced_subgraph(graf, vs)
  }, graf = graf)
  sgrafs
}
