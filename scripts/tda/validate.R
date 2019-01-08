library(data.table)
library(tidyverse)

subsample <- function(dt, r = 0.9) {
  nselect <- round(nrow(dt) * r)
  sample_n(dt, nselect)
}