library(data.table)

data.dir <- "../../data/"
gordon <- fread(paste0(data.dir, "cholera/gordon.otus"),
                col.names = c("sample", "otu", "count"))
# extract sample info
gordon[, c("subject", "state", "id") := tstrsplit(sample, "_")]
# mark 'end' time as last diarrhea pt + 1hr
gordon[, hour := {
  di.max <- max(as.numeric(id[state == "diarrhea" & id != "end"]))
  t <- sapply(id, function(id, di.max) {
    if (id == "end") {
      di.max + 1
    } else if (grepl("d", id)) {
      (24 * as.numeric(sub("d", "", id))) + di.max + 1
    } else {
      as.numeric(id)
    }
  }, di.max = di.max)
}, by = subject]

gordon[, sample := paste(subject, state, hour, sep = "_")]
