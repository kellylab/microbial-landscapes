library(data.table)

david <- fread(paste0("../../data/david/david.otus"),
               col.names = c("sample", "otu", "count")
               )
# parse subject and timepoints
david[, c("subject", "day") := tstrsplit(sample, "_")]