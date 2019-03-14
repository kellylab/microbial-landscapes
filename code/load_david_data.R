library(data.table)

david <- fread(paste0(data.dir, "david/david.otus"),
               col.names = c("sample", "otu", "count")
               )
# parse subject and timepoints
david[, c("subject", "day") := tstrsplit(sample, "_")]
david[, day := as.numeric(day)]
david.samples <- unique(david[, .(sample, subject, day)])
# event metadata
setkey(david.samples, subject)
events <- mapply(function(start, end, subject, event) {
  x <- data.table(day = seq(start, end, by = 1), subject = subject,
                  event = event)
  x
}, start = c(0, 71, 80, 104, 123, 0, 151, 160 ),
end = c(70, 122, 85, 113, david.samples["A", max(day)],
          150, 159, david.samples["B", max(day)]),
subject = c("A", "A", "A", "A", "A",
            "B", "B", "B"),
event = c("US (pre)", "travel", "diarrhea 1", "diarrhea 2", "US (post)",
          "pre-Salmonella", "Salmonella", "post-Salmonella"),
SIMPLIFY = FALSE) %>% rbindlist(use.names = TRUE)
# collapse event labels per day
events <- events[, .(event = paste(event, collapse = " + ")),
                 by = .(subject, day)]
david.samples <- merge(david.samples, events, by = c("subject", "day"))
is.healthy <- function(event) {
  if (event == "Salmonella" | grepl("diarrhea", event)) {
    FALSE
  } else {
    TRUE
  }
}
david.samples[, healthy := sapply(event, is.healthy)]