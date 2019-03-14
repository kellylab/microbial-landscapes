library(tidyverse)
library(data.table)

read_prochloro <- function(path) {
    dt <- fread(path)
    # replace all 'nd' with NA
    set(dt, i = NULL, names(dt), lapply(dt, function(col) {
      col[col == "nd"] <- NA
      col
    }))
    # TODO: don't throw away standard deviation data
    dt <- melt(dt,
               id.vars = c("cruiseid",
                           "lat",
                           "lon",
                           "date",
                           "year",
                           "month",
                           "day",
                           "yrday_gmt",
                           "depth",
                           "sigma_0",
                           "temp",
                           "sal"),
               measure.vars = patterns("^Abundance"),
               variable.name = "ecotype",
               value.name = "abundance")
    # discard rows with no time stamps
    dt <- dt[!is.na(date)]
    # format ordinal and calendar months
    dt[, c("month", "year") := lapply(list(month, year), as.numeric)]
    dt[, cal.month := month]
    dt[, month := month + 12 * (year - min(year))]
    # format abundances and ecotypes
    dt[, abundance := as.numeric(abundance)]
    # set NA abundance to minimum non-NA value
    min.abund <- min(dt$abundance, na.rm = TRUE)
    dt[is.na(abundance), abundance := min.abund]
    # format ecotype names
    dt[, ecotype := gsub("Abundance", "", ecotype)]
    dt[, ecotype := gsub("_", "", ecotype)]
    dt
}

prochlorococcus <- c(bats = "bats_orig.txt", hot = "hot_orig.txt") %>%
  sapply(function(fn, root) paste0(root, fn), root = data.dir) %>%
  lapply(read_prochloro) %>%
  rbindlist(idcol = "site")
prochlorococcus[, sample := paste(site, cruiseid, depth, sep = "-")]
prochlorococcus[, c("temp", "sal") := lapply(list(temp, sal), as.numeric)]
prochlorococcus.samples <- prochlorococcus[, -c("ecotype", "abundance")] %>%
  unique
