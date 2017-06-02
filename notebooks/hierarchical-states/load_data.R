library(data.table)
data.dir <- "../../data/"

# BATS
bats <- fread(paste0(data.dir, "/bats_orig.txt"))
bats[, cal.month := month]
bats[, month := month + 12 * (year - min(year))]
bats <- bats %>%
  melt(id.vars = c("cruiseid", "yrday_gmt", "year", "month", "depth", "cal.month",
                   "day"),
       measure.vars = patterns("Abundance"),
       variable.name = "ecotype",
       value.name = "abundance")
bats[, abundance := as.numeric(abundance)]
min.abund <- min(bats$abundance, na.rm = TRUE)
bats[is.na(abundance), abundance := min.abund]
bats[, ecotype := gsub("Abundance", "", ecotype)]
bats[, ecotype := gsub("_", "", ecotype)]
# extract sample annotations
bats.info <- unique(bats[, .(cruiseid, yrday_gmt, year, month, cal.month, day)])