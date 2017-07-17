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

# HOT
hot <- fread(paste0(data.dir, "/hot_orig.txt"))
# some data don't have timestamps. these are not useful
hot <- hot[!(year == "nd")]
# some samples didn't detect any ecotypes, leave these out
# manual sample exclusion is EXTREMELY UGLY
hot <- hot[Abundance_e9312 != "nd" &
             Abundance_eMED4 != "nd" &
             Abundance_NATL != "nd" &
             Abundancee_SS120 != "nd" &
             Abundancee_9313 != "nd"]
hot[, ":=" (year = as.numeric(year), month = as.numeric(month))]
hot[, cal.month := month]
hot[, month := month + 12 * (year - min(year))]
hot <- hot %>%
  melt(id.vars = c("cruiseid", "yrday_gmt", "year", "month", "depth", 
                   "cal.month", "day"),
       measure.vars = patterns("Abundance"),
       variable.name = "ecotype",
       value.name = "abundance")
hot[, abundance := as.numeric(abundance)]
min.abund <- min(hot$abundance, na.rm = TRUE)
hot[is.na(abundance), abundance := min.abund]
hot[, ecotype := gsub("Abundance", "", ecotype)]
hot[, ecotype := gsub("_", "", ecotype)]
# extract sample annotations
hot.info <- unique(hot[, .(cruiseid, yrday_gmt, year, month, cal.month, day)])