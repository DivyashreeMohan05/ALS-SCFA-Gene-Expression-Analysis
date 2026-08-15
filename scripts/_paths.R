#Project paths - sourced by all scripts

library(here)

DIR_RAW     <- here("data", "raw")
DIR_INTERIM <- here("data", "interim")
DIR_TABLES  <- here("results", "tables")
DIR_FIGURES <- here("results", "figures")

for (d in c(DIR_RAW, DIR_INTERIM, DIR_TABLES, DIR_FIGURES)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
