#Fetch and verify cached GEO series matrices
#Downloads each accession to DIR_RAW if absent; otherwise verifies the cached
#file's MD5 against the recorded checksum and aborts on mismatch.

library(GEOquery)

CHECKSUM_FILE <- here::here("results", "data_checksums.txt")

fetch_geo_cached <- function(acc) {
  pattern  <- paste0("^", acc, "_series_matrix\\.txt\\.gz$")
  existing <- list.files(DIR_RAW, pattern = pattern, full.names = TRUE)

  if (length(existing) == 0) {
    message("Downloading ", acc, " to ", DIR_RAW)
    getGEO(acc, GSEMatrix = TRUE, destdir = DIR_RAW)
    existing <- list.files(DIR_RAW, pattern = pattern, full.names = TRUE)
    if (length(existing) == 0) {
      stop("Expected series matrix for ", acc, " not found in ", DIR_RAW, " after download.")
    }
    md5    <- unname(tools::md5sum(existing[1]))
    record <- data.frame(file = basename(existing[1]), md5 = md5,
                          downloaded = format(Sys.Date()), stringsAsFactors = FALSE)
    write.table(record, CHECKSUM_FILE,
                append     = file.exists(CHECKSUM_FILE),
                col.names  = !file.exists(CHECKSUM_FILE),
                row.names  = FALSE, sep = "\t", quote = FALSE)
    message(basename(existing[1]), " downloaded and checksummed.")
  } else {
    if (!file.exists(CHECKSUM_FILE)) {
      stop("Cached file ", basename(existing[1]), " exists but ", CHECKSUM_FILE,
           " is missing. Delete data/raw/ and re-run to redownload with a fresh checksum record.")
    }
    recorded <- read.table(CHECKSUM_FILE, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    row      <- recorded[recorded$file == basename(existing[1]), ]
    if (nrow(row) == 0) {
      stop("No checksum recorded for ", basename(existing[1]), " in ", CHECKSUM_FILE,
           ". Delete data/raw/ and re-run to redownload with a fresh checksum record.")
    }
    actual <- unname(tools::md5sum(existing[1]))
    if (actual != row$md5[1]) {
      stop("MD5 mismatch for ", basename(existing[1]), ": expected ", row$md5[1],
           ", got ", actual, ". The cached file is corrupted or was replaced - delete it from ",
           DIR_RAW, " and re-run to redownload.")
    }
    message(basename(existing[1]), " checksum verified.")
  }
  invisible(existing[1])
}

for (.acc in c("GSE56500", "GSE68605")) fetch_geo_cached(.acc)
rm(.acc)
