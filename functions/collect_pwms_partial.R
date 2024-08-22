collect_pwms_partial <- function(df, output) {
  template <- c()
  mots <- c()
  for (i in 1:nrow(df)) {
    row <- df[i, ]
    motif <- getNodeSet(xmlParse(row$X2), glue::glue("//motif[@alt = '{row$X4}']"))
    values <- sprintf("%.5f", as.numeric(xmlValue(getNodeSet(motif[[1]], "./probabilities/alphabet_matrix/alphabet_array/value"))))
    motif <- t(matrix(values, 4))
    tf_ <- row$tf
    template <- c(template, glue::glue(">{tf_} {tf_}_motif{row$n}_{substr(digest::digest(paste0(basename(dirname(dirname(row$X2))), row$X3)), 0, 8)}"), apply(motif, 1, paste, collapse = " "), "\n")
    mots <- c(mots, paste0(c("A", "C", "G", "T")[apply(motif, 1, which.max)], collapse = ""))
  }
  writeLines(template, output)
  mots
}