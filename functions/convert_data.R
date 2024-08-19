merge_peak_files <- function(path, genome_path = GENOME_PATH) {
  all_files <- list.files(path, recursive = TRUE, full.names = TRUE)
  data <- lapply(all_files, function(f_name)
    read_delim(
      f_name,
      delim = "\t",
      escape_double = FALSE,
      trim_ws = TRUE
    ) %>%
      mutate(tf = basename(dirname(f_name)), method = basename(dirname(dirname(f_name))), f_name = basename(f_name)))
  data %<>% Reduce(rbind, .)
  if ("name" %in% names(data)) {
    data %<>% filter(name != "name") 
    data %<>% mutate_at(c(2:8), as.numeric)
    data %<>% group_by(`#CHROM`) %>% 
      mutate(
        seq = get_seqs(abs_summit, unique(`#CHROM`), genome_path)
      ) %>% 
      ungroup
  }
  data
}

get_seqs <- function(position, chr_name, genome_path, WINDOW = 60) {
  print(glue::glue("reading dna_seq {chr_name}"))
  dna_seq <- readDNAStringSet(glue::glue("{genome_path}/{chr_name}.fa.gz"), format = "fasta")
  print(glue::glue("processing dna_seq {chr_name}"))
  map_chr(position, ~ as.character(subseq(dna_seq, start = .x - WINDOW, end = .x + WINDOW)[[1]]))
}

merge_seqs <- function(path, all_files = NULL) {
  if(is.null(all_files)) all_files <- list.files(path, recursive = TRUE, full.names = TRUE, pattern = "*.gz")
  data <- lapply(all_files, function(f_name) {
    format_ <- if (tools::file_ext(f_name) == "fasta") "fasta" else "fastq"
    dna_seq <- as.character(readDNAStringSet(f_name, format = format_))
    frame <- data.frame(
      seq = dna_seq,
      tf = basename(dirname(f_name)),
      method = basename(dirname(dirname(f_name))),
      f_name = basename(f_name)
    )
  })
  data %<>% Reduce(rbind, .)
  data
}