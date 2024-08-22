get_train <- function(pt, exp = NULL) {
  if (!is.null(exp)) {
    return(readRDS(glue::glue("{TRAIN_PATH}{exp}_{pt}.RDS")))
  } else {
    return(readRDS(glue::glue("{TRAIN_PATH}_{pt}.RDS")))
  }
}

process_meme_exp <- function(data, tfs, exp, base_command = "docker exec memesuite meme", flags = NULL, nmotifs_ = 4, maxw = 30, pt = "lb") {
  for (tf_ in tfs) {
    dir_ <- glue::glue("meme_comparison/data/tmp_memedir_{pt}/{exp}/{tf_}")
    print(dir_)
    print(glue::glue("{dir_}.fasta"))
    if (file.exists(glue::glue("{dir_}/meme.xml")) & file.size(glue::glue("{dir_}/meme.xml")) != 0) next
    if (file.exists(glue::glue("{dir_}/meme.txt")) != 0) next
    if (!dir.exists(dir_)) {
      dir.create(dir_, recursive = TRUE)
    }
    writeLines(paste0(paste0(">", data$n[data$tf == tf_], "\n"), data$seq[data$tf == tf_], collapse = "\n"),
               glue::glue("{dir_}.fasta"))
    n_threads <- detectCores() - 1
    if (sum(data$tf == tf_) > 1500000) n_threads <- ceiling(n_threads/2)
    if (sum(data$tf == tf_) > 2500000) n_threads <- ceiling(n_threads/4)
    n_threads <- 1
    command <- paste(base_command,  glue::glue("{dir_}.fasta"), "-oc", dir_, "-nmotifs", nmotifs_, "-dna -maxw ", maxw, " -p", n_threads)
    if (!is.null(flags)) command <- paste(command, flags)
    system(command)
    file.remove(glue::glue("{dir_}.fasta"))
  }
}

process_row.artificial <- function(row, data, pt = "lb") {
  t <- Sys.time()
  print(row$n)
  exp <- paste0(row, collapse = "_")
  dir_ <- glue::glue("meme_comparison/data/tmp_memedir_{pt}/{exp}/")
  # if (length(list.files(dir_)) > 0) {
  #   print(paste(dir_, "exists"))
  #   return(list(param = row,
  #               time = Sys.time() - t))
  # }
  if (!row$SMS) data$SMS <- NULL
  if (!row$PBM) data$PBM <- NULL
  if (!row$HTS) data$HTS <- NULL
  if (row$PBM) {
    if (!is.na(row$PBM_subset)) {
      subset <- row$PBM_subset
      data$PBM %<>% filter(grepl(subset, f_name))
    }
    if (!is.na(row$PBM_filter)) {
      if (row$PBM_filter == "max") {
        data$PBM %<>% 
          mutate(method = gsub("_.*", "", f_name)) %>% 
          group_by(id_probe, pbm_sequence) %>% 
          filter(mean_signal_intensity == max(mean_signal_intensity)) %>%
          ungroup %>% 
          select(tf, seq = pbm_sequence) %>% 
          mutate(method = "PBM")
      }
      if (row$PBM_filter == "2sd") {
        data$PBM %<>% 
          mutate(method = gsub("_.*", "", f_name)) %>% 
          group_by(method) %>% 
          mutate(mean_sig = mean(mean_signal_intensity), sd_sig = sd(mean_signal_intensity), th = mean_sig + 2*sd_sig) %>% 
          ungroup %>% 
          filter(mean_signal_intensity > th) %>% 
          select(tf, seq = pbm_sequence) %>% 
          mutate(method = "PBM")
      }
      if (row$PBM_filter %in% c("5", "10")) {
        data$PBM %<>% 
          mutate(method = gsub("_.*", "", f_name)) %>% 
          group_by(method) %>% 
          mutate(mean_sig = mean(mean_signal_intensity), sd_sig = sd(mean_signal_intensity), th = mean_sig + 2*sd_sig) %>% 
          ungroup %>% 
          filter(mean_signal_intensity > th) %>% 
          select(tf, seq = pbm_sequence) %>% 
          mutate(method = "PBM")
        ind <- as.numeric(row$PBM_filter)
        data$PBM$seq <- substr(data$PBM$seq, nchar(data$PBM$seq)/2 - ind, nchar(data$PBM$seq)/2 + ind)
      }
    } else {
      data$PBM %<>% 
        mutate(method = gsub("_.*", "", f_name)) %>% 
        select(tf, seq = pbm_sequence) %>% 
        mutate(method = "PBM")
    }
  }
  if (row$HTS) {
    if (!is.na(row$HTS_subset)) {
      if (row$HTS_subset == "C3") {
        data$HTS %<>% filter(grepl("C3", f_name) | grepl("C4", f_name))
      }
      if (row$HTS_subset == "C4") {
        data$HTS %<>% filter(grepl("C4", f_name))
      }
      if (row$PBM_filter %in% c("5", "10")) {
        ind <- as.numeric(row$PBM_filter)
        data$HTS$seq <- substr(data$HTS$seq, nchar(data$HTS$seq)/2 - ind, nchar(data$HTS$seq)/2 + ind)
      }
    }
  }
  flags <- paste("-seed", row$seed)
  if (row$meme_rev) flags <- paste(flags, "-revcomp")
  data %<>% lapply(select, c("tf", "seq", "method"))
  data %<>% Reduce(rbind, .)
  data %<>% filter(!grepl("N", seq))
  data$n <- 1:nrow(data)
  tfs <- unique(data$tf)
  process_meme_exp(data, tfs, exp, base_command = "/home/ubuntu/bin/meme", nmotifs_ = 6, flags = flags, pt = pt)
  print(exp)
  print("done")
  list(
    param = row,
    time = Sys.time() - t
  )
}

process_row.genome <- function(row, data, pt = "lb") {
  t <- Sys.time()
  print(row$n)
  exp <- paste0(row, collapse = "_")
  dir_ <- glue::glue("meme_comparison/data/tmp_memedir_{pt}/{exp}/")
  # if (length(list.files(dir_)) > 0) {
  #   print(paste(dir_, "exists"))
  #   return(list(param = row,
  #               time = Sys.time() - t))
  # }
  tfs <- unique(data$tf)
  if (!is.na(row$length)) {
    ind <- row$length
    data$seq <- substr(data$seq, nchar(data$seq)/2 - ind, nchar(data$seq)/2 + ind)
  }
  if (!is.na(row$pileup)) {
    data %<>% group_by(method, tf) %>% filter(pileup > quantile(pileup, row$pileup/100)) %>% ungroup
  }
  if (!row$GHTS) data %<>% filter(method != "GHTS")
  if (!row$CHS) data %<>% filter(method != "CHS")
  if (length(tfs) != length(unique(data$tf))) return("not all tfs")
  flags <- paste("-seed", row$seed)
  if (row$meme_rev) flags <- paste(flags, "-revcomp")
  data %<>% filter(!grepl("N", seq))
  data$n <- 1:nrow(data)
  process_meme_exp(data, tfs, exp, base_command = "/home/ubuntu/bin/meme", nmotifs_ = 6, flags = flags, maxw = row$maxw, pt = pt)
  print(exp)
  print("done")
  list(
    param = row,
    time = Sys.time() - t
  )
}

process_ame <- function(row, dir_, base_command = "/home/ubuntu/bin/ame", suff = "_cut", type = "fin") {
  tf <- row$tf[1]
  num <- row$bucket[1]
  ame_file <- glue::glue("meme_comparison/data/ame{suff}_{type}/{tf}_{num}.RDS")
  if (!dir.exists(glue::glue("meme_comparison/data/ame{suff}_{type}"))) dir.create(glue::glue("meme_comparison/data/ame{suff}_{type}"), recursive = TRUE)
  ame_file <- glue::glue("meme_comparison/data/ame{suff}_{type}/{tf}_{num}.RDS")
  pos <- glue::glue("{dir_}{tf}{suff}_positive_{type}.fasta")
  if (!file.exists(pos)) pos <- glue::glue("{dir_}{tf}{suff}_ext_positive_{type}.fasta")
  neg <- glue::glue("{dir_}{tf}{suff}_negative_{type}.fasta")
  if (!file.exists(neg)) neg <- glue::glue("{dir_}{tf}{suff}_ext_negative_{type}.fasta")
  pwms <- glue::glue("{row$params}/{tf}/filtered_meme.xml")
  existing_pwm_paths <- pwms[file.exists(pwms)]
  if (length(existing_pwm_paths) > 0) {
    if (file.exists(ame_file)) {
      return("exists")
    }
    pwms <- paste(existing_pwm_paths, collapse = " ")
    command <- paste(base_command, "--oc", glue::glue("{getwd()}/ame{suff}_{type}_{tf}/"), " --noseq --control ", neg, pos, pwms)
    system(command)
    ame <- read_delim(glue::glue("{getwd()}/ame{suff}_{type}_{tf}/ame.tsv"), delim = "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE)
    ame %<>% filter(X1 != "rank", !is.na(X10))
    ame %<>% mutate_at(c("X15", "X17"), as.numeric)
    saveRDS(distinct(ame, X11, X12), glue::glue("meme_comparison/data/ame{suff}_{type}/counts_{tf}_{num}.RDS"))
    ame$X17 <- 100 - ame$X17
    ame$score <- (ame$X15 + ame$X17)/2
    ame %<>% select(X2, X3, X4, X15, X17, score)
    # ame %<>% arrange(desc(score)) %>% distinct(X2, .keep_all = TRUE)
    # ame %<>% head(10)
    saveRDS(ame, ame_file)
    return("success")
  } else {
    return("failure")
  }
}