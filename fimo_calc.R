require("Biostrings")
require("magrittr")
require("dplyr")
require("readr")
require("purrr")
require("parallel")
require("glue")
require("XML")
set.seed(53)
GENOME_PATH <- "data/genome/"
TRAIN_PATH <- "data/train_unified/"
DATA_PATHS <- list(
  "lb" = list(
    train = "data/raw_data/IBIS.train_data.Leaderboard.v2/",
    test = "data/raw_data/IBIS.test_data.Leaderboard.v3/",
    subm = "data/raw_data/Leaderboard_aaa_submissions.v3/",
    is_genome = c("GABPA", "PRDM5", "SP140", "ZNF362", "ZNF407"),
    is_artificial = c("LEF1", "NACC2", "RORB", "TIGD3", "NFKB1")
  ),
  "fin" = list(
    train = "data/raw_data/IBIS.train_data.Final.v1/train/",
    test = "data/raw_data/IBIS.test_data.Final.v1/",
    subm = "data/raw_data/Final_aaa_submissions.v1/",
    is_genome = c("CAMTA1", "LEUTX", "MYF6", "PRDM13", "SALL3", "USF3", "ZBED2", "ZBED5", "ZNF20", "ZNF251", "ZNF367", "ZNF395", "ZNF493", "ZNF518B", "ZNF648"),
    is_artificial = c("GCM1", "MKX", "MSANTD1", "MYPOP", "SP140L", "TPRX1", "ZFTA", "CREB3L3", "FIZ1", "ZBTB47", "ZNF286B", "ZNF500", "ZNF721", "ZNF780B", "ZNF831")
  )
)
DATA_TYPES <- names(DATA_PATHS)
sessionInfo()
sapply(list.files("functions/", full.names = TRUE), source)

base_command <- "/home/ubuntu/bin/fimo"
for (type in DATA_TYPES) {
  best_vars <- readRDS(glue::glue("meme_comparison/preds/best_vars_ext_selected_{type}12.RDS"))
  experiments <- list.files(DATA_PATHS[[type]]$subm, pattern = "random", full.names = TRUE)
  exp_name <- glue::glue("final_{type}")
  tryCatch({
    for (experiment in experiments[5]) {
      exp_ <- gsub("_aaa_random.tsv", "", gsub(paste0(DATA_PATHS[[type]]$subm, "/"), "", experiment))
      tfs_ <- readRDS(glue::glue("data/submissions/tfs_{exp_}_{type}.RDS"))
      for (tf_ in tfs_) {
        fimo_dir <- glue::glue("tmp/fimo_out_{type}/{exp_name}/{tf_}/{exp_}")
        if (!dir.exists(fimo_dir)) dir.create(fimo_dir, recursive = TRUE)
        best_vars_sub <- best_vars[best_vars$tf == tf_, ]
        for (i in 1:nrow(best_vars_sub)) {
          fimo_out_dir <- paste0("meme_comparison/data/fimo_out_", type, "/", exp_name, "/", exp_, "/")
          fimo_out_fname <- paste0(fimo_out_dir, basename(dirname(dirname(best_vars_sub$X2[i]))), "_", best_vars_sub$X4[i], "_", tf_, ".RDS")
          if (file.exists(fimo_out_fname)) next
          run_fimo <- glue::glue(
            "{base_command} --motif {best_vars_sub$X3[i]} --max-stored-scores 10000000 --oc {fimo_dir} {best_vars_sub$X2[i]} {DATA_PATHS[[type]]$test}/{exp_}_participants.fasta")
          system(run_fimo)
          for (file in list.files(fimo_dir, full.names = TRUE)[!grepl("tsv$", list.files(fimo_dir))]) {
            file.remove(file)
          }
          f_name <- glue::glue("{fimo_dir}/fimo.tsv")
          if (!dir.exists(fimo_out_dir)) dir.create(fimo_out_dir, recursive = TRUE)
          if (file.exists(f_name)) {
            dt <- read_delim(f_name)
            dt %<>% filter(!is.na(sequence_name))
            dt$prob_1 <- exp(dt$score) / (1 + exp(dt$score))
            dt$prob_2 <- 1 - dt$`q-value`
            dt %<>% filter(`p-value` < 0.001)
            dt$tf <- tf_
            if (nrow(dt) > 0) {
              if (nrow(dt) > 0) {
                saveRDS(dt, fimo_out_fname)
              }
            }
          }
        }
      }
    }
  },
  error = function(cond)
    print(cond))
}