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

dir_ <- "meme_comparison/references/"
names_g <- readRDS("data/params/genome_param_grid.RDS")
names_g <- unlist(lapply(split(names_g, seq(nrow(names_g))), function(row) paste(row, collapse = "_")))
names_a <- readRDS("data/params/artificial_param_grid.RDS")
names_a <- unlist(lapply(split(names_a, seq(nrow(names_a))), function(row) paste(row, collapse = "_")))

args <- commandArgs(trailingOnly=TRUE)
ind <- if (is.na(as.numeric(args[1]))) args[1] else as.numeric(args[1])

for (type in DATA_TYPES) {
  all_tfs <- c(DATA_PATHS[[type]]$is_genome, DATA_PATHS[[type]]$is_artificial)
  if (is.numeric(ind)) {
    tfs <- ind
  } else {
    tfs <- args[1]
  }
  if (length(tfs) == 0) next
  comparison_grid <- rbind(
    expand.grid(tf = DATA_PATHS[[type]]$is_genome, 
                params = names_g, 
                stringsAsFactors = FALSE),
    expand.grid(tf = DATA_PATHS[[type]]$is_artificial, 
                params = names_a, 
                stringsAsFactors = FALSE)
  )
  comparison_grid %<>% filter(tf %in% tfs)
  comparison_grid$params <- file.path(glue::glue("meme_comparison/data/tmp_memedir_{type}/"), comparison_grid$params)
  for (tf_ in unique(comparison_grid$tf)) {
    print(tf_)
    suff_ <- "_sep_cut"
    if (tf_ %in% DATA_PATHS[[type]]$is_genome) suff_ <- ""
    print(suff_)
    comparison_grid_pt <- comparison_grid %>% filter(tf == tf_)
    comparison_grid_pt$bucket <- ceiling(seq_along(comparison_grid_pt$tf) / 20)
    row_list <- split(comparison_grid_pt, comparison_grid_pt$bucket)
    results <- lapply(row_list, function(row) process_ame(row, dir_, suff = suff_, type = type))
  }
}