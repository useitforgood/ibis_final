require("Biostrings")
require("magrittr")
require("dplyr")
require("readr")
# require("purrr")
# require("tidyr")
# RAW_DATA_PATH <- "~/ibis_clean/final/data/raw_data/"
GENOME_PATH <- "~/ibis/raw_data/genome/"
TRAIN_PATH <- "~/ibis/clean/data/train_unified/"
DATA_PATHS <- list(
  "lb" = list(
    train = "raw_data/IBIS.train_data.Leaderboard.v2/",
    test = "raw_data/IBIS.test_data.Leaderboard.v2/",
    is_genome = c("GABPA", "PRDM5", "SP140", "ZNF362", "ZNF407")
  ),
  "fin" = list(
    train = "raw_data/IBIS.train_data.Final.v1/train/",
    test = "raw_data/IBIS.test_data.Final.v1/",
    is_genome = c("CAMTA1", "LEUTX", "MYF6", "PRDM13", "SALL3", "USF3", "ZBED2", "ZBED5", "ZNF20", "ZNF251", "ZNF367", "ZNF395", "ZNF493", "ZNF518B", "ZNF648")
  )
)
DATA_TYPES <- names(DATA_PATHS)