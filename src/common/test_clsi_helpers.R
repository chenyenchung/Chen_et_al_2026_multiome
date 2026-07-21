#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]))
source(file.path(dirname(script_path), "clsi_helpers.R"))

stopifnot(
  isTRUE(all.equal(
    simpsons_pair_concentration(c("A", "A", "B", "B")),
    1 / 3
  )),
  isTRUE(all.equal(
    simpsons_pair_concentration(c("A", "A", NA_character_)),
    1
  )),
  isTRUE(all.equal(
    simpsons_pair_concentration(c("A", "B", NA_character_)),
    0
  )),
  is.na(simpsons_pair_concentration(c("A", NA_character_)))
)

labels <- c("A", "A", "B", "B")
global_si <- simpsons_pair_concentration(labels)
all_pairs <- combn(seq_along(labels), 2, simplify = FALSE)
pair_clsi <- vapply(all_pairs, function(idx) {
  local_si <- simpsons_pair_concentration(labels[idx])
  (local_si - global_si) / (1 - global_si)
}, numeric(1))
stopifnot(isTRUE(all.equal(mean(pair_clsi), 0)))

nn_idx <- rbind(
  c(2, 3),
  c(1, 4),
  c(1, 4),
  c(2, 3)
)
scores <- compute_clsi(labels, nn_idx, label_name = "test labels")
stopifnot(length(scores) == length(labels), all(is.finite(scores)))

labels_with_na <- c("A", "A", "B", NA_character_)
sparse_nn <- rbind(
  c(2, 4),
  c(1, 4),
  c(1, 2),
  c(1, 3)
)
sparse_scores <- compute_clsi(
  labels_with_na,
  sparse_nn,
  global_labels = labels_with_na,
  label_name = "sparse labels"
)
stopifnot(is.na(sparse_scores[[1]]), is.na(sparse_scores[[2]]))

single_label_error <- tryCatch(
  {
    compute_clsi(rep("A", 4), nn_idx, label_name = "single label")
    FALSE
  },
  error = function(e) TRUE
)
stopifnot(single_label_error)

message("All cLSI helper tests passed.")
