simpsons_pair_concentration <- function(labels) {
  labels <- labels[!is.na(labels)]
  n_labels <- length(labels)
  if (n_labels < 2) {
    return(NA_real_)
  }

  label_freq <- table(labels)
  sum(label_freq * (label_freq - 1)) / (n_labels * (n_labels - 1))
}

compute_clsi <- function(
    labels,
    nn_idx,
    global_labels = labels,
    label_name = "labels"
) {
  if (!is.matrix(nn_idx)) {
    stop("nn_idx must be a matrix")
  }
  if (length(labels) != nrow(nn_idx)) {
    stop("labels and nn_idx must describe the same number of cells")
  }
  if (anyNA(nn_idx) || any(nn_idx < 1) || any(nn_idx > length(labels))) {
    stop("nn_idx contains invalid cell indices")
  }

  observed_global <- global_labels[!is.na(global_labels)]
  global_freq <- table(observed_global)
  if (length(observed_global) < 2 || sum(global_freq > 0) < 2) {
    stop(label_name, " has fewer than two observed labels globally")
  }

  global_si <- simpsons_pair_concentration(observed_global)
  if (!is.finite(global_si) || global_si >= 1) {
    stop("cannot correct local Simpson index when global Simpson index is 1")
  }

  vapply(seq_len(nrow(nn_idx)), function(i) {
    local_si <- simpsons_pair_concentration(labels[nn_idx[i, ]])
    if (is.na(local_si)) {
      return(NA_real_)
    }
    (local_si - global_si) / (1 - global_si)
  }, numeric(1))
}
