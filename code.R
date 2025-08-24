########################################################
# SBM vs HSBM (r-hyperedges)
########################################################

# Packages
library(igraph)
library(sbm)
library(SparseM)
library(igraphdata)
library(HyperSBM)
library(pheatmap)
library(gtools)      
library(caret)
library(blockmodels)
library(Matrix)
library(pROC)
library(aricode)

set.seed(231)

# ========================================================
# Stuctured randomness - SparseM plots
# ========================================================

.with_par <- function(newpar, code) {
  oldpar <- par(no.readonly = TRUE); on.exit(par(oldpar), add = TRUE)
  par(newpar); force(code)
}

plot_sparse_before_after <- function(M, memb = NULL,
                                     main_before = "Raw",
                                     main_after  = "Reordered",
                                     cols = c("green","red")) {
  .with_par(list(mfrow = c(1,2)), {
    SparseM::image(M, col = cols, main = main_before)
    if (!is.null(memb)) {
      ord <- order(memb)
      SparseM::image(M[ord, ord, drop = FALSE], col = cols, main = main_after)
    } else {
      SparseM::image(M, col = cols, main = paste(main_after, "(no memb)"))
    }
  })
}

# ========================================================
# Core
# ========================================================

fit_sbm_blockmodels <- function(A) {
  bm <- BM_bernoulli("SBM", A)
  bm$estimate()
  best <- which.max(bm$ICL)
  memb <- apply(bm$memberships[[best]]$Z, 1, which.max)
  params <- bm$model_parameters[[best]]
  list(bm = bm, q_best = best, memb = memb, params = params)
}

# ICL plotting 
plot_icl_curve <- function(bm, main = "ICL vs Q") {
  qs <- seq_along(bm$ICL)
  plot(qs, bm$ICL, type = "b", pch = 19, xlab = "Q", ylab = "ICL", main = main)
  abline(v = which.max(bm$ICL), lty = 2, col = "grey40")
}
fit_sbm_with_icl <- function(A, do_plots = TRUE, label = "SBM") {
  bm <- BM_bernoulli("SBM", A)
  bm$estimate()
  q_best <- which.max(bm$ICL)
  if (do_plots) plot_icl_curve(bm, sprintf("%s: ICL (Q* = %d)", label, q_best))
  memb <- apply(bm$memberships[[q_best]]$Z, 1, which.max)
  params <- bm$model_parameters[[q_best]]
  list(bm = bm, q_best = q_best, memb = memb, params = params)
}

# r-hyperedges - keep any r-tuple with >= min_edges among its pairs
open_hyperedges_from_graph <- function(g, r = 3, min_edges = 2, max_combos = Inf) {
  stopifnot(r >= 3)
  A <- as.matrix(as_adjacency_matrix(g, sparse = FALSE))
  A <- 1 * ((A + t(A)) > 0); diag(A) <- 0
  n <- vcount(g)
  
  cmb <- combinations(n, r)
  if (is.finite(max_combos) && nrow(cmb) > max_combos) {
    cmb <- cmb[sample(seq_len(nrow(cmb)), max_combos), , drop = FALSE]
  }
  
  keep <- logical(nrow(cmb))
  for (i in seq_len(nrow(cmb))) {
    idx <- cmb[i, ]
    ecount_induced <- sum(A[idx, idx]) / 2
    keep[i] <- (ecount_induced >= min_edges)
    keep[i] <- keep[i] && is_connected(induced_subgraph(g, idx))
  }
  if (!any(keep)) return(list())
  split(cmb[keep, , drop = FALSE], seq_len(sum(keep)))
}

# 2-section (motif projection) for any r>=3: W[i,j] = #hyperedges containing {i,j}
two_section_from_hyperedges <- function(n, hyperedges) {
  W <- Matrix(0, n, n, sparse = TRUE)
  for (h in hyperedges) {
    pairs <- t(combn(h, 2))
    W[pairs] <- W[pairs] + 1
    W[pairs[,2:1, drop = FALSE]] <- W[pairs[,2:1, drop = FALSE]] + 1
  }
  as.matrix(W)
}

# SBM edge link prediction AUC
auc_link_prediction_sbm <- function(A, test_frac = 0.1) {
  A <- 1 * ((A + t(A)) > 0); diag(A) <- 0
  pos <- which(A == 1, arr.ind = TRUE); pos <- pos[pos[,1] < pos[,2], , drop = FALSE]
  if (nrow(pos) < 5) return(NA_real_)
  n_test <- max(1, floor(test_frac * nrow(pos)))
  test_idx <- sample(seq_len(nrow(pos)), n_test)
  test_edges <- pos[test_idx, , drop = FALSE]
  
  A_train <- A
  for (k in seq_len(nrow(test_edges))) {
    i <- test_edges[k,1]; j <- test_edges[k,2]
    A_train[i,j] <- 0; A_train[j,i] <- 0
  }
  
  fit <- fit_sbm_blockmodels(A_train)
  Bhat <- fit$params$pi; zhat <- fit$memb
  score_pair <- function(i, j) Bhat[zhat[i], zhat[j]]
  
  s_pos <- apply(test_edges, 1, function(x) score_pair(x[1], x[2]))
  neg <- which(A == 0, arr.ind = TRUE); neg <- neg[neg[,1] < neg[,2], , drop = FALSE]
  neg <- neg[sample(seq_len(nrow(neg)), length(s_pos)), , drop = FALSE]
  s_neg <- apply(neg, 1, function(x) score_pair(x[1], x[2]))
  
  roc(c(rep(1, length(s_pos)), rep(0, length(s_neg))),
      c(s_pos, s_neg))$auc
}

# Convenience I/O + split for hyperedges
write_hyperedges <- function(H, path) {
  sink(fn)
  for (h in H) cat(paste(h, collapse = ","), "\n")
  sink()
  HG <- HyperSBM::import_Hypergraph(file_name = fn, method = "full")
  
}
split_hyperedges <- function(H, test_frac = 0.1) {
  n_test <- max(1, floor(length(H) * test_frac))
  idx <- sample(seq_along(H), n_test)
  list(train = H[-idx], test = H[idx])
}

# Motif‑SBM hyperedge AUC with split
auc_motif_hyperedge_with_split <- function(n, H_train, H_test) {
  if (length(H_train) < 5 || length(H_test) < 5) return(NA_real_)
  W <- two_section_from_hyperedges(n, H_train)
  Wb <- 1 * (W > 0); diag(Wb) <- 0
  fit <- fit_sbm_blockmodels(Wb)
  B <- fit$params$pi; z <- fit$memb
  pair_score <- function(i,j) B[z[i], z[j]]
  set_score  <- function(h) {
    pairs <- t(combn(h, 2))
    mean(apply(pairs, 1, function(p) pair_score(p[1], p[2])))
  }
  s_pos <- vapply(H_test, set_score, numeric(1))
  
  r <- length(H_train[[1]])
  key <- function(h) paste(sort(h), collapse="-")
  seen <- new.env(hash=TRUE, parent=emptyenv()); for (h in H_train) assign(key(h), TRUE, envir=seen)
  
  s_neg <- c()
  while (length(s_neg) < length(s_pos)) {
    cand <- sort(sample.int(n, r))
    if (!exists(key(cand), envir=seen)) { s_neg <- c(s_neg, set_score(cand)); assign(key(cand), TRUE, envir=seen) }
  }
  pROC::roc(c(rep(1,length(s_pos)), rep(0,length(s_neg))), c(s_pos, s_neg))$auc
}

# quick motif AUC - wrapper (auto-split)
auc_hyperedge_prediction <- function(n, hyperedges, test_frac = 0.1) {
  sp <- split_hyperedges(hyperedges, test_frac)
  auc_motif_hyperedge_with_split(n, sp$train, sp$test)
}

# True HSBM tensor + scorer
.get_hsbm_tensor <- function(res) {
  for (nm in c("Pi","B","Theta")) if (!is.null(res[[nm]]) && length(res[[nm]]) > 0) return(res[[nm]])
  stop("HSBM tensor not found (tried Pi/B/Theta).")
}
._grid_Qr <- function(Q, r) as.matrix(expand.grid(rep(list(1:Q), r)))
build_hsbm_predictor <- function(res, r) {
  tau <- res$tau
  Tns <- .get_hsbm_tensor(res)
  Q   <- ncol(tau)
  grid <- ._grid_Qr(Q, r)                           # (Q^r) x r
  lin_idx <- as.integer(1 + rowSums((grid - 1L) * cumprod(c(1L, rep(Q, r-1L)))[1:r]))
  Tvec <- as.numeric(Tns)[lin_idx]
  function(h) {
    mats <- mapply(function(v_i, z_i) tau[v_i, z_i],
                   v_i = rep(h, each = nrow(grid)),
                   z_i = as.vector(t(grid)))
    mats <- matrix(mats, ncol = r)
    w <- exp(rowSums(log(pmax(mats, 1e-15))))
    sum(w * Tvec)
  }
}

fit_true_hsbm <- function(H_train, Q = 3, start = 1, model = 2,
                          tol = 1e-2, maxit_VEM = 20, maxit_FP = 20, n_threads = 2, verbose = TRUE) {
  fn <- tempfile(fileext = ".txt")
  write_hyperedges(H_train, fn)
  HG <- HyperSBM::import_Hypergraph(fn, method = "full")
  t0 <- proc.time()[3]
  res <- HyperSBM::HSBM(Hypergraph = HG, Q = Q, start = start, model = model,
                        tol = tol, maxit_VEM = maxit_VEM, maxit_FP = maxit_FP,
                        n_threads = n_threads, print = verbose)
  t1 <- proc.time()[3]
  list(res = res, seconds = t1 - t0)
}

auc_true_hsbm_with_split <- function(n, H_train, H_test, fit_res) {
  if (length(H_train) < 5 || length(H_test) < 5) return(NA_real_)
  r <- length(H_train[[1]])
  score <- build_hsbm_predictor(fit_res$res, r)
  s_pos <- vapply(H_test, score, numeric(1))
  
  key <- function(h) paste(sort(h), collapse="-")
  seen <- new.env(hash=TRUE, parent=emptyenv()); for (h in H_train) assign(key(h), TRUE, envir=seen)
  
  s_neg <- c()
  while (length(s_neg) < length(s_pos)) {
    cand <- sort(sample.int(n, r))
    if (!exists(key(cand), envir=seen)) { s_neg <- c(s_neg, score(cand)); assign(key(cand), TRUE, envir=seen) }
  }
  pROC::roc(c(rep(1,length(s_pos)), rep(0,length(s_neg))), c(s_pos, s_neg))$auc
}

# ========================================================
# Comparison (ICL+plots inside)
# ========================================================

compare_three_models <- function(g, H, Q_hsbm = 3,
                                 hsbm_args = list(start=1, model=2, tol=1e-2,
                                                  maxit_VEM=20, maxit_FP=20, n_threads=2, verbose = TRUE),
                                 test_frac = 0.1,
                                 do_plots = TRUE) {
  A <- as.matrix(as_adjacency_matrix(g)); A <- 1 * ((A + t(A)) > 0); diag(A) <- 0
  n <- vcount(g)
  
  # --- SBM (ICL + SparseM)
  t0 <- proc.time()[3]
  fit_sbm <- if (do_plots) fit_sbm_with_icl(A, do_plots = TRUE, label = "Adjacency SBM") else fit_sbm_blockmodels(A)
  t1 <- proc.time()[3]
  auc_edge <- auc_link_prediction_sbm(A, test_frac); time_sbm <- t1 - t0
  
  if (do_plots) {
    plot_sparse_before_after(
      A, fit_sbm$memb,
      main_before = "Adjacency (raw)",
      main_after  = sprintf("Adjacency (SBM Q*=%d)", fit_sbm$q_best),
      cols = c("green","red")
    )
  }
  
  # --- Shared split for hyperedges
  sp <- split_hyperedges(H, test_frac); H_tr <- sp$train; H_te <- sp$test
  
  # --- Motif‑SBM (ICL + SparseM visual on full 2‑section)
  t2 <- proc.time()[3]
  auc_h_motif <- auc_motif_hyperedge_with_split(n, H_tr, H_te)
  t3 <- proc.time()[3]; time_motif <- t3 - t2
  
  if (do_plots) {
    W_full <- two_section_from_hyperedges(n, H)
    Wb_full <- 1 * (W_full > 0); diag(Wb_full) <- 0
    fit_motif_full <- fit_sbm_with_icl(Wb_full, do_plots = TRUE, label = "Motif/2‑section SBM")
    plot_sparse_before_after(
      Wb_full, fit_motif_full$memb,
      main_before = "2‑section (raw)",
      main_after  = sprintf("2‑section (SBM Q*=%d)", fit_motif_full$q_best),
      cols = c("white","black")
    )
  }
  
  # --- True HSBM (fit on H_train; predict H_test)
  fit_true <- do.call(fit_true_hsbm, c(list(H_train = H_tr, Q = Q_hsbm), hsbm_args))
  auc_h_true <- auc_true_hsbm_with_split(n, H_tr, H_te, fit_true)
  
  list(
    sbm   = list(auc_edge = auc_edge, seconds = time_sbm,   Q_best = fit_sbm$q_best),
    motif = list(auc_hyper = auc_h_motif, seconds = time_motif),
    hsbm  = list(auc_hyper = auc_h_true, seconds = fit_true$seconds, Q = Q_hsbm)
  )
}

# ========================================================
# Running functions (picks Q via ICL per dataset)
# ========================================================

prep_graph <- function(g, max_nodes = 200) {
  g <- as.undirected(simplify(g, remove.loops = TRUE, remove.multiple = TRUE))
  comp <- components(g)
  g <- induced_subgraph(g, V(g)[comp$membership == which.max(comp$csize)])
  if (vcount(g) > max_nodes) g <- induced_subgraph(g, sample(V(g), max_nodes))
  g
}

load_datasets <- function() {
  data(enron, package = "igraphdata"); data(karate, package = "igraphdata"); data(dolphins, package = "igraphdata")
  list(
    enron    = upgrade_graph(enron),
    karate   = upgrade_graph(karate),
    dolphins = upgrade_graph(dolphins)
  )
}

run_on_dataset <- function(g, r = 3, min_edges = 2, max_nodes = 200, max_combos = 10000,
                           Q_hsbm = NULL, test_frac = 0.1, do_plots = TRUE,
                           hsbm_args = list(start=1, model=2, tol=1e-2, maxit_VEM=20, maxit_FP=20, n_threads=2, verbose = TRUE)) {
  
  g <- prep_graph(g, max_nodes = max_nodes)
  n <- vcount(g)
  A <- as.matrix(as_adjacency_matrix(g)); A <- 1 * ((A + t(A)) > 0); diag(A) <- 0
  
  # 1) Q from adjacency ICL
  fit_edge_icm  <- fit_sbm_with_icl(A, do_plots = do_plots, label = "Adjacency SBM")
  
  # 2) build hyperedges & 2‑section, Q from motif ICL
  H <- open_hyperedges_from_graph(g, r = r, min_edges = min_edges, max_combos = max_combos)
  cat(sprintf("[dataset] n=%d, r=%d, min_edges=%d, |H|=%d\n", n, r, min_edges, length(H)))
  if (length(H) < 10) return(NULL)
  W <- two_section_from_hyperedges(n, H)
  Wb <- 1 * (W > 0); diag(Wb) <- 0
  fit_motif_icm <- fit_sbm_with_icl(Wb, do_plots = do_plots, label = "Motif/2‑section SBM")
  
  # 3) choose Q for true HSBM (prefer motif's Q*, else adjacency's)
  Q_auto <- if (!is.null(fit_motif_icm$q_best) && !is.na(fit_motif_icm$q_best)) fit_motif_icm$q_best else fit_edge_icm$q_best
  Q_use  <- if (is.null(Q_hsbm)) Q_auto else Q_hsbm
  
  # 4) unified comparison (plots shown here if do_plots=TRUE)
  compare_three_models(
    g, H, Q_hsbm = Q_use,
    hsbm_args = hsbm_args,
    test_frac = test_frac,
    do_plots = do_plots
  )
}

# ========================================================
# Runs
# ========================================================

# Enron
data(enron)
enron<-upgrade_graph(enron)
res_enron_demo <- run_on_dataset(enron, r = 3, min_edges = 2, do_plots = TRUE, Q_hsbm = NULL)
print(res_enron_demo)

# Loop over multiple datasets with auto-Q; keep plots off to avoid flooding
datasets <- load_datasets()
results_real <- list()
first_plot_done <- FALSE
for (ds_name in names(datasets)) {
  g0 <- datasets[[ds_name]]
  for (rv in c(3,4)) {
    min_list <- unique(pmax(1, c(2, 3, choose(rv,2)-1, choose(rv,2))))
    for (me in min_list) {
      key <- sprintf("%s_r%d_me%d", ds_name, rv, me)
      cat("\n=== Running", key, "===\n")
      res <- run_on_dataset(
        g0, r = rv, min_edges = me,
        Q_hsbm = NULL,      # let ICL choose per dataset
        test_frac = 0.1,
        do_plots = FALSE
      )
      results_real[[key]] <- res
      print(res)
    }
  }
}
# saveRDS(results_real, "results_real_three_models.rds")

# ========================================================
# Synthetic sweep: vary K, r, p2_in/out, p_r_in/out
# ========================================================

simulate_hybrid <- function(K=3, n_per_block=40, p2_in=0.12, p2_out=0.02,
                            r=3, p_r_in=0.06, p_r_out=0.002) {
  n <- K * n_per_block
  z <- rep(1:K, each = n_per_block)
  # Graph edges
  A <- matrix(0, n, n)
  for (i in 1:(n-1)) for (j in (i+1):n) {
    pij <- if (z[i] == z[j]) p2_in else p2_out
    A[i,j] <- rbinom(1, 1, pij); A[j,i] <- A[i,j]
  }
  # r-hyperedges
  H <- list()
  for (k in 1:K) {
    nodes <- which(z == k); if (length(nodes) >= r) {
      cmb <- combinations(length(nodes), r)
      keep <- rbinom(nrow(cmb), 1, p_r_in) == 1
      if (any(keep)) {
        h <- nodes[cmb[keep, , drop = FALSE]]
        H <- c(H, split(h, seq_len(nrow(h))))
      }
    }
  }
  total_mixed <- choose(n, r) - sum(choose(rep(n_per_block, K), r))
  M <- rbinom(1, total_mixed, p_r_out)
  if (M > 0) {
    for (t in 1:M) {
      repeat {
        cand <- sort(sample.int(n, r))
        if (length(unique(z[cand])) > 1) { H <- c(H, list(cand)); break }
      }
    }
  }
  g <- graph_from_adjacency_matrix(A, mode = "undirected", diag = FALSE)
  list(g = g, H = H, z = z)
}

evaluate_three <- function(sim, Q_hsbm = NULL, test_frac = 0.1) {
  g <- sim$g; H <- sim$H
  if (length(H) < 10) return(NULL)
  if (is.null(Q_hsbm)) Q_hsbm <- length(unique(sim$z))
  compare_three_models(
    g, H, Q_hsbm = Q_hsbm,
    hsbm_args = list(start=1, model=2, tol=1e-2, maxit_VEM=20, maxit_FP=20, n_threads=2, verbose=TRUE),
    test_frac = test_frac,
    do_plots = TRUE
  )
}

grid <- expand.grid(
  K = c(2,3,4),
  n_per_block = c(30, 40),
  p2_in = c(0.06, 0.10, 0.16),
  p2_out = c(0.02, 0.04),
  r = c(3,4),
  p_r_in = c(0.02, 0.05, 0.10),
  p_r_out = c(0.001, 0.005)
)

synth_results <- list()
for (i in 1:nrow(grid)) {
  gpar <- grid[i, ]
  cat(sprintf("\n[SYNTH] K=%d npb=%d p2_in=%.3f p2_out=%.3f r=%d p_r_in=%.3f p_r_out=%.3f\n",
              gpar$K, gpar$n_per_block, gpar$p2_in, gpar$p2_out, gpar$r, gpar$p_r_in, gpar$p_r_out))
  sim <- simulate_hybrid(K = gpar$K, n_per_block = gpar$n_per_block,
                         p2_in = gpar$p2_in, p2_out = gpar$p2_out,
                         r = gpar$r, p_r_in = gpar$p_r_in, p_r_out = gpar$p_r_out)
  res <- evaluate_three(sim, Q_hsbm = gpar$K, test_frac = 0.1)
  synth_results[[i]] <- list(params = gpar, metrics = res)
  print(res)
}

synth_df <- do.call(rbind, lapply(seq_along(synth_results), function(i) {
  if (is.null(synth_results[[i]]$metrics)) return(NULL)
  par <- synth_results[[i]]$params
  met <- synth_results[[i]]$metrics
  data.frame(
    K = par$K, n_per_block = par$n_per_block, r = par$r,
    p2_in = par$p2_in, p2_out = par$p2_out, p_r_in = par$p_r_in, p_r_out = par$p_r_out,
    delta2 = par$p2_in - par$p2_out,
    delta_r = par$p_r_in - par$p_r_out,
    auc_edge_SBM    = met$sbm$auc_edge,
    auc_hyper_motif = met$motif$auc_hyper,
    auc_hyper_true  = met$hsbm$auc_hyper
  )
}))
# write.csv(synth_df, "synth_thresholds_three_models.csv", row.names = FALSE)
