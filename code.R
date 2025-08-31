########################################################
# SBM vs HSBM (open r-hyperedges) on Enron
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
library(RcppAlgos)

set.seed(231)

# ========================================================

# Fit SBM (blockmodels package) and return best membership & params
fit_sbm_blockmodels <- function(A) {
  bm <- BM_bernoulli("SBM", A)
  bm$estimate()
  best <- which.max(bm$ICL)
  memb <- apply(bm$memberships[[best]]$Z, 1, which.max)
  params <- bm$model_parameters[[best]]
  list(bm = bm, q_best = best, memb = memb, params = params)
}

# Build OPEN r-hyperedges from a simple graph:
#     include any r-tuple whose induced subgraph has >= min_edges edges.
#     To avoid combinatorial blowup, you can cap 'max_combos'.
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
  # compute edges in each induced subgraph via adjacency submatrix sum / 2
  for (i in seq_len(nrow(cmb))) {
    idx <- cmb[i, ]
    ecount_induced <- sum(A[idx, idx]) / 2
    if (ecount_induced >= min_edges) {
      sub_gr <- induced_subgraph(g, idx)
      keep[i] <- is_connected(sub_gr)
    }
  }
  if (!any(keep)) return(list())
  split(cmb[keep, , drop = FALSE], seq_len(sum(keep)))  # list of integer vectors (size r)
}

# 2-section (motif projection) from hyperedges for any r>=3:
#     W[i,j] = # of hyperedges containing both i and j.
two_section_from_hyperedges <- function(n, hyperedges) {
  W <- Matrix(0, n, n, sparse = TRUE)
  for (h in hyperedges) {
    pairs <- t(combn(h, 2))
    W[pairs] <- W[pairs] + 1
    W[pairs[,2:1, drop = FALSE]] <- W[pairs[,2:1, drop = FALSE]] + 1
  }
  as.matrix(W)
}

# Edge link prediction AUC for SBM on adjacency
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

split_hyperedges <- function(H, test_frac = 0.1) {
  n_test <- max(1, floor(length(H) * test_frac))
  idx <- sample(seq_along(H), n_test)
  list(train = H[-idx], test = H[idx])
}
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

auc_hyperedge_prediction <- function(n, hyperedges, test_frac = 0.1) {
  sp <- split_hyperedges(hyperedges, test_frac)
  auc_motif_hyperedge_with_split(n, sp$train, sp$test)
}

auc_true_hsbm_with_split <- function(n, H_train, H_test, res,
                                     align_by_train = TRUE,
                                     auto_flip = TRUE) {
  if (length(H_train) < 2 || length(H_test) < 2) return(NA_real_)
  
  # --- labels from fit ---
  z <- res$Z
  r <- length(H_train[[1]])
  
  # --- optional alignment: map hyperedges to 1..length(z) via TRAIN vertex order ---
  if (align_by_train || n != length(z) || max(unlist(H_train)) > length(z)) {
    Vord <- sort(unique(unlist(H_train)))
    if (length(Vord) != length(z))
      stop(sprintf("res$Z length (%d) != |V(train)| (%d). Fit must be on H_train.",
                   length(z), length(Vord)))
    remap <- function(H) {
      out <- lapply(H, function(h) {
        m <- match(h, Vord); if (any(is.na(m))) return(NULL)
        sort(as.integer(m))
      })
      Filter(Negate(is.null), out)
    }
    H_train <- remap(H_train)
    H_test  <- remap(H_test)
    n <- length(z)
    if (length(H_test) < 2) return(NA_real_)
  }
  
  # --- pair_score + set_score (same structure as motif function) ---
  pair_score <- function(i, j) as.integer(z[i] == z[j])
  set_score  <- function(h) {
    P <- t(combn(h, 2))
    mean(pair_score(P[,1], P[,2]))
  }
  
  # --- auto-orient score direction using TRAIN (assortative vs heterophilous) ---
  if (auto_flip) {
    m <- min(length(H_train), 2000L)
    idx <- sample.int(length(H_train), m)
    s_tr <- vapply(H_train[idx], set_score, numeric(1))
    s_rn <- replicate(m, { h <- sort(sample.int(length(z), r)); set_score(h) })
    if (mean(s_tr) < mean(s_rn)) {
      old_set_score <- set_score
      set_score <- function(h) 1 - old_set_score(h)
    }
  }
  
  # --- positives (H_test) ---
  s_pos <- vapply(H_test, set_score, numeric(1))
  
  # --- negatives: random unseen r-sets (avoid train & test) ---
  key  <- function(h) paste(sort(h), collapse = "-")
  seen <- new.env(hash = TRUE, parent = emptyenv())
  for (h in H_train) assign(key(h), TRUE, envir = seen)
  for (h in H_test)  assign(key(h), TRUE, envir = seen)
  
  s_neg <- numeric(0)
  while (length(s_neg) < length(s_pos)) {
    cand <- sort(sample.int(length(z), r)); kc <- key(cand)
    if (!exists(kc, envir = seen, inherits = FALSE)) {
      s_neg <- c(s_neg, set_score(cand))
      assign(kc, TRUE, envir = seen)
    }
  }
  
  pROC::roc(c(rep(1, length(s_pos)), rep(0, length(s_neg))),
            c(s_pos, s_neg))$auc
}

plot_icl_curve <- function(bm, main = "ICL vs Q") {
  qs <- seq_along(bm$ICL)
  plot(qs, bm$ICL, type = "b", pch = 19, xlab = "Q", ylab = "ICL", main = main)
  abline(v = which.max(bm$ICL), lty = 2, col = "grey40")
}

# ========================================================
data(enron)
enron <- upgrade_graph(enron)
graph <- enron
# ========================================================
# 1) SBM 
# ========================================================

adj <- as_adjacency_matrix(graph, sparse = FALSE)
adj <- 1 * (((adj + t(adj)) > 0))  # undirected simple
diag(adj) <- 0

# Visualise matrix structure
SparseM::image(adj, col = c("green", "red"))

# Fit SBM (blockmodels)
bm <- BM_bernoulli("SBM", adj_enron)
bm$estimate()
best_q <- which.max(bm$ICL)
plot_icl_curve(bm = bm, main = sprintf("ICL Q=%d", best_q))
cat("Best Q (SBM):", best_q, "\n")

model <- bm$memberships[[best_q]]
membership <- apply(model$Z, 1, which.max)
ordered_adj <- adj[order(membership), order(membership)]
SparseM::image(ordered_adj, col = c("green", "red"))


# Edge link-pred for SBM
auc_edges <- auc_link_prediction_sbm(adj, test_frac = 0.2)
cat("SBM edge AUC:", round(auc_edges, 3), "\n")

# ========================================================
# 2) Build open r-hyperedges / Motif SBM
# ========================================================
g <- as.undirected(simplify(graph, remove.multiple = TRUE, remove.loops = TRUE))
cl <- components(g)
vert_ids <- V(g)[cl$membership == which.max(cl$csize)]
# keep a manageable subgraph for hyperedge mining:
g_main <- induced_subgraph(g, vert_ids[1:min(200, length(vert_ids))])

# ---- Choose r and "openness" threshold ----
r <- 3                 # you can set r > 3
min_edges <- 2         # OPEN triangles: >=2 edges among the r nodes
max_combos <- 10000    # cap combos for speed

hyperedges <- open_hyperedges_from_graph(g_main, r = r, min_edges = min_edges, max_combos = max_combos)
cat(sprintf("r=%d, min_edges=%d -> %d hyperedges kept\n", r, min_edges, length(hyperedges)))

# 2-section (motif projection) & SBM fit
n_main <- vcount(g_main)
W <- two_section_from_hyperedges(n_main, hyperedges)
W_bin <- 1 * (W > 0); diag(W_bin) <- 0
SparseM::image(W_bin, col = c("white","black"))

fit_motif <- fit_sbm_blockmodels(W_bin)
memb_motif <- fit_motif$memb
plot_icl_curve(bm = fit_motif$bm, main = sprintf("ICL Q=%d", fit_motif$q_best))
cat("Best Q (H(motif)SBM):", fit_motif$q_best, "\n")
W_ord <- W_bin[order(memb_motif), order(memb_motif)]
SparseM::image(W_ord, col = c("white","black"))

# Hyperedge (r-set) prediction AUC using motif SBM
auc_hyper <- auc_hyperedge_prediction(n_main, hyperedges, test_frac = 0.2)
cat(sprintf("Motif-SBM hyperedge AUC (r=%d, open): %.3f\n", r, auc_hyper))

# ========================================================
# 3) True HSBM
# ========================================================

RUN_TRUE_HSBM <- TRUE
if (RUN_TRUE_HSBM && length(hyperedges) > 0) {
  # Write hyperedges to file as 1-based indices per line
  #sp <- split_hyperedges(hyperedges, test_frac = 0.2)
  sink("./HG_open.txt")
  for (h in sp$train) cat(paste(h, collapse = ","), "\n")
  sink()
  HG <- HyperSBM::import_Hypergraph(file_name = "./HG_open.txt", method = "full")
  
  # Speed: small n, r fixed, affiliation model, spectral start, looser tol
  res_hsbm <- HyperSBM::HSBM(
    Hypergraph = HG, Q = 6, start = 1, model = 2,
    tol = 1e-2, maxit_VEM = 15, maxit_FP = 20, n_threads = 2, print = TRUE
  )
  save(res_hsbm, file = "res_hsbm_enron_open_small.RData")
}


auc_h_true <- auc_true_hsbm_with_split(n = length(res_hsbm$Z), H_train = sp$train, H_test = sp$test, res = res_hsbm)
cat("HSBM (true model) AUC:", auc_h_true, "\n")






               
