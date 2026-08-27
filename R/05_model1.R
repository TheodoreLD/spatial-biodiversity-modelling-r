# =============================================================================
# 05_model1.R 
#
# GOAL: characterise the environmental conditions where spiny-trunk richness is high,
# via ALE response shapes.
#
# Sections
#   0  inputs, collinearity and concurvity
#   1  block size from the residual correlogram (averaged over allocations + subsamples)
#   2  hyperparameter tuning by the 1-SE rule
#   3  block-size re-check with the tuned model
#   4  performance (pooled within each allocation, then averaged)
#   5  residual diagnostics
#   6  block-bootstrap ensemble -> ALE curves and CIs
#   7  resolution sweep         (design justification)
#   8  transform test           (design justification)
#   9  species-drop sensitivity (robustness)
#  10  figures
#
#
# Sections 6-9 are skippable, since 9 alone takes roughly 50 minutes:
#   M1_RUN="all"                    (default)
#   M1_RUN="core"                   sections 0-5 and 10 only
#   M1_RUN="core,ale,resolution"    any subset of ale / resolution / transform / drop
# Repetition counts (env vars): M1_CV_REPS (10) performance allocations, M1_BLOCK_REPS (5)
# correlogram allocations/subsamples, M1_CV_SEEDS (5) resolution + transform allocations,
# M1_BOOT (500) ALE bootstrap resamples, M1_DROP_REPS (300) species-drop replicates,
# M1_DROP_FRACS (0.3) deletion fractions, M1_MAX_TREES (3000) boosting-iteration ceiling.
# =============================================================================
suppressPackageStartupMessages({ library(terra); library(catboost); library(dplyr); library(readr)
  library(ggplot2); library(patchwork) })
set.seed(123)

source("R/config.R")
STACK  <- paths$predictor_stack   # outputs/raster_model/predictor_stack.tif (built by 04)
outdir <- paths$model1; dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# SETTINGS 
# =============================================================================

# ---- which sections run (see the header for the M1_RUN vocabulary) ----
.M1_RUN <- trimws(strsplit(Sys.getenv("M1_RUN", unset = "all"), "[ ,]+")[[1]])
run_sec <- function(x) "all" %in% .M1_RUN || x %in% .M1_RUN
message("sections requested: ", paste(.M1_RUN, collapse = ", "))

# ---- analysis unit: aggregate 0.1 deg -> 1 deg (AGG = 10) ----
# Chosen by the resolution sweep (resolution_sweep.csv)
AGG <- as.integer(Sys.getenv("MODEL1_AGG", unset = "10"))

# ---- repetition counts ----
# Every estimate that depends on a random draw is averaged over several of them, so these
# are the runtime levers.
NBLK  <- as.integer(Sys.getenv("M1_BLOCK_REPS", unset = "5"))    # correlogram: allocations x subsamples
NCV   <- as.integer(Sys.getenv("M1_CV_REPS",    unset = "10"))   # performance: fold allocations
NSEED <- as.integer(Sys.getenv("M1_CV_SEEDS",   unset = "5"))    # resolution + transform: allocations
# 500 gives ~12 draws in each 2.5% tail. Each replicate costs two
# fits (early stopping, then refit).
NB    <- as.integer(Sys.getenv("M1_BOOT",       unset = "500"))
# 300: the envelope is a 2.5/97.5 percentile, gives ~7.5 draws in each tail.
NREP  <- as.integer(Sys.getenv("M1_DROP_REPS",  unset = "300"))
# 30% of the fraction removed.
FRACS    <- as.numeric(strsplit(Sys.getenv("M1_DROP_FRACS", unset = "0.3"), ",")[[1]])
ALE_FRAC <- .3            

# ---- what is modelled ----
ENV <- c("vpi_mean","cold_mean","dry_mean","height")
HERB <- c("mammal_richness","consumption")
PREDS <- c(ENV, HERB); TARGET <- "trunkrichness"
NAME <- c(vpi_mean="Productivity (VPI)", cold_mean="Cold-season deciduousness",
          dry_mean="Dry-season deciduousness", height="Vegetation height",
          mammal_richness="Megafauna richness (>5 kg)", consumption="Megafauna consumption")

# ---- CatBoost ----
# shared CatBoost settings: RMSE loss, 0.8 row/feature subsampling, silent
# random_seed and thread_count pinned.
BASE <- list(loss_function = "RMSE", rsm = 0.8, subsample = 0.8,
             random_seed = 1, thread_count = 4,
             allow_writing_files = FALSE, logging_level = "Silent")
PRELIM <- list(depth=6, learning_rate=.05, l2_leaf_reg=3)  

# Ceiling on boosting iterations, so early stopping selects the tree count rather than
# hitting a cap. 
MAX_TREES <- as.integer(Sys.getenv("M1_MAX_TREES", unset = "3000"))


# =============================================================================
# INPUTS
# =============================================================================
ps <- load_stack(STACK, c(TARGET, PREDS), "predictor stack")


# =============================================================================
# SHARED MACHINERY
# =============================================================================

# ---- generic ----
# R2 of predictions yh against observed y
r2 <- function(y, yh) 1 - sum((y - yh)^2)/sum((y - mean(y))^2)

# assign each cell a spatial-block ID (deg-degree lon/lat tiles)
blk <- function(lon, lat, deg) floor((lon + 180)/deg)*100000 + floor((lat + 90)/deg)

# randomly assign each block to one of k CV folds (reproducible via seed).
folds_of <- function(b, k = 5, seed = 0){
  set.seed(seed); u <- sample(unique(b))
  f <- setNames(rep(0:(k - 1), length.out = length(u)), as.character(u))
  as.integer(f[as.character(b)])
}

# great-circle (haversine) distance in km
hav <- function(lo1, la1, lo2, la2){
  R <- 6371
  a <- sin((la2 - la1)*pi/360)^2 + cos(la1*pi/180)*cos(la2*pi/180)*sin((lo2 - lo1)*pi/360)^2
  2*R*asin(sqrt(a))
}

# ---- CatBoost ----
# create the data format for CatBoost
pool <- function(d, lab = NULL, preds = PREDS) {
  X <- as.matrix(d[, preds, drop = FALSE])
  if (is.null(lab)) catboost::catboost.load_pool(X)
  else catboost::catboost.load_pool(X, label = lab)
}
# predict with model m on data d
pr <- function(m, d) catboost::catboost.predict(m, pool(d))

# ---- ALE ----
.empty_bins <- 0L
ale_curve <- function(m, d, p, q){
  bi  <- pmin(pmax(findInterval(d[[p]], q[-c(1,length(q))]), 0), length(q)-2)
  eff <- numeric(length(q)-1)
  for (k in seq_len(length(q)-1)){ mm <- bi==(k-1)
    if (!any(mm)) { .empty_bins <<- .empty_bins + 1L; next }
    dl <- d[mm,]; du <- d[mm,]; dl[[p]] <- q[k]; du[[p]] <- q[k+1]
    eff[k] <- mean(pr(m,du) - pr(m,dl)) }
  a <- c(0, cumsum(eff))
  a - mean(a[bi + 2])
}

# create a validation set (20%) to detect when CatBoost is overfitting and should
# stop training. 
eval_mask <- function(b, seed = 0){
  set.seed(seed); u <- unique(b)
  b %in% sample(u, max(1, floor(.2*length(u))))
}
# train the model to determine the boosting iterations
fit_es <- function(d, b, par, seed = 0, target = TARGET){
  m <- eval_mask(b, seed)
  if (sum(m) < 50 || sum(!m) < 50)
    return(catboost::catboost.train(pool(d, d[[target]]), NULL,
                                    params = c(BASE, par, list(iterations = 400))))
  catboost::catboost.train(
    pool(d[!m, ], d[[target]][!m]),
    pool(d[m, ], d[[target]][m]),
    params = c(BASE, par, list(iterations = MAX_TREES, od_type = "Iter",
                               od_wait = 50, use_best_model = TRUE)))
}

fit_final <- function(d, b, par, seed = 0){
  mes <- fit_es(d, b, par, seed)
  nt  <- tryCatch(mes$tree_count, error = function(e) NULL)
  if (is.null(nt) || !is.finite(nt) || nt < 1) nt <- 400L
  list(model = catboost.train(pool(d, d[[TARGET]]), NULL,
                              params = c(BASE, par, list(iterations = as.integer(nt)))),
       n_trees = as.integer(nt))
}

# ---- resampling ----
# Bootstrap the data by resampling whole blocks with replacement: create new dataset
block_boot_rows <- function(bid){
  by <- split(seq_along(bid), bid)
  ub <- names(by)
  unlist(by[sample(ub, length(ub), replace = TRUE)], use.names = FALSE)
}

# =============================================================================
# 0. INPUTS + COLLINEARITY
# =============================================================================
message("0. INPUTS + COLLINEARITY")
df0 <- as.data.frame(ps, xy=TRUE, na.rm=TRUE) |> rename(longitude=x, latitude=y)
# ---- analysis unit: aggregate 0.1 deg -> 1 deg (AGG set in SETTINGS) --------
if (AGG > 1) {
  df <- df0 |>
    mutate(R = floor((90 - latitude)/(0.1*AGG)), C = floor((longitude + 180)/(0.1*AGG))) |>
    group_by(R, C) |>
    summarise(across(all_of(c(TARGET, PREDS, "longitude", "latitude")), \(x) mean(x, na.rm = TRUE)),
              n = n(), .groups = "drop") |>
    filter(n >= max(1, floor(0.10*AGG*AGG))) |>
    select(-n) |> as.data.frame()
} else df <- df0
message(sprintf("  analysis unit %.1f deg -> %d cells (from %d at 0.1 deg; still only 74 mapped ranges)",
                0.1*AGG, nrow(df), nrow(df0)))
write_csv(as.data.frame(t(sapply(df[c(TARGET,PREDS)], summary))) |> tibble::rownames_to_column("var"),
          file.path(outdir,"input_summary.csv"))

# Spearman rank-correlation matrix among predictors
cc <- round(cor(df[PREDS], method="spearman"),3); rownames(cc) <- NAME[PREDS]; colnames(cc) <- NAME[PREDS]
write_csv(as.data.frame(cc) |> tibble::rownames_to_column("var"), file.path(outdir,"predictor_correlation.csv"))
# VIF: regress each predictor linearly on all the others → 1/(1−R²), the degree of linear multicollinearity
vif <- sapply(PREDS, function(p) 1/(1-summary(lm(as.formula(paste(p,"~",paste(setdiff(PREDS,p),collapse="+"))), df))$r.squared))
# concurvity: 3-fold CV R² of predicting each predictor from the others with CatBoost — the nonlinear analogue of VIF
Z <- scale(df[PREDS]); set.seed(0); fc <- sample(rep(0:2, length.out=nrow(df)))
concurv <- sapply(PREDS, function(p){
  oth <- setdiff(PREDS,p); oo <- rep(NA_real_, nrow(df))
  for (f in 0:2){ tr <- fc!=f; te <- fc==f
    m <- catboost.train(catboost.load_pool(Z[tr,oth,drop=FALSE], label=Z[tr,p]), NULL,
                        params=list(loss_function="RMSE",depth=4,l2_leaf_reg=5,learning_rate=.05,
                                    iterations=200,allow_writing_files=FALSE,logging_level="Silent"))
    oo[te] <- catboost.predict(m, catboost.load_pool(Z[te,oth,drop=FALSE])) }
  max(0, r2(Z[,p], oo)) })
write_csv(tibble::tibble(predictor=as.character(NAME[PREDS]), VIF=round(as.numeric(vif),2),
          concurvity_R2=round(as.numeric(concurv),3)), file.path(outdir,"predictor_concurvity.csv"))
message("  VIF: ", paste(PREDS, round(vif,2), collapse="  "))
message("  concurvity R2: ", paste(PREDS, round(concurv,2), collapse="  "))

# =============================================================================
# 1. BLOCK SIZE FROM RESIDUAL RANGE
# =============================================================================
message("1. BLOCK SIZE")

b5 <- blk(df$longitude, df$latitude, 5)
NSUB <- min(2500, nrow(df))
bins <- seq(0,3000,150); mids <- (bins[-1]+bins[-length(bins)])/2 # define the distance classes

# calculate resiudal spatial autocorrelation (Moran-like): compare the model residuals of cells sperated
# by a given distance. Do resiudals within bin tend to be sililar?
correlogram <- function(z, D, bins){
  P <- outer(z, z)
  sapply(seq_len(length(bins)-1), function(k){
    m <- D>=bins[k] & D<bins[k+1] & D>0
    if (!any(m)) NA_real_ else mean(P[m], na.rm=TRUE) }) }

# builds a pariwise distance matrix for one random subsample of cells
sub_D <- function(idx) outer(seq_along(idx), seq_along(idx),
  function(a,b) hav(df$longitude[idx][a], df$latitude[idx][a], df$longitude[idx][b], df$latitude[idx][b]))

# how far does the spatial pattern left in the model resiudal extend.
resid_correlograms <- function(par, seeds){
  sapply(seeds, function(s){
    f5 <- folds_of(b5, 5, s); o <- rep(NA_real_, nrow(df)) # splits the spatial blocks into 5 folds
    for (f in 0:4){ tr <- f5!=f; o[f5==f] <- pr(fit_es(df[tr,], b5[tr], par), df[f5==f,]) } # predict each fold using a model trained on the others
    set.seed(1000+s); idx <- sample(nrow(df), NSUB)
    correlogram(scale((df[[TARGET]] - o)[idx])[,1], sub_D(idx), bins) }) } # calculates how similar those residuals are at different distance ranges
cross_km <- function(mv){ h <- which(mv < 0.1)[1]
  if (is.na(h)) { warning("residual correlogram never drops below 0.1 within ", max(bins),
                          " km; falling back to the widest band", call.=FALSE); max(mids) } else mids[h] }

M0 <- resid_correlograms(PRELIM, seq_len(NBLK))
mI <- rowMeans(M0, na.rm=TRUE)
range_km <- cross_km(mI)

# convert km into degree to define block size
LAT90 <- as.numeric(quantile(abs(df$latitude), .90))
KM_PER_DEG <- 111.32 * cos(LAT90 * pi/180)
deg_from_km <- function(km) max(5, ceiling(km / KM_PER_DEG))
BLOCK <- deg_from_km(range_km)
per_alloc0 <- apply(M0, 2, cross_km)
message(sprintf("  residual range ~%.0f km | %.0f km per deg lon at |lat| p90 = %.1f -> block = %d deg (per-allocation: %s)",
                range_km, KM_PER_DEG, LAT90, BLOCK, paste(sprintf("%.0f", per_alloc0), collapse=", ")))
bid <- blk(df$longitude, df$latitude, BLOCK)

# =============================================================================
# 2. TUNING (spatial CV + early stopping)
# =============================================================================
message("2. TUNING")
# 1-SE rule: among configurations whose spatial-CV RMSE is within one standard error
# of the best, take the SIMPLEST (shallowest, most regularised).
GRID <- unlist(lapply(c(3,4,6,8), function(dp) lapply(c(3,10,30), function(l2)
          list(depth=dp, learning_rate=.05, l2_leaf_reg=l2))), recursive=FALSE)
fT <- folds_of(bid,5,1)
stat <- lapply(GRID, function(g){
  e <- sapply(0:4, function(f){
    tr <- fT!=f; te <- fT==f
    if (sum(te) < 10) return(NA_real_)
    sqrt(mean((df[[TARGET]][te]-pr(fit_es(df[tr,], bid[tr], g), df[te,]))^2)) })
  e <- e[!is.na(e)]; list(par=g, mean=mean(e), se=sd(e)/sqrt(length(e))) })
for (s in stat) message(sprintf("  depth=%d l2=%2d  RMSE=%.5f", s$par$depth, s$par$l2_leaf_reg, s$mean))
mu   <- sapply(stat, `[[`, "mean"); ib <- which.min(mu)
thr  <- stat[[ib]]$mean + stat[[ib]]$se
elig <- stat[sapply(stat, function(s) s$mean <= thr)]
ord  <- order(sapply(elig, function(s) s$par$depth), -sapply(elig, function(s) s$par$l2_leaf_reg))
best <- elig[[ord[1]]]$par
message(sprintf("  best RMSE=%.5f (depth=%d); 1-SE pick -> depth=%d l2=%d",
                stat[[ib]]$mean, stat[[ib]]$par$depth, best$depth, best$l2_leaf_reg))

write_csv(data.frame(
  depth        = sapply(stat, function(s) s$par$depth),
  l2_leaf_reg  = sapply(stat, function(s) s$par$l2_leaf_reg),
  cv_rmse      = mu,
  cv_rmse_se   = sapply(stat, `[[`, "se"),
  within_1se   = mu <= thr,
  selected     = sapply(stat, function(s) s$par$depth == best$depth &&
                                          s$par$l2_leaf_reg == best$l2_leaf_reg)),
  file.path(outdir, "tuning_grid.csv"))

# =============================================================================
# 3. BLOCK-SIZE RE-CHECK WITH THE TUNED MODEL
# =============================================================================
message("3. BLOCK RE-CHECK")
# Does the tuned model give the same block size as the preliminary fit used to derive it?

# `recheck_block()` repeats the whole averaged derivation
# with the tuned hyperparameters, on the same subsample and bins, once tuning has run.
recheck_block <- function(par){
  M <- resid_correlograms(par, seq_len(NBLK))
  rk <- cross_km(rowMeans(M, na.rm=TRUE))
  list(range_km = rk, block = deg_from_km(rk),
       per_alloc = apply(M, 2, cross_km))
}

rb <- recheck_block(best)
message(sprintf("  block re-check with tuned model: range ~%.0f km -> %d deg (preliminary gave %d deg; per-allocation: %s)",
                rb$range_km, rb$block, BLOCK, paste(sprintf("%.0f", rb$per_alloc), collapse=", ")))
if (rb$block != BLOCK)
  warning(sprintf("block size depends on the preliminary hyperparameters: %d deg from the tuned model vs %d deg used. Downstream CV and bootstrap intervals assume %d deg.",
                  rb$block, BLOCK, BLOCK), call.=FALSE)
write_csv(data.frame(source=c("preliminary (depth 6, l2 3)", "tuned"),
                     depth=c(6L, best$depth), l2_leaf_reg=c(3L, best$l2_leaf_reg),
                     n_allocations=c(NBLK, NBLK),
                     range_km_mean_correlogram=c(range_km, rb$range_km),
                     block_deg=c(BLOCK, rb$block),
                     range_km_per_allocation=c(paste(sprintf("%.0f", per_alloc0), collapse="; "),
                                               paste(sprintf("%.0f", rb$per_alloc), collapse="; "))),
          file.path(outdir,"blocksize_recheck.csv"))

# =============================================================================
# 4. PERFORMANCE (repeated spatial CV + overfit gap)
# =============================================================================
message("4. PERFORMANCE")
fr <- c(); trr <- c(); pooled_v <- numeric(NCV); rmse_v <- numeric(NCV)
oof <- rep(NA_real_, nrow(df))          # allocation 1, for the plotted diagnostic panels
OOF <- matrix(NA_real_, nrow(df), NCV)  # all allocations, for the averaged diagnostics
for (rep in 1:NCV){ fo <- folds_of(bid,5,100+rep); oo <- rep(NA_real_, nrow(df))
  for (f in 0:4){ tr <- fo!=f; te <- fo==f
    m <- fit_final(df[tr,], bid[tr], best, seed=rep)$model; p <- pr(m, df[te,])
    oo[te] <- p
    fr  <- c(fr,  r2(df[[TARGET]][te], p))
    # fit_final refits on ALL of df[tr,], so every training cell really is a training cell
    trr <- c(trr, r2(df[[TARGET]][tr], pr(m, df[tr,])))
    if (rep==1) oof[te] <- p }
  OOF[, rep] <- oo
  ok <- is.finite(oo)
  pooled_v[rep] <- r2(df[[TARGET]][ok], oo[ok])
  rmse_v[rep]   <- sqrt(mean((df[[TARGET]][ok]-oo[ok])^2)) }
pooled <- mean(pooled_v)

write_csv(data.frame(pooled_oof_r2=pooled, pooled_oof_r2_sd=sd(pooled_v), # how well does the model predict spatially unseen cells overall
                     pooled_oof_r2_min=min(pooled_v), pooled_oof_r2_max=max(pooled_v),
                     n_cv_reps=NCV, train_r2=mean(trr), # is the model strongly overfitting
                     spatial_cv_r2=mean(fr), spatial_cv_sd=sd(fr), # is predictive performance reasonably consistent across different geographic test regions
                     fold_r2_p2.5=unname(quantile(fr,.025)), fold_r2_p97.5=unname(quantile(fr,.975)),
                     rmse=mean(rmse_v)),
          file.path(outdir,"performance_summary.csv"))
message(sprintf("  REPORT: pooled OOF R2=%.3f (sd %.3f over %d allocations, range %.3f-%.3f) | train R2=%.3f (gap %.3f)",
                pooled, sd(pooled_v), NCV, min(pooled_v), max(pooled_v), mean(trr), mean(trr)-pooled))
message(sprintf("  archive: mean-of-folds %.3f (SD %.3f)", mean(fr), sd(fr)))

# =============================================================================
# 5. RESIDUAL DIAGNOSTICS
# =============================================================================
message("5. DIAGNOSTICS")
# check final out-of-fold residuals to determine if it is still spatially autocorrelated
resid <- df[[TARGET]] - oof          # allocation 1, for the plotted panels

# ---- WHAT THE RESIDUALS STILL CARRY: spatial structure, and calibration -----------------
# Two questions about the out-of-fold errors. Moran's I by distance: is there spatial signal
# the model failed to capture, i.e. do nearby cells still err together? Calibration slope:
# are the predictions on the right scale, or systematically compressed towards the mean?
# (1.0 = correct scale, below 1 = under-dispersed.)
#
# Both are averaged over ALL allocations rather than taken from allocation 1: on this
# dataset allocation 1 is the least favourable of the ten (its pooled R2 equals
# pooled_oof_r2_min exactly), so its diagnostics alone would understate the model.
Mr <- do.call(cbind, lapply(seq_len(NCV), function(rp){
  rr <- df[[TARGET]] - OOF[, rp]
  sapply(seq_len(NBLK), function(s){ set.seed(2000 + 100*rp + s); idx <- sample(nrow(df), NSUB)
    correlogram(scale(rr[idx])[,1], sub_D(idx), bins) }) }))
mIr <- rowMeans(Mr, na.rm=TRUE)
cal_v <- vapply(seq_len(NCV), function(rp){
  ok <- is.finite(OOF[, rp]); unname(coef(lm(df[[TARGET]][ok] ~ OOF[ok, rp]))[2]) }, numeric(1))
write_csv(tibble::tibble(dist_km=mids, moranI_residual=mIr,
                         moranI_sd=apply(Mr, 1, sd, na.rm=TRUE),
                         n_allocations=NCV, n_subsamples=NBLK),
          file.path(outdir,"residual_correlogram.csv"))
write_csv(data.frame(calibration_slope=mean(cal_v), calibration_slope_sd=sd(cal_v),
                     calibration_slope_min=min(cal_v), calibration_slope_max=max(cal_v),
                     n_allocations=NCV),
          file.path(outdir,"calibration.csv"))
message(sprintf("  calibration slope = %.3f (sd %.3f over %d allocations)",
                mean(cal_v), sd(cal_v), NCV))
message(sprintf("  residual Moran I (75km) = %.2f -> residual spatial dependence remains; uncertainty uses spatial resampling", mIr[1]))

# ---- DOES THE BLOCK SIZE CHANGE THE ANSWER? --------------------------------------------
# BLOCK was derived in section 1 from residuals produced under 5 deg CV -- 5 being the floor
# deg_from_km can return, so the derivation starts at the smallest block it will ever use.
# That start biases it one way: narrow blocks put training cells close to test cells, the
# model absorbs some of the local structure, and its residuals therefore look LESS
# autocorrelated than they really are. So section 1 can only under-shoot.
#
# It returned 8 deg. The final model's own residuals (mIr, computed above under 8 deg CV)
# imply 9 deg -- the predicted direction, and the size of the gap. Refit once at 9 deg to
# see whether the reported R2 actually depends on the choice.
#
# Caveat when reading the output: this is ONE fold allocation, while `pooled` above is the
# mean of NCV. The allocation spread at 8 deg is roughly 0.33-0.39, so a single draw at
# 9 deg landing near 0.385 is not by itself evidence that the larger block scores better.
rk2 <- cross_km(mIr); blk2 <- deg_from_km(rk2)
if (blk2 > BLOCK) {
  b2 <- blk(df$longitude, df$latitude, blk2); fo2 <- folds_of(b2, 5, 101)
  oo2 <- rep(NA_real_, nrow(df))
  for (f in 0:4){ tr <- fo2!=f; te <- fo2==f
    if (sum(te) < 10 || sum(tr) < 50) next
    oo2[te] <- pr(fit_final(df[tr,], b2[tr], best, seed=1)$model, df[te,]) }
  ok2 <- is.finite(oo2); r2_big <- r2(df[[TARGET]][ok2], oo2[ok2])
  write_csv(data.frame(block_deg=c(BLOCK, blk2),
                       source=c("derived from 5 deg preliminary residuals",
                                "implied by the final model's own residuals"),
                       range_km=c(range_km, rk2),
                       n_blocks=c(length(unique(bid)), length(unique(b2))),
                       pooled_oof_r2=c(pooled, r2_big)),
            file.path(outdir,"blocksize_confirmatory.csv"))
  message(sprintf("  confirmatory CV at %d deg (%d blocks): pooled OOF R2 = %.3f vs %.3f at %d deg",
                  blk2, length(unique(b2)), r2_big, pooled, BLOCK))
}

if (run_sec("ale")) {
# =============================================================================
# 6. BLOCK-BOOTSTRAP ENSEMBLE -> ALE CIs
# =============================================================================
message("6. ENSEMBLE (ALE CIs)")

grids <- lapply(PREDS, function(p) unique(quantile(df[[p]], seq(0,1,length.out=41))))
names(grids) <- PREDS

# calculate how model predictions change across the observed range of one predictor
# use local changes within binns and accumulating them into an ALE response curve

# Determine how many CatBoost trees to use
ff      <- fit_final(df, bid, best, seed = 0)
n_trees <- ff$n_trees
# train final model on all the data
m_full  <- ff$model
message(sprintf("  final model: early stopping chose %d trees, refitted on all %d cells",
                as.integer(n_trees), nrow(df)))
# calcul the main ALE curve
ale_point <- setNames(lapply(PREDS, function(p) ale_curve(m_full, df, p, grids[[p]])), PREDS)

# create a bootstrap versions of the data and fit model to each
ale <- setNames(lapply(PREDS, function(p) matrix(NA, length(grids[[p]]), NB)), PREDS)
for (b in 1:NB){
  set.seed(10000 + b)
  idx <- block_boot_rows(bid)
  d <- df[idx, , drop=FALSE]
  m <- fit_final(d, bid[idx], best, seed = 0)$model
  for (p in PREDS) ale[[p]][,b] <- ale_curve(m, d, p, grids[[p]])
  if (b %% 10 == 0) message("  ", b, "/", NB)
}
  
# use the bootstrap cruves to make the confidence interval
# (NB bootstrap curves + 1 point estimate) per predictor, each over that predictor's own
# number of intervals -- not a hardcoded 40, and not excluding the point-estimate calls.
.n_bins <- (NB + 1L) * sum(vapply(grids, function(q) length(q) - 1L, integer(1)))
message(sprintf("  empty ALE intervals: %d of %d (%.3f%%)",
                .empty_bins, .n_bins, 100*.empty_bins/.n_bins))
# An empty interval leaves eff = 0, which reads as "no effect" when it means "no data", and
# because the curve is a cumsum that error propagates to every later point. It has never
# fired here -- the grid is a quantile grid, so each interval holds ~120 cells drawn from
# 11-53 distinct blocks, and the chance a bootstrap draw misses all of them is ~e^-11 at
# worst. Warn rather than only counting, so a future grid or a coarser block cannot make
# this silent.
if (.empty_bins > 0)
  warning(sprintf(paste("%d ALE intervals were empty in some bootstrap replicate and were",
                        "accumulated as zero effect. The grid is too fine for the block",
                        "structure: check unique blocks per interval before trusting the band."),
                  .empty_bins), call. = FALSE)
ale_tab <- bind_rows(lapply(PREDS, function(p){
  M <- ale[[p]]
  tibble::tibble(predictor  = as.character(NAME[[p]]),
                 x          = as.numeric(grids[[p]]),
                 ale        = as.numeric(ale_point[[p]]),
                 lo         = as.numeric(apply(M, 1, function(v) quantile(v, 0.025, na.rm=TRUE))),
                 hi         = as.numeric(apply(M, 1, function(v) quantile(v, 0.975, na.rm=TRUE))),
                 ale_boot_mean = as.numeric(rowMeans(M, na.rm=TRUE)),
                 n_boot     = NB) }))
write_csv(ale_tab, file.path(outdir,"ale_curves.csv"))

}

# =============================================================================
# 7-8. DESIGN JUSTIFICATIONS (resolution, transform)
# =============================================================================
# Fit and test the model under different resolution and select the best pooled out-of-fold R²
if (run_sec("resolution") || run_sec("transform")) local({
# fit_es, not fit_final: this sweep compares resolutions against each other, so one
# consistent estimator matters more than matching section 4's absolute R2 (see fit_final).
cvfit <- function(d, b) fit_es(d, b, best)
agg <- function(A){ if (A==1) return(df0)
  df0 |> mutate(R=floor((90-latitude)/(0.1*A)), C=floor((longitude+180)/(0.1*A))) |>
    group_by(R,C) |> summarise(across(all_of(c(TARGET,PREDS,"longitude","latitude")), \(x) mean(x,na.rm=TRUE)), n=n(), .groups="drop") |>
    filter(n >= max(1, floor(0.10*A*A))) |> select(-n) |> as.data.frame() }
cv_r2 <- function(d, deg=BLOCK, seed=0){ b<-blk(d$longitude,d$latitude,deg); fo<-folds_of(b,5,seed); oof<-rep(NA_real_,nrow(d))
  for (f in 0:4){ tr<-fo!=f; te<-fo==f; if(sum(te)<10) next; oof[te]<-pr(cvfit(d[tr,],b[tr]), d[te,]) }
  ok<-is.finite(oof); r2(d[[TARGET]][ok], oof[ok]) }

# 7. RESOLUTION ----------------------------------------------------------------
if (run_sec("resolution")) {
message("7. RESOLUTION SWEEP (", NSEED, " fold allocations each)")
resrows <- lapply(c(1,5,10,20), function(A){ d<-agg(A)
  v <- sapply(seq_len(NSEED), function(s) cv_r2(d, seed=s))
  data.frame(deg=0.1*A, n_cells=nrow(d), mean_cv_r2=round(mean(v),3), sd_cv_r2=round(sd(v),3),
             min_cv_r2=round(min(v),3), max_cv_r2=round(max(v),3), n_seeds=NSEED) })
resolution <- bind_rows(resrows); write_csv(resolution, file.path(outdir,"resolution_sweep.csv"))
print(resolution)
}

# 8. TRANSFORM (at 1 deg) -------------------------------------------------------
if (run_sec("transform")) {
message("8. TRANSFORM TEST")
d1 <- agg(AGG); y <- d1[[TARGET]]
if (any(y <= 0) || any(y >= 1))
  warning(sum(y<=0), " values <=0 and ", sum(y>=1), " values >=1 are clamped before log/logit", call.=FALSE)
EPS <- 1e-9
tf  <- list(raw=function(y)y, log=function(y)log(pmax(y,EPS)),
            logit=function(y)qlogis(pmin(pmax(y,EPS),1-EPS)), sqrt=function(y)sqrt(pmax(y,0)))
itf <- list(raw=function(z)z, log=function(z)exp(z), logit=function(z)plogis(z),
            sqrt=function(z)pmax(z,0)^2)

trows <- lapply(names(tf), function(nm){ d<-d1; d[[TARGET]]<-tf[[nm]](y); b<-blk(d$longitude,d$latitude,BLOCK)
  st <- sapply(seq_len(NSEED), function(s){ fo<-folds_of(b,5,s); oof<-rep(NA_real_,nrow(d))
    for (f in 0:4){ tr<-fo!=f; te<-fo==f; if(sum(te)<10) next; oof[te]<-itf[[nm]](pr(cvfit(d[tr,],b[tr]), d[te,])) }
    ok<-is.finite(oof); res<-y[ok]-oof[ok]
    c(r2=r2(y[ok],oof[ok]), skew=mean((res-mean(res))^3)/sd(res)^3) })
  data.frame(transform=nm, mean_cv_r2_on_raw=round(mean(st["r2",]),3), sd_cv_r2=round(sd(st["r2",]),3),
             min_cv_r2=round(min(st["r2",]),3), max_cv_r2=round(max(st["r2",]),3),
             resid_skew=round(mean(st["skew",]),2), n_seeds=NSEED) })
transform <- bind_rows(trows); write_csv(transform, file.path(outdir,"transform_test.csv")); print(transform)
}


})

# =============================================================================
# 9. SPECIES-DROP SENSITIVITY
# =============================================================================
if (run_sec("drop")) local({
message("9. SPECIES-DROP SENSITIVITY")
ABS <- paths$abs_richness; REL <- paths$rel_richness; RANGES <- paths$spiny_ranges
NAME <- c(vpi_mean="Productivity (VPI)", cold_mean="Cold-season decid.", dry_mean="Dry-season decid.",
          height="Vegetation height", mammal_richness="Megafauna richness (>5 kg)", consumption="Megafauna consumption")
ALE_ORDER <- c("mammal_richness","consumption","vpi_mean","height","cold_mean","dry_mean")

# makes sure all rasters match, reconstruct richness, get cell coordiantes and precitors, and converts range into presence raster cells
rr <- rast(REL); tmpl <- rr
compareGeom(tmpl, rast(ABS), crs=TRUE, ext=TRUE, rowcol=TRUE, res=TRUE, stopOnError=TRUE)
compareGeom(tmpl, ps, crs=TRUE, ext=TRUE, rowcol=TRUE, res=TRUE, stopOnError=TRUE)
Av <- values(rast(ABS))[,1]; Rv <- values(rr)[,1]
W <- ifelse(is.finite(Rv) & Rv>0, Av/Rv, NA)
Pv <- values(ps[[PREDS]]); predfin <- rowSums(is.finite(Pv))==length(PREDS)

n <- ncell(tmpl); nc <- ncol(tmpl); ev <- as.vector(ext(tmpl)); rs <- res(tmpl)
col0 <- (seq_len(n)-1)%%nc; row0 <- (seq_len(n)-1)%/%nc
lon <- ev[["xmin"]]+(col0+0.5)*rs[1]; lat <- ev[["ymax"]]-(row0+0.5)*rs[2]
# NOT re-derived: the global AGG comes from MODEL1_AGG and defines the unit the model was
# tuned at. Shadowing it here would pin this section to 1 deg while the model ran elsewhere.
if (round(1/rs[1]) != AGG)
  stop(sprintf("species-drop aggregation (%d) disagrees with the model's AGG (%d)",
               round(1/rs[1]), AGG), call.=FALSE)
Rb <- row0%/%AGG; Cb <- col0%/%AGG
MINSUB <- max(1, floor(0.10*AGG*AGG))

gpkgs <- list.files(RANGES, "\\.gpkg$", full.names=TRUE)
message(length(gpkgs)," species; rasterising presences")
pres <- lapply(gpkgs, function(g){ r <- rasterize(vect(g), tmpl, background=0, touches=TRUE); which(values(r)[,1]>0) })
nsp <- length(pres)

# Confirm the re-rasterised ranges reproduce the richness product Model 1 was fitted to.
Sfull <- integer(n); for (k in seq_len(nsp)) Sfull[pres[[k]]] <- Sfull[pres[[k]]] + 1L
okA <- is.finite(Av)
if (!all(Sfull[okA] == Av[okA]))
  stop(sprintf("re-rasterised richness does not reproduce ABS (%.4f%% of cells match, max diff %d)",
               100*mean(Sfull[okA]==Av[okA]), max(abs(Sfull[okA]-Av[okA]))), call.=FALSE)
message("  reconstruction check: re-rasterised richness reproduces ABS exactly")

build <- function(keep){
  S <- integer(n); for (k in keep) S[pres[[k]]] <- S[pres[[k]]]+1L
  retention <- length(keep)/nsp
  m <- is.finite(W) & predfin
  rel <- rep(NA_real_,n); rel[m] <- S[m]/(W[m]*retention)
  d <- data.frame(trunkrichness=rel[m], lon=lon[m], lat=lat[m], R=Rb[m], C=Cb[m]); d <- cbind(d, Pv[m,,drop=FALSE])
  d |> group_by(R,C) |> summarise(across(everything(), \(x) mean(x,na.rm=TRUE)), nn=n(), .groups="drop") |>
    filter(nn>=MINSUB) |> select(-nn) |> as.data.frame() }
# fit_es, not fit_final: full pool and dropped pools are compared against each other, so
# these R2 values are internally consistent but not on the pooled_oof_r2 scale.
fit <- function(d) fit_es(d, blk(d$lon, d$lat, BLOCK), best)
cvr2 <- function(d){ b<-blk(d$lon,d$lat,BLOCK); fo<-folds_of(b); oof<-rep(NA_real_,nrow(d))
  for (f in 0:4){ tr<-fo!=f; te<-fo==f; if(sum(te)<10) next
    oof[te]<-pr(fit_es(d[tr,], b[tr], best), d[te,]) }
  ok<-is.finite(oof); r2(d[[TARGET]][ok], oof[ok]) }

full <- build(seq_len(nsp))
grids <- lapply(PREDS, function(p) unique(quantile(full[[p]], seq(0,1,length.out=31)))); names(grids)<-PREDS
mfull <- fit(full); ale_full <- lapply(PREDS, function(p) ale_curve(mfull, full, p, grids[[p]])); names(ale_full)<-PREDS
cv_full <- cvr2(full)
message(sprintf("FULL (%d sp): n=%d cv_r2=%.3f", nsp, nrow(full), cv_full))

rows <- list(); shape <- list(); alest <- setNames(lapply(PREDS, function(p) list()), PREDS)
for (fr in FRACS){ for (rep in 1:NREP){ set.seed(as.integer(1e5*fr + rep))
    keep <- sort(sample(nsp, round((1-fr)*nsp))); d <- build(keep)
    rows[[length(rows)+1]] <- data.frame(frac=fr, rep=rep, n_species=length(keep), n_cells=nrow(d), cv_r2=cvr2(d))
    m <- fit(d)
    for (p in PREDS){ a <- ale_curve(m, d, p, grids[[p]])
      flat <- sd(a) == 0 || !is.finite(sd(a))
      shape[[length(shape)+1]] <- data.frame(frac=fr, rep=rep, predictor=p,
        ale_cor=if (flat) NA_real_ else cor(a, ale_full[[p]]),
        ale_rmse=sqrt(mean((a-ale_full[[p]])^2)),
        ale_amp_ratio=sd(a)/sd(ale_full[[p]]))
      if (isTRUE(all.equal(fr, ALE_FRAC))) alest[[p]][[length(alest[[p]])+1]] <- a } }
  message(sprintf("  drop %d%% done (%d reps)", as.integer(fr*100), NREP)) }

res <- bind_rows(data.frame(frac=0, rep=1L, n_species=nsp, n_cells=nrow(full), cv_r2=cv_full),
                 bind_rows(rows))
write_csv(res, file.path(outdir,"jackknife.csv"))

shp <- bind_rows(shape)
write_csv(shp |> group_by(frac, predictor) |>
            summarise(n_rep=n(), n_flat=sum(is.na(ale_cor)),
                      across(c(ale_cor, ale_rmse, ale_amp_ratio),
                             list(mean=\(x) mean(x, na.rm=TRUE), sd=\(x) sd(x, na.rm=TRUE)),
                             .names="{.col}_{.fn}"), .groups="drop"),
          file.path(outdir,"species_drop_ale_stability.csv"))

dd <- bind_rows(lapply(ALE_ORDER, function(p){ M <- do.call(cbind, alest[[p]])
  data.frame(predictor=NAME[[p]], x=grids[[p]], full=ale_full[[p]],
             lo=apply(M,1,function(v)quantile(v,.025)), hi=apply(M,1,function(v)quantile(v,.975))) }))
write_csv(dd, file.path(outdir,"species_drop_ale_envelope.csv"))
dd$predictor <- factor(dd$predictor, levels=NAME[ALE_ORDER])
ggsave(file.path(outdir,"M1_jackknife_ALE.png"),
  ggplot(dd, aes(x)) + geom_ribbon(aes(ymin=lo,ymax=hi), fill="#e8a13a", alpha=.35) +
    geom_line(aes(y=full), colour="#2c6fbb", linewidth=.8) + geom_hline(yintercept=0, colour="grey60") +
    facet_wrap(~predictor, scales="free", ncol=3) + labs(x=NULL, y="ALE",
    title=sprintf("Model 1 species-drop sensitivity: ALE stability when %d%% of species are removed (blue = all %d species)",
                  as.integer(ALE_FRAC*100), nsp),
    subtitle=sprintf("Envelope = 2.5-97.5%% of %d random-deletion replicates; response rescaled by the retention fraction so amplitudes are comparable", NREP)) +
    theme_minimal(base_size=10) + theme(plot.title=element_text(face="bold",size=11),
    plot.subtitle=element_text(size=8.5, colour="grey30"), axis.line=element_line(colour="grey30")),
  width=13, height=7, dpi=175)

})

# =============================================================================
# 10. FIGURES
# =============================================================================
message("10. FIGURES")
theme_set(theme_minimal(base_size=11) + theme(panel.grid=element_blank(),
          plot.title=element_text(face="bold", size=11), plot.title.position="plot",
          axis.line=element_line(colour="grey30")))
BLUE <- "#2c6fbb"; ok <- is.finite(oof)
ord <- as.character(NAME[c("mammal_richness","consumption","vpi_mean","height","cold_mean","dry_mean")])
write_csv(tibble::tibble(longitude=df$longitude, latitude=df$latitude, observed=df[[TARGET]], oof=oof, resid=resid),
          file.path(outdir,"model1_oof.csv"))

if (run_sec("ale")) {
at <- ale_tab; at$predictor <- factor(at$predictor, levels=ord)
ggsave(file.path(outdir,"M1_ALE.png"),
  ggplot(at, aes(x)) + geom_ribbon(aes(ymin=lo,ymax=hi),fill=BLUE,alpha=.25) +
    geom_line(aes(y=ale),colour=BLUE,linewidth=.8) + geom_hline(yintercept=0,colour="grey60",linewidth=.3) +
    # Free y as well as x. All six ALEs share response units, so a common y would let the
    # panels be compared for effect size -- but this model is read for response SHAPE, not
    # for which predictor matters most (see the README method note), and a common scale
    # flattens the weaker curves until their shape is unreadable. The caption says the
    # scales are independent so nobody compares heights across panels.
    facet_wrap(~predictor, scales="free", ncol=2) +
    scale_y_continuous(labels=scales::label_number(scale=1e3)) +
    labs(x=NULL, y=expression("ALE ("%*%10^-3*")"),
         title="Model 1 ALE with 95% pointwise block-bootstrap percentile interval",
         subtitle=paste0("Band = 2.5-97.5% across ", NB, " block-bootstrap refits, taken separately at each x ",
                         "(pointwise, not a simultaneous band).\nThe line is the single full-data fit, so it is ",
                         "not the centre of the band by construction.\nEach panel has its own y scale: read the ",
                         "shapes, not the relative heights. Rows: megafauna, vegetation, seasonality.")) +
    theme(plot.subtitle=element_text(size=9, colour="grey35"),
          strip.text=element_text(face="bold", size=10)),
  width=7.5, height=9, dpi=300)
}
# collinearity: correlation heatmap + VIF + concurvity
cmL <- as.data.frame(cc) |> tibble::rownames_to_column("var") |> tidyr::pivot_longer(-var, names_to="col", values_to="rho")
lev <- as.character(NAME[PREDS]); cmL$var <- factor(cmL$var, levels=lev); cmL$col <- factor(cmL$col, levels=lev)
ph <- ggplot(cmL, aes(col,var,fill=rho)) + geom_tile() + geom_text(aes(label=sprintf("%.2f",rho)),size=2.5) +
  scale_fill_gradient2(low="#3A7CA5",mid="white",high="#c0392b",limits=c(-1,1),name=expression(rho)) +
  labs(x=NULL,y=NULL,title="Spearman correlation") + theme(axis.text.x=element_text(angle=40,hjust=1), axis.line=element_blank())
vt <- tibble::tibble(pred=as.character(NAME[PREDS]), VIF=as.numeric(vif), C=as.numeric(concurv))
pv <- ggplot(vt, aes(VIF, reorder(pred,VIF))) + geom_col(fill="#8a9a5b") + geom_vline(xintercept=5,linetype=2,colour="grey") + labs(y=NULL,title="VIF (linear collinearity)")
pc <- ggplot(vt, aes(C, reorder(pred,C))) + geom_col(fill="#cc7777") + xlim(0,1) + labs(x=expression(R^2),y=NULL,title="Concurvity (R² from other predictors)")
ggsave(file.path(outdir,"M1_collinearity.png"), ph + pv + pc + plot_layout(widths=c(1.35,1,1)), width=17, height=5, dpi=170)
# diagnostics: obs-vs-pred, calibration, QQ, residual correlogram
dobs <- tibble::tibble(o=df[[TARGET]][ok], p=oof[ok])

p1 <- ggplot(dobs, aes(o,p)) + geom_point(alpha=.1,size=.5,colour=BLUE) + geom_abline(linetype=2) +
  labs(x="Observed", y="Predicted (OOF)", title=sprintf("Obs vs pred (allocation 1: pooled OOF R²=%.3f)", pooled_v[1]))
cal <- dobs; cal$b <- dplyr::ntile(cal$p, 10); cal <- aggregate(cbind(o,p)~b, cal, mean)  # ntile cannot fail on tied predictions
p2 <- ggplot(cal, aes(p,o)) + geom_line(colour=BLUE) + geom_point(colour=BLUE) + geom_abline(linetype=2) + labs(x="Predicted (decile)",y="Observed",title="Calibration")
p3 <- ggplot(tibble::tibble(s=resid[ok]), aes(sample=s)) + stat_qq(size=.4,colour=BLUE) + stat_qq_line(linetype=2) + labs(x="Theoretical",y="Ordered",title="Residual QQ")
cg <- tibble::tibble(dist_km=mids, moranI=mIr)
p4 <- ggplot(cg, aes(dist_km,moranI)) + geom_line(colour="#c44") + geom_point(colour="#c44",size=1) + geom_hline(yintercept=.1,linetype=3,colour="grey") + coord_cartesian(xlim=c(0,2000)) + labs(x="Distance (km)",y="Moran's I",title="Residual autocorrelation")
ggsave(file.path(outdir,"M1_diagnostics.png"), p1+p2+p3+p4+plot_layout(nrow=1), width=19, height=4.3, dpi=170)

message("DONE -> ", outdir)
