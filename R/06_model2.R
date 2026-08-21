# =============================================================================
# 06_model2.R — Model 2 (mammal-clade associations), native CatBoost.
# Per continent x condition x clade, TWO covariate settings are fitted:
#   (1) clade alone                (`slope_uni`, `r2_uni_*`, model2_curves*.csv)
#   (2) clade + 4 env covariates   (`slope_env`, `r2_env_*`, model2_curves_env.csv)
# and within each, TWO models:
#   * unconstrained -> response-curve shape, out-of-bag R2, and the direction of the
#     response (the constrained model cannot express a negative one)
#   * positive-monotonic -> `slope_std`, the clade-ranking metric, + OOB R2
# The block resamples are seeded per bootstrap iteration and are therefore IDENTICAL
# across the two settings, so the with/without-environment comparison is paired and any
# difference is attributable to the environment rather than to the resampling.
#
# NB the monotonic constraint is one-sided (non-decreasing), so slope_std is >= 0 by
# construction: it ranks the STRENGTH of a positive association and is not a test of
# direction. Direction is not reported as a statistic; the unconstrained response curves
# (model2_curves*.csv, FIG_M2_shape_unconstrained_*) are where the shape can be read.
#
# Also: variance partitioning (unique clade vs shared with environment).
# Honest spatial block-bootstrap CIs (7 deg blocks); R2 out-of-bag.
# Inputs (paths in R/config.R): the predictor stack, the per-clade richness stack,
# and the Biogeography rasters folder (Relative_Richness_*.tif).
# =============================================================================
suppressPackageStartupMessages({ library(terra); library(catboost); library(dplyr); library(readr) })
source("R/config.R")
set.seed(0)
PRED  <- paths$predictor_stack
CLADE <- paths$clade_stack
RELD  <- paths$biogeo_rasters
outdir <- paths$model2; dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

CLADES <- c("Pilosa","Cetartiodactyla","Diprotodontia","Primates","Rodentia",
            "Perissodactyla","Proboscidea","Cingulata","Litopterna","Notoungulata")
# PHYLACINE labels even-toed ungulates "Cetartiodactyla"; cetaceans are excluded
# (Terrestrial==1), so the analysed group is reported as Artiodactyla.
DISP <- function(x) ifelse(x=="Cetartiodactyla", "Artiodactyla", x)
ENV <- c("vpi_mean","cold_mean","dry_mean","height")
CONTS <- c("Americas","Africa","Asia")
B <- 200   # block-bootstrap resamples; 200 so the confidence intervals are not limited
           # by resample count (at 50 the 2.5/97.5 percentiles sat near the extremes)
cbp <- list(loss_function="RMSE", depth=4, l2_leaf_reg=10, learning_rate=0.05,
            iterations=300, allow_writing_files=FALSE, logging_level="Silent")
r2 <- function(a,p) 1 - sum((a-p)^2)/sum((a-mean(a))^2)
P1 <- function(x) catboost.load_pool(matrix(x, ncol=1))
Pm <- function(M) catboost.load_pool(as.matrix(M))
# Continent assignment. The Africa and Asia boxes overlap in lon [25,55] x
# lat [-12,40], so the order they are applied in decides who wins there. Africa
# is therefore applied LAST, restricted to the west of the Afro-Asian boundary
# (Suez -> Red Sea -> Bab-el-Mandeb): that keeps East Africa African without
# pulling Arabia, the Levant, Cyprus or Turkey out of Asia.
afro_asia_lon <- function(la)
  ifelse(la >= 31.5, 25,                                    # Mediterranean: Cyprus/Levant/Turkey are Asia
  ifelse(la >= 12.5, 43 - (la - 12.5) * (43 - 32) / 19,     # Bab-el-Mandeb (43E) -> Suez (32E)
                     60))                                    # Gulf of Aden: the Horn is all African
cont_of <- function(lo,la){ out <- rep("Other", length(lo))
  out[lo>=-170 & lo<=-30 & la>=-60 & la<=80] <- "Americas"
  out[lo>= 25  & lo<=180 & la>=-12 & la<=80] <- "Asia"
  out[lo>=112  & lo<=154 & la>=-44 & la<=-10] <- "Other"     # Cape York: Australia is not analysed as a continent
  out[lo>=-20  & lo<= 55 & la>=-35 & la<=40 & lo < afro_asia_lon(la)] <- "Africa"
  out }
blocks <- function(sub, deg=7) floor((sub$lon+180)/deg)*100000 + floor((sub$lat+90)/deg)

# ---- build 1-degree tables ---------------------------------------------------
message("build 1-degree tables")
ps <- load_stack(PRED, ENV, "predictor stack (environment layers)")
cs <- rast(require_file(CLADE, "clade richness stack"))
missing_clades <- setdiff(CLADES, names(cs))
if (length(missing_clades))
  stop("Clade layers absent from ", CLADE, ": ", paste(missing_clades, collapse=", "),
       "\n  (04_data_preparation.R builds one layer per PHYLACINE order passing the >5 kg",
       " herbivore filter; CLADES above lists the ten analysed here.)", call.=FALSE)
env_v <- values(ps); envfin <- rowSums(is.finite(env_v))==length(ENV)
# 04 writes clade absence as 0; any remaining NA is outside the modelled domain and is
# masked out below anyway, so treat it as absence rather than dropping the cell.
cl_v <- values(cs[[CLADES]]); cl_v[is.na(cl_v)] <- 0
# Cell -> lon/lat and the 0.1 deg -> 1 deg aggregation index, taken from the raster
# itself rather than assuming a 3600-column global 0.1 deg grid.
n <- ncell(ps); nc <- ncol(ps); ev <- as.vector(ext(ps)); rs <- res(ps)
col0 <- (seq_len(n)-1) %% nc; row0 <- (seq_len(n)-1) %/% nc
lon <- ev[["xmin"]] + (col0+0.5)*rs[1]; lat <- ev[["ymax"]] - (row0+0.5)*rs[2]
AGG <- round(1/rs[1])                      # source cells per degree (10 at 0.1 deg)
MINSUB <- max(1, floor(0.10*AGG*AGG))      # same 10% coverage rule as Model 1
Rb <- row0 %/% AGG; Cb <- col0 %/% AGG
rel_files <- c(total="Relative_Richness_Tot.tif", prickly="Relative_Richness_Prickly.tif", thorny="Relative_Richness_Thorny.tif")
build <- function(rf){
  rel <- values(rast(file.path(RELD, rf)))[,1]
  m <- is.finite(rel) & rel>0 & envfin
  d <- data.frame(rel=rel[m], R=Rb[m], C=Cb[m], lon=lon[m], lat=lat[m])
  d <- cbind(d, env_v[m,,drop=FALSE]); d <- cbind(d, cl_v[m,,drop=FALSE])
  g <- d |> group_by(R,C) |> summarise(across(everything(), \(z) mean(z, na.rm=TRUE)), n=n(), .groups="drop") |>
    filter(n>=MINSUB) |> as.data.frame()
  g$continent <- cont_of(g$lon, g$lat); g
}
tables <- lapply(rel_files, build); names(tables) <- names(rel_files)
for (cond in names(tables)) message(sprintf("  %s: %d cells | %s", cond, nrow(tables[[cond]]),
  paste(names(table(tables[[cond]]$continent)), table(tables[[cond]]$continent), collapse=" ")))

# Accumulated local effect of column 1 of X over the grid q, centred on zero. For a
# univariate model this equals the prediction curve up to a constant, so the same
# estimator serves both the with- and without-environment conditions and any difference
# between them is the environment rather than the method.
ale_at <- function(m, X, q){
  bi <- pmin(pmax(findInterval(X[,1], q[-c(1,length(q))]),0), length(q)-2); eff <- numeric(length(q)-1)
  for (k in seq_len(length(q)-1)){ mm <- bi==(k-1); if(!any(mm)) next
    Xl <- X[mm,,drop=FALSE]; Xu <- X[mm,,drop=FALSE]; Xl[,1] <- q[k]; Xu[,1] <- q[k+1]
    eff[k] <- mean(catboost.predict(m, Pm(Xu)) - catboost.predict(m, Pm(Xl))) }
  a <- c(0, cumsum(eff)); a - mean(a) }
# THE RANKING METRIC (`slope_std`). Fully standardised slope of the constrained response:
# the least-squares slope of the fitted curve against z-scored clade richness, divided by
# the standard deviation of the response. Its units are therefore
#   SDs of relative spiny-trunk richness per SD of clade richness,
# i.e. standardised on BOTH axes, which is what makes it comparable across clades (whose
# richness scales differ) and across continents (whose response scales differ).
# The slope is measured at the DISTINCT values of a 21-point quantile grid, so it spans the
# range the data cover; unique() collapses repeated quantiles, so each distinct richness
# level contributes equally to the fit rather than each cell.
# sd(y) is the FULL-SAMPLE value, not each resample's own. That matters: it makes sd(y) a
# per-continent constant, so ranking by the raw slope and by slope_std are identical and
# cannot disagree, and the bootstrap CI reflects uncertainty in the slope alone.
# NB >= 0 by construction: it ranks the STRENGTH of a positive association, never its
# direction.
zslope <- function(yv, xv, mn, sdx){ z <- (xv-mn)/sdx; unname(coef(lm(yv ~ z))[2]) }
std_slope <- function(m, q, mn, sdx, sdy){
  if (!is.finite(sdy) || sdy <= 0 || length(q) < 3) return(NA_real_)
  zslope(catboost.predict(m, P1(q)), q, mn, sdx) / sdy }
# Same metric for the environment-adjusted model, where the clade effect is its ALE
# (the model also sees four covariates). Identical grid and scaling, so the two conditions
# live on one scale, so the two are directly comparable.
ale_std_slope <- function(m, X, mn, sdx, sdy, nq=21){
  q <- unique(quantile(X[,1], seq(0,1,length.out=nq)))
  if (length(q) < 3 || !is.finite(sdy) || sdy <= 0) return(NA_real_)
  zslope(ale_at(m, X, q), q, mn, sdx) / sdy }

# ---- core models: slope_std + R2 (mono & unconstrained) + curves --------------
message("core models (constrained slope + unconstrained shape + R2)")
rows <- list(); curves <- list(); ccur <- list(); draws <- list()
ri <- 0; ci <- 0; cci <- 0; rsi <- 0
for (cond in names(tables)){ tab <- tables[[cond]]
  for (cont in CONTS){ sub <- tab[tab$continent==cont,]; if (!nrow(sub)) next
    bid <- blocks(sub); ublk <- unique(bid); y <- sub$rel
    for (cl in CLADES){ x <- sub[[cl]]
      # Only clades ABSENT from the continent are skipped. A clade that is present but fits
      # flat is a reported result, not an exclusion -- Diprotodontia occupies 10% of Asian
      # cells and returns a zero slope; `n_present` travels with it so sparseness is visible.
      if (sum(x>0)==0) next
      grid <- seq(min(x), max(x), length.out=50)   # uniform, for plotting the curves
      # the slope is measured on a QUANTILE grid, so it is weighted by cell density
      qgrid <- unique(quantile(x, seq(0,1,length.out=21)))
      mn <- mean(x); sdx <- sd(x)
      sdy_all <- sd(y); sl_b <- rep(NA_real_, B)   # aligned to b so clades stay comparable
      r2u <- c(); r2m <- c(); pmat <- matrix(NA_real_,50,0); pmm <- matrix(NA_real_,50,0)
      # Seeded per iteration (not once per clade) so the environment-adjusted loop below
      # draws exactly the same blocks: the two conditions are paired.
      for (b in 1:B){ set.seed(b); samp <- sample(ublk, length(ublk), TRUE); mtr <- bid %in% samp; moob <- !mtr
        if (length(unique(x[mtr]))<3) next
        po <- catboost.load_pool(matrix(x[mtr],ncol=1), label=y[mtr])
        mono <- catboost.train(po, NULL, params=c(cbp, list(monotone_constraints="(1)")))
        pmm <- cbind(pmm, catboost.predict(mono, P1(grid)))
        sl_b[b] <- std_slope(mono, qgrid, mn, sdx, sdy_all)
        unc <- catboost.train(po, NULL, params=cbp)
        pmat <- cbind(pmat, catboost.predict(unc, P1(grid)))
        if (sum(moob)>=20 && length(unique(y[moob]))>1){ Xo <- P1(x[moob])
          r2u <- c(r2u, r2(y[moob], catboost.predict(unc, Xo)))
          r2m <- c(r2m, r2(y[moob], catboost.predict(mono, Xo))) } }
      if (all(!is.finite(sl_b))) next
      draws[[length(draws)+1]] <- data.frame(condition=cond, continent=cont, clade=DISP(cl),
                                             b=seq_len(B), slope_std=sl_b)
      if (!length(r2u)) r2u <- NA_real_; if (!length(r2m)) r2m <- NA_real_
      slv <- sl_b[is.finite(sl_b)]; if (!length(slv)) slv <- NA_real_
      ri <- ri+1; rows[[ri]] <- data.frame(condition=cond, continent=cont, clade=DISP(cl), n=nrow(sub),
        n_present=sum(x>0), sd_response=sd(y),   # the denominator of slope_std
        slope_std=mean(slv,na.rm=TRUE), slope_std_lo=quantile(slv,.025,na.rm=TRUE),
        slope_std_hi=quantile(slv,.975,na.rm=TRUE),
        x_max=max(x),
        r2=mean(r2u,na.rm=TRUE), r2_lo=quantile(r2u,.025,na.rm=TRUE), r2_hi=quantile(r2u,.975,na.rm=TRUE),
        r2_mono=mean(r2m,na.rm=TRUE), r2_mono_lo=quantile(r2m,.025,na.rm=TRUE), r2_mono_hi=quantile(r2m,.975,na.rm=TRUE))
      ci <- ci+1; curves[[ci]] <- data.frame(condition=cond, continent=cont, clade=DISP(cl), x=grid,
        y=rowMeans(pmat), y_lo=apply(pmat,1,function(v) quantile(v,.025)), y_hi=apply(pmat,1,function(v) quantile(v,.975)))
      cci <- cci+1; ccur[[cci]] <- data.frame(condition=cond, continent=cont, clade=DISP(cl),
        x=grid, y=rowMeans(pmm),
        y_lo=apply(pmm,1,function(v) quantile(v,.025)), y_hi=apply(pmm,1,function(v) quantile(v,.975))) }
    message(sprintf("  %s/%s done", cond, cont)) } }
res <- bind_rows(rows); write_csv(res, file.path(outdir,"model2_results.csv"))
write_csv(bind_rows(curves), file.path(outdir,"model2_curves.csv"))
write_csv(bind_rows(ccur), file.path(outdir,"model2_curves_constrained.csv"))
write_csv(bind_rows(draws), file.path(outdir,"model2_bootstrap_draws.csv"))
message(sprintf("  %d clade-continent-condition models", nrow(res)))

# ---- condition 2: the same models WITH the environmental covariates ------------
# Mirrors the loop above exactly — unconstrained model for shape + OOB R2, monotonic
# model for the ranking metric + OOB R2 — but each model now also sees the four
# environmental predictors. The block resamples are identical (same per-iteration
# seeds), so the with/without comparison is paired.
message("environment-adjusted models (shape + R2 + slope)")
mono5 <- "(1,0,0,0,0)"; erows <- list(); ecur <- list(); ecc <- list()
ei <- 0; eci <- 0; ecci <- 0
for (cond in names(tables)){ tab <- tables[[cond]]
  for (cont in CONTS){ sub <- tab[tab$continent==cont,]; if (!nrow(sub)) next
    bid <- blocks(sub); ublk <- unique(bid); y <- sub$rel
    sdy_all <- sd(y)
    for (cl in CLADES){ xc <- sub[[cl]]
      if (sum(xc>0)==0) next   # see core loop
      X <- as.matrix(sub[,c(cl,ENV)]); mn <- mean(xc); sdx <- sd(xc)
      grid <- seq(min(xc), max(xc), length.out=50)
      r2u <- c(); r2m <- c(); amat <- matrix(NA_real_,50,0)
      # constrained ALE, kept so the constrained SHAPE exists under both conditions --
      # previously only the unconstrained curve was saved for the environment setting
      cmat <- matrix(NA_real_,50,0)
      sl_e <- rep(NA_real_, B)   # aligned to b so clades are re-rankable per resample
      for (b in 1:B){ set.seed(b); mtr <- bid %in% sample(ublk,length(ublk),TRUE); moob <- !mtr
        if (length(unique(X[mtr,1]))<3) next
        Xtr <- X[mtr,,drop=FALSE]; po <- catboost.load_pool(Xtr, label=y[mtr])
        mono <- catboost.train(po, NULL, params=c(cbp, list(monotone_constraints=mono5)))
        sl_e[b] <- ale_std_slope(mono, Xtr, mn, sdx, sdy_all)
        cmat <- cbind(cmat, ale_at(mono, Xtr, grid))
        unc <- catboost.train(po, NULL, params=cbp)
        amat <- cbind(amat, ale_at(unc, Xtr, grid))
        if (sum(moob)>=20 && length(unique(y[moob]))>1){ Xo <- Pm(X[moob,,drop=FALSE])
          r2u <- c(r2u, r2(y[moob], catboost.predict(unc, Xo)))
          r2m <- c(r2m, r2(y[moob], catboost.predict(mono, Xo))) } }
      if (!any(is.finite(sl_e))) next
      if (!length(r2u)) r2u <- NA_real_; if (!length(r2m)) r2m <- NA_real_
      sle <- sl_e[is.finite(sl_e)]
      if (ncol(cmat)){ ecci <- ecci+1; ecc[[ecci]] <- data.frame(condition=cond, continent=cont, clade=DISP(cl),
        x=grid, y=rowMeans(cmat), y_lo=apply(cmat,1,function(v) quantile(v,.025)),
        y_hi=apply(cmat,1,function(v) quantile(v,.975))) }
      ei <- ei+1; erows[[ei]] <- data.frame(condition=cond, continent=cont, clade=DISP(cl),
        slope_env=mean(sle,na.rm=TRUE), slope_env_lo=quantile(sle,.025,na.rm=TRUE),
        slope_env_hi=quantile(sle,.975,na.rm=TRUE),
        r2_env_unc=mean(r2u,na.rm=TRUE), r2_env_unc_lo=quantile(r2u,.025,na.rm=TRUE), r2_env_unc_hi=quantile(r2u,.975,na.rm=TRUE),
        r2_env_mono=mean(r2m,na.rm=TRUE), r2_env_mono_lo=quantile(r2m,.025,na.rm=TRUE), r2_env_mono_hi=quantile(r2m,.975,na.rm=TRUE))
      if (ncol(amat)){ eci <- eci+1; ecur[[eci]] <- data.frame(condition=cond, continent=cont, clade=DISP(cl),
        x=grid, y=rowMeans(amat), y_lo=apply(amat,1,function(v) quantile(v,.025)),
        y_hi=apply(amat,1,function(v) quantile(v,.975))) } }
    message(sprintf("  %s/%s done", cond, cont)) } }
write_csv(bind_rows(ecur), file.path(outdir,"model2_curves_env.csv"))
write_csv(bind_rows(ecc),  file.path(outdir,"model2_curves_env_constrained.csv"))
# One row per clade x continent x condition, carrying BOTH conditions side by side.
env <- bind_rows(erows) |>
  left_join(res[,c("condition","continent","clade","r2","r2_mono",
                   "slope_std","slope_std_lo","slope_std_hi")],
            by=c("condition","continent","clade")) |>
  rename(r2_uni_unc=r2, r2_uni_mono=r2_mono,
         slope_uni=slope_std, slope_uni_lo=slope_std_lo, slope_uni_hi=slope_std_hi)
write_csv(env, file.path(outdir,"model2_env_results.csv"))

# ---- diagnostics: variance partitioning -------------------------------------
# These fits are UNCONSTRAINED on purpose. Partitioning asks how much variance a clade
# explains uniquely versus jointly with environment; a monotone constraint deliberately
# restricts the fit, so partitioning a constrained R2 would measure the constraint as much
# as the data.
message("diagnostics (variance partitioning)")
# Averaged over NSEED fold assignments. One assignment is not enough: with 44-65 blocks
# per continent the raw r2_clade of a weak predictor moves by up to 0.21 between seeds
# (Americas thorny Notoungulata: -0.178 to +0.031), larger than the quantity itself.
# The components are computed PER SEED and then averaged, so each is a paired difference
# within its own fold assignment and the fold noise largely cancels -- `unique` varies by
# ~0.03 across seeds where the raw r2_clade varies by 0.21. The *_sd columns carry that
# spread so a reader can see which parts are resolvable.
# The identity holds on the means as well as per seed:
#   unique_env + shared + unique_clade = r2_full.
NSEED <- 10
# balanced fold assignment, identical to 05_model1.R: shuffle the blocks, deal round-robin,
# so each fold holds the same number of blocks +/-1. Matters more here than in Model 1
# because the per-continent subsets have far fewer blocks (44-65) than the global 199.
foldsf <- function(sub, deg=7, k=5, seed=1){ b <- blocks(sub, deg); set.seed(seed)
  u <- sample(unique(b)); f <- setNames(rep(0:(k-1), length.out=length(u)), as.character(u))
  as.integer(f[as.character(b)]) }
cvfit <- function(sub, feats, fo){ X <- as.matrix(sub[,feats,drop=FALSE]); y <- sub$rel; oof <- rep(NA_real_, nrow(sub))
  for (f in unique(fo)){ tr <- fo!=f; te <- fo==f; if (sum(te)<10) next
    m <- catboost.train(catboost.load_pool(X[tr,,drop=FALSE], label=y[tr]), NULL, params=cbp)
    oof[te] <- catboost.predict(m, Pm(X[te,,drop=FALSE])) }
  ok <- is.finite(oof); r2(y[ok], oof[ok]) }
drows <- list(); di <- 0
for (cond in names(tables)){ tab <- tables[[cond]]
  for (cont in CONTS){ sub <- tab[tab$continent==cont,]; if (!nrow(sub)) next
    act <- CLADES[sapply(CLADES, function(c) sum(sub[[c]]>0)>0)]
    RE <- numeric(NSEED)
    RC <- matrix(NA_real_, NSEED, length(act)); RF <- matrix(NA_real_, NSEED, length(act))
    for (sd_i in seq_len(NSEED)){ fo <- foldsf(sub, seed=sd_i)
      RE[sd_i] <- cvfit(sub, ENV, fo)
      for (j in seq_along(act)){ RC[sd_i,j] <- cvfit(sub, act[j], fo)
                                 RF[sd_i,j] <- cvfit(sub, c(ENV,act[j]), fo) } }
    for (j in seq_along(act)){
      uq <- RF[,j] - RE            # what the clade adds on top of environment
      ue <- RF[,j] - RC[,j]        # what environment adds on top of the clade
      sh <- RE + RC[,j] - RF[,j]   # what either could have supplied (the confounded part)
      di <- di+1; drows[[di]] <- data.frame(condition=cond, continent=cont, clade=DISP(act[j]),
        r2_env=mean(RE),       r2_env_sd=sd(RE),
        r2_clade=mean(RC[,j]), r2_clade_sd=sd(RC[,j]),
        r2_full=mean(RF[,j]),  r2_full_sd=sd(RF[,j]),
        unique_env=mean(ue),   unique_env_sd=sd(ue),
        shared=mean(sh),       shared_sd=sd(sh),
        unique_clade=mean(uq), unique_clade_sd=sd(uq), n_seed=NSEED) }
    message(sprintf("  %s/%s: %d clades x %d fold seeds", cond, cont, length(act), NSEED)) } }
write_csv(bind_rows(drows), file.path(outdir,"model2_diagnostics.csv"))

# =============================================================================
# figures (ggplot): slope ranking (both settings), constrained + unconstrained shapes,
# and the variance partition — per condition, continents faceted
# =============================================================================
message("figures")
suppressPackageStartupMessages({ library(ggplot2); library(tidyr) })
theme_set(theme_minimal(base_size=10) + theme(panel.grid.minor=element_blank(),
          plot.title=element_text(face="bold",size=11), axis.line=element_line(colour="grey30"),
          strip.text=element_text(face="bold")))
COL <- c(Americas="#E1812C",Africa="#3A7CA5",Asia="#4C9A4C")
ED <- c("Pilosa","Cingulata","Litopterna","Notoungulata","Diprotodontia")
ORDER <- c("Artiodactyla","Perissodactyla","Rodentia","Primates","Proboscidea","Pilosa","Cingulata","Litopterna","Notoungulata","Diprotodontia")
lab_ed <- function(c) ifelse(c %in% ED, paste0(c," (ED)"), c)
res <- read_csv(file.path(outdir,"model2_results.csv"), show_col_types=FALSE)
unc <- read_csv(file.path(outdir,"model2_curves.csv"), show_col_types=FALSE)
con <- read_csv(file.path(outdir,"model2_curves_constrained.csv"), show_col_types=FALSE)
cenv <- read_csv(file.path(outdir,"model2_curves_env.csv"), show_col_types=FALSE)
cenvc <- read_csv(file.path(outdir,"model2_curves_env_constrained.csv"), show_col_types=FALSE)
dia <- read_csv(file.path(outdir,"model2_diagnostics.csv"), show_col_types=FALSE)
env <- read_csv(file.path(outdir,"model2_env_results.csv"), show_col_types=FALSE)
for (cond in c("total","prickly","thorny")){
  # ---- ranking, both settings, all three conditions --------------------------
  # The species-drop replicates exist only for `total` (06_model2_jackknife.R does not
  # jackknife the subsets, because removing species from prickly/thorny changes what the
  # subset means), so this figure carries what IS defined for every condition: the slope
  # under both settings with its block-bootstrap interval, plus the share of explained
  # variance that cannot be attributed to either predictor. FIG_M2_jackknife_slope adds
  # the species-drop points on top of this, for `total` only.
  STG <- c("clade alone", "clade + environment")
  ee  <- env |> filter(condition==cond) |> mutate(continent=factor(as.character(continent),CONTS))
  ord <- ee |> arrange(continent, slope_uni)
  LEV <- unique(paste(as.integer(ord$continent), ord$clade))
  mk  <- function(ct, cl) factor(paste(as.integer(ct), cl), levels=LEV)
  rk <- bind_rows(
    ee |> transmute(continent, clade, setting=STG[1], slope=slope_uni,
                    lo=slope_uni_lo, hi=slope_uni_hi),
    ee |> transmute(continent, clade, setting=STG[2], slope=slope_env,
                    lo=slope_env_lo, hi=slope_env_hi)) |>
    mutate(setting=factor(setting, STG), key=mk(continent, clade))
  cf <- dia |> filter(condition==cond) |>
    mutate(continent=factor(as.character(continent),CONTS),
           txt=ifelse(shared>0 & shared>2*shared_sd,
                      sprintf("%.0f%%", 100*shared/r2_full), "0%"),
           setting=factor(STG[2], STG), key=mk(continent, clade))
  xr <- max(rk$hi, na.rm=TRUE) * 1.14
  ggsave(file.path(outdir,paste0("FIG_M2_rank_",cond,".png")),
    ggplot(rk, aes(slope, key, colour=continent)) +
      geom_vline(xintercept=0, linetype=2, colour="grey40") +
      geom_errorbarh(aes(xmin=lo, xmax=hi), height=0, linewidth=.5, alpha=.85) +
      geom_point(size=2.4) +
      geom_text(data=cf, aes(x=xr, y=key, label=txt), inherit.aes=FALSE,
                hjust=1, size=2.6, colour="grey25") +
      expand_limits(x=xr) +
      scale_colour_manual(values=COL, guide="none") +
      scale_y_discrete(labels=function(x) sub("^[0-9]+ ", "", x)) +
      facet_grid(continent ~ setting, scales="free_y", space="free_y") +
      labs(x="Standardised slope (SDs of relative richness per SD of clade richness)", y=NULL,
      title=sprintf("Model 2 — clade ranking with and without environmental covariates, relative %s-trunk richness",cond)),
    width=13, height=7, dpi=170)
  su <- unc |> filter(condition==cond); sc <- con |> filter(condition==cond)
  cl <- intersect(ORDER, unique(su$clade))
  su$lab <- factor(lab_ed(su$clade), levels=lab_ed(cl)); sc$lab <- factor(lab_ed(sc$clade), levels=lab_ed(cl))
  ggsave(file.path(outdir,paste0("FIG_M2_shape_unconstrained_",cond,".png")),
    ggplot() + geom_ribbon(data=su, aes(x,ymin=y_lo,ymax=y_hi,fill=continent),alpha=.13) +
      geom_line(data=su, aes(x,y,colour=continent),linewidth=.7) +
      scale_colour_manual(values=COL,name=NULL) + scale_fill_manual(values=COL,guide="none") +
      facet_wrap(~lab, scales="free", ncol=5) + labs(x="clade richness", y="pred. rel. richness",
      title=sprintf("Model 2 — UNCONSTRAINED clade response, relative %s-trunk richness",cond),
      ) +
      theme(legend.position="bottom"),
    width=15, height=3.1*ceiling(length(cl)/5), dpi=160)
  # The constrained (non-decreasing) curve that the ranking metric is actually measured
  # on. Without this panel the constrained fit
  # is invisible: model2_curves_constrained.csv$y/y_lo/y_hi were written but never
  # plotted, so there was no way to judge how well the single number summarises the curve.
  ggsave(file.path(outdir,paste0("FIG_M2_shape_constrained_",cond,".png")),
    ggplot(sc) +
      geom_ribbon(aes(x,ymin=y_lo,ymax=y_hi,fill=continent),alpha=.13) +
      geom_line(aes(x,y,colour=continent),linewidth=.7) +
      scale_colour_manual(values=COL,name=NULL) + scale_fill_manual(values=COL,guide="none") +
      facet_wrap(~lab, scales="free", ncol=5) + labs(x="clade richness", y="pred. rel. richness",
      title=sprintf("Model 2 — monotonic-constrained clade response, relative %s-trunk richness",cond),
      ) +
      theme(legend.position="bottom"),
    width=15, height=3.1*ceiling(length(cl)/5), dpi=160)
  # same unconstrained shape under both conditions, centred so they are comparable:
  # solid = clade alone, dashed = clade + environment (ALE, already centred)
  cmp <- bind_rows(
    su |> group_by(continent, clade) |> mutate(yc = y - mean(y)) |> ungroup() |>
      transmute(continent, clade, x, yc, condition_lab = "clade alone"),
    cenv |> filter(condition==cond) |> transmute(continent, clade, x, yc = y,
      condition_lab = "clade + environment"))
  cmp$lab <- factor(lab_ed(cmp$clade), levels=lab_ed(cl))
  ggsave(file.path(outdir,paste0("FIG_M2_shape_env_",cond,".png")),
    ggplot(cmp, aes(x, yc, colour=continent, linetype=condition_lab)) +
      geom_hline(yintercept=0, colour="grey70", linewidth=.3) + geom_line(linewidth=.7) +
      scale_colour_manual(values=COL, name=NULL) +
      scale_linetype_manual(values=c("clade alone"="solid","clade + environment"="22"), name=NULL) +
      facet_wrap(~lab, scales="free", ncol=5) +
      labs(x="clade richness", y="centred response",
      title=sprintf("Model 2 — UNCONSTRAINED clade response with and without environmental covariates, relative %s-trunk richness",cond)) +
      theme(legend.position="bottom"),
    width=15, height=3.1*ceiling(length(cl)/5), dpi=160)
  # the same comparison for the CONSTRAINED model -- this is the shape the slope is
  # measured on, so it is the one that explains why a ranking moves under adjustment
  cmpc <- bind_rows(
    sc |> group_by(continent, clade) |> mutate(yc = y - mean(y)) |> ungroup() |>
      transmute(continent, clade, x, yc, condition_lab = "clade alone"),
    cenvc |> filter(condition==cond) |> transmute(continent, clade, x, yc = y,
      condition_lab = "clade + environment"))
  cmpc$lab <- factor(lab_ed(cmpc$clade), levels=lab_ed(cl))
  ggsave(file.path(outdir,paste0("FIG_M2_shape_env_constrained_",cond,".png")),
    ggplot(cmpc, aes(x, yc, colour=continent, linetype=condition_lab)) +
      geom_hline(yintercept=0, colour="grey70", linewidth=.3) + geom_line(linewidth=.7) +
      scale_colour_manual(values=COL, name=NULL) +
      scale_linetype_manual(values=c("clade alone"="solid","clade + environment"="22"), name=NULL) +
      facet_wrap(~lab, scales="free", ncol=5) +
      labs(x="clade richness", y="centred response",
      title=sprintf("Model 2 — CONSTRAINED clade response with and without environmental covariates, relative %s-trunk richness",cond),
      ) +
      theme(legend.position="bottom"),
    width=15, height=3.1*ceiling(length(cl)/5), dpi=160)
  # Variance partitioning: what the clade explains ON ITS OWN, what it adds ON TOP OF
  # environment, and what the two share. Comparing clade-alone R2 with clade+environment R2
  # does not isolate the clade -- most of the latter is environment. These fits are
  # unconstrained (see the diagnostics block).
  vp0 <- dia |> filter(condition==cond)
  # Three-way partition: what only environment supplies, what only the clade supplies, and
  # what either could have (the confounded part). They sum exactly to r2_full.
  # Drawn as three bars from zero rather than a stack: `shared` is NEGATIVE in 36 of 60
  # cells (suppression -- the two together beat the sum of their parts) and exceeds twice
  # its own seed-spread in only 31, so stacking it would give the least resolvable
  # component the most visual weight and send negative segments to the wrong side of zero.
  # Error bars are +/-1 SD across the NSEED fold assignments: a bar shorter than its own
  # error bar is not resolvable, which is the honest reading for much of `shared`.
  PARTS <- c("environment only","confounded (either could supply)","clade only")
  vp <- vp0 |>
    mutate(continent=factor(as.character(continent),CONTS)) |>
    # order clades WITHIN each panel: a plain reorder() pools across continents, which for
    # a partition puts Americas-only clades above widespread ones with larger contributions
    arrange(continent, unique_clade) |>
    mutate(key=factor(paste(as.integer(continent), clade),
                      levels=unique(paste(as.integer(continent), clade))))
  pl <- bind_rows(
    vp |> transmute(continent, key, part=PARTS[1], v=unique_env,   sd=unique_env_sd),
    vp |> transmute(continent, key, part=PARTS[2], v=shared,       sd=shared_sd),
    vp |> transmute(continent, key, part=PARTS[3], v=unique_clade, sd=unique_clade_sd)) |>
    mutate(part=factor(part, PARTS))
  # The confounded share as a percentage is NOT annotated here -- this figure already
  # shows `shared` directly as a bar with its own error bar, so a percentage would repeat
  # it. The percentage is carried on the ranking figure (FIG_M2_jackknife_slope), where
  # `shared` is otherwise invisible.
  dg <- position_dodge(width=.72)
  ggsave(file.path(outdir,paste0("FIG_M2_r2_partition_",cond,".png")),
    ggplot(pl, aes(v, key, fill=part)) +
      geom_vline(xintercept=0, colour="grey45", linewidth=.4) +
      geom_col(width=.68, position=dg) +
      geom_errorbarh(aes(xmin=v-sd, xmax=v+sd), height=0, linewidth=.4,
                     colour="grey25", position=dg) +
      geom_point(data=vp, aes(x=r2_full, y=key), inherit.aes=FALSE,
                 shape=23, size=2.1, fill="white", colour="grey15", stroke=.6) +
      scale_fill_manual(values=setNames(c("#c6d3dc","#e0b23c","#b2182b"), PARTS), name=NULL) +
      scale_y_discrete(labels=function(x) sub("^[0-9]+ ", "", x)) +
      facet_wrap(~continent, scales="free_y") +
      labs(x="Out-of-fold R²", y=NULL,
      title=sprintf("Model 2 — how the explained variance splits between clade and environment, relative %s-trunk richness",cond)) +
      theme(legend.position="bottom"),
    width=14, height=6.5, dpi=170)

}
message("DONE -> ", outdir)
