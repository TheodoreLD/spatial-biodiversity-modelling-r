# =============================================================================
# 06_model2.R — Model 2 (mammal-clade associations), native CatBoost.
# Per continent x condition x clade, TWO covariate settings are fitted:
#   (1) clade alone                (`slope_uni`, `r2_uni_*`, model2_curves*.csv)
#   (2) clade + 4 env covariates   (`slope_env`, `r2_env_*`, model2_curves_env.csv)
# and within each, TWO models:
#   * unconstrained -> response-curve shape + out-of-bag R2
#   * positive-monotonic -> z-scored linear slope = the clade-ranking metric, + OOB R2
# The block resamples are seeded per bootstrap iteration and are therefore IDENTICAL
# across the two settings, so the with/without-environment comparison is paired and any
# difference is attributable to the environment rather than to the resampling.
#
# NB the monotonic constraint is one-sided (non-decreasing), so its slope is >= 0 by
# construction. It ranks the STRENGTH of positive association; `slope_lo > 0` is not a
# test of direction. Read the unconstrained shape for whether a response is positive.
#
# Also: variance partitioning (unique clade vs shared with environment) + residual Moran.
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
B <- 50
# Smallest standardized slope treated as a real (non-zero) association. Slopes are
# ~1e-4; anything below this is the monotonic model's zero floor.
SLOPE_MIN <- 5e-5
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
# The ranking metric: slope of a linear fit of a response curve against z-scored x.
zslope <- function(yv, xv, mn, sd){ z <- (xv-mn)/sd; unname(coef(lm(yv ~ z))[2]) }
std_slope <- function(m, grid, mn, sd) zslope(catboost.predict(m, P1(grid)), grid, mn, sd)
ale_slope <- function(m, X, mn, sd, nq=21){
  q <- unique(quantile(X[,1], seq(0,1,length.out=nq))); if (length(q)<3) return(NA_real_)
  zslope(ale_at(m, X, q), q, mn, sd) }

# ---- core models: slope + R2 (mono & unconstrained) + curves -----------------
message("core models (monotonic slope + unconstrained shape + R2)")
rows <- list(); curves <- list(); ccur <- list(); rstab <- list(); ri <- 0; ci <- 0; cci <- 0; rsi <- 0
for (cond in names(tables)){ tab <- tables[[cond]]
  for (cont in CONTS){ sub <- tab[tab$continent==cont,]; if (nrow(sub)<80) next
    bid <- blocks(sub); ublk <- unique(bid); y <- sub$rel
    slope_by_b <- list()   # clade -> slope in each bootstrap resample, for rank stability
    for (cl in CLADES){ x <- sub[[cl]]
      if (sum(x>0)<30 || length(unique(x))<5) next
      mn <- mean(x); sd <- sd(x); grid <- seq(min(x), max(x), length.out=50); zg <- (grid-mn)/sd
      # `grid` (uniform) is for plotting the response curves; the reported slope is
      # measured on a quantile grid so that it weights the x-range by data density
      # exactly as the environment-adjusted `ale_slope` does — the two are divided
      # by each other downstream (`retained`), so they must share a grid convention.
      qgrid <- unique(quantile(x, seq(0,1,length.out=21)))
      slopes_b <- rep(NA_real_, B)   # kept aligned to b so clades stay comparable per resample
      # `slope` is in response units per +1 SD of clade richness, and BOTH of those
      # scales are continent-specific -- e.g. 1 SD of Proboscidea richness is 0.23
      # species in Africa but 0.71 in Asia. Dividing by the resample's own sd(y) gives
      # a unit-free effect size (SDs of response per SD of clade) that IS comparable
      # across continents. Rank clades within a continent with `slope`; compare the
      # same clade between continents with `slope_std`.
      slopes_std_b <- rep(NA_real_, B)
      r2u <- c(); r2m <- c(); pmat <- matrix(NA_real_,50,0); pmm <- matrix(NA_real_,50,0)
      # Seeded per iteration (not once per clade) so the environment-adjusted loop below
      # draws exactly the same blocks: the two conditions are paired.
      for (b in 1:B){ set.seed(b); samp <- sample(ublk, length(ublk), TRUE); mtr <- bid %in% samp; moob <- !mtr
        if (length(unique(x[mtr]))<3) next
        po <- catboost.load_pool(matrix(x[mtr],ncol=1), label=y[mtr])
        mono <- catboost.train(po, NULL, params=c(cbp, list(monotone_constraints="(1)")))
        slopes_b[b] <- std_slope(mono, qgrid, mn, sd); pmm <- cbind(pmm, catboost.predict(mono, P1(grid)))
        sdy <- sd(y[mtr]); if (is.finite(sdy) && sdy > 0) slopes_std_b[b] <- slopes_b[b]/sdy
        unc <- catboost.train(po, NULL, params=cbp)
        pmat <- cbind(pmat, catboost.predict(unc, P1(grid)))
        if (sum(moob)>=20 && length(unique(y[moob]))>1){ Xo <- P1(x[moob])
          r2u <- c(r2u, r2(y[moob], catboost.predict(unc, Xo)))
          r2m <- c(r2m, r2(y[moob], catboost.predict(mono, Xo))) } }
      slopes <- slopes_b[is.finite(slopes_b)]
      if (!length(slopes)) next
      slope_by_b[[DISP(cl)]] <- slopes_b
      if (!length(r2u)) r2u <- NA_real_; if (!length(r2m)) r2m <- NA_real_
      sstd <- slopes_std_b[is.finite(slopes_std_b)]; if (!length(sstd)) sstd <- NA_real_
      ri <- ri+1; rows[[ri]] <- data.frame(condition=cond, continent=cont, clade=DISP(cl), n=nrow(sub),
        n_present=sum(x>0), slope=mean(slopes), slope_lo=quantile(slopes,.025), slope_hi=quantile(slopes,.975),
        # slope_std divides each resample's slope by THAT resample's sd(y), so the CI
        # carries uncertainty in both. sd_response is the full-sample sd(y), recorded for
        # reference only -- slope/sd_response is close to slope_std but not identical.
        slope_std=mean(sstd,na.rm=TRUE), slope_std_lo=quantile(sstd,.025,na.rm=TRUE),
        slope_std_hi=quantile(sstd,.975,na.rm=TRUE), sd_response=sd(y),
        r2=mean(r2u,na.rm=TRUE), r2_lo=quantile(r2u,.025,na.rm=TRUE), r2_hi=quantile(r2u,.975,na.rm=TRUE),
        r2_mono=mean(r2m,na.rm=TRUE), r2_mono_lo=quantile(r2m,.025,na.rm=TRUE), r2_mono_hi=quantile(r2m,.975,na.rm=TRUE))
      ci <- ci+1; curves[[ci]] <- data.frame(condition=cond, continent=cont, clade=DISP(cl), x=grid,
        y=rowMeans(pmat), y_lo=apply(pmat,1,function(v) quantile(v,.025)), y_hi=apply(pmat,1,function(v) quantile(v,.975)))
      # constrained curve + the reported slope drawn through its centre (so the dashed
      # overlay in the figures is the number in model2_results.csv, not a refit)
      ym <- rowMeans(pmm); b1 <- mean(slopes); b0 <- mean(ym) - b1*mean(zg)
      cci <- cci+1; ccur[[cci]] <- data.frame(condition=cond, continent=cont, clade=DISP(cl), x=grid, y=ym,
        y_lo=apply(pmm,1,function(v) quantile(v,.025)), y_hi=apply(pmm,1,function(v) quantile(v,.975)), y_fit=b0+b1*zg) }
    # ---- rank stability -------------------------------------------------------
    # Every clade saw the SAME blocks in resample b (per-iteration seed), so the
    # clades can be re-ranked within each resample. top1/top3_freq = share of
    # resamples in which the clade came 1st / in the top 3. A firm ranking gives
    # values near 1 for the leaders; a soft one spreads them out.
    if (length(slope_by_b) >= 2) {
      M <- do.call(cbind, slope_by_b); nm <- colnames(M)
      t1 <- setNames(numeric(ncol(M)), nm); t3 <- t1; nb <- 0
      for (b in seq_len(nrow(M))) {
        v <- M[b, ]; if (all(!is.finite(v))) next
        o <- order(v, decreasing = TRUE, na.last = NA); nb <- nb + 1
        t1[nm[o[1]]] <- t1[nm[o[1]]] + 1
        for (k in seq_len(min(3, length(o)))) t3[nm[o[k]]] <- t3[nm[o[k]]] + 1
      }
      rsi <- rsi + 1
      rstab[[rsi]] <- data.frame(condition = cond, continent = cont, clade = nm,
                                 top1_freq = as.numeric(t1)/nb, top3_freq = as.numeric(t3)/nb,
                                 n_resamples = nb, row.names = NULL)
    }
    message(sprintf("  %s/%s done", cond, cont)) } }
res <- bind_rows(rows); write_csv(res, file.path(outdir,"model2_results.csv"))
write_csv(bind_rows(curves), file.path(outdir,"model2_curves.csv"))
write_csv(bind_rows(ccur), file.path(outdir,"model2_curves_constrained.csv"))
rank_stab <- bind_rows(rstab); write_csv(rank_stab, file.path(outdir,"model2_rank_stability.csv"))
message(sprintf("  %d clade-continent-condition models", nrow(res)))

# ---- condition 2: the same models WITH the environmental covariates ------------
# Mirrors the loop above exactly — unconstrained model for shape + OOB R2, monotonic
# model for the ranking slope + OOB R2 — but each model now also sees the four
# environmental predictors. The block resamples are identical (same per-iteration
# seeds), so the with/without comparison is paired.
message("environment-adjusted models (shape + R2 + slope)")
mono5 <- "(1,0,0,0,0)"; erows <- list(); ecur <- list(); ei <- 0; eci <- 0
for (cond in names(tables)){ tab <- tables[[cond]]
  for (cont in CONTS){ sub <- tab[tab$continent==cont,]; if (nrow(sub)<80) next
    bid <- blocks(sub); ublk <- unique(bid); y <- sub$rel
    for (cl in CLADES){ xc <- sub[[cl]]; if (sum(xc>0)<30 || length(unique(xc))<5) next
      X <- as.matrix(sub[,c(cl,ENV)]); mn <- mean(xc); sd <- sd(xc)
      grid <- seq(min(xc), max(xc), length.out=50)
      slp <- c(); slp_std <- c(); r2u <- c(); r2m <- c(); amat <- matrix(NA_real_,50,0)
      for (b in 1:B){ set.seed(b); mtr <- bid %in% sample(ublk,length(ublk),TRUE); moob <- !mtr
        if (length(unique(X[mtr,1]))<3) next
        Xtr <- X[mtr,,drop=FALSE]; po <- catboost.load_pool(Xtr, label=y[mtr])
        mono <- catboost.train(po, NULL, params=c(cbp, list(monotone_constraints=mono5)))
        s <- ale_slope(mono, Xtr, mn, sd)
        if (is.finite(s)) { slp <- c(slp, s)
          sdy <- sd(y[mtr]); if (is.finite(sdy) && sdy > 0) slp_std <- c(slp_std, s/sdy) }
        unc <- catboost.train(po, NULL, params=cbp)
        amat <- cbind(amat, ale_at(unc, Xtr, grid))
        if (sum(moob)>=20 && length(unique(y[moob]))>1){ Xo <- Pm(X[moob,,drop=FALSE])
          r2u <- c(r2u, r2(y[moob], catboost.predict(unc, Xo)))
          r2m <- c(r2m, r2(y[moob], catboost.predict(mono, Xo))) } }
      if (!length(slp)) next
      if (!length(r2u)) r2u <- NA_real_; if (!length(r2m)) r2m <- NA_real_
      if (!length(slp_std)) slp_std <- NA_real_
      ei <- ei+1; erows[[ei]] <- data.frame(condition=cond, continent=cont, clade=DISP(cl),
        slope_env=mean(slp), slope_env_lo=quantile(slp,.025), slope_env_hi=quantile(slp,.975),
        slope_env_std=mean(slp_std,na.rm=TRUE), slope_env_std_lo=quantile(slp_std,.025,na.rm=TRUE),
        slope_env_std_hi=quantile(slp_std,.975,na.rm=TRUE),
        r2_env_unc=mean(r2u,na.rm=TRUE), r2_env_unc_lo=quantile(r2u,.025,na.rm=TRUE), r2_env_unc_hi=quantile(r2u,.975,na.rm=TRUE),
        r2_env_mono=mean(r2m,na.rm=TRUE), r2_env_mono_lo=quantile(r2m,.025,na.rm=TRUE), r2_env_mono_hi=quantile(r2m,.975,na.rm=TRUE))
      if (ncol(amat)){ eci <- eci+1; ecur[[eci]] <- data.frame(condition=cond, continent=cont, clade=DISP(cl),
        x=grid, y=rowMeans(amat), y_lo=apply(amat,1,function(v) quantile(v,.025)),
        y_hi=apply(amat,1,function(v) quantile(v,.975))) } }
    message(sprintf("  %s/%s done", cond, cont)) } }
write_csv(bind_rows(ecur), file.path(outdir,"model2_curves_env.csv"))
# One row per clade x continent x condition, carrying BOTH conditions side by side.
env <- bind_rows(erows) |>
  left_join(res[,c("condition","continent","clade","slope","slope_lo","slope_hi","r2","r2_mono")],
            by=c("condition","continent","clade")) |>
  rename(slope_uni=slope, slope_uni_lo=slope_lo, slope_uni_hi=slope_hi,
         r2_uni_unc=r2, r2_uni_mono=r2_mono)
# `retained` = how much of the raw slope survives adjustment for environment. It is a
# ratio, so it is only meaningful when the denominator is a real slope: below SLOPE_MIN
# the raw slope is numerically indistinguishable from the monotonic model's zero floor
# and the ratio explodes (e.g. 82x off a raw slope of 8e-7). Report NA there instead.
env$retained <- ifelse(env$slope_uni >= SLOPE_MIN, env$slope_env/env$slope_uni, NA_real_)
write_csv(env, file.path(outdir,"model2_env_results.csv"))

# ---- diagnostics: variance partitioning + residual Moran ---------------------
message("diagnostics (variance partitioning + residual autocorrelation)")
hav <- function(lo1,la1,lo2,la2){ R<-6371; 2*R*asin(sqrt(sin((la2-la1)*pi/360)^2 +
  cos(la1*pi/180)*cos(la2*pi/180)*sin((lo2-lo1)*pi/360)^2)) }
foldsf <- function(sub, deg=7, k=5, seed=1){ b <- blocks(sub, deg); set.seed(seed)
  u <- unique(b); f <- setNames(sample(0:(k-1), length(u), TRUE), as.character(u)); as.integer(f[as.character(b)]) }
cvfit <- function(sub, feats, fo){ X <- as.matrix(sub[,feats,drop=FALSE]); y <- sub$rel; oof <- rep(NA_real_, nrow(sub))
  for (f in unique(fo)){ tr <- fo!=f; te <- fo==f; if (sum(te)<10) next
    m <- catboost.train(catboost.load_pool(X[tr,,drop=FALSE], label=y[tr]), NULL, params=cbp)
    oof[te] <- catboost.predict(m, Pm(X[te,,drop=FALSE])) }
  ok <- is.finite(oof); list(r2=r2(y[ok], oof[ok]), resid=y-oof) }
moran <- function(sub, v, seed=3){ ok <- is.finite(v); sub <- sub[ok,]; v <- v[ok]
  set.seed(seed); i <- sample(nrow(sub), min(2000, nrow(sub))); z <- scale(v[i])[,1]
  D <- outer(seq_along(i), seq_along(i), function(a,b) hav(sub$lon[i][a],sub$lat[i][a],sub$lon[i][b],sub$lat[i][b]))
  W <- (D>0 & D<300)*1; sum(W*outer(z,z))/sum(W) }
drows <- list(); di <- 0
for (cond in names(tables)){ tab <- tables[[cond]]
  for (cont in CONTS){ sub <- tab[tab$continent==cont,]; if (nrow(sub)<80) next
    fo <- foldsf(sub); act <- CLADES[sapply(CLADES, function(c) sum(sub[[c]]>0)>=30 && length(unique(sub[[c]]))>=5)]
    r2e <- cvfit(sub, ENV, fo)$r2
    for (cl in act){ rc <- cvfit(sub, cl, fo); rf <- cvfit(sub, c(ENV,cl), fo)
      di <- di+1; drows[[di]] <- data.frame(condition=cond, continent=cont, clade=DISP(cl),
        r2_env=r2e, r2_clade=rc$r2, r2_full=rf$r2, unique_clade=rf$r2-r2e, shared=r2e+rc$r2-rf$r2,
        max_corr_env=max(abs(cor(sub[[cl]], sub[,ENV], method="spearman"))), resid_moran=moran(sub, rf$resid)) }
    message(sprintf("  %s/%s: %d clades", cond, cont, length(act))) } }
write_csv(bind_rows(drows), file.path(outdir,"model2_diagnostics.csv"))

# =============================================================================
# figures (ggplot): slope ranking, shape+constrained-slope, R2 (mono vs unc),
# confounding synthesis — per condition, continents faceted
# =============================================================================
message("figures")
suppressPackageStartupMessages(library(ggplot2))
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
env <- read_csv(file.path(outdir,"model2_env_results.csv"), show_col_types=FALSE)
cenv <- read_csv(file.path(outdir,"model2_curves_env.csv"), show_col_types=FALSE)
dia <- read_csv(file.path(outdir,"model2_diagnostics.csv"), show_col_types=FALSE)
rst <- read_csv(file.path(outdir,"model2_rank_stability.csv"), show_col_types=FALSE)
for (cond in c("total","prickly","thorny")){
  d <- res |> filter(condition==cond) |>
    left_join(rst[,c("condition","continent","clade","top3_freq")], by=c("condition","continent","clade")) |>
    mutate(continent=factor(continent,CONTS), s=slope*1e4,
           lab=paste0(lab_ed(clade), "  (", round(top3_freq*100), "%)"))
  ggsave(file.path(outdir,paste0("FIG_M2_slope_",cond,".png")),
    ggplot(d, aes(s, reorder(lab,s), colour=continent)) + geom_vline(xintercept=0,linetype=2,colour="grey40") +
      geom_errorbarh(aes(xmin=slope_lo*1e4,xmax=slope_hi*1e4),height=0,colour="grey60") + geom_point(size=2.6) +
      scale_colour_manual(values=COL,guide="none") + facet_wrap(~continent, scales="free_y") +
      labs(x="Std. monotonic slope (×10⁻⁴ per +1 SD)", y=NULL,
      title=sprintf("Model 2 — mammal-clade association with relative %s-trunk richness (ED = endemic)",cond),
      caption="(x%) = rank stability: share of spatial block-bootstrap resamples placing the clade in the top 3") +
      theme(plot.caption=element_text(hjust=0, size=7.5, colour="grey30")),
    width=14, height=5, dpi=170)
  # ---- is the SAME clade important on every continent? ------------------------
  # slope_std is unit-free, so unlike `slope` it can be compared BETWEEN continents.
  # Restricted to clades present in >= 2 continents, which is where the question means
  # anything (Pilosa/Cingulata/Litopterna/Notoungulata are Americas-only, Diprotodontia
  # Asia-only). Overlapping CIs across continents = no evidence the clade matters more
  # in one than another.
  cc <- d |> filter(is.finite(slope_std)) |> group_by(clade) |> filter(n() >= 2) |> ungroup()
  if (nrow(cc)) ggsave(file.path(outdir,paste0("FIG_M2_cross_continent_",cond,".png")),
    ggplot(cc, aes(reorder(clade, slope_std), slope_std, colour=continent)) +
      geom_hline(yintercept=0, linetype=2, colour="grey40") +
      geom_errorbar(aes(ymin=slope_std_lo, ymax=slope_std_hi), width=0,
                    position=position_dodge(width=.6), alpha=.65) +
      geom_point(size=3, position=position_dodge(width=.6)) + coord_flip() +
      scale_colour_manual(values=COL, name=NULL) +
      labs(x=NULL, y="Standardised effect (SD of relative richness per +1 SD of clade richness)",
      title=sprintf("Model 2 — the same clade across continents, relative %s-trunk richness",cond),
      caption="Unit-free, so comparable BETWEEN continents (the raw slope is not: 1 SD of clade richness means a different number of species per continent).") +
      theme(legend.position="bottom", plot.caption=element_text(hjust=0, size=7, colour="grey30")),
    width=10, height=5.5, dpi=170)
  su <- unc |> filter(condition==cond); sc <- con |> filter(condition==cond)
  cl <- intersect(ORDER, unique(su$clade))
  su$lab <- factor(lab_ed(su$clade), levels=lab_ed(cl)); sc$lab <- factor(lab_ed(sc$clade), levels=lab_ed(cl))
  ggsave(file.path(outdir,paste0("FIG_M2_shape_plus_slope_",cond,".png")),
    ggplot() + geom_ribbon(data=su, aes(x,ymin=y_lo,ymax=y_hi,fill=continent),alpha=.13) +
      geom_line(data=su, aes(x,y,colour=continent),linewidth=.7) +
      geom_line(data=sc, aes(x,y_fit,colour=continent),linetype=2,linewidth=.5) +
      scale_colour_manual(values=COL,name=NULL) + scale_fill_manual(values=COL,guide="none") +
      facet_wrap(~lab, scales="free", ncol=5) + labs(x="clade richness", y="pred. rel. richness",
      title=sprintf("Model 2 — unconstrained shape + constrained slope (dashed), relative %s-trunk richness",cond),
      caption="The dashed line is the slope of the MONOTONIC model, not a fit to the solid unconstrained curve; see FIG_M2_shape_constrained_* for it drawn on the curve it summarises.") +
      theme(legend.position="bottom", plot.caption=element_text(hjust=0,size=7,colour="grey30")),
    width=15, height=3.1*ceiling(length(cl)/5), dpi=160)
  # The constrained (non-decreasing) curve that the ranking metric is actually measured
  # on, with the reported slope drawn through it. Without this panel the constrained fit
  # is invisible: model2_curves_constrained.csv$y/y_lo/y_hi were written but never
  # plotted, so there was no way to judge how well the single number summarises the curve.
  ggsave(file.path(outdir,paste0("FIG_M2_shape_constrained_",cond,".png")),
    ggplot(sc) + geom_ribbon(aes(x,ymin=y_lo,ymax=y_hi,fill=continent),alpha=.13) +
      geom_line(aes(x,y,colour=continent),linewidth=.7) +
      geom_line(aes(x,y_fit,colour=continent),linetype=2,linewidth=.5) +
      scale_colour_manual(values=COL,name=NULL) + scale_fill_manual(values=COL,guide="none") +
      facet_wrap(~lab, scales="free", ncol=5) + labs(x="clade richness", y="pred. rel. richness",
      title=sprintf("Model 2 — monotonic-constrained shape (solid) + reported slope (dashed), relative %s-trunk richness",cond),
      caption="Constraint is one-sided (non-decreasing), so the curve cannot fall and the slope is >= 0 by construction: it ranks the STRENGTH of positive association, not its direction.") +
      theme(legend.position="bottom", plot.caption=element_text(hjust=0,size=7,colour="grey30")),
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
      title=sprintf("Model 2 — clade response shape with and without environmental covariates, relative %s-trunk richness",cond)) +
      theme(legend.position="bottom"),
    width=15, height=3.1*ceiling(length(cl)/5), dpi=160)
  ggsave(file.path(outdir,paste0("FIG_M2_r2_mono_",cond,".png")),
    ggplot(d, aes(y=reorder(clade,r2_mono))) + geom_vline(xintercept=0,linetype=2,colour="grey40") +
      geom_errorbarh(aes(xmin=r2_mono_lo,xmax=r2_mono_hi),height=0,colour="grey60") +
      geom_point(aes(x=r2_mono,colour=continent),size=2.4) + geom_point(aes(x=r2),shape=1,size=2.4,colour="grey30") +
      scale_colour_manual(values=COL,guide="none") + facet_wrap(~continent,scales="free_y") +
      labs(x="Out-of-bag R² (filled = monotonic, open = unconstrained)", y=NULL,
      title=sprintf("Model 2 — constrained vs unconstrained clade-model R², relative %s-trunk richness",cond)),
    width=14, height=5, dpi=170)
  m <- env |> filter(condition==cond) |> left_join(dia[,c("condition","continent","clade","r2_clade","shared")],
        by=c("condition","continent","clade")) |>
    mutate(continent=factor(continent,CONTS), conf=ifelse(r2_clade>0.01, pmin(pmax(shared/r2_clade,0),1), NA_real_),
           stable=slope_uni>=SLOPE_MIN, retlab=ifelse(stable & !is.na(retained), sprintf("%.0f%%",retained*100),"n/a"),
           lab=paste0(clade, ifelse(slope_env_lo>0,"*",""), "  (ret ",retlab,")"))
  ggsave(file.path(outdir,paste0("FIG_M2_confounding_synthesis_",cond,".png")),
    ggplot(m, aes(y=reorder(lab,slope_env))) + geom_vline(xintercept=0,linetype=2,colour="grey40") +
      geom_segment(aes(x=slope_uni*1e4,xend=slope_env*1e4,yend=reorder(lab,slope_env)),colour="grey75",linewidth=1) +
      geom_point(aes(x=slope_uni*1e4),shape=1,size=2.3,colour="grey50") +
      geom_errorbarh(aes(xmin=slope_env_lo*1e4,xmax=slope_env_hi*1e4),height=0,colour="grey55") +
      geom_point(aes(x=slope_env*1e4,fill=conf),shape=21,size=3.2,colour="black") +
      scale_fill_gradient(low="#fee8c8",high="#b30000",limits=c(0,1),name="confounded frac.",na.value="grey80") +
      facet_wrap(~continent,scales="free_y") + labs(x="Std. slope ×10⁻⁴ (open = raw → filled = env-adjusted)", y=NULL,
      title=sprintf("Model 2 — clade ranking robustness to environmental confounding, relative %s-trunk richness (* = env-adj CI>0)",cond)),
    width=15, height=5, dpi=165)
}
message("DONE -> ", outdir)
