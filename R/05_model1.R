# =============================================================================
# 05_model1.R — Model 1, native CatBoost.
# Response = relative spiny-trunk richness (bounded proportion).
#
# GOAL: characterise the environmental conditions where spiny-trunk richness is high,
# via ALE response shapes. Deliberately NOT a predictor-importance ranking — the two
# megafauna metrics are near-collinear (see predictor_concurvity.csv) and are read as
# one axis, with the redundancy disclosed rather than resolved by ranking.
#
# Design decisions, each backed by a companion output:
#   * analysis unit 1 deg, aggregated from 0.1 deg  -> resolution_sweep.csv
#   * no response transform                          -> transform_test.csv
#   * no spatial term                                -> spatial_term_test.csv
#   * tuning by the 1-SE rule (simplest adequate model wins)
#   * block size derived from the residual correlogram, not assumed
#   * ALE, which stays valid under correlated predictors
#   * one block-bootstrap ensemble supplies the ALE CIs
#   * early stopping on a spatially separated subset of the training blocks
# (The first three are produced by 05_model1_sensitivity.R; species-drop robustness
#  by 05_model1_jackknife.R.)
# =============================================================================
suppressPackageStartupMessages({ library(terra); library(catboost); library(dplyr); library(readr) })
set.seed(123)

# ---- paths (from R/config.R) ------------------------------------------------
source("R/config.R")
STACK  <- paths$predictor_stack   # outputs/raster_model/predictor_stack.tif (built by 04)
outdir <- paths$model1; dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

ENV <- c("vpi_mean","cold_mean","dry_mean","height")
HERB <- c("mammal_richness","consumption")
PREDS <- c(ENV, HERB); TARGET <- "trunkrichness"
NAME <- c(vpi_mean="Productivity (VPI)", cold_mean="Cold-season deciduousness",
          dry_mean="Dry-season deciduousness", height="Vegetation height",
          mammal_richness="Megafauna richness (>5 kg)", consumption="Megafauna consumption")

ps <- load_stack(STACK, c(TARGET, PREDS), "predictor stack")

# ---------------- helpers ----------------
# R² of predictions yh against observed y
r2 <- function(y,yh) 1 - sum((y-yh)^2)/sum((y-mean(y))^2)
# assign each cell a spatial-block ID (deg-degree lon/lat tiles)
blk <- function(lon,lat,deg) floor((lon+180)/deg)*100000 + floor((lat+90)/deg)
# randomly assign each block to one of k CV folds (reproducible via seed)
folds_of <- function(b,k=5,seed=0){ set.seed(seed); u<-unique(b); f<-setNames(sample(0:(k-1),length(u),TRUE),u); as.integer(f[as.character(b)]) }
# great-circle (haversine) distance in km between two lon/lat points
hav <- function(lo1,la1,lo2,la2){ R<-6371; a<-sin((la2-la1)*pi/360)^2+cos(la1*pi/180)*cos(la2*pi/180)*sin((lo2-lo1)*pi/360)^2; 2*R*asin(sqrt(a)) }
# shared CatBoost settings: RMSE loss, 0.8 row/feature subsampling, silent
BASE <- list(loss_function="RMSE", rsm=0.8, subsample=0.8, allow_writing_files=FALSE, logging_level="Silent")
# build a CatBoost data pool from the predictor columns, with or without a label
pool <- function(d, lab=NULL) if (is.null(lab)) catboost.load_pool(as.matrix(d[,PREDS])) else
  catboost.load_pool(as.matrix(d[,PREDS]), label=lab)
# train CatBoost with early stopping on a spatially-separated ~20% of blocks held out as the eval set
fit_es <- function(d, b, par, seed=0){
  set.seed(seed); u <- unique(b); ev <- sample(u, max(1, floor(.2*length(u)))); m <- b %in% ev
  if (sum(m) < 50 || sum(!m) < 50)
    return(catboost.train(pool(d, d[[TARGET]]), NULL, params=c(BASE, par, list(iterations=400))))
  catboost.train(pool(d[!m,], d[[TARGET]][!m]), pool(d[m,], d[[TARGET]][m]),
                 params=c(BASE, par, list(iterations=1000, od_type="Iter", od_wait=50, use_best_model=TRUE)))
}
# predict with model m on data d
pr <- function(m,d) catboost.predict(m, pool(d))

# =============================================================================
# 0. INPUTS + COLLINEARITY
# =============================================================================
message("0. INPUTS + COLLINEARITY")
df0 <- as.data.frame(ps, xy=TRUE, na.rm=TRUE) |> rename(longitude=x, latitude=y)
# ---- analysis unit: aggregate 0.1 deg -> 1 deg (AGG = 10) -------------------
# Chosen by the resolution sweep (resolution_sweep.csv): 1 deg maximises spatial-CV R2
# and cuts n from ~404k pseudo-replicates to ~4.9k.
AGG <- as.integer(Sys.getenv("MODEL1_AGG", unset = "10"))
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
Z <- scale(df[PREDS]); set.seed(0); fc <- sample(0:2, nrow(df), TRUE)
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
b5 <- blk(df$longitude, df$latitude, 5); f5 <- folds_of(b5, 5, 1); oof0 <- rep(NA_real_, nrow(df))
for (f in 0:4){ tr <- f5!=f; oof0[f5==f] <- pr(fit_es(df[tr,], b5[tr], list(depth=6, learning_rate=.05, l2_leaf_reg=3)), df[f5==f,]) }
res0 <- df[[TARGET]] - oof0
i <- sample(nrow(df), min(2500, nrow(df))); z <- scale(res0[i])[,1]
D <- outer(seq_along(i), seq_along(i), function(a,b) hav(df$longitude[i][a], df$latitude[i][a], df$longitude[i][b], df$latitude[i][b]))
bins <- seq(0,3000,150); mids <- (bins[-1]+bins[-length(bins)])/2

correlogram <- function(z, D, bins) sapply(seq_len(length(bins)-1), function(k){
  m <- D>=bins[k] & D<bins[k+1] & D>0
  if (!any(m)) NA_real_ else mean(outer(z,z)[m], na.rm=TRUE) })
mI <- correlogram(z, D, bins)
hit <- which(mI < 0.1)[1]
if (is.na(hit)) { warning("residual correlogram never drops below 0.1 within ", max(bins),
                          " km; falling back to the widest band", call.=FALSE)
                  range_km <- max(mids) } else range_km <- mids[hit]
BLOCK <- max(5, ceiling(range_km/111))
message(sprintf("  residual range ~%.0f km -> block = %d deg", range_km, BLOCK))
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
# The grid is written out so the 1-SE selection can be audited without re-running: every
# other design decision in Model 1 has a CSV behind it, and this one previously did not.
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
# 3. PERFORMANCE (repeated spatial CV + block-size sensitivity + overfit gap)
# =============================================================================
message("3. PERFORMANCE")
fr <- c(); oof <- rep(NA_real_, nrow(df)); trr <- c()
for (rep in 1:3){ fo <- folds_of(bid,5,100+rep)
  for (f in 0:4){ tr <- fo!=f; te <- fo==f
    m <- fit_es(df[tr,], bid[tr], best, seed=rep); p <- pr(m, df[te,])
    fr <- c(fr, r2(df[[TARGET]][te], p))
    if (rep==1){ oof[te] <- p; trr <- c(trr, r2(df[[TARGET]][tr], pr(m, df[tr,]))) } } }
sens <- do.call(rbind, lapply(c(max(3,BLOCK-3), BLOCK, BLOCK+3), function(deg){
  bb <- blk(df$longitude, df$latitude, deg); fo <- folds_of(bb,5,7)
  v <- sapply(0:4, function(f){ tr<-fo!=f; te<-fo==f; r2(df[[TARGET]][te], pr(fit_es(df[tr,], bb[tr], best), df[te,])) })
  data.frame(block_deg=deg, r2=round(mean(v),3), sd=round(sd(v),3)) }))
write_csv(sens, file.path(outdir,"blocksize_sensitivity.csv"))

ok_oof <- is.finite(oof); pooled <- r2(df[[TARGET]][ok_oof], oof[ok_oof])

# Two numbers are reported: pooled_oof_r2 (every out-of-fold prediction scored once
# against the global mean) and train_r2 (the overfitting check). 
write_csv(data.frame(pooled_oof_r2=pooled, train_r2=mean(trr),
                     spatial_cv_r2=mean(fr), spatial_cv_sd=sd(fr), spatial_cv_se=sd(fr)/sqrt(length(fr)),
                     fold_r2_p2.5=unname(quantile(fr,.025)), fold_r2_p97.5=unname(quantile(fr,.975)),
                     rmse=sqrt(mean((df[[TARGET]]-oof)^2, na.rm=TRUE))),
          file.path(outdir,"performance_summary.csv"))
message(sprintf("  REPORT: pooled OOF R2=%.3f | train R2=%.3f (gap %.3f)",
                pooled, mean(trr), mean(trr)-pooled))
message(sprintf("  archive: mean-of-folds %.3f (SD %.3f)", mean(fr), sd(fr)))

# =============================================================================
# 4. RESIDUAL DIAGNOSTICS
# =============================================================================
message("4. DIAGNOSTICS")
resid <- df[[TARGET]] - oof
z <- scale(resid[i])[,1]
mIr <- correlogram(z, D, bins)
write_csv(tibble::tibble(dist_km=mids, moranI_residual=mIr),
          file.path(outdir,"residual_correlogram.csv"))
message(sprintf("  residual Moran I (75km) = %.2f -> CIs must come from spatial folds", mIr[1]))

# =============================================================================
# 5. BLOCK-BOOTSTRAP ENSEMBLE -> ALE CIs
# =============================================================================
message("5. ENSEMBLE (ALE CIs)")
set.seed(5)

NB <- 200; ub <- unique(bid)
grids <- lapply(PREDS, function(p) unique(quantile(df[[p]], seq(0,1,length.out=41))))
names(grids) <- PREDS
ale <- setNames(lapply(PREDS, function(p) matrix(NA, length(grids[[p]]), NB)), PREDS)
for (b in 1:NB){
  d <- df[bid %in% sample(ub, length(ub), TRUE), ]
  m <- catboost.train(pool(d, d[[TARGET]]), NULL, params=c(BASE, best, list(iterations=300)))
  for (p in PREDS){ q <- grids[[p]]; bi <- pmin(pmax(findInterval(df[[p]], q[-c(1,length(q))]),0), length(q)-2)
    eff <- numeric(length(q)-1)
    for (k in seq_len(length(q)-1)){ mm <- bi==(k-1)
      if (!any(mm)) next
      dl <- df[mm,]; du <- df[mm,]; dl[[p]] <- q[k]; du[[p]] <- q[k+1]
      eff[k] <- mean(pr(m,du) - pr(m,dl)) }
    a <- c(0, cumsum(eff)); ale[[p]][,b] <- a - mean(a) }
  if (b %% 5 == 0) message("  ", b, "/", NB)
}
ale_tab <- bind_rows(lapply(PREDS, function(p){
  M <- ale[[p]]
  tibble::tibble(predictor = as.character(NAME[[p]]),
                 x   = as.numeric(grids[[p]]),
                 ale = as.numeric(rowMeans(M, na.rm=TRUE)),
                 lo  = as.numeric(apply(M, 1, function(v) quantile(v, 0.025, na.rm=TRUE))),
                 hi  = as.numeric(apply(M, 1, function(v) quantile(v, 0.975, na.rm=TRUE)))) }))
write_csv(ale_tab, file.path(outdir,"ale_curves.csv"))

# =============================================================================
# 6. SPATIAL THINNING CHECK
# =============================================================================
message("6. THINNING")
k2 <- blk(df$longitude, df$latitude, 2)
thin <- sapply(1:10, function(b){ set.seed(b); o <- sample(nrow(df)); ti <- o[!duplicated(k2[o])]
  d <- df[ti,]; bb <- bid[ti]; fo <- folds_of(bb,5,b)
  mean(sapply(unique(fo), function(f){ tr <- fo!=f; te <- fo==f
    if (sum(te)<5 || sum(tr)<20) return(NA)
    r2(d[[TARGET]][te], pr(catboost.train(pool(d[tr,], d[[TARGET]][tr]), NULL,
       params=c(BASE,best,list(iterations=300))), d[te,])) }), na.rm=TRUE) })
message(sprintf("  thinned spatial-CV R2 = %.3f", mean(thin)))
write_csv(data.frame(thinned_cv_r2=mean(thin)), file.path(outdir,"thinning.csv"))

# =============================================================================
# 7. FIGURES
# =============================================================================
message("7. FIGURES")
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
theme_set(theme_minimal(base_size=11) + theme(panel.grid=element_blank(),
          plot.title=element_text(face="bold", size=11), plot.title.position="plot",
          axis.line=element_line(colour="grey30")))
BLUE <- "#2c6fbb"; ok <- is.finite(oof)
ord <- as.character(NAME[c("mammal_richness","consumption","vpi_mean","height","cold_mean","dry_mean")])
write_csv(tibble::tibble(longitude=df$longitude, latitude=df$latitude, observed=df[[TARGET]], oof=oof, resid=resid),
          file.path(outdir,"model1_oof.csv"))

at <- ale_tab; at$predictor <- factor(at$predictor, levels=ord)
ggsave(file.path(outdir,"M1_ALE.png"),
  ggplot(at, aes(x)) + geom_ribbon(aes(ymin=lo,ymax=hi),fill=BLUE,alpha=.25) +
    geom_line(aes(y=ale),colour=BLUE,linewidth=.8) + geom_hline(yintercept=0,colour="grey60",linewidth=.3) +
    facet_wrap(~predictor, scales="free", ncol=2) +
    scale_y_continuous(labels=scales::label_number(scale=1e3)) +
    labs(x=NULL, y=expression("ALE ("%*%10^-3*")"),
         title="Model 1 ALE with 95% block-bootstrap CI",
         subtitle="rows: megafauna, vegetation, seasonality") +
    theme(plot.subtitle=element_text(size=9, colour="grey35"),
          strip.text=element_text(face="bold", size=10)),
  width=7.5, height=9, dpi=300)
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
  labs(x="Observed", y="Predicted (OOF)", title=sprintf("Obs vs pred (pooled OOF R²=%.3f)", pooled))
cal <- dobs; cal$b <- cut(cal$p, quantile(cal$p, seq(0,1,.1)), include.lowest=TRUE); cal <- aggregate(cbind(o,p)~b, cal, mean)
p2 <- ggplot(cal, aes(p,o)) + geom_line(colour=BLUE) + geom_point(colour=BLUE) + geom_abline(linetype=2) + labs(x="Predicted (decile)",y="Observed",title="Calibration")
p3 <- ggplot(tibble::tibble(s=resid[ok]), aes(sample=s)) + stat_qq(size=.4,colour=BLUE) + stat_qq_line(linetype=2) + labs(x="Theoretical",y="Ordered",title="Residual QQ")
cg <- tibble::tibble(dist_km=mids, moranI=mIr)
p4 <- ggplot(cg, aes(dist_km,moranI)) + geom_line(colour="#c44") + geom_point(colour="#c44",size=1) + geom_hline(yintercept=.1,linetype=3,colour="grey") + coord_cartesian(xlim=c(0,2000)) + labs(x="Distance (km)",y="Moran's I",title="Residual autocorrelation")
ggsave(file.path(outdir,"M1_diagnostics.png"), p1+p2+p3+p4+plot_layout(nrow=1), width=19, height=4.3, dpi=170)
message("DONE -> ", outdir)
