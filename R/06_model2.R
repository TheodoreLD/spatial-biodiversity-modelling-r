# =============================================================================
# 06_model2.R - Model 2: are spiny-trunk trees relatively richer where particular
# large-herbivore clades are richer?
#
# Per condition x continent x clade, four fits: clade alone and clade + 4 environmental
# covariates, each unconstrained and monotone-constrained. Plus one environment-only model per
# condition x continent, which carries no clade and so serves every clade in that stratum.
#   * unconstrained      -> the shape of a response and its direction 
#   * positive-monotonic -> the reported range, and the R2 that supports it
#
#
# Where each quantity lives:
#   model2_results.csv      clade alone -- range_std (+_lo/_hi), r2_clade_unc, r2_clade_mono
#   model2_env_results.csv  clade + environment -- range_env, r2_cladeenv_*, r2_env, and the
#                           decomposition unique_clade_mono / unique_env_mono / shared_mono.
#                           The clade-alone columns are joined in here under the names
#                           range_uni and r2_clade_*, so the two settings sit side by side.
#   model2_curves.csv                  clade alone, unconstrained
#   model2_curves_constrained.csv      clade alone, monotone
#   model2_curves_env.csv              clade within clade + environment, unconstrained
#   model2_curves_env_constrained.csv  clade within clade + environment, monotone
#   model2_jackknife_*.csv             species-drop robustness (section 6a)
#
#
# Sections
#   1  SETTINGS     every number that changes the analysis, in one place
#   2  APPEARANCE   colours, clade order, theme (used only by section 7)
#   3  HELPERS      continent rules, spatial blocks, ALE, the ranking metric
#   4  INPUTS       rasters -> one 1-degree table per condition (+ species ranges)
#   5  MODELS       clade alone, then clade + environment: CV estimates, bootstrap intervals
#   6  ROBUSTNESS   species-drop jackknife -- resamples the species pool, not space
#   7  FIGURES      shapes, R2, the range, the decomposition, robustness, coverage,
#                   and the manuscript composite; reads back the CSVs 5-6 write
#
#
# Inputs (paths in R/config.R): the predictor stack (environment layers), the per-clade
# richness stack, the Biogeography rasters folder (Relative_Richness_*.tif, absolute
# richness), and the per-species .gpkg ranges (section 6a only).
# Requires 04_data_preparation.R to have been run. Nothing here reads Model 1.
#
# =============================================================================
suppressPackageStartupMessages({ library(terra); library(catboost); library(dplyr)
  library(readr); library(ggplot2); library(patchwork) })
for (.p in c("ggrepel", "scales", "tidyr", "tibble"))
  if (!requireNamespace(.p, quietly = TRUE)) stop("package '", .p, "' is required")
source("R/config.R")
set.seed(0)
# =============================================================================
# 1. SETTINGS
# =============================================================================
PRED  <- paths$predictor_stack
CLADE <- paths$clade_stack
RELD  <- paths$biogeo_rasters
outdir <- paths$model2; dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

# ---- what is analysed --------------------------------------------------------
CLADES <- c("Pilosa","Cetartiodactyla","Diprotodontia","Primates","Rodentia",
            "Perissodactyla","Proboscidea","Cingulata","Litopterna","Notoungulata")

DISP <- function(x) ifelse(x=="Cetartiodactyla", "Artiodactyla", x)
ENV   <- c("vpi_mean","cold_mean","dry_mean","height")
CONTS <- c("Americas","Africa","Asia")
rel_files <- c(total="Relative_Richness_Tot.tif",
               prickly="Relative_Richness_Prickly.tif",
               thorny="Relative_Richness_Thorny.tif")

# ---- resampling --------------------------------------------------------------
# Spatial block size, in degrees
KFOLD <- 5    # folds per cross-validation assignment
NSEED <- 10   # independent fold assignments, averaged for the point estimates
BB    <- 250  # block-bootstrap resamples behind every confidence interval
NGRID <- 50   # One quantile grid 

# ---- robustness (section 6) ---------------------------------------------------
FRACS <- 0.3   # drop 30% of the mapped species per replicate
NREP  <- 90    # replicates per fraction

# ---- the learner --------------------------------------------------------------
# random_seed and thread_count are pinned
MONO1 <- "(1)"  # clade alone: non-decreasing in the clade
# clade + env: non-decreasing in the clade, environment unconstrained. 
MONO5 <- sprintf("(1,%s)", paste(rep(0, length(ENV)), collapse=","))
cbp <- list(loss_function="RMSE", depth=4, l2_leaf_reg=10, learning_rate=0.05,
            iterations=300, allow_writing_files=FALSE, logging_level="Silent",
            random_seed=1, thread_count=4)
# =============================================================================
# 2. APPEARANCE
# =============================================================================

COL   <- c(Americas="#E1812C", Africa="#3A7CA5", Asia="#4C9A4C")
ED    <- c("Pilosa","Cingulata","Litopterna","Notoungulata","Diprotodontia")  
ORDER <- c("Artiodactyla","Perissodactyla","Rodentia","Primates","Proboscidea",
           "Pilosa","Cingulata","Litopterna","Notoungulata","Diprotodontia")
lab_ed <- function(c) ifelse(c %in% ED, paste0(c," (ED)"), c)
theme_set(theme_minimal(base_size=10) +
  theme(panel.grid=element_blank(), plot.title=element_text(face="bold", size=11),
        plot.title.position="plot", axis.line=element_line(colour="grey30"),
        strip.text=element_text(face="bold")))
# =============================================================================
# 3. HELPERS
# =============================================================================
# ---- 3a. scoring and pools ---------------------------------------------------
# Out-of-fold R2, and the two CatBoost pool constructors (one column / a matrix).

Pm <- function(M) catboost.load_pool(as.matrix(M))

# ---- 3b. continents ----------------------------------------------------------

afro_asia_lon <- function(la)
  ifelse(la >= 31.5, 25,                                    
  ifelse(la >= 12.5, 43 - (la - 12.5) * (43 - 32) / 19,     
                     60))                                    
cont_of <- function(lo,la){ out <- rep("Other", length(lo))
  out[lo>=-170 & lo<=-30 & la>=-60 & la<=80] <- "Americas"
  out[lo>= 25  & lo<=180 & la>=-12 & la<=80] <- "Asia"
  out[lo>=112  & lo<=154 & la>=-44 & la<=-10] <- "Other"     
  out[lo>=128  & la>=-12 & la<=0] <- "Other"
  out[lo>=-20  & lo<= 55 & la>=-35 & la<=40 & lo < afro_asia_lon(la)] <- "Africa"
  out }

# WWF biogeographic realm
realm_at <- function(lo, la, eco)
  terra::extract(eco[, "REALM"], terra::vect(data.frame(lon=lo, lat=la),
                 geom=c("lon","lat"), crs=terra::crs(eco)))$REALM

assign_continent <- function(lo, la, eco) {
  out <- cont_of(lo, la)
  r <- realm_at(lo, la, eco)
  out[!is.na(r) & r %in% c("Australasia","Oceania")] <- "Other"
  out }

# Deal whole blocks into k folds; needed both to size the blocks (4d)
# and to score the models (5).
foldsf <- function(bid, k = KFOLD, seed = 1){ set.seed(seed)
  u <- sample(unique(bid))
  f <- setNames(rep(seq_len(k), length.out = length(u)), as.character(u))
  unname(f[as.character(bid)]) }

# ---- 3c. spatial blocks ------------------------------------------------------
# A deg x deg cell of longitude/latitude, encoded as one integer.
blocks <- function(sub, deg)
  floor((sub$lon+180)/deg)*100000 + floor((sub$lat+90)/deg)
# Great-circle distance, and the residual correlogram Model 1 uses to size its blocks.
hav <- function(lo1, la1, lo2, la2){
  a <- sin((la2-la1)*pi/360)^2 + cos(la1*pi/180)*cos(la2*pi/180)*sin((lo2-lo1)*pi/360)^2
  2*6371*asin(pmin(1, sqrt(a))) }
BINS <- seq(0, 3000, 150); MIDS <- (BINS[-1] + BINS[-length(BINS)])/2
correlogram <- function(z, D){ Pp <- outer(z, z)
  sapply(seq_len(length(BINS)-1), function(k){
    mk <- D >= BINS[k] & D < BINS[k+1] & D > 0
    if (!any(mk)) NA_real_ else mean(Pp[mk], na.rm = TRUE) }) }
# first bin midpoint at which the correlogram drops below 0.1, as in Model 1
cross_km <- function(v){ h <- which(v < 0.1)[1]; if (is.na(h)) max(MIDS) else MIDS[h] }
# km -> degrees of longitude at the 90th percentile of |latitude|, rounded up
deg_from_km <- function(km, lat) max(5, ceiling(km / (111.32*cos(
  as.numeric(quantile(abs(lat), .90))*pi/180))))

# ---- 3d. accumulated local effects and the ranking metric --------------------

ale_at <- function(m, X, q){
  bi <- pmin(pmax(findInterval(X[,1], q[-c(1,length(q))]),0), length(q)-2); eff <- numeric(length(q)-1)
  for (k in seq_len(length(q)-1)){ mm <- bi==(k-1); if(!any(mm)) next
    Xl <- X[mm,,drop=FALSE]; Xu <- X[mm,,drop=FALSE]; Xl[,1] <- q[k]; Xu[,1] <- q[k+1]
    eff[k] <- mean(catboost.predict(m, Pm(Xu)) - catboost.predict(m, Pm(Xl))) }
  a <- c(0, cumsum(eff)); a - mean(a) }

# The response reduced to one number: the total change across the clade's richness range,
# in SDs of relative spiny-trunk richness. 
resp_range <- function(yv, sdy){
  if (!is.finite(sdy) || sdy <= 0 || length(yv) < 2) return(NA_real_)
  rg <- max(yv) - min(yv)
  if (!is.finite(rg) || rg <= 0) return(0)
  rg / sdy }

# =============================================================================
# 4. INPUTS
# =============================================================================
# ---- 4a. rasters and the 1-degree geometry -----------------------------------
# Open the two stacks and derive lon/lat plus the 0.1 -> 1 degree aggregation index.
message("build 1-degree tables")
ps <- load_stack(PRED, ENV, "predictor stack (environment layers)")
cs <- rast(require_file(CLADE, "clade richness stack"))
missing_clades <- setdiff(CLADES, names(cs))
if (length(missing_clades))
  stop("Clade layers absent from ", CLADE, ": ", paste(missing_clades, collapse=", "),
       "\n  (04_data_preparation.R builds one layer per PHYLACINE order passing the >5 kg",
       " herbivore filter; CLADES above lists the ten analysed here.)", call.=FALSE)

terra::compareGeom(ps, cs,
                   rast(require_file(paths$abs_richness, "absolute richness")),
                   rast(require_file(paths$rel_richness, "relative richness")),
                   stopOnError=TRUE)
for (.rf in unname(rel_files))         
  terra::compareGeom(ps, rast(file.path(RELD, .rf)), stopOnError=TRUE)
env_v <- values(ps); envfin <- rowSums(is.finite(env_v))==length(ENV)

# -- clade richness
cl_v <- values(cs[[CLADES]]); cl_v[is.na(cl_v)] <- 0
eco <- vect(require_file(paths$ecoregions, "WWF ecoregions (realm lookup)"))

# -- cell -> lon/lat and the 0.1 -> 1 degree index
n <- ncell(ps); nc <- ncol(ps); ev <- as.vector(ext(ps)); rs <- res(ps)
col0 <- (seq_len(n)-1) %% nc; row0 <- (seq_len(n)-1) %/% nc
lon <- ev[["xmin"]] + (col0+0.5)*rs[1]; lat <- ev[["ymax"]] - (row0+0.5)*rs[2]
AGG <- round(1/rs[1])   # source cells per degree (10 at 0.1 deg)

stopifnot(AGG >= 1, isTRUE(all.equal(rs[1], rs[2])),
          isTRUE(all.equal(AGG * rs[1], 1)))
MINSUB <- max(1, floor(0.10*AGG*AGG))      # same 10% coverage rule as Model 1
Rb <- row0 %/% AGG; Cb <- col0 %/% AGG

# ---- 4b. one 1-degree table per richness condition ----------------------------
# The 0.1 deg RATIO is averaged to 1 deg, rather than the ratio being recomputed from
# aggregated numerator and denominator.
build_condition <- function(rf){
  rel <- values(rast(file.path(RELD, rf)))[,1]
  m <- is.finite(rel) & rel>0 & envfin
  d <- data.frame(rel=rel[m], R=Rb[m], C=Cb[m], lon=lon[m], lat=lat[m])
  d <- cbind(d, env_v[m,,drop=FALSE]); d <- cbind(d, cl_v[m,,drop=FALSE])
  g <- d |> group_by(R,C) |> summarise(across(everything(), \(z) mean(z, na.rm=TRUE)), n=n(), .groups="drop") |>
    filter(n>=MINSUB) |> as.data.frame()
  g$continent <- assign_continent(g$lon, g$lat, eco); g
}
tables <- lapply(rel_files, build_condition); names(tables) <- names(rel_files)
for (cond in names(tables)) message(sprintf("  %s: %d cells | %s", cond, nrow(tables[[cond]]),
  paste(names(table(tables[[cond]]$continent)), table(tables[[cond]]$continent), collapse=" ")))

# ---- 4c. per-species ranges (section 6a only) --------------------------------

Av <- values(rast(require_file(paths$abs_richness, "absolute richness")))[,1]
Rv <- values(rast(require_file(paths$rel_richness, "relative richness")))[,1]
W  <- ifelse(is.finite(Rv) & Rv > 0, Av/Rv, NA)
tmpl  <- rast(ps)   # same geometry as everything else, rather than assuming a square res
gpkgs <- list.files(paths$spiny_ranges, "\\.gpkg$", full.names=TRUE)
message(length(gpkgs), " species ranges; rasterising presences")
pres <- lapply(gpkgs, function(g)
  which(values(rasterize(vect(g), tmpl, background=0, touches=TRUE))[,1] > 0))
nsp <- length(pres)

# =============================================================================
# ---- 4d. spatial block size, measured per continent --------------------------
message("block size: residual autocorrelation per continent")
NBLK <- 5   # fold allocations averaged for the correlogram
.tot <- tables[["total"]]
AC <- bind_rows(lapply(CONTS, function(ct){
  sub <- .tot[.tot$continent == ct, ]
  b5 <- blocks(sub, 5); X <- as.matrix(sub[, ENV]); y <- sub$rel
  M <- sapply(seq_len(NBLK), function(s){
    fo <- foldsf(b5, KFOLD, s); oof <- rep(NA_real_, nrow(sub))
    for (f in unique(fo)){ tr <- fo != f; te <- fo == f
      if (sum(te) < 5 || sum(tr) < 20) next
      mo <- catboost.train(catboost.load_pool(X[tr, , drop=FALSE], label = y[tr]),
                           NULL, params = cbp)
      oof[te] <- catboost.predict(mo, Pm(X[te, , drop=FALSE])) }
    ok <- is.finite(oof)
    set.seed(1000 + s); idx <- sample(which(ok), min(2500, sum(ok)))
    D <- outer(seq_along(idx), seq_along(idx), function(a, b)
      hav(sub$lon[idx][a], sub$lat[idx][a], sub$lon[idx][b], sub$lat[idx][b]))
    correlogram(scale((y - oof)[idx])[, 1], D) })
  mI <- rowMeans(M, na.rm = TRUE)
  data.frame(continent = ct, n_cells = nrow(sub), range_km = cross_km(mI),
             lat_p90 = as.numeric(quantile(abs(sub$lat), .90)),
             km_per_degree = 111.32*cos(as.numeric(quantile(abs(sub$lat), .90))*pi/180),
             block_deg = deg_from_km(cross_km(mI), sub$lat),
             range_km_per_allocation = paste(sprintf("%.0f", apply(M, 2, cross_km)),
                                             collapse = "; ")) }))
write_csv(AC, file.path(outdir, "TAB_S_model2_autocorrelation.csv"))
BLOCK_DEG <- setNames(AC$block_deg, AC$continent)
for (i in seq_len(nrow(AC)))
  message(sprintf("  %-9s %d cells | range ~%.0f km | %.0f km/deg at |lat| p90 %.1f -> %d deg",
          AC$continent[i], AC$n_cells[i], AC$range_km[i], AC$km_per_degree[i],
          AC$lat_p90[i], AC$block_deg[i]))

# =============================================================================
# 5. MODELS
# =============================================================================


# -- pooled out-of-fold R2 against a fixed baseline mu
cv_r2 <- function(X, y, fo, mu, mono = NULL){
  pr <- if (is.null(mono)) cbp else c(cbp, list(monotone_constraints = mono))
  oof <- rep(NA_real_, length(y))
  for (f in unique(fo)){ tr <- fo != f; te <- fo == f
    if (sum(te) < 10 || sum(tr) < 20) next
    m <- catboost.train(catboost.load_pool(X[tr, , drop = FALSE], label = y[tr]),
                        NULL, params = pr)
    oof[te] <- catboost.predict(m, Pm(X[te, , drop = FALSE])) }
  ok <- is.finite(oof)
  if (sum(ok) < 20) return(NA_real_)
  ss <- sum((y[ok] - mu)^2); if (!is.finite(ss) || ss <= 0) return(NA_real_)
  1 - sum((y[ok] - oof[ok])^2) / ss }

# -- the same estimate under block-bootstrap resamples.
cv_boot <- function(X, y, bid, mu, mono = NULL){
  pr <- if (is.null(mono)) cbp else c(cbp, list(monotone_constraints = mono))
  ublk <- unique(bid); out <- rep(NA_real_, BB)
  for (b in seq_len(BB)){
    set.seed(b)
    samp <- sample(ublk, length(ublk), TRUE)
    idx  <- unlist(lapply(samp, function(u) which(bid == u)), use.names = FALSE)
    uf <- setNames(sample(rep(seq_len(KFOLD), length.out = length(ublk))), as.character(ublk))
    fo <- unname(uf[as.character(bid[idx])])
    Xb <- X[idx, , drop = FALSE]; yb <- y[idx]
    oof <- rep(NA_real_, length(yb))
    for (f in unique(fo)){ tr <- fo != f; te <- fo == f
      if (sum(te) < 10 || sum(tr) < 20) next
      m <- catboost.train(catboost.load_pool(Xb[tr, , drop = FALSE], label = yb[tr]),
                          NULL, params = pr)
      oof[te] <- catboost.predict(m, Pm(Xb[te, , drop = FALSE])) }
    ok <- is.finite(oof)
    if (sum(ok) >= 20){ ss <- sum((yb[ok] - mu)^2)
      if (is.finite(ss) && ss > 0) out[b] <- 1 - sum((yb[ok] - oof[ok])^2)/ss } }
  out }

# -- the response under the same block resamples: 
boot_curve <- function(X, y, bid, g, mono = NULL){
  ublk <- unique(bid); M <- matrix(NA_real_, length(g), 0)
  for (b in seq_len(BB)){
    set.seed(b)
    samp <- sample(ublk, length(ublk), TRUE)
    idx  <- unlist(lapply(samp, function(u) which(bid == u)), use.names = FALSE)
    Xb <- X[idx, , drop = FALSE]
    if (length(unique(Xb[,1])) < 3) next
    m <- catboost.train(catboost.load_pool(Xb, label = y[idx]), NULL,
                        params = if (is.null(mono)) cbp else
                                 c(cbp, list(monotone_constraints = mono)))
    a <- ale_at(m, Xb, g); M <- cbind(M, a - a[1]) }
  M }

# -- mean and 95% percentile interval
mlh <- function(v){ v <- v[is.finite(v)]
  if (!length(v)) return(c(NA_real_, NA_real_, NA_real_))
  unname(c(mean(v), quantile(v, .025), quantile(v, .975))) }

# ---- 5a. clade alone: the response range + R2, and the response curves ---------
message("clade-alone models (CV estimates + bootstrap intervals)")
rows <- list(); curves <- list(); ccur <- list(); cvst <- list()
ri <- 0; ci <- 0; cci <- 0
for (cond in names(tables)){ tab <- tables[[cond]]
  for (cont in CONTS){ sub <- tab[tab$continent==cont,]; if (!nrow(sub)) next
    bid <- blocks(sub, BLOCK_DEG[[cont]]); ublk <- unique(bid)
    y <- sub$rel; mu <- mean(y); sdy_all <- sd(y)
    FO <- lapply(seq_len(NSEED), function(s) foldsf(bid, seed = s))
    for (cl in CLADES){ x <- sub[[cl]]
      if (sum(x>0)==0) next
      cgrid <- unique(quantile(x, seq(0,1,length.out=NGRID)))
      if (length(cgrid) < 3){
        message(sprintf("  skip %s/%s/%s: present in %d cells, only %d distinct quantile levels",
                        cond, cont, DISP(cl), sum(x>0), length(cgrid))); next }
      Xu <- matrix(x, ncol=1)
      r2u <- sapply(FO, function(fo) cv_r2(Xu, y, fo, mu))
      r2m <- sapply(FO, function(fo) cv_r2(Xu, y, fo, mu, MONO1))
      bm  <- cv_boot(Xu, y, bid, mu, MONO1)
      .bm <- mlh(bm)
      mf <- catboost.train(catboost.load_pool(Xu, label=y), NULL,
                           params=c(cbp, list(monotone_constraints=MONO1)))
      av <- ale_at(mf, Xu, cgrid)
      stopifnot(which.min(av) == 1)
      CM <- boot_curve(Xu, y, bid, cgrid, MONO1)
      rb <- mlh(apply(CM, 2, function(z) resp_range(z, sdy_all)))
      cq <- if (ncol(CM)) t(apply(CM, 1, quantile, c(.025,.975), na.rm=TRUE)) else
            matrix(NA_real_, length(cgrid), 2)
      uf  <- catboost.train(catboost.load_pool(Xu, label=y), NULL, params=cbp)
      avu <- ale_at(uf, Xu, cgrid)
      UM  <- boot_curve(Xu, y, bid, cgrid)
      uq  <- if (ncol(UM)) t(apply(UM, 1, quantile, c(.025,.975), na.rm=TRUE)) else
             matrix(NA_real_, length(cgrid), 2)
      ri <- ri+1; rows[[ri]] <- data.frame(condition=cond, continent=cont, clade=DISP(cl),
        n=nrow(sub), block_deg=BLOCK_DEG[[cont]], n_blocks=length(ublk),
        n_present=sum(x>0), sd_response=sdy_all, n_grid=length(cgrid), x_max=max(x),
        range_std=resp_range(av, sdy_all),
        range_std_boot=rb[1], range_std_lo=rb[2], range_std_hi=rb[3],
        n_boot_curve=ncol(CM),
        n_seed_ok=sum(is.finite(r2m)),
        r2_clade_unc=mean(r2u,na.rm=TRUE),  r2_clade_unc_sd=sd(r2u,na.rm=TRUE),
        r2_clade_mono=mean(r2m,na.rm=TRUE), r2_clade_mono_sd=sd(r2m,na.rm=TRUE),
        r2_clade_mono_boot=.bm[1], r2_clade_mono_lo=.bm[2], r2_clade_mono_hi=.bm[3],
        n_boot=sum(is.finite(bm)))
      cvst[[paste(cond,cont,cl)]] <- list(unc=r2u, mono=r2m, boot=bm)
      ci <- ci+1; curves[[ci]] <- data.frame(condition=cond, continent=cont, clade=DISP(cl),
        x=cgrid, y=avu - avu[1], y_lo=uq[,1], y_hi=uq[,2])
      cci <- cci+1; ccur[[cci]] <- data.frame(condition=cond, continent=cont, clade=DISP(cl),
        x=cgrid, y=av - av[1], y_lo=cq[,1], y_hi=cq[,2]) }
    message(sprintf("  %s/%s done", cond, cont)) } }
res <- bind_rows(rows); write_csv(res, file.path(outdir,"model2_results.csv"))
write_csv(bind_rows(curves), file.path(outdir,"model2_curves.csv"))
write_csv(bind_rows(ccur),   file.path(outdir,"model2_curves_constrained.csv"))
message(sprintf("  %d clade-continent-condition models", nrow(res)))

# ---- 5b. the same, with the four environmental covariates
message("environment-adjusted models (CV estimates + bootstrap intervals)")
erows <- list(); ecur <- list(); ecc <- list()
ei <- 0; eci <- 0; ecci <- 0
for (cond in names(tables)){ tab <- tables[[cond]]
  for (cont in CONTS){ sub <- tab[tab$continent==cont,]; if (!nrow(sub)) next
    bid <- blocks(sub, BLOCK_DEG[[cont]]); ublk <- unique(bid)
    y <- sub$rel; mu <- mean(y); sdy_all <- sd(y)
    FO <- lapply(seq_len(NSEED), function(s) foldsf(bid, seed = s))
    XE <- as.matrix(sub[, ENV])
    r2e <- sapply(FO, function(fo) cv_r2(XE, y, fo, mu))
    be  <- cv_boot(XE, y, bid, mu)
    .be <- mlh(be)
    for (cl in CLADES){ xc <- sub[[cl]]
      if (sum(xc>0)==0) next
      X <- as.matrix(sub[,c(cl,ENV)])
      cgrid <- unique(quantile(xc, seq(0,1,length.out=NGRID)))
      if (length(cgrid) < 3) next
      ca <- cvst[[paste(cond,cont,cl)]]; if (is.null(ca)) next
      r2u <- sapply(FO, function(fo) cv_r2(X, y, fo, mu))
      r2m <- sapply(FO, function(fo) cv_r2(X, y, fo, mu, MONO5))
      bj  <- cv_boot(X, y, bid, mu, MONO5)
      .bj <- mlh(bj)
      Suc <- mlh(bj - be); Sue <- mlh(bj - ca$boot); Ssh <- mlh(ca$boot + be - bj)
      .jm <- mean(r2m, na.rm=TRUE); .em <- mean(r2e, na.rm=TRUE); .cm <- mean(ca$mono, na.rm=TRUE)
      mfe <- catboost.train(catboost.load_pool(X, label=y), NULL,
                            params=c(cbp, list(monotone_constraints=MONO5)))
      ave <- ale_at(mfe, X, cgrid)
      stopifnot(which.min(ave) == 1)
      EM  <- boot_curve(X, y, bid, cgrid, MONO5)
      rbe <- mlh(apply(EM, 2, function(z) resp_range(z, sdy_all)))
      eq  <- if (ncol(EM)) t(apply(EM, 1, quantile, c(.025,.975), na.rm=TRUE)) else
             matrix(NA_real_, length(cgrid), 2)
      ufe <- catboost.train(catboost.load_pool(X, label=y), NULL, params=cbp)
      avu <- ale_at(ufe, X, cgrid)
      UME <- boot_curve(X, y, bid, cgrid)
      uqe <- if (ncol(UME)) t(apply(UME, 1, quantile, c(.025,.975), na.rm=TRUE)) else
             matrix(NA_real_, length(cgrid), 2)
      ei <- ei+1; erows[[ei]] <- data.frame(condition=cond, continent=cont, clade=DISP(cl),
        r2_cladeenv_unc=mean(r2u,na.rm=TRUE),  r2_cladeenv_unc_sd=sd(r2u,na.rm=TRUE),
        r2_cladeenv_mono=mean(r2m,na.rm=TRUE), r2_cladeenv_mono_sd=sd(r2m,na.rm=TRUE),
        r2_cladeenv_mono_boot=.bj[1], r2_cladeenv_mono_lo=.bj[2], r2_cladeenv_mono_hi=.bj[3],
        r2_env=mean(r2e,na.rm=TRUE), r2_env_sd=sd(r2e,na.rm=TRUE),
        r2_env_boot=.be[1], r2_env_lo=.be[2], r2_env_hi=.be[3],
        unique_clade_mono = .jm - .em,
        unique_clade_boot = Suc[1], unique_clade_mono_lo = Suc[2], unique_clade_mono_hi = Suc[3],
        unique_env_mono   = .jm - .cm,
        unique_env_boot   = Sue[1], unique_env_mono_lo = Sue[2], unique_env_mono_hi = Sue[3],
        shared_mono       = .cm + .em - .jm,
        shared_boot       = Ssh[1], shared_mono_lo = Ssh[2], shared_mono_hi = Ssh[3],
        range_env=resp_range(ave, sdy_all),
        range_env_boot=rbe[1], range_env_lo=rbe[2], range_env_hi=rbe[3],
        n_boot_curve=ncol(EM), n_seed_ok=sum(is.finite(r2m)),
        n_boot=sum(is.finite(bj)))
      ecci <- ecci+1; ecc[[ecci]] <- data.frame(condition=cond, continent=cont, clade=DISP(cl),
        x=cgrid, y=ave - ave[1], y_lo=eq[,1], y_hi=eq[,2])
      eci <- eci+1; ecur[[eci]] <- data.frame(condition=cond, continent=cont, clade=DISP(cl),
        x=cgrid, y=avu - avu[1], y_lo=uqe[,1], y_hi=uqe[,2]) }
    message(sprintf("  %s/%s done", cond, cont)) } }
write_csv(bind_rows(ecur), file.path(outdir,"model2_curves_env.csv"))
write_csv(bind_rows(ecc),  file.path(outdir,"model2_curves_env_constrained.csv"))
env <- bind_rows(erows) |>
  left_join(res[,c("condition","continent","clade","r2_clade_unc","r2_clade_mono",
                   "r2_clade_mono_boot","r2_clade_mono_lo","r2_clade_mono_hi",
                   "range_std","range_std_lo","range_std_hi")],
            by=c("condition","continent","clade")) |>
  rename(range_uni=range_std, range_uni_lo=range_std_lo, range_uni_hi=range_std_hi)
write_csv(env, file.path(outdir,"model2_env_results.csv"))

# =============================================================================
# 6. ROBUSTNESS
# =============================================================================

message("species-drop jackknife")
build_from_species <- function(keep) {
  S <- integer(n); for (k in keep) S[pres[[k]]] <- S[pres[[k]]] + 1L
  m <- S > 0 & is.finite(W) & envfin
  d <- data.frame(rel = S[m]/W[m], R = Rb[m], C = Cb[m], lon = lon[m], lat = lat[m])
  d <- cbind(d, env_v[m, , drop = FALSE]); d <- cbind(d, cl_v[m, , drop = FALSE])
  g <- d |> group_by(R, C) |>
    summarise(across(everything(), \(z) mean(z, na.rm = TRUE)), nn = n(), .groups = "drop") |>
    filter(nn >= MINSUB) |> select(-nn) |> as.data.frame()
  g$continent <- assign_continent(g$lon, g$lat, eco); g
}

# ---- the range, both settings, for one table ---------------------------------

scores_for <- function(g) {
  out <- list()
  for (cont in CONTS) {
    sub <- g[g$continent == cont, ]; if (!nrow(sub)) next
    y <- sub$rel; sdy <- sd(y); if (!is.finite(sdy) || sdy <= 0) next
    for (cl in CLADES) {
      x <- sub[[cl]]
      if (sum(x > 0) == 0) next          
      q <- unique(quantile(x, seq(0, 1, length.out = NGRID)))
      if (length(q) < 3) next
      mono <- catboost.train(catboost.load_pool(matrix(x, ncol = 1), label = y), NULL,
                             params = c(cbp, list(monotone_constraints = MONO1)))
      Xu <- matrix(x, ncol = 1)         
      au <- ale_at(mono, Xu, q)
      X <- as.matrix(sub[, c(cl, ENV)])
      menv <- catboost.train(catboost.load_pool(X, label = y), NULL,
                             params = c(cbp, list(monotone_constraints = MONO5)))
      ae <- ale_at(menv, X, q)
      cu <- resp_range(au, sdy); ce <- resp_range(ae, sdy)
      if (!all(is.finite(c(cu, ce)))) next
      out[[length(out)+1]] <- data.frame(continent = cont, clade = DISP(cl),
                                         range_uni = cu, range_env = ce)
    }
  }
  bind_rows(out)
}

# -- the full-species reference ranking, used to order the rows

{ .chk <- build_from_species(seq_len(nsp))
  .a <- tables[["total"]][order(tables[["total"]]$R, tables[["total"]]$C),]
  .b <- .chk[order(.chk$R, .chk$C),]
  stopifnot(nrow(.a) == nrow(.b), all(.a$R == .b$R), all(.a$C == .b$C),
            max(abs(.a$rel - .b$rel)) < 1e-9)
  message(sprintf("  full-species reconstruction matches the total table (%d cells)", nrow(.b))) }
full <- scores_for(build_from_species(seq_len(nsp))) |>
  group_by(continent) |> mutate(rank_full = rank(-range_env, ties.method = "min")) |> ungroup()
message(sprintf("full %d-species reference ranking computed", nsp))

# -- the replicates: NREP random species subsets per fraction
rows <- list()
for (fr in FRACS) {
  for (rep in seq_len(NREP)) {
    set.seed(as.integer(1000*fr*100 + rep))
    keep <- sort(sample(nsp, round((1 - fr)*nsp)))
    tb <- build_from_species(keep)
    s  <- scores_for(tb)
    if (!nrow(s)) next
    s <- s |> group_by(continent) |> mutate(rank = rank(-range_env, ties.method = "min")) |> ungroup()
    rows[[length(rows)+1]] <- s |> mutate(frac = fr, rep = rep, n_species = length(keep),
                                          n_cells = nrow(tb))
  }
  message(sprintf("  drop %d%% done (%d reps)", as.integer(fr*100), NREP))
}
jk <- bind_rows(rows)
write_csv(jk, file.path(outdir, "model2_jackknife_ranks.csv"))

# -- how far each clade's range moves when species are dropped
summ <- jk |> group_by(continent, clade, frac) |>
  summarise(n_fits = n(),
            range_env_mean = mean(range_env), range_env_sd = sd(range_env),
            range_uni_mean = mean(range_uni), range_uni_sd = sd(range_uni),
            .groups = "drop") |>
  left_join(full |> select(continent, clade, range_env_full = range_env,
                           range_uni_full = range_uni, rank_full),
            by = c("continent","clade")) |>
  arrange(continent, rank_full)
write_csv(summ, file.path(outdir, "model2_jackknife_summary.csv"))
print(as.data.frame(summ))

# =============================================================================
# 7. FIGURES AND TABLES
# =============================================================================
# Five elements, each a tidy table plus its figures:
#   7a shapes        unconstrained + constrained, both settings, range marked on the latter
#   7b R2            the three models and the clade's unique part, both fits
#   7c scoring       the response range against the R2 that supports it, both settings
#   7d confounding   the decomposition, both fits
#   7e species-drop  does the scoring survive a thinner species pool
#   7f coverage      which clade x continent combinations were estimable, and which were not

CONDS <- names(rel_files)
message("figures and tables")

rd  <- function(f) read_csv(file.path(outdir, f), show_col_types = FALSE)
sav <- function(p, nm, w, h){ ggsave(file.path(outdir, nm), p, width=w, height=h,
                                     dpi=170, bg="white"); message("  ", nm) }
# every join below drops the continent factor's level order, so re-apply it after each one
cf   <- function(d) mutate(d, continent = factor(continent, CONTS))
cl_f <- function(d) mutate(d, lab = factor(lab_ed(clade),
                     levels = lab_ed(intersect(ORDER, unique(d$clade)))))
# one row key per clade x continent, ordered within continent, for the interval plots
keyf <- function(d, by) d |> arrange(continent, .data[[by]]) |>
  mutate(key = factor(paste(as.integer(continent), clade),
                      levels = unique(paste(as.integer(continent), clade))))
ystrip <- list(scale_y_discrete(labels = function(z) sub("^[0-9]+ ", "", z)),
               facet_grid(continent ~ ., scales = "free_y", space = "free_y"),
               labs(y = NULL))
# the reported setting first; "+" sorts before "a", so left alone this order is an accident
SETL  <- c("clade + environment", "clade alone")
setf  <- function(d) mutate(d, setting = factor(setting, SETL))

RES <- rd("model2_results.csv"); ENVR <- rd("model2_env_results.csv")

# ---- 7a. the response shapes -------------------------------------------------

SDR <- RES[, c("condition", "continent", "clade", "sd_response")]
z0 <- function(d) d |> left_join(SDR, by = c("condition","continent","clade")) |>
  mutate(y = y/sd_response, y_lo = y_lo/sd_response, y_hi = y_hi/sd_response) |>
  select(-sd_response)

CURVES <- bind_rows(
  z0(rd("model2_curves.csv"))                 |> mutate(setting="clade alone",     fit="unconstrained"),
  z0(rd("model2_curves_constrained.csv"))     |> mutate(setting="clade alone",     fit="constrained"),
  z0(rd("model2_curves_env.csv"))             |> mutate(setting="clade + environment", fit="unconstrained"),
  z0(rd("model2_curves_env_constrained.csv")) |> mutate(setting="clade + environment", fit="constrained"))
write_csv(CURVES, file.path(outdir, "TAB_M2_shapes.csv"))

# The level the response reaches, and where it first reaches it.
topline <- function(z){ z <- z[order(z$x),]
  ym <- max(z$y)
  if (!is.finite(ym) || ym <= 0) return(data.frame(x0=z$x[1], xe=z$x[1], ym=z$y[1]))
  data.frame(x0=z$x[1], xe=z$x[which(z$y >= ym - 1e-12)[1]], ym=ym) }

fig_shape <- function(cond, set, ft){
  d <- CURVES |> filter(condition == cond, setting == set, fit == ft) |> cf() |> cl_f()
  p <- ggplot(d, aes(x, y, colour = continent)) +
    geom_ribbon(data = function(z) z[is.finite(z$y_hi), ],
                aes(ymin = y_lo, ymax = y_hi, fill = continent), colour = NA, alpha = .15,
                show.legend = FALSE) +
    geom_line(linewidth = .6) +
    scale_colour_manual(values = COL, name = NULL) +
    scale_fill_manual(values = COL, guide = "none") +
    scale_y_continuous(breaks = scales::pretty_breaks(3)) +
    facet_wrap(~ lab, scales = "free", ncol = 5) +
    theme(legend.position = "inside", legend.position.inside = c(.90, .22),
          legend.background = element_blank(), legend.key = element_blank(),
          strip.text = element_text(face = "bold", size = 8.5)) +
    labs(x = "Mammal clade richness",
         y = "Change in relative spiny-trunk richness (SD)",
         title = sprintf("Model 2 \u2014 %s response shape, %s (%s richness)", ft, set, cond))
  if (ft == "constrained"){
    ch <- d |> group_by(continent, lab) |> group_modify(~ topline(.x)) |> ungroup() |> cf()
    p <- p + geom_segment(data = ch, aes(x=x0, xend=xe, y=ym, yend=ym), linewidth = .45,
                          linetype = 2, show.legend = FALSE) +
             geom_point(data = ch, aes(x=xe, y=ym), size = 1.5, show.legend = FALSE) }
  sav(p, sprintf("FIG_M2_shape_%s_%s_%s.png",
                 if (set == "clade alone") "alone" else "env",
                 if (ft == "constrained") "mono" else "unc", cond), 15, 8.5) }

# ---- 7b. how much each model actually predicts -------------------------------
# Four quantities on one axis, all cross-validated the same way. 
R2TAB <- bind_rows(
  RES  |> transmute(condition, continent, clade, quantity = "clade alone", est = r2_clade_mono,
                    lo = r2_clade_mono_lo,
                    hi = r2_clade_mono_hi),
  ENVR |> transmute(condition, continent, clade, quantity = "clade + environment", est = r2_cladeenv_mono,
                    lo = r2_cladeenv_mono_lo,
                    hi = r2_cladeenv_mono_hi),
  ENVR |> transmute(condition, continent, clade, quantity = "environment alone", est = r2_env,
                    lo = r2_env_lo,
                    hi = r2_env_hi),
  ENVR |> transmute(condition, continent, clade, quantity = "clade only (unique)",
                    est = unique_clade_mono,
                    lo = unique_clade_mono_lo,
                    hi = unique_clade_mono_hi))
write_csv(R2TAB, file.path(outdir, "TAB_M2_r2.csv"))

# environment alone is identical for every clade in a continent
QL   <- c("clade alone", "clade + environment", "clade only (unique)")
QCOL <- c(`clade alone`="#7B4FA8", `clade + environment`="#1F78B4",
          `clade only (unique)`="#D94801")
fig_r2 <- function(cond){
  ord <- ENVR |> filter(condition == cond) |> cf() |> keyf("unique_clade_mono") |>
    select(continent, clade, key)
  ref <- ENVR |> filter(condition == cond) |> cf() |>
    group_by(continent) |> summarise(r2_env = first(r2_env), .groups = "drop")
  d <- R2TAB |> filter(condition == cond, quantity != "environment alone") |>
    inner_join(ord, by = c("continent","clade")) |> cf() |>
    mutate(quantity = factor(quantity, QL))
  p <- ggplot(d, aes(est, key, colour = quantity)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
    geom_vline(data = ref, aes(xintercept = r2_env), colour = "grey55", linewidth = .5) +
    geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0, linewidth = .45,
                  alpha = .55, position = position_dodge(width = .6), show.legend = FALSE) +
    geom_point(size = 2.1, position = position_dodge(width = .6)) +
    scale_colour_manual(values = QCOL, name = NULL) +
    ystrip +
    labs(x = "Cross-validated R\u00b2",
         title = sprintf("Model 2 \u2014 predictive skill, constrained fit, %s richness", cond)) +
    theme(legend.position = "bottom")
  sav(p, sprintf("FIG_M2_r2_%s.png", cond), 11.5, 9) }

# ---- 7c. the response range, against the R2 that supports it -----------------
# The range on one axis and the R2 that decides whether to believe it on the other. 
SCORE <- bind_rows(
  RES  |> transmute(condition, continent, clade, setting = "clade alone", range = range_std,
                    range_lo = range_std_lo,
                    range_hi = range_std_hi,
                    r2 = r2_clade_mono,
                    r2_lo = r2_clade_mono_lo,
                    r2_hi = r2_clade_mono_hi),
  ENVR |> transmute(condition, continent, clade, setting = "clade + environment", range = range_env,
                    range_lo = range_env_lo,
                    range_hi = range_env_hi,
                    r2 = unique_clade_mono,
                    r2_lo = unique_clade_mono_lo,
                    r2_hi = unique_clade_mono_hi))
SCORE <- setf(SCORE) |>
  left_join(RES |> select(condition, continent, clade, n_grid, x_max),
            by = c("condition","continent","clade"))
write_csv(SCORE, file.path(outdir, "TAB_M2_scoring.csv"))

fig_score <- function(cond, set){
  d <- SCORE |> filter(condition == cond, setting == set) |> cf()
  xl <- if (set == "clade alone") "Cross-validated R\u00b2  \u2014  clade richness on its own"
        else "Cross-validated R\u00b2 that only the clade supplies, over the environment"
  p <- ggplot(d, aes(r2, range, colour = continent)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
    geom_errorbar(aes(ymin=range_lo, ymax=range_hi), width=0, linewidth=.35, alpha=.3) +
    geom_errorbar(aes(xmin=r2_lo, xmax=r2_hi), orientation="y", width=0,
                  linewidth=.35, alpha=.3) +
    geom_point(aes(size = n_grid), alpha = .85) +
    ggrepel::geom_text_repel(aes(label = clade), size = 2.6, max.overlaps = Inf,
      box.padding = .4, point.padding = .25, min.segment.length = .4, segment.size = .2,
      segment.colour = "grey75", show.legend = FALSE) +
    scale_size_continuous(name = "distinct richness levels", range = c(1.1, 5)) +
    scale_colour_manual(values = COL, guide = "none") +
    facet_wrap(~ continent, ncol = 3, scales = "free_x") +
    labs(x = xl, y = "Response magnitude (SD of relative spiny-trunk richness)",
         title = sprintf("Model 2 \u2014 strength against support, %s (%s richness)", set, cond))
  sav(p, sprintf("FIG_M2_scoring_%s_%s.png",
                 if (set == "clade alone") "alone" else "env", cond), 12.5, 5.8) }

# ---- 7d. what the clade and environment share ------------------------------------
# The three parts sum to R2(clade + environment) within every draw. 
CONF <- bind_rows(lapply(c("unique_clade", "unique_env", "shared"), function(p)
  ENVR |> transmute(condition, continent, clade, part = p,
    est   = .data[[paste0(p, "_mono")]],
    lo    = .data[[paste0(p, "_mono_lo")]],
    hi    = .data[[paste0(p, "_mono_hi")]],
    joint = r2_cladeenv_mono)))
write_csv(CONF, file.path(outdir, "TAB_M2_confounding.csv"))

PL   <- c("unique_clade", "shared", "unique_env")
PCOL <- c(unique_clade="#6A3D9A", shared="#8C8C8C", unique_env="#01665A")
PLAB <- c(unique_clade="only the clade", shared="shared", unique_env="only the environment")
fig_conf <- function(cond){
  ord <- ENVR |> filter(condition == cond) |> cf() |> keyf("shared_mono") |>
    select(continent, clade, key)
  d <- CONF |> filter(condition == cond) |>
    inner_join(ord, by = c("continent","clade")) |> cf() |>
    mutate(part = factor(part, PL))
  p <- ggplot(d, aes(est, key, colour = part)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
    geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0, linewidth = .45,
                  alpha = .55, position = position_dodge(width = .6), show.legend = FALSE) +
    geom_point(size = 2.1, position = position_dodge(width = .6)) +
    scale_colour_manual(values = PCOL, labels = PLAB, name = NULL) +
    ystrip +
    labs(x = "Cross-validated R\u00b2",
         title = sprintf("Model 2 \u2014 clade, environment, and what they share, %s richness", cond)) +
    theme(legend.position = "bottom")
  sav(p, sprintf("FIG_M2_confounding_%s.png", cond), 11.5, 9) }

# ---- 7e. species-drop robustness ---------------------------------------------
# Does the score survive not having mapped every species? Drops 30% of the pool, 90 times,
# and refits. 
JK  <- rd("model2_jackknife_ranks.csv")   |> cf()
JKS <- rd("model2_jackknife_summary.csv") |> cf()
write_csv(JKS, file.path(outdir, "TAB_M2_jackknife.csv"))

jk_long <- bind_rows(
  JK |> transmute(continent, clade, setting = "clade alone",     rng = range_uni),
  JK |> transmute(continent, clade, setting = "clade + environment", rng = range_env)) |> cf() |> setf()
fl_long <- bind_rows(
  JKS |> transmute(continent, clade, setting = "clade alone",     full = range_uni_full),
  JKS |> transmute(continent, clade, setting = "clade + environment", full = range_env_full)) |> cf() |> setf()
ordj <- JKS |> keyf("range_env_full") |> select(continent, clade, key)
jk_long <- jk_long |> inner_join(ordj, by = c("continent","clade")) |> cf()
fl_long <- fl_long |> inner_join(ordj, by = c("continent","clade")) |> cf()
med <- jk_long |> group_by(continent, setting, key) |> summarise(m = median(rng), .groups = "drop")
KEYS <- c("one species-drop replicate", "median of the replicates", "full-pool estimate")
sav(ggplot(mapping = aes(y = key, colour = continent)) +
      geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
      geom_point(data = jk_long, aes(x = rng, shape = KEYS[1]), alpha = .3, size = 1.3,
                 position = position_jitter(height = .16, width = 0)) +
      geom_point(data = med, aes(x = m, shape = KEYS[2]), size = 4, stroke = 1.1) +
      geom_point(data = fl_long, aes(x = full, shape = KEYS[3]), size = 2.5,
                 fill = "white", stroke = .9) +
      scale_shape_manual(values = setNames(c(16, 124, 23), KEYS), breaks = KEYS, name = NULL) +
      scale_colour_manual(values = COL, guide = "none") +
      scale_y_discrete(labels = function(z) sub("^[0-9]+ ", "", z)) +
      guides(shape = guide_legend(override.aes = list(alpha = 1, size = c(1.8, 4, 2.5)))) +
      facet_grid(continent ~ setting, scales = "free_y", space = "free_y") +
      labs(x = "Response magnitude (SD of relative spiny-trunk richness)", y = NULL,
           title = "Model 2 \u2014 does the range survive dropping 30% of the mapped species?") +
      theme(legend.position = "bottom"),
    "FIG_M2_jackknife.png", 12, 8)

# ---- 7f. coverage ------------------------------------------------------------
# Which clade x continent x condition combinations were estimable. 
COV <- expand.grid(condition = CONDS, continent = CONTS, clade = sort(unique(DISP(CLADES))),
                   stringsAsFactors = FALSE) |>
  left_join(RES |> select(condition, continent, clade, n, n_present, n_grid, sd_response, x_max),
            by = c("condition", "continent", "clade")) |>
  mutate(estimated = !is.na(n_grid)) |> arrange(condition, continent, clade)
write_csv(COV, file.path(outdir, "TAB_M2_coverage.csv"))
message(sprintf("  coverage: %d of %d combinations estimable", sum(COV$estimated), nrow(COV)))

# ---- 7f-bis. the printable supplementary table -------------------------------
fmt <- function(e, l, h) ifelse(is.na(e), "", sprintf("%+.3f [%+.3f, %+.3f]", e, l, h))
r2w <- R2TAB |> select(condition, continent, clade, quantity, est, lo, hi) |>
  mutate(v = fmt(est, lo, hi)) |> select(-est, -lo, -hi) |>
  tidyr::pivot_wider(names_from = quantity, values_from = v)
cfw <- CONF |> select(condition, continent, clade, part, est, lo, hi) |>
  mutate(v = fmt(est, lo, hi)) |> select(-est, -lo, -hi) |>
  tidyr::pivot_wider(names_from = part, values_from = v)
scw <- SCORE |> filter(setting == "clade + environment") |>
  transmute(condition, continent, clade,
            magnitude = fmt(range, range_lo, range_hi))
SUP <- COV |>
  left_join(r2w, by = c("condition","continent","clade")) |>
  left_join(cfw, by = c("condition","continent","clade")) |>
  left_join(scw, by = c("condition","continent","clade")) |>
  transmute(condition, continent, clade,
            cells = n, cells_occupied = n_present, richness_levels = n_grid,
            max_richness = x_max, response_sd = signif(sd_response, 3),
            R2_clade_alone   = `clade alone`,
            R2_environment_alone = `environment alone`,
            R2_clade_environment = `clade + environment`,
            contribution     = `clade only (unique)`,
            magnitude_SD     = magnitude,
            unique_environment = unique_env, shared) |>
  arrange(condition, continent, clade) |>
  mutate(across(where(is.character), \(z) ifelse(is.na(z), "", z)),
         status = ifelse(is.na(richness_levels), "not estimable", "estimated"))
write_csv(SUP, file.path(outdir, "TAB_S_model2_summary.csv"))
message(sprintf("  supplementary summary: %d rows (%d estimated)",
                nrow(SUP), sum(SUP$status == "estimated")))

# ---- 7f-ter. sampling design per stratum -------------------------------------
message("sampling design per stratum")
.ps <- load_stack(PRED, ENV, "predictor stack (environment layers)")
.ev <- as.vector(ext(.ps)); .rs <- res(.ps); .nc <- ncol(.ps); .n <- ncell(.ps)
.c0 <- (seq_len(.n)-1) %% .nc; .r0 <- (seq_len(.n)-1) %/% .nc
.lon <- .ev[["xmin"]] + (.c0+0.5)*.rs[1]; .lat <- .ev[["ymax"]] - (.r0+0.5)*.rs[2]
.AGG <- round(1/.rs[1]); .MIN <- max(1, floor(0.10*.AGG*.AGG))
.Rb <- .r0 %/% .AGG; .Cb <- .c0 %/% .AGG
.envfin <- rowSums(is.finite(values(.ps))) == length(ENV)
.eco <- vect(require_file(paths$ecoregions, "WWF ecoregions (realm lookup)"))
.abs_files <- c(total="Absolute_Richness_Tot.tif", prickly="Absolute_Richness_Prickly.tif",
                thorny="Absolute_Richness_Thorny.tif")
STRATA <- bind_rows(lapply(names(rel_files), function(cn){
  rl <- values(rast(file.path(RELD, rel_files[[cn]])))[,1]
  ab <- values(rast(file.path(RELD, .abs_files[[cn]])))[,1]
  m  <- is.finite(rl) & rl > 0 & .envfin
  g <- data.frame(rel=rl[m], abs=ab[m], R=.Rb[m], C=.Cb[m], lon=.lon[m], lat=.lat[m]) |>
    group_by(R, C) |>
    summarise(rel=mean(rel), abs=mean(abs), lon=mean(lon), lat=mean(lat), k=n(), .groups="drop") |>
    filter(k >= .MIN) |> as.data.frame()
  g$continent <- assign_continent(g$lon, g$lat, .eco)
  bind_rows(lapply(CONTS, function(ct){ z <- g[g$continent == ct, ]
    data.frame(condition=cn, continent=ct, cells=nrow(z),
               rel_mean=mean(z$rel), rel_sd=sd(z$rel), rel_max=max(z$rel),
               abs_mean=mean(z$abs), abs_max=max(z$abs)) })) }))
STRATA <- STRATA |>
  left_join(RES |> distinct(condition, continent, n, n_blocks, block_deg),
            by = c("condition","continent"))
# the recomputation must land on the same cells the models were fitted to
stopifnot(all(STRATA$cells == STRATA$n),
          max(abs(STRATA$rel_sd - RES$sd_response[match(paste(STRATA$condition, STRATA$continent),
                paste(RES$condition, RES$continent))])) < 1e-9)
STRATA <- STRATA |> select(condition, continent, cells, n_blocks, block_deg,
                          rel_mean, rel_sd, rel_max, abs_mean, abs_max) |>
  arrange(condition, continent)
write_csv(STRATA, file.path(outdir, "TAB_S_model2_strata.csv"))
message(sprintf("  %d strata, cell counts and response SDs match the fitted models",
                nrow(STRATA)))

# ---- 7f-quater. the settings that fix a Model 2 result -----------------------
.ac <- read_csv(file.path(outdir, "TAB_S_model2_autocorrelation.csv"),
                show_col_types = FALSE)
SETTINGS <- tibble::tibble(
  setting = c("Learner", "Loss", "Tree depth", "L2 regularisation", "Learning rate",
              "Boosting iterations", "Row/feature subsampling", "Hyperparameter tuning",
              "Monotone constraint (clade alone)", "Monotone constraint (clade + environment)",
              "Environmental covariates", "Cross-validation folds",
              "Independent fold assignments", "Block-bootstrap resamples",
              "Spatial block size", "Response grid", "Estimability rule",
              "Species-drop fraction", "Species-drop replicates", "Species-drop conditions"),
  value = c("CatBoost gradient boosting", cbp$loss_function, as.character(cbp$depth),
            as.character(cbp$l2_leaf_reg), as.character(cbp$learning_rate),
            as.character(cbp$iterations), "none",
            "none; settings fixed so scores are comparable across clades",
            MONO1, MONO5, paste(ENV, collapse=", "), as.character(KFOLD),
            as.character(NSEED), as.character(BB),
            paste(sprintf("%s %d deg", .ac$continent, .ac$block_deg), collapse=", "),
            sprintf("up to %d quantiles of clade richness, deduplicated", NGRID),
            "at least three distinct richness values",
            sprintf("%.0f%%", 100*FRACS), as.character(NREP), "total richness only"))
write_csv(SETTINGS, file.path(outdir, "TAB_S_model2_settings.csv"))
message(sprintf("  settings table: %d rows", nrow(SETTINGS)))

# ---- 7g. the manuscript composite --------------------------------------------

composite <- function(cond){
  sc <- SCORE |> filter(condition == cond, setting == "clade + environment") |> cf()
  pA <- ggplot(sc, aes(r2, range, colour = continent)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
    geom_errorbar(aes(ymin=range_lo, ymax=range_hi), width=0, linewidth=.4, alpha=.35) +
    geom_errorbar(aes(xmin=r2_lo, xmax=r2_hi), orientation="y", width=0,
                  linewidth=.4, alpha=.35) +
    geom_point(aes(size = n_grid), alpha=.85) +
    ggrepel::geom_text_repel(aes(label=clade), size=2.7, max.overlaps=Inf, box.padding=.4,
      point.padding=.25, min.segment.length=.4, segment.size=.2, segment.colour="grey75",
      show.legend=FALSE) +
    scale_size_continuous(name="distinct richness\nlevels", range=c(1.1, 5),
                          breaks=c(5, 15, 30, 45)) +
    scale_colour_manual(values=COL, guide="none") +
    facet_wrap(~ continent, ncol=3, scales="free_x") +
    labs(x="Mammal clade contribution (cross-validated R\u00b2)",
         y="Response magnitude (SD of relative spiny-trunk richness)",
         title="A)") +
    theme(legend.position="bottom", legend.key.height=unit(.8, "lines"),
          legend.title=element_text(size=8))


  ord <- ENVR |> filter(condition == cond) |> cf() |> keyf("unique_clade_mono") |>
    select(continent, clade, key)
  bb <- CONF |> filter(condition == cond) |>
    inner_join(ord, by = c("continent", "clade")) |> cf() |>
    mutate(part = factor(part, PL))
  pB <- ggplot(bb, aes(est, key, colour = part)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
    geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0, linewidth = .45,
                  alpha = .55, position = position_dodge(width = .7), show.legend = FALSE) +
    geom_point(size = 2.1, position = position_dodge(width = .7)) +
    scale_colour_manual(values = PCOL, labels = PLAB, name = NULL) +
    ystrip +
    labs(x = "Cross-validated R\u00b2", title = "B)") +
    theme(legend.position = "bottom")

  cc <- CURVES |> filter(condition == cond, setting == "clade + environment",
                         fit == "unconstrained") |> cf() |> cl_f()
  pC <- ggplot(cc, aes(x, y, colour = continent)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey70", linewidth = .3) +
    geom_ribbon(aes(ymin=y_lo, ymax=y_hi, fill=continent), colour=NA, alpha=.15,
                show.legend=FALSE) +
    geom_line(linewidth=.6) +
    scale_fill_manual(values=COL, guide="none") +
    scale_colour_manual(values=COL) +
    scale_y_continuous(breaks=scales::pretty_breaks(3)) +
    facet_wrap(~ lab, scales="free", ncol=5) +
    theme(strip.text=element_text(face="bold", size=8.5),
          legend.position="inside", legend.position.inside=c(0.90, 0.24),
          legend.title=element_blank(), legend.background=element_blank(),
          legend.key=element_blank(), legend.text=element_text(size=10),
          legend.key.size=unit(1.1, "lines")) +
    labs(x="Mammal clade richness",
         y="Change in relative spiny-trunk richness (SD)",
         title="C)")

  sav(((pA | pB) + plot_layout(widths=c(1.3, 1))) / pC + plot_layout(heights=c(1, 1.15)),
      sprintf("FIG_M2_composite%s.png", if (cond == "total") "" else paste0("_", cond)),
      17.5, 12) }

# ---- draw everything ---------------------------------------------------------
FITS <- c("constrained", "unconstrained")
for (cond in CONDS){
  for (set in SETL){ for (ft in FITS) fig_shape(cond, set, ft); fig_score(cond, set) }
  fig_r2(cond); fig_conf(cond); composite(cond) }
message("DONE -> ", outdir)
