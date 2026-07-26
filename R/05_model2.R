# =============================================================================
# 05_model2.R — Model 2 (mammal-clade associations), native CatBoost. Mirrors the
# Python model2 pipeline. Univariate per continent x condition x clade:
#   * monotonic model  -> standardized slope (the clade-ranking metric)
#   * unconstrained    -> response-curve shape + OOB R2
#   * BOTH models' OOB R2 (constrained + unconstrained)
#   * environment-adjusted ALE slope (clade + 4 env covariates)
#   * variance partitioning (unique clade vs shared with environment) + residual Moran
# Honest spatial block-bootstrap CIs (7 deg blocks); R2 out-of-bag.
# Inputs via env vars: TRUNK_STACK (predictor_stack.tif), TRUNK_CLADE
# (Clade_richness_stack.tif), TRUNK_REL_DIR (folder with Relative_Richness_*.tif).
# =============================================================================
suppressPackageStartupMessages({ library(terra); library(catboost); library(dplyr); library(readr) })
set.seed(0)
PRED  <- Sys.getenv("TRUNK_STACK", "predictor_stack.tif")
CLADE <- Sys.getenv("TRUNK_CLADE", "Clade_richness_stack.tif")
RELD  <- Sys.getenv("TRUNK_REL_DIR", ".")
outdir <- file.path("outputs","model2"); dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

CLADES <- c("Pilosa","Cetartiodactyla","Diprotodontia","Primates","Rodentia",
            "Perissodactyla","Proboscidea","Cingulata","Litopterna","Notoungulata")
ENV <- c("vpi_mean","cold_mean","dry_mean","height")
CONTS <- c("Americas","Africa","Asia")
B <- 50
cbp <- list(loss_function="RMSE", depth=4, l2_leaf_reg=10, learning_rate=0.05,
            iterations=300, allow_writing_files=FALSE, logging_level="Silent")
r2 <- function(a,p) 1 - sum((a-p)^2)/sum((a-mean(a))^2)
P1 <- function(x) catboost.load_pool(matrix(x, ncol=1))
Pm <- function(M) catboost.load_pool(as.matrix(M))
cont_of <- function(lo,la){ out <- rep("Other", length(lo))
  out[lo>=-170 & lo<=-30 & la>=-60 & la<=80] <- "Americas"
  out[lo>=-20  & lo<= 55 & la>=-35 & la<=40] <- "Africa"
  out[lo>= 25  & lo<=180 & la>=-12 & la<=80] <- "Asia"; out }
blocks <- function(sub, deg=7) floor((sub$lon+180)/deg)*100000 + floor((sub$lat+90)/deg)

# ---- build 1-degree tables (mirror model2_build.py) -------------------------
message("build 1-degree tables")
ps <- rast(PRED); cs <- rast(CLADE)
env_v <- values(ps[[ENV]]); envfin <- rowSums(is.finite(env_v))==length(ENV)
cl_v <- values(cs[[CLADES]]); cl_v[is.na(cl_v)] <- 0
n <- ncell(ps); col0 <- (seq_len(n)-1) %% 3600; row0 <- (seq_len(n)-1) %/% 3600
lon <- -180 + (col0+0.5)*0.1; lat <- 90 - (row0+0.5)*0.1
Rb <- row0 %/% 10; Cb <- col0 %/% 10
rel_files <- c(total="Relative_Richness_Tot.tif", prickly="Relative_Richness_Prickly.tif", thorny="Relative_Richness_Thorny.tif")
build <- function(rf){
  rel <- values(rast(file.path(RELD, rf)))[,1]
  m <- is.finite(rel) & rel>0 & envfin
  d <- data.frame(rel=rel[m], R=Rb[m], C=Cb[m], lon=lon[m], lat=lat[m])
  d <- cbind(d, env_v[m,,drop=FALSE]); d <- cbind(d, cl_v[m,,drop=FALSE])
  g <- d |> group_by(R,C) |> summarise(across(everything(), \(z) mean(z, na.rm=TRUE)), n=n(), .groups="drop") |>
    filter(n>=10) |> as.data.frame()
  g$continent <- cont_of(g$lon, g$lat); g
}
tables <- lapply(rel_files, build); names(tables) <- names(rel_files)
for (cond in names(tables)) message(sprintf("  %s: %d cells | %s", cond, nrow(tables[[cond]]),
  paste(names(table(tables[[cond]]$continent)), table(tables[[cond]]$continent), collapse=" ")))

std_slope <- function(m, grid, mn, sd){ p <- catboost.predict(m, P1(grid)); z <- (grid-mn)/sd; unname(coef(lm(p ~ z))[2]) }
ale_slope <- function(m, X, mn, sd, nq=21){
  q <- unique(quantile(X[,1], seq(0,1,length.out=nq))); if (length(q)<3) return(NA_real_)
  bi <- pmin(pmax(findInterval(X[,1], q[-c(1,length(q))]),0), length(q)-2); eff <- numeric(length(q)-1)
  for (k in seq_len(length(q)-1)){ mm <- bi==(k-1); if(!any(mm)) next
    Xl <- X[mm,,drop=FALSE]; Xu <- X[mm,,drop=FALSE]; Xl[,1] <- q[k]; Xu[,1] <- q[k+1]
    eff[k] <- mean(catboost.predict(m, Pm(Xu)) - catboost.predict(m, Pm(Xl))) }
  a <- c(0, cumsum(eff)); a <- a - mean(a); z <- (q-mn)/sd; unname(coef(lm(a ~ z))[2]) }

# ---- core models: slope + R2 (mono & unconstrained) + curves -----------------
message("core models (monotonic slope + unconstrained shape + R2)")
rows <- list(); curves <- list(); ccur <- list(); ri <- 0; ci <- 0; cci <- 0
for (cond in names(tables)){ tab <- tables[[cond]]
  for (cont in CONTS){ sub <- tab[tab$continent==cont,]; if (nrow(sub)<80) next
    bid <- blocks(sub); ublk <- unique(bid); y <- sub$rel
    for (cl in CLADES){ x <- sub[[cl]]
      if (sum(x>0)<30 || length(unique(x))<5) next
      mn <- mean(x); sd <- sd(x); grid <- seq(min(x), max(x), length.out=50); zg <- (grid-mn)/sd
      slopes <- c(); r2u <- c(); r2m <- c(); pmat <- matrix(NA_real_,50,0); pmm <- matrix(NA_real_,50,0)
      set.seed(0)
      for (b in 1:B){ samp <- sample(ublk, length(ublk), TRUE); mtr <- bid %in% samp; moob <- !(bid %in% samp)
        if (length(unique(x[mtr]))<3) next
        po <- catboost.load_pool(matrix(x[mtr],ncol=1), label=y[mtr])
        mono <- catboost.train(po, NULL, params=c(cbp, list(monotone_constraints="(1)")))
        slopes <- c(slopes, std_slope(mono, grid, mn, sd)); pmm <- cbind(pmm, catboost.predict(mono, P1(grid)))
        unc <- catboost.train(po, NULL, params=cbp)
        pmat <- cbind(pmat, catboost.predict(unc, P1(grid)))
        if (sum(moob)>=20 && length(unique(y[moob]))>1){ Xo <- P1(x[moob])
          r2u <- c(r2u, r2(y[moob], catboost.predict(unc, Xo)))
          r2m <- c(r2m, r2(y[moob], catboost.predict(mono, Xo))) } }
      if (!length(slopes)) next
      if (!length(r2u)) r2u <- NA_real_; if (!length(r2m)) r2m <- NA_real_
      ri <- ri+1; rows[[ri]] <- data.frame(condition=cond, continent=cont, clade=cl, n=nrow(sub),
        n_present=sum(x>0), slope=mean(slopes), slope_lo=quantile(slopes,.025), slope_hi=quantile(slopes,.975),
        r2=mean(r2u,na.rm=TRUE), r2_lo=quantile(r2u,.025,na.rm=TRUE), r2_hi=quantile(r2u,.975,na.rm=TRUE),
        r2_mono=mean(r2m,na.rm=TRUE), r2_mono_lo=quantile(r2m,.025,na.rm=TRUE), r2_mono_hi=quantile(r2m,.975,na.rm=TRUE))
      ci <- ci+1; curves[[ci]] <- data.frame(condition=cond, continent=cont, clade=cl, x=grid,
        y=rowMeans(pmat), y_lo=apply(pmat,1,function(v) quantile(v,.025)), y_hi=apply(pmat,1,function(v) quantile(v,.975)))
      ym <- rowMeans(pmm); bf <- coef(lm(ym ~ zg))     # constrained curve + standardized linear fit
      cci <- cci+1; ccur[[cci]] <- data.frame(condition=cond, continent=cont, clade=cl, x=grid, y=ym,
        y_lo=apply(pmm,1,function(v) quantile(v,.025)), y_hi=apply(pmm,1,function(v) quantile(v,.975)), y_fit=bf[1]+bf[2]*zg) }
    message(sprintf("  %s/%s done", cond, cont)) } }
res <- bind_rows(rows); write_csv(res, file.path(outdir,"model2_results.csv"))
write_csv(bind_rows(curves), file.path(outdir,"model2_curves.csv"))
write_csv(bind_rows(ccur), file.path(outdir,"model2_curves_constrained.csv"))
message(sprintf("  %d clade-continent-condition models", nrow(res)))

# ---- environment-adjusted ALE slope ------------------------------------------
message("environment-adjusted ALE slope")
mono5 <- "(1,0,0,0,0)"; erows <- list(); ei <- 0
for (cond in names(tables)){ tab <- tables[[cond]]
  for (cont in CONTS){ sub <- tab[tab$continent==cont,]; if (nrow(sub)<80) next
    bid <- blocks(sub); ublk <- unique(bid); y <- sub$rel
    for (cl in CLADES){ xc <- sub[[cl]]; if (sum(xc>0)<30 || length(unique(xc))<5) next
      X <- as.matrix(sub[,c(cl,ENV)]); mn <- mean(xc); sd <- sd(xc); slp <- c(); set.seed(0)
      for (b in 1:B){ m <- bid %in% sample(ublk,length(ublk),TRUE); if (length(unique(X[m,1]))<3) next
        mod <- catboost.train(catboost.load_pool(X[m,,drop=FALSE], label=y[m]), NULL,
                             params=c(cbp, list(monotone_constraints=mono5)))
        s <- ale_slope(mod, X[m,,drop=FALSE], mn, sd); if (is.finite(s)) slp <- c(slp, s) }
      if (!length(slp)) next
      ei <- ei+1; erows[[ei]] <- data.frame(condition=cond, continent=cont, clade=cl,
        slope_env=mean(slp), slope_env_lo=quantile(slp,.025), slope_env_hi=quantile(slp,.975)) }
    message(sprintf("  %s/%s done", cond, cont)) } }
env <- bind_rows(erows) |> left_join(res[,c("condition","continent","clade","slope","slope_lo","slope_hi")],
  by=c("condition","continent","clade")) |> rename(slope_uni=slope, slope_uni_lo=slope_lo, slope_uni_hi=slope_hi)
env$retained <- ifelse(env$slope_uni>1e-9, env$slope_env/env$slope_uni, NA_real_)
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
      di <- di+1; drows[[di]] <- data.frame(condition=cond, continent=cont, clade=cl,
        r2_env=r2e, r2_clade=rc$r2, r2_full=rf$r2, unique_clade=rf$r2-r2e, shared=r2e+rc$r2-rf$r2,
        max_corr_env=max(abs(cor(sub[[cl]], sub[,ENV], method="spearman"))), resid_moran=moran(sub, rf$resid)) }
    message(sprintf("  %s/%s: %d clades", cond, cont, length(act))) } }
write_csv(bind_rows(drows), file.path(outdir,"model2_diagnostics.csv"))
message("DONE -> ", outdir)
