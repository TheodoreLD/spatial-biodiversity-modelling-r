# =============================================================================
# 05_model1_jackknife.R — species-drop jackknife. Drop a fraction of the 74 mapped species,
# re-derive relative richness (S_subset / W), aggregate to 1 deg, refit CatBoost, and
# record whether the ALE shapes and the spatial-CV R2 stay stable.
# Reads (paths in R/config.R): the predictor stack, the absolute + relative Tot
# richness rasters, and the folder of the 74 per-species .gpkg ranges.
# =============================================================================
suppressPackageStartupMessages({ library(terra); library(catboost); library(dplyr); library(readr); library(ggplot2) })
set.seed(1)
source("R/config.R")
STACK  <- paths$predictor_stack
ABS    <- paths$abs_richness
REL    <- paths$rel_richness
RANGES <- paths$spiny_ranges
outdir <- paths$model1; dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
ENV <- c("vpi_mean","cold_mean","dry_mean","height"); HERB <- c("mammal_richness","consumption")
PREDS <- c(ENV,HERB); TARGET <- "trunkrichness"
NAME <- c(vpi_mean="Productivity (VPI)", cold_mean="Cold-season decid.", dry_mean="Dry-season decid.",
          height="Vegetation height", mammal_richness="Megafauna richness (>5 kg)", consumption="Megafauna consumption")
ALE_ORDER <- c("mammal_richness","consumption","vpi_mean","height","cold_mean","dry_mean")
r2 <- function(y,yh) 1 - sum((y-yh)^2)/sum((y-mean(y))^2)
blk <- function(lon,lat,deg) floor((lon+180)/deg)*100000 + floor((lat+90)/deg)
folds_of <- function(b,k=5,seed=1){ set.seed(seed); u<-unique(b); f<-setNames(sample(0:(k-1),length(u),TRUE),as.character(u)); as.integer(f[as.character(b)]) }
cbp <- list(loss_function="RMSE", depth=3, l2_leaf_reg=30, learning_rate=.05, iterations=400, allow_writing_files=FALSE, logging_level="Silent")

tmpl <- rast(ext(-180,180,-90,90), res=0.1, crs="EPSG:4326")
ps <- load_stack(STACK, c(TARGET, PREDS), "predictor stack")
Av <- values(rast(ABS))[,1]; Rv <- values(rast(REL))[,1]
W <- ifelse(is.finite(Rv) & Rv>0, Av/Rv, NA)                 
Pv <- values(ps[[PREDS]]); predfin <- rowSums(is.finite(Pv))==length(PREDS)

n <- ncell(tmpl); nc <- ncol(tmpl); ev <- as.vector(ext(tmpl)); rs <- res(tmpl)
col0 <- (seq_len(n)-1)%%nc; row0 <- (seq_len(n)-1)%/%nc
lon <- ev[["xmin"]]+(col0+0.5)*rs[1]; lat <- ev[["ymax"]]-(row0+0.5)*rs[2]
AGG <- round(1/rs[1]); Rb <- row0%/%AGG; Cb <- col0%/%AGG
MINSUB <- max(1, floor(0.10*AGG*AGG))

gpkgs <- list.files(RANGES, "\\.gpkg$", full.names=TRUE)
message(length(gpkgs)," species; rasterising presences")
pres <- lapply(gpkgs, function(g){ r <- rasterize(vect(g), tmpl, background=0, touches=TRUE); which(values(r)[,1]>0) })
nsp <- length(pres)

build <- function(keep){ S <- integer(n); for (k in keep) S[pres[[k]]] <- S[pres[[k]]]+1L
  m <- S>0 & is.finite(W) & predfin; rel <- rep(NA_real_,n); rel[m] <- S[m]/W[m]
  d <- data.frame(trunkrichness=rel[m], lon=lon[m], lat=lat[m], R=Rb[m], C=Cb[m]); d <- cbind(d, Pv[m,,drop=FALSE])
  d |> group_by(R,C) |> summarise(across(everything(), \(x) mean(x,na.rm=TRUE)), nn=n(), .groups="drop") |>
    filter(nn>=MINSUB) |> select(-nn) |> as.data.frame() }
fit <- function(d) catboost.train(catboost.load_pool(as.matrix(d[,PREDS]), label=d[[TARGET]]), NULL, params=cbp)
ale1 <- function(m, X, j, q){ bi <- pmin(pmax(findInterval(X[,j], q[-c(1,length(q))]),0), length(q)-2); eff <- numeric(length(q)-1)
  for (k in seq_len(length(q)-1)){ mm <- bi==(k-1); if(!any(mm)) next; Xl<-X[mm,,drop=FALSE]; Xu<-X[mm,,drop=FALSE]; Xl[,j]<-q[k]; Xu[,j]<-q[k+1]
    eff[k] <- mean(catboost.predict(m, catboost.load_pool(Xu)) - catboost.predict(m, catboost.load_pool(Xl))) }
  a <- c(0,cumsum(eff)); a-mean(a) }
cvr2 <- function(d){ b<-blk(d$lon,d$lat,7); fo<-folds_of(b); oof<-rep(NA_real_,nrow(d))
  for (f in 0:4){ tr<-fo!=f; te<-fo==f; if(sum(te)<10) next; oof[te]<-catboost.predict(fit(d[tr,]), catboost.load_pool(as.matrix(d[te,PREDS]))) }
  ok<-is.finite(oof); r2(d[[TARGET]][ok], oof[ok]) }

full <- build(seq_len(nsp)); Xf <- as.matrix(full[,PREDS])
grids <- lapply(PREDS, function(p) unique(quantile(Xf[,match(p,PREDS)], seq(0,1,length.out=31)))); names(grids)<-PREDS
mfull <- fit(full); ale_full <- lapply(PREDS, function(p) ale1(mfull, Xf, match(p,PREDS), grids[[p]])); names(ale_full)<-PREDS
message(sprintf("FULL (%d sp): n=%d cv_r2=%.3f", nsp, nrow(full), cvr2(full)))

rows <- list(); alest <- setNames(lapply(PREDS, function(p) list()), PREDS)
for (fr in c(.1,.2,.3)){ for (rep in 1:20){ set.seed(as.integer(1000*fr*100+rep))
    keep <- sort(sample(nsp, round((1-fr)*nsp))); d <- build(keep); X <- as.matrix(d[,PREDS])
    rows[[length(rows)+1]] <- data.frame(frac=fr, rep=rep, n_species=length(keep), n_cells=nrow(d), cv_r2=cvr2(d))
    if (fr==.3){ m <- fit(d); for (p in PREDS) alest[[p]][[length(alest[[p]])+1]] <- ale1(m, X, match(p,PREDS), grids[[p]]) } }
  message(sprintf("  drop %d%% done", as.integer(fr*100))) }
res <- bind_rows(rows); write_csv(res, file.path(outdir,"jackknife.csv"))

dd <- bind_rows(lapply(ALE_ORDER, function(p){ M <- do.call(cbind, alest[[p]])
  data.frame(predictor=NAME[[p]], x=grids[[p]], full=ale_full[[p]],
             lo=apply(M,1,function(v)quantile(v,.025)), hi=apply(M,1,function(v)quantile(v,.975))) }))
dd$predictor <- factor(dd$predictor, levels=NAME[ALE_ORDER])
ggsave(file.path(outdir,"M1_jackknife_ALE.png"),
  ggplot(dd, aes(x)) + geom_ribbon(aes(ymin=lo,ymax=hi), fill="#e8a13a", alpha=.35) +
    geom_line(aes(y=full), colour="#2c6fbb", linewidth=.8) + geom_hline(yintercept=0, colour="grey60") +
    facet_wrap(~predictor, scales="free", ncol=3) + labs(x=NULL, y="ALE",
    title="Model 1 species-drop jackknife: ALE stability when 30% of species are removed (blue = all 74 species)") +
    theme_minimal(base_size=10) + theme(plot.title=element_text(face="bold",size=11), axis.line=element_line(colour="grey30")),
  width=13, height=7, dpi=175)
message("DONE -> ", outdir)
