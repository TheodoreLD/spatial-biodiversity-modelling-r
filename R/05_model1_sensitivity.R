# =============================================================================
# 05_model1_sensitivity.R — Model 1 design-justification diagnostics (companion to
# 05_model1.R). Three sweeps that justify the modelling choices reviewers questioned:
#   1. RESOLUTION   — spatial-CV R2 vs aggregation (0.1 / 0.5 / 1 / 2 deg) -> why 1 deg
#   2. TRANSFORM    — raw vs log / logit / sqrt response -> why no transform
#   3. SPATIAL TERM — adding lon/lat: it absorbs the ecological signal being measured
# Reads the predictor stack (path in R/config.R). Writes CSVs to outputs/model1/.
# =============================================================================
suppressPackageStartupMessages({ library(terra); library(catboost); library(dplyr); library(readr) })
source("R/config.R")
STACK <- paths$predictor_stack
outdir <- paths$model1; dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
ENV <- c("vpi_mean","cold_mean","dry_mean","height"); HERB <- c("mammal_richness","consumption")
PREDS <- c(ENV,HERB); TARGET <- "trunkrichness"
ps <- load_stack(STACK, c(TARGET, PREDS), "predictor stack")
df0 <- as.data.frame(ps, xy=TRUE, na.rm=TRUE) |> rename(longitude=x, latitude=y)
r2 <- function(y,yh) 1 - sum((y-yh)^2)/sum((y-mean(y))^2)
blk <- function(lon,lat,deg) floor((lon+180)/deg)*100000 + floor((lat+90)/deg)
folds_of <- function(b,k=5,seed=1){ set.seed(seed); u<-unique(b); f<-setNames(sample(0:(k-1),length(u),TRUE),as.character(u)); as.integer(f[as.character(b)]) }
cbp <- list(loss_function="RMSE", depth=3, l2_leaf_reg=30, learning_rate=.05, iterations=400,
            allow_writing_files=FALSE, logging_level="Silent")
fit <- function(d, feats){ catboost.train(catboost.load_pool(as.matrix(d[,feats,drop=FALSE]), label=d[[TARGET]]), NULL, params=cbp) }
prd <- function(m,d,feats) catboost.predict(m, catboost.load_pool(as.matrix(d[,feats,drop=FALSE])))
agg <- function(A){ if (A==1) return(df0)
  df0 |> mutate(R=floor((90-latitude)/(0.1*A)), C=floor((longitude+180)/(0.1*A))) |>
    group_by(R,C) |> summarise(across(all_of(c(TARGET,PREDS,"longitude","latitude")), \(x) mean(x,na.rm=TRUE)), n=n(), .groups="drop") |>
    filter(n >= max(1, floor(0.10*A*A))) |> select(-n) |> as.data.frame() }
cv_r2 <- function(d, feats, deg=7){ b<-blk(d$longitude,d$latitude,deg); fo<-folds_of(b); oof<-rep(NA_real_,nrow(d))
  for (f in 0:4){ tr<-fo!=f; te<-fo==f; if(sum(te)<10) next; oof[te]<-prd(fit(d[tr,],feats),d[te,],feats) }
  ok<-is.finite(oof); r2(d[[TARGET]][ok], oof[ok]) }

# 1. RESOLUTION -----------------------------------------------------------------
message("1. resolution sweep")
resrows <- lapply(c(1,5,10,20), function(A){ d<-agg(A); data.frame(deg=0.1*A, n_cells=nrow(d), spatial_cv_r2=round(cv_r2(d,PREDS),3)) })
resolution <- bind_rows(resrows); write_csv(resolution, file.path(outdir,"resolution_sweep.csv"))
print(resolution)

# 2. TRANSFORM (at 1 deg) -------------------------------------------------------
message("2. transform test")
d1 <- agg(10); y <- d1[[TARGET]]
tf <- list(raw=function(y)y, log=function(y)log(y), logit=function(y)log(y/(1-y)), sqrt=function(y)sqrt(y))
itf <- list(raw=function(z)z, log=function(z)exp(z), logit=function(z)1/(1+exp(-z)), sqrt=function(z)z^2)
trows <- lapply(names(tf), function(nm){ d<-d1; d[[TARGET]]<-tf[[nm]](y); b<-blk(d$longitude,d$latitude,7); fo<-folds_of(b); oof<-rep(NA_real_,nrow(d))
  for (f in 0:4){ tr<-fo!=f; te<-fo==f; if(sum(te)<10) next; oof[te]<-itf[[nm]](prd(fit(d[tr,],PREDS),d[te,],PREDS)) }
  ok<-is.finite(oof); res<-y[ok]-oof[ok]; data.frame(transform=nm, cv_r2_on_raw=round(r2(y[ok],oof[ok]),3), resid_skew=round(mean((res-mean(res))^3)/sd(res)^3,2)) })
transform <- bind_rows(trows); write_csv(transform, file.path(outdir,"transform_test.csv")); print(transform)

# 3. SPATIAL TERM (at 1 deg): what do coordinates do to the ecological signal? -------
# Reported as a trade-off, not as a win for either model: adding lon/lat RAISES CV R2
# (the coordinates are excellent predictors of a spatially structured response) while
# roughly halving the permutation importance of both predictor groups. That is the
# reason Model 1 omits them — they absorb the very signal being characterised — but the
# improvement in fit is real and is reported alongside.
message("3. spatial-term test")
grp_imp <- function(feats){ b<-blk(d1$longitude,d1$latitude,7); fo<-folds_of(b,seed=2); out<-c(megafauna=0,environment=0); nn<-0
  for (f in 0:4){ tr<-fo!=f; te<-fo==f; if(sum(te)<10) next; m<-fit(d1[tr,],feats); base<-r2(d1[[TARGET]][te],prd(m,d1[te,],feats)); nn<-nn+1
    for (g in c("megafauna","environment")){ gg<-if(g=="megafauna")HERB else ENV; tp<-d1[te,]; k<-sample(nrow(tp)); for(v in gg) tp[[v]]<-d1[te,][[v]][k]
      out[g]<-out[g]+(base-r2(d1[[TARGET]][te],prd(m,tp,feats))) } }
  out/nn }
srows <- list(
  data.frame(model="predictors only", t(round(grp_imp(PREDS),3)), cv_r2=round(cv_r2(d1,PREDS),3)),
  data.frame(model="+ lon/lat",       t(round(grp_imp(c(PREDS,"longitude","latitude")),3)), cv_r2=round(cv_r2(d1,c(PREDS,"longitude","latitude")),3)))
spatial <- bind_rows(srows); write_csv(spatial, file.path(outdir,"spatial_term_test.csv")); print(spatial)
message("DONE -> ", outdir, " (resolution_sweep / transform_test / spatial_term_test).csv")
