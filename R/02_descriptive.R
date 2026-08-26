# =============================================================================
# 02_descriptive.R — richness maps, biome & ecoregion distributions (block-bootstrap
# CIs), and the Fig-3 composite.
# Inputs (paths from R/config.R): the Biogeography rasters folder (Absolute_/Relative_
# Richness_*.tif), the Ecoregions2017 shapefile, and data/assets/world_borders.geojson.
#
# LAYOUT. Everything is computed before anything is drawn, so the script reads in the
# order the analysis actually depends on:
#   1  settings      every tunable constant, in one place
#   2  appearance    theme, biome palette, label shortening
#   3  helpers       weighted statistics and great-circle distance (pure functions)
#   4  overlay       assign each 0.1 deg cell to its ecoregion, biome and area weight
#   5  block size    within-ecoregion autocorrelation range -> bootstrap block in degrees
#   6  summaries     biome densities + ecoregion means with block-bootstrap CIs
#   7  figures       maps, biome panel, ecoregion panel, Fig-3 composite
#
# DOMAIN. Relative richness is defined only where spiny-trunk species occur: 01 sets an
# absolute richness of 0 to NA before dividing, so a cell holding trees but no spiny-trunk
# species carries no value rather than a zero. Every statistic below is therefore
# conditional on presence -- "where these trees occur, what share of the woody flora are
# they" -- over 403,587 cells for the total condition: 27% of the retained land by cell
# count, 36% by area. The two differ because cells are cos-latitude weighted while the
# unoccupied ones skew polar; the run-time message reports the area figure, which is the
# one consistent with the weighting used everywhere else.
#
# Biomes holding no spiny-trunk species at all therefore cannot appear in the figures
# (Boreal Forests/Taiga and Tundra for every condition, plus Mediterranean for thorny).
# summarise_cond() names them at run time rather than letting them vanish silently.
# =============================================================================
suppressPackageStartupMessages({ library(terra); library(sf); library(dplyr); library(readr)
  library(ggplot2); library(patchwork); library(scales) })
source("R/config.R")
create_project_dirs()
set.seed(0)
RELD   <- paths$biogeo_rasters   # folder with Relative_Richness_*.tif
ECOSHP <- paths$ecoregions       # Ecoregions2017.shp
BORD   <- paths$borders          # world_borders.geojson
outdir <- paths$descriptive


# =============================================================================
# 1. SETTINGS — every knob the analysis exposes
# =============================================================================

# The three richness conditions, used by the correlogram, the summaries and the figures.
RELF <- c(total   = "Relative_Richness_Tot.tif",
          prickly = "Relative_Richness_Prickly.tif",
          thorny  = "Relative_Richness_Thorny.tif")
CLAB <- c(total = "spiny", prickly = "prickly", thorny = "thorny")   # for figure titles

# Biomes that are not vegetated land, dropped everywhere.
EXCLUDE <- c("N/A","Rock and Ice","Rock & Ice","Inland Water","Lakes","Lake")

# --- ecoregion bootstrap ---
NB_ECO        <- 2000   # Resample count for the ecoregion bootstrap.
MIN_EFF_CI    <- 3      # fewer effective blocks than this and no interval is reported
MIN_ECO_CELLS <- 10     # an ecoregion below this many cells is not summarised at all

# --- biome density curves (section 6 computes them, section 7 draws them) ---
# BW_FRAC is a fraction of each biome's OWN robust SD, not bw.nrd0: that rule's n^(-1/5)
# term would smooth the 155,544-cell biome about six times harder than the 421-cell one,
# so apparent smoothness would track cell count rather than spread. Under bw.nrd0 all
# twelve biomes look multimodal and none of those modes survives wider smoothing, so
# structure finer than this width is not interpretable.
BW_FRAC  <- 0.30  # how much the biome density curves are smoothed: 30% of that biome's own spread
RIDGE_H  <- 0.92  # how tall each density curve is drawn

# --- within-ecoregion correlogram, which the block size is derived from ---
CORR_MIN_CELLS <- 200               # ecoregions smaller than this cannot support a curve
CORR_PER_UNIT  <- 150               # cells sampled per ecoregion (cost is quadratic in this)
CORR_BINS      <- seq(0, 1200, 50)  # distance bins, km
CORR_THRESH    <- 0.1               # correlation below which cells count as independent
SUMM_MAX       <- 475               # longest distance shown in the correlogram summary table


# =============================================================================
# 2. APPEARANCE — theme, biome palette, label shortening
# =============================================================================
theme_set(theme_minimal(base_size=10) + theme(panel.grid=element_blank(),
          plot.title=element_text(face="bold",size=11), plot.title.position="plot",
          axis.line=element_line(colour="grey30")))

BCOL <- c(
  "Tropical & Subtropical Moist Broadleaf Forests"           = "#2e8b57",
  "Tropical & Subtropical Dry Broadleaf Forests"             = "#c0392b",
  "Tropical & Subtropical Coniferous Forests"                = "#5fa8a0",
  "Temperate Broadleaf & Mixed Forests"                      = "#7fb069",
  "Temperate Conifer Forests"                                = "#35748c",
  "Boreal Forests/Taiga"                                     = "#4a6f8a",
  "Tropical & Subtropical Grasslands, Savannas & Shrublands" = "#8c6d31",
  "Temperate Grasslands, Savannas & Shrublands"              = "#c98c5a",
  "Flooded Grasslands & Savannas"                            = "#7b3fa0",
  "Montane Grasslands & Shrublands"                          = "#d1738f",
  "Tundra"                                                   = "#9db4c0",
  "Mediterranean Forests, Woodlands & Scrub"                 = "#e08a3c",
  "Deserts & Xeric Shrublands"                               = "#d4b106",
  "Mangroves"                                                = "#2c6fbb")
GREY <- "#9e9e9e"
bcol <- function(b) ifelse(b %in% names(BCOL), BCOL[b], GREY)

# Biomes absent from the top-25 ecoregion panel are greyed, so panels B and C read together.
bcol_keyed <- function(b, keyed) ifelse(b %in% keyed, bcol(b), GREY)
blab <- function(b) gsub("Xeric Shrublands","Xeric Shrub.", gsub("Grasslands, Savannas & Shrublands","Grasslands & Savannas",
        gsub("Tropical & Subtropical","Trop./Subtrop.", b)))


# =============================================================================
# 3. HELPERS — pure functions, no side effects
# =============================================================================

# randomly resamples all elements of x with replacement, preserving duplicates for bootstrap analysis.
boot_resample <- function(x) x[sample.int(length(x), length(x), replace=TRUE)]

# calculates quantiles while accounting for unequal raster-cell areas
# preventing smaller high-latitude cells from being overrepresented
wquantile <- function(x, w, p){ o <- order(x); x <- x[o]; w <- w[o]
  x[which(cumsum(w)/sum(w) >= p)[1]] }

# estimates how many equally weighted independent units the unequal weights effectively represent.
# converts unequal block or ecoregion weights into an “equivalent” number of equally sized units,
# so one very large unit is not misleadingly counted as many independent units when judging whether there is enough information for a confidence interval
kish <- function(W) sum(W)^2 / sum(W^2) # Kish effective sample size

# repeatedly resamples entire spatial blocks and recalculates the area-weighted mean, producing
# a distribution used to estimate uncertainty and confidence intervals while respecting spatial dependence
boot_wmean <- function(v, w, b, nb=1000, seed=0){ set.seed(seed)
  svw <- tapply(v*w, b, sum); sw <- tapply(w, b, sum); ub <- names(sw); o <- numeric(nb)
  for (i in 1:nb){ s <- boot_resample(ub); o[i] <- sum(svw[s])/sum(sw[s]) }
  o }

# great-circle distance in kilometers for coordinates (Haversine formula)
# coordinates into distance km to measure the range of spatial autocorrelation
hav <- function(a,b,c,e){ p <- pi/180; dl <- (c-a)*p; dp <- (e-b)*p
  6371*2*asin(pmin(1, sqrt(sin(dp/2)^2 + cos(b*p)*cos(e*p)*sin(dl/2)^2))) }


# =============================================================================
# 4. OVERLAY — assign each 0.1° cell to its ecoregion, biome and area weight
# =============================================================================
message("overlay (rasterize ecoregions)")
# prepare ecoregion rasters and lookup table
eco <- st_read(ECOSHP, quiet=TRUE) |> st_make_valid() |> st_transform(4326)
tmpl <- rast(ext(-180,180,-90,90), resolution=0.1, crs="EPSG:4326")
eco_id <- rasterize(vect(eco), tmpl, field="ECO_ID")
lut <- eco |> st_drop_geometry() |> distinct(ECO_ID, ECO_NAME, BIOME_NAME)
if (anyDuplicated(lut$ECO_ID))
  stop("ECO_ID maps to more than one ECO_NAME/BIOME_NAME in ", ECOSHP, call.=FALSE)
# store cell coordinates, block IDs, ecoregion IDs, and area weights
n <- ncell(tmpl); nc <- ncol(tmpl); ev <- as.vector(ext(tmpl)); rs <- res(tmpl)
col0 <- (seq_len(n)-1) %% nc; row0 <- (seq_len(n)-1) %/% nc
lon <- ev[["xmin"]] + (col0+0.5)*rs[1]; lat <- ev[["ymax"]] - (row0+0.5)*rs[2] # calculates the centre of every raster cell
blk_at <- function(deg){ ix <- floor((lon+180)/deg); iy <- floor((lat+90)/deg); ix*(max(iy)+1) + iy } # group cells into blocks
eid <- values(eco_id)[,1]
message(sprintf("  %d of %d ecoregions captured a 0.1 deg cell centre",
                sum(!is.na(unique(eid))), nrow(lut)))
wgt <- cos(lat * pi/180) # account for change cells size with poles

# assigns each raster cell its biome name by matching the cell’s ecoregion ID to the ecoregion lookup table
bio_of  <- lut$BIOME_NAME[match(eid, lut$ECO_ID)]
keepbase <- is.finite(eid) & !is.na(bio_of) & !bio_of %in% EXCLUDE

# validate relative richness raster
read_rel <- function(relfile){
  r <- rast(file.path(RELD, relfile))
  if (nlyr(r) != 1) stop(relfile, " must have exactly one layer", call.=FALSE)
  compareGeom(r, tmpl, crs=TRUE, ext=TRUE, rowcol=TRUE, res=TRUE, stopOnError=TRUE)
  values(r)[,1]
}


# =============================================================================
# 5. BLOCK SIZE — within-ecoregion autocorrelation range, converted to degrees
# =============================================================================

# estimates how far spatial autocorrelation extends within ecoregions or biomes, so an appropriate spatial-bootstrap block size can be chosen
range_within <- function(relfile, tag, level, min_cells, per_unit, bins,
                         thresh=CORR_THRESH, seed=1){
  # Read relative-richness values and retain valid occupied land cells
  rel  <- read_rel(relfile)
  keep <- keepbase & is.finite(rel) & rel > 0
  # Assign each retained cell to an ecoregion or biome
  unit <- if (level == "ecoregion") as.character(eid)[keep] else bio_of[keep]
  # Richness values and cell coordinates
  v <- rel[keep]; x <- lon[keep]; y <- lat[keep]
  # Retain units containing enough cells for the calculation
  tb <- table(unit); big <- names(tb)[tb >= min_cells]
  set.seed(seed)
  # Centre of each distance bin, in kilometres and prepare objects
  mid <- (bins[-1] + bins[-length(bins)])/2
  num <- numeric(length(mid)); den <- numeric(length(mid))
  nub <- numeric(length(mid)); used <- 0L
  # Analyse each sufficiently large ecoregion or biome separately
  for (uu in big) {
    i <- which(unit == uu); if (length(i) > per_unit) i <- sample(i, per_unit)
    z <- scale(v[i])[,1]; if (!all(is.finite(z))) next
    used <- used + 1L
    D <- outer(seq_along(i), seq_along(i), function(a,b) hav(x[i][a], y[i][a], x[i][b], y[i][b]))
    P <- outer(z, z); k <- findInterval(D, bins); ok <- D > 0 & k >= 1 & k <= length(mid)
    sm <- tapply(P[ok], factor(k[ok], levels=seq_along(mid)), sum); sm[is.na(sm)] <- 0  # Sum similarity values within each distance bin
    cnt <- as.numeric(table(factor(k[ok], levels=seq_along(mid)))) # Count cell pairs within each distance bin
    num <- num + as.numeric(sm)
    den <- den + cnt
    nub <- nub + (cnt > 0) }
  # Mean spatial correlation in each distance bin
  I <- ifelse(den > 0, num/den, NA_real_)
  below <- is.finite(I) & I < thresh
  hit <- which(below & c(below[-1], FALSE))[1]
  if (is.na(hit)) warning(sprintf("%s/%s: correlation never falls below %.2f within %.0f km",
                                  tag, level, thresh, max(mid)), call.=FALSE)
  # Return the correlation curve and estimated autocorrelation range
  data.frame(condition=tag, level=level, dist_km=mid, spatial_corr=round(I,4),
             n_ordered_pairs=den, n_units_bin=nub,
             range_km=if (is.na(hit)) NA_real_ else mid[hit], n_units=used)
}

message("within-ecoregion autocorrelation range")

# Estimate within-ecoregion spatial autocorrelation for every richness condition
cor_all <- do.call(rbind, lapply(names(RELF), function(k)
  range_within(RELF[[k]], k, "ecoregion", min_cells=CORR_MIN_CELLS,
               per_unit=CORR_PER_UNIT, bins=CORR_BINS)))
write_csv(cor_all, file.path(outdir, "within_ecoregion_correlogram.csv"))

# Summary table
mid_s <- sort(unique(cor_all$dist_km[cor_all$dist_km <= SUMM_MAX]))
pick  <- function(k, f) f(cor_all[cor_all$condition == k, ])
cor_summ <- data.frame(metric=c(sprintf("Spatial correlation at %g km", mid_s), "Range (km)",
                                "Ecoregions", "Ecoregions contributing (min over rows above)"))
for (k in names(RELF))
  cor_summ[[k]] <- c(pick(k, function(x) x$spatial_corr[match(mid_s, x$dist_km)]),
                     pick(k, function(x) x$range_km[1]),
                     pick(k, function(x) x$n_units[1]),
                     pick(k, function(x) min(x$n_units_bin[x$dist_km <= SUMM_MAX])))
names(cor_summ) <- c("metric", "spiny", "prickly", "thorny")
write_csv(cor_summ, file.path(outdir, "within_ecoregion_correlogram_summary.csv"))

# converts the largest estimated spatial-autocorrelation range from kilometres into a conservative block size in degrees for the spatial bootstrap
LAT90 <- max(vapply(RELF, function(f){
  r <- read_rel(f)
  as.numeric(quantile(abs(lat[keepbase & is.finite(r) & r > 0]), .90, na.rm=TRUE))
}, numeric(1)))
KM_PER_DEG <- 111.32 * cos(LAT90 * pi/180)

.rk <- tapply(cor_all$range_km, cor_all$condition, `[`, 1)
if (anyNA(.rk))
  stop("within-ecoregion autocorrelation range is NA for: ",
       paste(names(.rk)[is.na(.rk)], collapse=", "),
       " -- the correlogram never fell below the threshold, so no block size can be derived",
       call.=FALSE)
RANGE_KM  <- max(.rk)
BLOCK_ECO <- ceiling(RANGE_KM / KM_PER_DEG)
message(sprintf("  1 deg lon = %.0f km at |lat| p90 = %.1f", KM_PER_DEG, LAT90))
message(sprintf("  within-ecoregion range %.0f km -> ecoregion block %d deg; biomes resample ecoregions",
                RANGE_KM, BLOCK_ECO))


# =============================================================================
# 6. SUMMARIES — biome densities, ecoregion means, block-bootstrap CIs
# =============================================================================

# Build the analysis table for one relative-richness condition
build_d <- function(relfile){
  rel <- read_rel(relfile)
  keep <- is.finite(rel) & rel>0 & is.finite(eid)
  data.frame(rel=rel[keep], w=wgt[keep], ECO_ID=eid[keep],
             blk=blk_at(BLOCK_ECO)[keep]) |>
    left_join(lut, by="ECO_ID") |> filter(!BIOME_NAME %in% EXCLUDE, !is.na(BIOME_NAME))
}

# produces biome-level descriptive summaries and ecoregion-level means with spatial-bootstrap confidence intervals for one richness condition
summarise_cond <- function(relfile, tag){
  d <- build_d(relfile)
  # Area-weighted mean of relative richness
  wmean <- function(x) sum(x$rel * x$w)/sum(x$w)
  # Biome summaries
  bio <- d |> group_by(BIOME_NAME) |> group_modify(~{
      tibble(mean=wmean(.x),
             q25=wquantile(.x$rel, .x$w, .25), q75=wquantile(.x$rel, .x$w, .75), # Area-weighted cell-level interquartile
             cell_median=wquantile(.x$rel, .x$w, .50),
             hi=wquantile(.x$rel, .x$w, .995),   # right end of the density for this biome
             n_cells=nrow(.x), n_eco=n_distinct(.x$ECO_ID),
             eff_eco=kish(tapply(.x$w, .x$ECO_ID, sum))) }) |> # Effective number of ecoregions after accounting for unequal areas
    ungroup() |> arrange(desc(mean)) |> rename(biome=BIOME_NAME) |>
    mutate(rsd=(q75-q25)/1.349)   # robust SD: the scale the kernel width is set from

  # Area-weighted density of the CELL distribution, one curve per biome. Each is truncated
  # at its own 99.5th percentile so no curve is cut off by a limit set by another biome,
  # and normalised to equal PEAK height: cell counts span 369x, so equal area would leave
  # the smallest biome invisible. The count travels with each row instead.
  dens <- do.call(rbind, lapply(seq_len(nrow(bio)), function(i){
    z <- d[d$BIOME_NAME == bio$biome[i], ]
    dn <- density(z$rel, weights = z$w/sum(z$w), bw = bio$rsd[i]*BW_FRAC,
                  from = 0, to = bio$hi[i], n = 512)
    data.frame(biome = bio$biome[i], x = dn$x, y = dn$y/max(dn$y)) }))
  # Ecoregion summaries
  ec <- d |> group_by(ECO_ID) |> filter(n()>=MIN_ECO_CELLS) |> group_modify(~{
      bb <- boot_wmean(.x$rel, .x$w, .x$blk, nb=NB_ECO, seed=.y$ECO_ID%%9999) # Spatial block-bootstrap distribution of the area-weighted mean
      tibble(ECO_NAME=.x$ECO_NAME[1], BIOME_NAME=.x$BIOME_NAME[1], mean=wmean(.x),
             lo=unname(quantile(bb,.025,na.rm=TRUE)), hi=unname(quantile(bb,.975,na.rm=TRUE)), # Bootstrap 95% confidence interval
             n_cells=nrow(.x), n_blocks=n_distinct(.x$blk),
             eff_blocks=kish(tapply(.x$w, .x$blk, sum))) }) |> # Effective block count after accounting for unequal block areas
    ungroup() |> arrange(desc(mean))

  # Remove unreliable confidence intervals
  thin_ec <- ec$eff_blocks < MIN_EFF_CI
  ec$lo[thin_ec] <- NA_real_; ec$hi[thin_ec] <- NA_real_
  message(sprintf("  %s: %d of %d ecoregions have fewer than %.0f effective blocks -> CI not reported (median effective %.1f, nominal %.1f)",
                  tag, sum(thin_ec), nrow(ec), MIN_EFF_CI,
                  median(ec$eff_blocks), median(ec$n_blocks)))
  message(sprintf("  %s: biomes hold %d-%d ecoregions (effective %.1f-%.1f)",
                  tag, min(bio$n_eco), max(bio$n_eco), min(bio$eff_eco), max(bio$eff_eco)))
  # A biome with no spiny-trunk species anywhere contributes no cells, so it drops out of
  # every table and figure. Report it: an absent bar is a result, not missing data.
  absent <- setdiff(na.omit(unique(lut$BIOME_NAME[!lut$BIOME_NAME %in% EXCLUDE])), bio$biome)
  if (length(absent))
    message(sprintf("  %s: no spiny-trunk species anywhere in %s -> absent from the figures by construction",
                    tag, paste(absent, collapse="; ")))
  message(sprintf("  %s: %d cells carry a value (%.1f%% of the retained land area)",
                  tag, sum(bio$n_cells), 100*sum(d$w)/sum(wgt[keepbase])))

  # Check for strongly skewed biomde distribution
  skewed <- bio$biome[which(bio$mean < bio$q25 | bio$mean > bio$q75)] # Find biomes where the weighted mean lies outside the weighted IQR
  if (length(skewed))
    message("  ", tag, ": area-weighted mean outside its own weighted IQR (heavy-tailed): ",
            paste(skewed, collapse="; "))

  write_csv(bio, file.path(outdir, sprintf("biome_summary_%s.csv", tag)))
  write_csv(ec, file.path(outdir, sprintf("ecoregion_summary_%s.csv", tag)))
  list(bio=bio, ec=ec, dens=dens)
}

message("biome/ecoregion distributions (block bootstrap)")
S <- setNames(lapply(names(RELF), function(k) summarise_cond(RELF[[k]], k)), names(RELF))


# =============================================================================
# 7. FIGURES — maps, biome panel, ecoregion panel, Fig-3 composite
# =============================================================================

# ---- 7a. maps ----------------------------------------------------------------
message("maps")
borders <- st_read(BORD, quiet=TRUE) |> st_make_valid() |> st_wrap_dateline(options=c("WRAPDATELINE=YES"))

mk_map <- function(tif, ttl, clab, pal, scale_1e3=FALSE, inset=FALSE){
  r <- rast(file.path(RELD,tif)); r[r<=0] <- NA; d <- as.data.frame(r, xy=TRUE); names(d) <- c("x","y","v")
  vmax <- as.numeric(quantile(d$v, .99, na.rm=TRUE))
  g <- ggplot() + geom_raster(data=d, aes(x,y,fill=v)) +
    geom_sf(data=borders, fill=NA, colour="grey25", linewidth=.12) +
    scale_fill_viridis_c(option=pal, limits=c(min(d$v,na.rm=TRUE),vmax), oob=squish, name=clab,
      labels=if (scale_1e3) label_number(scale=1e3) else waiver(),
      guide=guide_colourbar(barwidth=unit(.38,"cm"), barheight=unit(if (inset) 4.4 else 5.4,"cm"))) +
    coord_sf(xlim=c(-180,180), ylim=c(-58,84), expand=FALSE) + labs(title=ttl) +
    theme_void(base_size=11) +
    theme(plot.title=element_text(face="bold", size=12, hjust=0, margin=margin(t=3,b=10)),
          plot.title.position="plot",
          legend.title=element_text(size=8.5, hjust=0), legend.text=element_text(size=7.5))
  if (inset)
    g + theme(legend.position="inside", legend.position.inside=c(0.075, 0.33),
              legend.justification=c(0.5, 0.5), legend.direction="vertical",
              legend.background=element_rect(fill=alpha("white", .75), colour=NA),
              legend.margin=margin(3,4,3,4), plot.margin=margin(2,2,2,2))
  else
    g + theme(legend.position="right", legend.justification="center",
              legend.margin=margin(l=4), plot.margin=margin(2,2,2,2))
}
maps <- list(
  c("Absolute_Richness_Tot.tif","Absolute richness of spiny-trunk species (0.1°)","Species richness","D","MAP_absolute_richness"),
  c("Absolute_Richness_Prickly.tif","Absolute richness of prickly-trunk species (0.1°)","Species richness","D","MAP_absolute_prickly"),
  c("Absolute_Richness_Thorny.tif","Absolute richness of thorny-trunk species (0.1°)","Species richness","D","MAP_absolute_thorny"),
  c("Relative_Richness_Tot.tif","Relative richness of spiny-trunk species (0.1°)","Relative richness","D","MAP_relative_richness"),
  c("Relative_Richness_Prickly.tif","Relative richness of prickly-trunk species (0.1°)","Relative richness","D","MAP_relative_prickly"),
  c("Relative_Richness_Thorny.tif","Relative richness of thorny-trunk species (0.1°)","Relative richness","D","MAP_relative_thorny"))
for (m in maps) ggsave(file.path(outdir,paste0(m[5],".png")), mk_map(m[1],m[2],m[3],m[4]), width=12, height=5.2, dpi=170)

# ---- 7b. axis labels and captions --------------------------------------------
XLAB_B <- expression("Relative richness ("%*%10^-3*")")
XLAB_E <- expression("Mean relative richness ("%*%10^-3*", 95% CI)")
CAP_B  <- paste0("Restricted to cells holding at least one spiny-trunk species. Area-weighted density of\n",
                 "the cell distribution, truncated at the 99.5th percentile of each biome and scaled to equal\n",
                 "peak height; cell count at right.  Diamond = area-weighted mean;  tick = area-weighted median.\n",
                 "Bandwidth = ", BW_FRAC, " x IQR/1.349 within each biome, so all biomes are smoothed at the same\n",
                 "fraction of their own spread rather than by cell count. The curves show shape, not fine structure.")
# %.0f, not %d: ceiling() returns a double and sprintf("%d", 2.0) is an error in R.
CAP_E  <- sprintf(paste("Restricted to cells holding at least one spiny-trunk species.",
                        "point = area-weighted mean;\nbar = 95%% block-bootstrap CI over %.0f-deg blocks;",
                        "(n) = effective independent blocks.\nOpen symbol: fewer than %.0f effective",
                        "blocks, too few for a percentile bootstrap, so no interval is reported."),
                  BLOCK_ECO, MIN_EFF_CI)

# ---- 7c. panel builders ------------------------------------------------------
# Ridgeline of the weighted cell densities. Rows keep the descending-mean order of the
# biome table, and biomes absent from panel C stay grey, so B and C read as one figure.
mk_biome <- function(bio, dens, keyed) {
  bio$pos  <- match(bio$biome, rev(bio$biome))
  dens$pos <- bio$pos[match(dens$biome, bio$biome)]
  xmax <- max(bio$hi)*1.02
  ggplot() +
    geom_ribbon(data=dens, aes(x=x, ymin=pos, ymax=pos + y*RIDGE_H, fill=biome, group=pos),
                colour="grey25", linewidth=.22, alpha=.82) +
    geom_segment(data=bio, aes(x=cell_median, xend=cell_median, y=pos, yend=pos+.33),
                 colour="grey15", linewidth=.65) +
    geom_point(data=bio, aes(x=mean, y=pos), shape=18, size=2.9, colour="grey10") +
    geom_text(data=bio, aes(x=xmax, y=pos+.30, label=format(n_cells, big.mark=",")),
              hjust=1, size=2.5, colour="grey45") +
    scale_fill_manual(values=sapply(bio$biome, bcol_keyed, keyed=keyed), guide="none") +
    scale_y_continuous(breaks=bio$pos, labels=blab(bio$biome), expand=expansion(add=c(.1,1.05))) +
    scale_x_continuous(labels=label_number(scale=1e3), limits=c(0,xmax),
                       expand=expansion(mult=c(.01,.02))) +
    labs(x=XLAB_B, y=NULL, caption=CAP_B) +
    theme(axis.text.y=element_text(size=9),
          panel.grid.major.x=element_line(colour="grey92", linewidth=.3),
          axis.line.x=element_line(colour="grey30"),
          plot.caption=element_text(size=7, colour="grey35", hjust=0)) }

mk_eco <- function(ec, legend=TRUE) { top <- head(ec,25)
  # the effective block count travels with the name: an interval the reader can see is backed by
  # 34 blocks should not look like one backed by 4
  labs_txt <- sprintf("%s  (%.0f)", top$ECO_NAME, top$eff_blocks)
  top$lab  <- factor(labs_txt, levels=rev(labs_txt))
  top$no_ci <- is.na(top$lo)
  ggplot(top, aes(mean, lab)) + geom_errorbar(aes(xmin=lo,xmax=hi), orientation="y", width=0, colour="grey65", linewidth=.5) +
    geom_point(aes(colour=BIOME_NAME, shape=no_ci), size=3.1, fill="white", stroke=1) +
    scale_shape_manual(values=c(`FALSE`=16, `TRUE`=21), guide="none") +
    scale_colour_manual(values=sapply(sort(unique(top$BIOME_NAME)),bcol),
      labels=function(b) blab(b), name="Biome", guide=if (legend) "legend" else "none") +
    scale_x_continuous(labels=label_number(scale=1e3)) +
    labs(x=XLAB_E, y=NULL, caption=CAP_E) +
    theme(axis.text.y=element_text(size=9), legend.position="right",
          legend.text=element_text(size=7), legend.title=element_text(size=8),
          plot.caption=element_text(size=7, colour="grey35", hjust=0)) }

# ---- 7d. per-condition figures and the composite -----------------------------
message("figures + composite")
for (cond in names(S)){
  bio <- S[[cond]]$bio; ec <- S[[cond]]$ec; dens <- S[[cond]]$dens
  keyed <- unique(head(ec,25)$BIOME_NAME)
  miss <- setdiff(unique(c(bio$biome, keyed)), names(BCOL))
  if (length(miss)) warning(cond, ": no colour assigned for ", paste(miss, collapse="; "), call.=FALSE)
  ggsave(file.path(outdir,sprintf("FIG_biome_%s.png",cond)),
         mk_biome(bio, dens, keyed) + labs(title=sprintf("Relative %s-trunk richness by biome", CLAB[[cond]])),
         width=9.4, height=6.8, dpi=170)
  ggsave(file.path(outdir,sprintf("FIG_ecoregion_%s.png",cond)),
         mk_eco(ec, legend=TRUE) + labs(title=sprintf("Top-25 ecoregions, relative %s-trunk richness", CLAB[[cond]])),
         width=11, height=9, dpi=170)
  pa <- mk_map(RELF[[cond]], sprintf("A)  Relative richness of %s-trunk species (0.1°)", CLAB[[cond]]),
               expression(atop("Relative richness","("%*%10^-3*")")), "D", scale_1e3=TRUE, inset=TRUE)
  pb <- mk_biome(bio, dens, keyed) + labs(title="B)  By biome")
  pc <- mk_eco(ec, legend=FALSE) + labs(title="C)  Top 25 ecoregions")

  # The ridgeline needs more vertical room than the point-and-bar panel it replaced, so the
  # second row takes more of the height and the whole figure is taller.
  comp <- free(pa, side="l") / (pb + pc + plot_layout(widths=c(1,1.12))) +
          plot_layout(heights=c(0.95,1.32))
  ggsave(file.path(outdir,sprintf("FIG_composite_%s.png",cond)), comp,
         width=15.5, height=15.0, dpi=170, bg="white")
}
message("DONE -> ", outdir)
