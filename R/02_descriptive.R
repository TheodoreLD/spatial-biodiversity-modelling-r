# =============================================================================
# 02_descriptive.R — richness maps, biome & ecoregion distributions (block-bootstrap
# CIs + rank-stability), and the Fig-3 composite. All R (ggplot).
# Inputs via env vars: TRUNK_REL_DIR (folder with Absolute_/Relative_
# Richness_*.tif), TRUNK_ECO (Ecoregions2017.shp), TRUNK_BORDERS (world_borders.geojson).
# =============================================================================
suppressPackageStartupMessages({ library(terra); library(sf); library(dplyr); library(readr)
  library(ggplot2); library(patchwork); library(scales) })
set.seed(0)
RELD <- Sys.getenv("TRUNK_REL_DIR", ".")
ECOSHP <- Sys.getenv("TRUNK_ECO", "Ecoregions2017.shp")
BORD <- Sys.getenv("TRUNK_BORDERS", "")
outdir <- file.path("outputs","descriptive"); dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
theme_set(theme_minimal(base_size=10) + theme(panel.grid=element_blank(),
          plot.title=element_text(face="bold",size=11), axis.line=element_line(colour="grey30")))
BCOL <- c("Tropical & Subtropical Grasslands, Savannas & Shrublands"="#8c6d31","Flooded Grasslands & Savannas"="#7b3fa0",
          "Tropical & Subtropical Dry Broadleaf Forests"="#c0392b","Tropical & Subtropical Moist Broadleaf Forests"="#2e8b57",
          "Deserts & Xeric Shrublands"="#d4b106","Mangroves"="#2c6fbb")
bcol <- function(b) ifelse(b %in% names(BCOL), BCOL[b], "#888888")
blab <- function(b) gsub("Xeric Shrublands","Xeric Shrub.", gsub("Grasslands, Savannas & Shrublands","Grasslands & Savannas",
        gsub("Tropical & Subtropical","Trop./Subtrop.", b)))
EXCLUDE <- c("N/A","Rock and Ice","Rock & Ice","Inland Water","Lakes","Lake")

# ---- borders + maps ----------------------------------------------------------
message("maps")
borders <- st_read(BORD, quiet=TRUE) |> st_make_valid() |> st_wrap_dateline(options=c("WRAPDATELINE=YES"))
mk_map <- function(tif, ttl, clab, pal){
  r <- rast(file.path(RELD,tif)); r[r<=0] <- NA; d <- as.data.frame(r, xy=TRUE); names(d) <- c("x","y","v")
  vmax <- as.numeric(quantile(d$v, .99, na.rm=TRUE))
  ggplot() + geom_raster(data=d, aes(x,y,fill=v)) + geom_sf(data=borders, fill=NA, colour="black", linewidth=.12) +
    scale_fill_viridis_c(option=pal, limits=c(min(d$v,na.rm=TRUE),vmax), oob=squish, name=clab) +
    coord_sf(xlim=c(-180,180), ylim=c(-58,84), expand=FALSE) + labs(title=ttl) +
    theme_void(base_size=11) + theme(plot.title=element_text(face="bold"),
      legend.key.width=unit(.28,"cm"), legend.key.height=unit(1.1,"cm"))
}
maps <- list(
  c("Absolute_Richness_Tot.tif","Absolute richness of spiny-trunk species (0.1°)","Species richness","D","MAP_absolute_richness"),
  c("Absolute_Richness_Prickly.tif","Absolute richness of prickly-trunk species (0.1°)","Species richness","D","MAP_absolute_prickly"),
  c("Absolute_Richness_Thorny.tif","Absolute richness of thorny-trunk species (0.1°)","Species richness","D","MAP_absolute_thorny"),
  c("Relative_Richness_Tot.tif","Relative richness of spiny-trunk species (0.1°)","Relative richness","D","MAP_relative_richness"),
  c("Relative_Richness_Prickly.tif","Relative richness of prickly-trunk species (0.1°)","Relative richness","D","MAP_relative_prickly"),
  c("Relative_Richness_Thorny.tif","Relative richness of thorny-trunk species (0.1°)","Relative richness","D","MAP_relative_thorny"))
for (m in maps) ggsave(file.path(outdir,paste0(m[5],".png")), mk_map(m[1],m[2],m[3],m[4]), width=12, height=5.2, dpi=170)

# ---- overlay: assign each 0.1° cell to its ecoregion + biome -----------------
message("overlay (rasterize ecoregions)")
eco <- st_read(ECOSHP, quiet=TRUE) |> st_make_valid()
tmpl <- rast(ext(-180,180,-90,90), resolution=0.1, crs="EPSG:4326")
eco_id <- rasterize(vect(eco), tmpl, field="ECO_ID")
lut <- eco |> st_drop_geometry() |> distinct(ECO_ID, ECO_NAME, BIOME_NAME)
n <- ncell(tmpl); col0 <- (seq_len(n)-1) %% 3600; row0 <- (seq_len(n)-1) %/% 3600
lon <- -180 + (col0+0.5)*0.1; lat <- 90 - (row0+0.5)*0.1
blk5 <- floor((lon+180)/5)*1000 + floor((lat+90)/5)
eid <- values(eco_id)[,1]

boot <- function(v, b, nb=1000, seed=0){ set.seed(seed); ub <- unique(b); by <- split(v, b)
  o <- numeric(nb); for (i in 1:nb) o[i] <- mean(unlist(by[as.character(sample(ub,length(ub),TRUE))])); o }

summarise_cond <- function(relfile, tag){
  rel <- values(rast(file.path(RELD,relfile)))[,1]
  keep <- is.finite(rel) & rel>0 & is.finite(eid)
  d <- data.frame(rel=rel[keep], ECO_ID=eid[keep], blk=blk5[keep]) |>
       left_join(lut, by="ECO_ID") |> filter(!BIOME_NAME %in% EXCLUDE, !is.na(BIOME_NAME))
  # biome
  bio <- d |> group_by(BIOME_NAME) |> group_modify(~{ bb <- boot(.x$rel,.x$blk)
      tibble(mean=mean(.x$rel), lo=quantile(bb,.025,na.rm=TRUE), hi=quantile(bb,.975,na.rm=TRUE), n_cells=nrow(.x), n_eco=n_distinct(.x$ECO_ID)) }) |>
    ungroup() |> arrange(desc(mean)) |> rename(biome=BIOME_NAME)
  write_csv(bio, file.path(outdir, sprintf("biome_summary_%s.csv", tag)))
  # ecoregion (>=10 cells) + rank stability (top-15 across block resamples)
  ec <- d |> group_by(ECO_ID) |> filter(n()>=10) |> group_modify(~{ bb <- boot(.x$rel,.x$blk, nb=400, seed=.y$ECO_ID%%9999)
      tibble(ECO_NAME=.x$ECO_NAME[1], BIOME_NAME=.x$BIOME_NAME[1], mean=mean(.x$rel), lo=quantile(bb,.025,na.rm=TRUE), hi=quantile(bb,.975,na.rm=TRUE), n_cells=nrow(.x)) }) |>
    ungroup() |> arrange(desc(mean))
  ub <- unique(d$blk); cnt <- setNames(numeric(nrow(ec)), ec$ECO_ID); NR <- 250
  for (i in 1:NR){ s <- d[d$blk %in% sample(ub,length(ub),TRUE),]
    mm <- s |> group_by(ECO_ID) |> summarise(m=mean(rel), k=n(), .groups="drop") |> filter(k>=10) |> arrange(desc(m)) |> head(15)
    for (e in as.character(mm$ECO_ID)) if (e %in% names(cnt)) cnt[e] <- cnt[e]+1 }
  ec$top15_freq <- as.numeric(cnt[as.character(ec$ECO_ID)])/NR
  write_csv(ec, file.path(outdir, sprintf("ecoregion_summary_%s.csv", tag)))
  list(bio=bio, ec=ec)
}
message("biome/ecoregion distributions (block bootstrap)")
S <- list(total=summarise_cond("Relative_Richness_Tot.tif","total"),
          prickly=summarise_cond("Relative_Richness_Prickly.tif","prickly"),
          thorny=summarise_cond("Relative_Richness_Thorny.tif","thorny"))

# ---- biome & ecoregion figures + composite -----------------------------------
message("figures + composite")
mk_biome <- function(bio) { bio$biome <- factor(bio$biome, levels=rev(bio$biome))
  ggplot(bio, aes(mean, biome)) + geom_errorbarh(aes(xmin=lo,xmax=hi),height=0,colour="grey60") +
    geom_point(aes(colour=biome), size=3) + scale_colour_manual(values=sapply(levels(bio$biome),bcol),guide="none") +
    scale_y_discrete(labels=function(b) blab(as.character(b))) + labs(x="Mean relative richness (± 95% CI)", y=NULL) }
mk_eco <- function(ec) { top <- head(ec,25); top$ECO_NAME <- factor(top$ECO_NAME, levels=rev(top$ECO_NAME))
  ggplot(top, aes(mean, ECO_NAME)) + geom_errorbarh(aes(xmin=lo,xmax=hi),height=0,colour="grey60") +
    geom_point(aes(colour=BIOME_NAME), size=2.6) + scale_colour_manual(values=sapply(unique(top$BIOME_NAME),bcol),
      labels=function(b) blab(b), name="Biome") + labs(x="Mean relative richness (± 95% CI)", y=NULL) +
    theme(legend.position="right", legend.text=element_text(size=6.5), legend.title=element_text(size=7.5)) }
for (cond in names(S)){
  ggsave(file.path(outdir,sprintf("FIG_biome_%s.png",cond)), mk_biome(S[[cond]]$bio) + labs(title=sprintf("Relative %s-trunk richness by biome",cond)), width=10, height=6, dpi=170)
  ggsave(file.path(outdir,sprintf("FIG_ecoregion_%s.png",cond)), mk_eco(S[[cond]]$ec) + labs(title=sprintf("Top-25 ecoregions, relative %s-trunk richness",cond)), width=10.5, height=9, dpi=170)
  # composite: A map, B biome, C ecoregions
  relmap <- switch(cond, total="Relative_Richness_Tot.tif", prickly="Relative_Richness_Prickly.tif", thorny="Relative_Richness_Thorny.tif")
  pa <- mk_map(relmap, sprintf("A)  Relative richness of %s-trunk species (0.1°)",cond), "Relative richness", "D")
  pb <- mk_biome(S[[cond]]$bio) + labs(title="B)  By biome") + theme(plot.title=element_text(hjust=0))
  pc <- mk_eco(S[[cond]]$ec) + labs(title="C)  Top 25 ecoregions") + theme(plot.title=element_text(hjust=0))
  comp <- pa / (pb | pc) + plot_layout(heights=c(1.25,1.5))
  ggsave(file.path(outdir,sprintf("FIG_composite_%s.png",cond)), comp, width=16, height=13.4, dpi=160)
}
message("DONE -> ", outdir)
