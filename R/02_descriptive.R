# =============================================================================
# 02_descriptive.R — richness maps, biome & ecoregion distributions (block-bootstrap
# CIs + rank-stability), and the Fig-3 composite.
# Inputs (paths from R/config.R): the Biogeography rasters folder (Absolute_/Relative_
# Richness_*.tif), the Ecoregions2017 shapefile, and data/assets/world_borders.geojson.
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
theme_set(theme_minimal(base_size=10) + theme(panel.grid=element_blank(),
          plot.title=element_text(face="bold",size=11), plot.title.position="plot",
          axis.line=element_line(colour="grey30")))
# Full WWF biome palette so nothing that reaches a figure falls back to grey.
# Warm = dry/open, green = forest, purple = flooded, blue = mangrove; temperate
# tones are lighter than their tropical counterpart.
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
# In the biome panel only biomes that also appear among the plotted ecoregions are
# coloured; the rest stay grey. The colours then function purely as the key for the
# ecoregion panel, which carries no legend of its own in the composite.
bcol_keyed <- function(b, keyed) ifelse(b %in% keyed, bcol(b), GREY)
blab <- function(b) gsub("Xeric Shrublands","Xeric Shrub.", gsub("Grasslands, Savannas & Shrublands","Grasslands & Savannas",
        gsub("Tropical & Subtropical","Trop./Subtrop.", b)))
EXCLUDE <- c("N/A","Rock and Ice","Rock & Ice","Inland Water","Lakes","Lake")

# ---- borders + maps ----------------------------------------------------------
message("maps")
borders <- st_read(BORD, quiet=TRUE) |> st_make_valid() |> st_wrap_dateline(options=c("WRAPDATELINE=YES"))
# Full land extent (-58/84, -180/180) so no country is sliced at the poles or the
# dateline. `inset` puts the colourbar inside the map, in the empty South Pacific:
# outside on the right it consumes panel width, which shifts the map frame off centre
# relative to panels B/C below. Used for the composite; the standalone maps, which
# have nothing to align with, keep the colourbar on the right.
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
    # theme_void leaves no axis strip, so the map panel starts immediately below the
    # title: without a bottom margin here the Arctic coastlines run into the title.
    theme(plot.title=element_text(face="bold", size=12, hjust=0, margin=margin(t=3,b=14)),
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

# ---- overlay: assign each 0.1° cell to its ecoregion + biome -----------------
message("overlay (rasterize ecoregions)")
eco <- st_read(ECOSHP, quiet=TRUE) |> st_make_valid()
tmpl <- rast(ext(-180,180,-90,90), resolution=0.1, crs="EPSG:4326")
eco_id <- rasterize(vect(eco), tmpl, field="ECO_ID")
lut <- eco |> st_drop_geometry() |> distinct(ECO_ID, ECO_NAME, BIOME_NAME)
# Cell -> lon/lat from the template itself, not from an assumed 3600-column grid.
n <- ncell(tmpl); nc <- ncol(tmpl); ev <- as.vector(ext(tmpl)); rs <- res(tmpl)
col0 <- (seq_len(n)-1) %% nc; row0 <- (seq_len(n)-1) %/% nc
lon <- ev[["xmin"]] + (col0+0.5)*rs[1]; lat <- ev[["ymax"]] - (row0+0.5)*rs[2]
blk5 <- floor((lon+180)/5)*1000 + floor((lat+90)/5)   # 5 deg block id for the bootstrap
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
XLAB <- expression("Mean relative richness ("%*%10^-3*", ± 95% CI)")
# `keyed` = the biomes represented among the plotted ecoregions; only those are
# coloured, so the biome panel reads as the key for the ecoregion panel.
mk_biome <- function(bio, keyed) { bio$biome <- factor(bio$biome, levels=rev(bio$biome))
  ggplot(bio, aes(mean, biome)) + geom_errorbarh(aes(xmin=lo,xmax=hi),height=0,colour="grey65",linewidth=.5) +
    geom_point(aes(colour=biome), size=3.1) +
    scale_colour_manual(values=sapply(levels(bio$biome), bcol_keyed, keyed=keyed), guide="none") +
    scale_y_discrete(labels=function(b) blab(as.character(b))) +
    scale_x_continuous(labels=label_number(scale=1e3)) +
    labs(x=XLAB, y=NULL) + theme(axis.text.y=element_text(size=9)) }
# legend=TRUE for the standalone figure (nothing else to key off); FALSE in the
# composite, where panel B supplies the key.
mk_eco <- function(ec, legend=TRUE) { top <- head(ec,25)
  lab <- paste0(top$ECO_NAME, "  (", round(top$top15_freq*100), "%)")   # append rank stability
  top$lab <- factor(lab, levels=rev(lab))
  cap <- paste0(if (legend) "" else "Colours as in panel B.   ",
                "(x%) = rank stability: share of spatial block-bootstrap resamples ranking the ecoregion in the global top 15")
  ggplot(top, aes(mean, lab)) + geom_errorbarh(aes(xmin=lo,xmax=hi),height=0,colour="grey65",linewidth=.5) +
    geom_point(aes(colour=BIOME_NAME), size=3.1) +
    scale_colour_manual(values=sapply(sort(unique(top$BIOME_NAME)),bcol),
      labels=function(b) blab(b), name="Biome", guide=if (legend) "legend" else "none") +
    scale_x_continuous(labels=label_number(scale=1e3)) +
    labs(x=XLAB, y=NULL, caption=cap) +
    theme(axis.text.y=element_text(size=9), legend.position="right",
          legend.text=element_text(size=7), legend.title=element_text(size=8),
          plot.caption=element_text(hjust=0, size=7.5, colour="grey35")) }
CLAB <- c(total="spiny", prickly="prickly", thorny="thorny")   # "total-trunk species" was wrong
RELF <- c(total="Relative_Richness_Tot.tif", prickly="Relative_Richness_Prickly.tif",
          thorny="Relative_Richness_Thorny.tif")
for (cond in names(S)){
  bio <- S[[cond]]$bio; ec <- S[[cond]]$ec
  keyed <- unique(head(ec,25)$BIOME_NAME)
  miss <- setdiff(unique(c(bio$biome, keyed)), names(BCOL))
  if (length(miss)) warning(cond, ": no colour assigned for ", paste(miss, collapse="; "), call.=FALSE)
  ggsave(file.path(outdir,sprintf("FIG_biome_%s.png",cond)),
         mk_biome(bio, keyed) + labs(title=sprintf("Relative %s-trunk richness by biome", CLAB[[cond]])),
         width=10, height=6, dpi=170)
  ggsave(file.path(outdir,sprintf("FIG_ecoregion_%s.png",cond)),
         mk_eco(ec, legend=TRUE) + labs(title=sprintf("Top-25 ecoregions, relative %s-trunk richness", CLAB[[cond]])),
         width=11, height=9, dpi=170)
  # composite: A map, B biome (the key), C ecoregions (no legend of its own)
  pa <- mk_map(RELF[[cond]], sprintf("A)  Relative richness of %s-trunk species (0.1°)", CLAB[[cond]]),
               expression(atop("Relative richness","("%*%10^-3*")")), "D", scale_1e3=TRUE, inset=TRUE)
  pb <- mk_biome(bio, keyed) + labs(title="B)  By biome")
  pc <- mk_eco(ec, legend=FALSE) + labs(title="C)  Top 25 ecoregions")
  # free(): patchwork aligns plot panels down a column, so without it the map is
  # indented to start where panel B's panel starts -- i.e. after B's long biome labels
  # (~450px), leaving a gap under the "A)" title and pushing the frame off centre.
  # Panel A's height is set just above what the width-limited map needs, so the map
  # spans the full width with little vertical slack left for it to centre in.
  comp <- free(pa) / (pb + pc + plot_layout(widths=c(1,1.15))) + plot_layout(heights=c(1.08,1.03))
  ggsave(file.path(outdir,sprintf("FIG_composite_%s.png",cond)), comp,
         width=15, height=12.9, dpi=170, bg="white")
}
message("DONE -> ", outdir)
