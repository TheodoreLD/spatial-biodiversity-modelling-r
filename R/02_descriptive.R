# =============================================================================
# 02_descriptive.R — richness maps, biome & ecoregion distributions (block-bootstrap
# CIs), and the Fig-3 composite.
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

bcol_keyed <- function(b, keyed) ifelse(b %in% keyed, bcol(b), GREY)
blab <- function(b) gsub("Xeric Shrublands","Xeric Shrub.", gsub("Grasslands, Savannas & Shrublands","Grasslands & Savannas",
        gsub("Tropical & Subtropical","Trop./Subtrop.", b)))
EXCLUDE <- c("N/A","Rock and Ice","Rock & Ice","Inland Water","Lakes","Lake")

# ---- borders + maps ----------------------------------------------------------
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

# ---- overlay: assign each 0.1° cell to its ecoregion + biome -----------------
message("overlay (rasterize ecoregions)")
eco <- st_read(ECOSHP, quiet=TRUE) |> st_make_valid()
tmpl <- rast(ext(-180,180,-90,90), resolution=0.1, crs="EPSG:4326")
eco_id <- rasterize(vect(eco), tmpl, field="ECO_ID")
lut <- eco |> st_drop_geometry() |> distinct(ECO_ID, ECO_NAME, BIOME_NAME)
n <- ncell(tmpl); nc <- ncol(tmpl); ev <- as.vector(ext(tmpl)); rs <- res(tmpl)
col0 <- (seq_len(n)-1) %% nc; row0 <- (seq_len(n)-1) %/% nc
lon <- ev[["xmin"]] + (col0+0.5)*rs[1]; lat <- ev[["ymax"]] - (row0+0.5)*rs[2]
blk5 <- floor((lon+180)/5)*1000 + floor((lat+90)/5)   
eid <- values(eco_id)[,1]

# 5 deg blocks: at or above the measured range of within-ecoregion autocorrelation
# (375-525 km across all three conditions), which is the dependence a block bootstrap
# on an ecoregion mean has to absorb. 

resample <- function(x) x[sample.int(length(x), length(x), replace=TRUE)]
boot <- function(v, b, nb=1000, seed=0){ set.seed(seed); ub <- unique(b); by <- split(v, b)
  o <- numeric(nb); for (i in 1:nb) o[i] <- mean(unlist(by[as.character(resample(ub))])); o }

summarise_cond <- function(relfile, tag){
  rel <- values(rast(file.path(RELD,relfile)))[,1]
  keep <- is.finite(rel) & rel>0 & is.finite(eid)
  d <- data.frame(rel=rel[keep], ECO_ID=eid[keep], blk=blk5[keep]) |>
       left_join(lut, by="ECO_ID") |> filter(!BIOME_NAME %in% EXCLUDE, !is.na(BIOME_NAME))

  bio <- d |> group_by(BIOME_NAME) |> group_modify(~{ bb <- boot(.x$rel,.x$blk)
      tibble(mean=mean(.x$rel), lo=quantile(bb,.025,na.rm=TRUE), hi=quantile(bb,.975,na.rm=TRUE),
             n_cells=nrow(.x), n_eco=n_distinct(.x$ECO_ID), n_blocks=n_distinct(.x$blk)) }) |>
    ungroup() |> arrange(desc(mean)) |> rename(biome=BIOME_NAME)

  ec <- d |> group_by(ECO_ID) |> filter(n()>=10) |> group_modify(~{ bb <- boot(.x$rel,.x$blk, nb=400, seed=.y$ECO_ID%%9999)
      tibble(ECO_NAME=.x$ECO_NAME[1], BIOME_NAME=.x$BIOME_NAME[1], mean=mean(.x$rel),
             lo=quantile(bb,.025,na.rm=TRUE), hi=quantile(bb,.975,na.rm=TRUE),
             n_cells=nrow(.x), n_blocks=n_distinct(.x$blk)) }) |>
    ungroup() |> arrange(desc(mean))
  # A unit contained in ONE 5 deg block has nothing to resample: every bootstrap draw
  # returns that block, so lo == hi and the interval reads as perfect precision when the
  # method could not estimate it at all. Report the mean and no interval.
  ec$lo[ec$n_blocks == 1] <- NA_real_; ec$hi[ec$n_blocks == 1] <- NA_real_
  bio$lo[bio$n_blocks == 1] <- NA_real_; bio$hi[bio$n_blocks == 1] <- NA_real_
  message(sprintf("  %s: %d of %d ecoregions occupy a single block -> CI not reported",
                  tag, sum(ec$n_blocks == 1), nrow(ec)))

  write_csv(bio, file.path(outdir, sprintf("biome_summary_%s.csv", tag)))
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
mk_biome <- function(bio, keyed) { bio$biome <- factor(bio$biome, levels=rev(bio$biome))
  ggplot(bio, aes(mean, biome)) + geom_errorbarh(aes(xmin=lo,xmax=hi),height=0,colour="grey65",linewidth=.5) +
    geom_point(aes(colour=biome), size=3.1) +
    scale_colour_manual(values=sapply(levels(bio$biome), bcol_keyed, keyed=keyed), guide="none") +
    scale_y_discrete(labels=function(b) blab(as.character(b))) +
    scale_x_continuous(labels=label_number(scale=1e3)) +
    labs(x=XLAB, y=NULL) + theme(axis.text.y=element_text(size=9)) }

mk_eco <- function(ec, legend=TRUE) { top <- head(ec,25)
  top$lab <- factor(top$ECO_NAME, levels=rev(top$ECO_NAME))
  top$no_ci <- is.na(top$lo)                       # single-block: interval not estimable
  cap <- paste0(if (legend) "" else "Colours as in panel B.  ",
                if (any(top$no_ci))
                  "Open symbols: the ecoregion lies within one 5\u00b0 block, so the block bootstrap cannot estimate an interval; the mean is shown alone."
                else "")
  ggplot(top, aes(mean, lab)) + geom_errorbarh(aes(xmin=lo,xmax=hi),height=0,colour="grey65",linewidth=.5) +
    geom_point(aes(colour=BIOME_NAME, shape=no_ci), size=3.1, fill="white", stroke=1) +
    scale_shape_manual(values=c(`FALSE`=16, `TRUE`=21), guide="none") +
    scale_colour_manual(values=sapply(sort(unique(top$BIOME_NAME)),bcol),
      labels=function(b) blab(b), name="Biome", guide=if (legend) "legend" else "none") +
    scale_x_continuous(labels=label_number(scale=1e3)) +
    labs(x=XLAB, y=NULL, caption=cap) +
    theme(axis.text.y=element_text(size=9), legend.position="right",
          legend.text=element_text(size=7), legend.title=element_text(size=8),
          plot.caption=element_text(hjust=0, size=7.5, colour="grey35")) }
CLAB <- c(total="spiny", prickly="prickly", thorny="thorny")   
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
  pa <- mk_map(RELF[[cond]], sprintf("A)  Relative richness of %s-trunk species (0.1°)", CLAB[[cond]]),
               expression(atop("Relative richness","("%*%10^-3*")")), "D", scale_1e3=TRUE, inset=TRUE)
  pb <- mk_biome(bio, keyed) + labs(title="B)  By biome")
  pc <- mk_eco(ec, legend=FALSE) + labs(title="C)  Top 25 ecoregions")
 
  comp <- free(pa, side="l") / (pb + pc + plot_layout(widths=c(1,1.15))) +
          plot_layout(heights=c(1.12,1.03))
  ggsave(file.path(outdir,sprintf("FIG_composite_%s.png",cond)), comp,
         width=15, height=13.0, dpi=170, bg="white")
}
message("DONE -> ", outdir)
