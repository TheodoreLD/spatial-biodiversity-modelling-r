# =============================================================================
# 06_model2_jackknife.R — species-drop robustness of the Model-2 clade RANKING
# (companion to 06_model2.R; the Model-1 analogue is 05_model1_jackknife.R).
#
# Reviewer 1's concern is that the response is built from only 74 mapped species,
# so any ordering of clades may hinge on a handful of them. Here we drop a fraction
# of the 74 species, re-derive relative richness (S_subset / W), re-aggregate to the
# same 1 deg unit, refit the positive-monotonic clade models per continent, and ask
# whether the RANKING survives.
#
# The analysis unit stays at 1 deg throughout: this varies the species pool, not the
# resolution, so it does not reopen the pseudo-replication question that the
# resolution sweep settled. Only the `total` condition is jackknifed — dropping
# species from the prickly/thorny subsets would change what the subset means.
#
# Reads (paths in R/config.R): predictor stack (env layers), clade stack, absolute +
# relative Tot richness rasters, and the 74 per-species .gpkg ranges.
# Writes: outputs/model2/model2_jackknife_ranks.csv (+ _summary.csv, + figure)
# =============================================================================
suppressPackageStartupMessages({ library(terra); library(catboost); library(dplyr)
  library(readr); library(ggplot2) })
source("R/config.R")
outdir <- paths$model2; dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

CLADES <- c("Pilosa","Cetartiodactyla","Diprotodontia","Primates","Rodentia",
            "Perissodactyla","Proboscidea","Cingulata","Litopterna","Notoungulata")
DISP <- function(x) ifelse(x == "Cetartiodactyla", "Artiodactyla", x)
ENV   <- c("vpi_mean","cold_mean","dry_mean","height")
CONTS <- c("Americas","Africa","Asia")
FRACS <- c(0.1, 0.2, 0.3); NREP <- 10
cbp <- list(loss_function="RMSE", depth=4, l2_leaf_reg=10, learning_rate=0.05,
            iterations=300, allow_writing_files=FALSE, logging_level="Silent")

# --- identical geometry / continent rules to 06_model2.R ----------------------
afro_asia_lon <- function(la)
  ifelse(la >= 31.5, 25, ifelse(la >= 12.5, 43 - (la - 12.5) * (43 - 32) / 19, 60))
cont_of <- function(lo, la) { out <- rep("Other", length(lo))
  out[lo>=-170 & lo<=-30 & la>=-60 & la<= 80] <- "Americas"
  out[lo>=  25 & lo<=180 & la>=-12 & la<= 80] <- "Asia"
  out[lo>= 112 & lo<=154 & la>=-44 & la<=-10] <- "Other"
  out[lo>= -20 & lo<= 55 & la>=-35 & la<= 40 & lo < afro_asia_lon(la)] <- "Africa"; out }
# Same ranking metric as 06_model2.R: the standardised slope of the constrained response,
# in SDs of relative richness per SD of clade richness. Kept identical here so the
# jackknife tests the statistic the paper reports rather than a different one.
zslope <- function(yv, xv, mn, sdx){ z <- (xv-mn)/sdx; unname(coef(lm(yv ~ z))[2]) }

# --- inputs -------------------------------------------------------------------
ps <- load_stack(paths$predictor_stack, ENV, "predictor stack (environment layers)")
cs <- rast(require_file(paths$clade_stack, "clade richness stack"))
miss <- setdiff(CLADES, names(cs))
if (length(miss)) stop("clade layers absent: ", paste(miss, collapse = ", "), call. = FALSE)

Av <- values(rast(require_file(paths$abs_richness, "absolute richness")))[,1]
Rv <- values(rast(require_file(paths$rel_richness, "relative richness")))[,1]
W  <- ifelse(is.finite(Rv) & Rv > 0, Av/Rv, NA)          # woody-plant richness denominator

env_v  <- values(ps); envfin <- rowSums(is.finite(env_v)) == length(ENV)
cl_v   <- values(cs[[CLADES]]); cl_v[is.na(cl_v)] <- 0

n <- ncell(ps); nc <- ncol(ps); ev <- as.vector(ext(ps)); rs <- res(ps)
col0 <- (seq_len(n)-1) %% nc; row0 <- (seq_len(n)-1) %/% nc
lon <- ev[["xmin"]] + (col0+0.5)*rs[1]; lat <- ev[["ymax"]] - (row0+0.5)*rs[2]
AGG <- round(1/rs[1]); Rb <- row0 %/% AGG; Cb <- col0 %/% AGG
MINSUB <- max(1, floor(0.10*AGG*AGG))

tmpl  <- rast(ext(-180,180,-90,90), res = rs[1], crs = "EPSG:4326")
gpkgs <- list.files(paths$spiny_ranges, "\\.gpkg$", full.names = TRUE)
message(length(gpkgs), " species; rasterising presences")
pres <- lapply(gpkgs, function(g) which(values(rasterize(vect(g), tmpl, background = 0, touches = TRUE))[,1] > 0))
nsp <- length(pres)

# --- rebuild the 1 deg table from an arbitrary species subset -----------------
build <- function(keep) {
  S <- integer(n); for (k in keep) S[pres[[k]]] <- S[pres[[k]]] + 1L
  m <- S > 0 & is.finite(W) & envfin
  d <- data.frame(rel = S[m]/W[m], R = Rb[m], C = Cb[m], lon = lon[m], lat = lat[m])
  d <- cbind(d, env_v[m, , drop = FALSE]); d <- cbind(d, cl_v[m, , drop = FALSE])
  g <- d |> group_by(R, C) |>
    summarise(across(everything(), \(z) mean(z, na.rm = TRUE)), nn = n(), .groups = "drop") |>
    filter(nn >= MINSUB) |> select(-nn) |> as.data.frame()
  g$continent <- cont_of(g$lon, g$lat); g
}

# --- slope_std for one table --------------------------------------------------
# Only the ranking metric is recomputed per species subset. The jackknife's question is
# whether the ORDER depends on which of the 74 species were mapped, so the spread of
# slope_std across replicates answers it; refitting shapes and R2 as well multiplies the
# runtime without adding to that answer.
scores_for <- function(g) {
  out <- list()
  for (cont in CONTS) {
    sub <- g[g$continent == cont, ]; if (nrow(sub) < 80) next
    y <- sub$rel; sdy <- sd(y); if (!is.finite(sdy) || sdy <= 0) next
    for (cl in CLADES) {
      x <- sub[[cl]]
      if (sum(x > 0) < 30 || length(unique(x)) < 5) next
      mono <- catboost.train(catboost.load_pool(matrix(x, ncol = 1), label = y), NULL,
                             params = c(cbp, list(monotone_constraints = "(1)")))
      q <- unique(quantile(x, seq(0, 1, length.out = 21)))
      out[[length(out)+1]] <- data.frame(continent = cont, clade = DISP(cl),
        slope_std = zslope(catboost.predict(mono, catboost.load_pool(matrix(q, ncol = 1))),
                           q, mean(x), sd(x)) / sdy)
    }
  }
  bind_rows(out)
}

# --- full model as the reference ordering -------------------------------------
full <- scores_for(build(seq_len(nsp))) |>
  group_by(continent) |> mutate(rank_full = rank(-slope_std, ties.method = "min")) |> ungroup()
message("full 74-species reference ranking computed")

# --- jackknife ----------------------------------------------------------------
rows <- list()
for (fr in FRACS) {
  for (rep in seq_len(NREP)) {
    set.seed(as.integer(1000*fr*100 + rep))
    keep <- sort(sample(nsp, round((1 - fr)*nsp)))
    s <- scores_for(build(keep))
    if (!nrow(s)) next
    s <- s |> group_by(continent) |> mutate(rank = rank(-slope_std, ties.method = "min")) |> ungroup()
    rows[[length(rows)+1]] <- s |> mutate(frac = fr, rep = rep, n_species = length(keep))
  }
  message(sprintf("  drop %d%% done (%d reps)", as.integer(fr*100), NREP))
}
jk <- bind_rows(rows)
write_csv(jk, file.path(outdir, "model2_jackknife_ranks.csv"))
# --- summary: how often does each clade stay in the top 3? --------------------
summ <- jk |> group_by(continent, clade) |>
  summarise(n_fits = n(), top1_freq = mean(rank == 1), top3_freq = mean(rank <= 3),
            median_rank = median(rank), slope_std_mean = mean(slope_std),
            slope_std_sd = sd(slope_std), .groups = "drop") |>
  left_join(full |> select(continent, clade, slope_std_full = slope_std, rank_full),
            by = c("continent","clade")) |>
  arrange(continent, rank_full)
write_csv(summ, file.path(outdir, "model2_jackknife_summary.csv"))
print(as.data.frame(summ))

# --- figure -------------------------------------------------------------------
COL <- c(Americas = "#E1812C", Africa = "#3A7CA5", Asia = "#4C9A4C")
pd <- jk |> mutate(continent = factor(continent, CONTS))
ord <- full |> arrange(continent, rank_full)
pd$clade <- factor(pd$clade, levels = rev(unique(ord$clade)))
ggsave(file.path(outdir, "FIG_M2_jackknife_rank.png"),
  ggplot(pd, aes(rank, clade, colour = continent)) +
    geom_jitter(width = 0.18, height = 0.18, alpha = .55, size = 1.5) +
    scale_colour_manual(values = COL, guide = "none") +
    scale_x_continuous(breaks = 1:10) +
    facet_wrap(~continent, scales = "free_y") +
    labs(x = "rank by the standardised slope of the constrained response", y = NULL,
         title = "Model 2 — clade ranking under species-drop jackknife (10/20/30% of the 74 species removed)",
         caption = "Each point is one jackknife replicate. A firm ranking gives tight columns; a spread means the order is not identifiable.") +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 11),
          axis.line = element_line(colour = "grey30"), panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold"),
          plot.caption = element_text(hjust = 0, size = 7.5, colour = "grey30")),
  width = 14, height = 5, dpi = 170)
message("DONE -> ", outdir)
