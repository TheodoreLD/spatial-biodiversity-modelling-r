# =============================================================================
# 06_model2_paper_figure.R — the combined Model-2 figure for the manuscript.
#
#   (A) clade slopes after dropping 30% of the mapped species, with and without the
#       environmental covariates -- how strong each association is and how much it
#       depends on which species happen to be mapped.
#   (B) the UNCONSTRAINED clade response curves, clade alone -- what the association
#       actually looks like. This panel is the only place in Model 2 where a negative
#       or humped response can appear: (A) is measured on the monotone-constrained fit
#       and is >= 0 by construction, so it ranks strength and cannot show direction.
#
# Formatting follows the other manuscript figures (02_descriptive.R): 15 x 13 in at
# 170 dpi on white, theme_minimal(10), bold "A)" / "B)" tags in the title slot,
# no panel grid. Each panel keeps its key below it; no plot captions -- the prose
# explanation belongs in the manuscript caption, not burned into the image.
#
# Run AFTER 06_model2.R and 06_model2_jackknife.R -- it only assembles their outputs
# and fits nothing, so it takes seconds and can be re-run freely while styling.
# Writes: outputs/model2/FIG_M2_paper.png
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(readr); library(ggplot2)
  library(tidyr); library(patchwork) })
source("R/config.R")
outdir <- paths$model2
CONTS <- c("Americas","Africa","Asia")
COL   <- c(Americas = "#E1812C", Africa = "#3A7CA5", Asia = "#4C9A4C")
ED    <- c("Pilosa","Cingulata","Litopterna","Notoungulata","Diprotodontia")
ORDER <- c("Artiodactyla","Perissodactyla","Rodentia","Primates","Proboscidea",
           "Pilosa","Cingulata","Litopterna","Notoungulata","Diprotodontia")
lab_ed <- function(c) ifelse(c %in% ED, paste0(c, " (ED)"), c)
# identical to the theme in 02_descriptive.R so the figure sits with the others
theme_set(theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(), plot.title = element_text(face = "bold", size = 11),
        plot.title.position = "plot", axis.line = element_line(colour = "grey30"),
        strip.text = element_text(face = "bold"), legend.position = "bottom",
        legend.text = element_text(size = 8.5)))
KEYS <- c("one species-drop replicate", "median of the 30 replicates",
          "bootstrap mean ± 95% CI")

# ---- (A) species-drop jackknife, both settings -------------------------------
# The diamond is the MEAN of the 200 block-bootstrap draws (= model2_results.csv$slope_std),
# not a single full-data fit, and its bar is that bootstrap's 95% interval. NB comparing
# two clades by whether these bars overlap UNDERSTATES the difference: the resamples are
# shared across clades, so the paired difference has much lower variance than the two
# marginals suggest (21 of 36 Americas pairs separable paired, vs 11 by overlap). Use
# model2_bootstrap_draws.csv for any clade-vs-clade claim.
jk  <- read_csv(file.path(outdir, "model2_jackknife_ranks.csv"), show_col_types = FALSE)
ev  <- read_csv(file.path(outdir, "model2_env_results.csv"), show_col_types = FALSE) |>
  filter(condition == "total") |> mutate(continent = factor(continent, CONTS))
STG  <- c("clade alone", "clade + environment")
ord  <- ev |> arrange(continent, slope_uni)
LEV  <- unique(paste(as.integer(ord$continent), ord$clade))
mkey <- function(cont, clade) factor(paste(as.integer(cont), clade), levels = LEV)
ref <- bind_rows(
  ev |> transmute(continent, clade, setting = STG[1],
                  slope = slope_uni, lo = slope_uni_lo, hi = slope_uni_hi),
  ev |> transmute(continent, clade, setting = STG[2],
                  slope = slope_env, lo = slope_env_lo, hi = slope_env_hi)) |>
  mutate(setting = factor(setting, STG), key = mkey(continent, clade))
rep30 <- jk |> filter(frac == 0.3) |> mutate(continent = factor(continent, CONTS)) |>
  tidyr::pivot_longer(c(slope_uni, slope_env), names_to = "setting", values_to = "slope") |>
  mutate(setting = factor(ifelse(setting == "slope_uni", STG[1], STG[2]), STG),
         key = mkey(continent, clade))
med <- rep30 |> group_by(continent, key, setting) |>
  summarise(slope = median(slope), .groups = "drop")
# Share of explained variance that cannot be attributed to clade or environment
# individually (shared / r2_full from the variance partition). This is a split of
# VARIANCE, not the gap between the two diamonds on a row, which is a change in the slope
# ESTIMATE -- they answer different questions and may disagree. 0% means nothing is
# unattributable: the overlap is either within the fold-seed spread or negative, i.e.
# the two predictors are complementary rather than overlapping.
cf <- read_csv(file.path(outdir, "model2_diagnostics.csv"), show_col_types = FALSE) |>
  filter(condition == "total") |>
  mutate(continent = factor(continent, CONTS),
         txt = ifelse(shared > 0 & shared > 2*shared_sd,
                      sprintf("%.0f%%", 100*shared/r2_full), "0%"),
         setting = factor(STG[2], STG), key = mkey(continent, clade))
xcf <- max(c(ref$slope, rep30$slope), na.rm = TRUE) * 1.16
pA <- ggplot(mapping = aes(y = key, colour = continent)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
  geom_point(data = rep30, aes(x = slope, shape = KEYS[1]), alpha = .35, size = 1.3,
             position = position_jitter(height = .16, width = 0)) +
  geom_point(data = med, aes(x = slope, shape = KEYS[2]), size = 4, stroke = 1.1) +
  geom_errorbarh(data = ref, aes(xmin = lo, xmax = hi), height = 0,
                 linewidth = .45, alpha = .9) +
  geom_point(data = ref, aes(x = slope, shape = KEYS[3]), size = 2.5, fill = "white", stroke = .8) +
  geom_text(data = cf, aes(x = xcf, y = key, label = txt), hjust = 1, size = 2.5,
            colour = "grey25", inherit.aes = FALSE, show.legend = FALSE) +
  expand_limits(x = xcf) +
  scale_shape_manual(values = setNames(c(16, 124, 23), KEYS), breaks = KEYS, name = NULL) +
  scale_colour_manual(values = COL, guide = "none") +
  scale_y_discrete(labels = function(x) sub("^[0-9]+ ", "", x)) +
  guides(shape = guide_legend(override.aes = list(alpha = 1, size = c(1.8, 4, 2.5)))) +
  facet_grid(continent ~ setting, scales = "free_y", space = "free_y") +
  labs(x = "Standardised slope (SDs of relative richness per SD of clade richness)",
       y = NULL, title = "A)")

# ---- (B) unconstrained response curves, clade alone --------------------------
# Single-predictor models, so the curve is the model's prediction over the observed
# range -- not an ALE (ALE is used only for the 5-feature environment setting) and not
# a partial dependence (there is nothing to marginalise over). Values are absolute
# predicted relative richness, NOT centred on zero as an ALE would be.
su <- read_csv(file.path(outdir, "model2_curves.csv"), show_col_types = FALSE) |>
  filter(condition == "total")
cl <- intersect(ORDER, unique(su$clade))
su$lab <- factor(lab_ed(su$clade), levels = lab_ed(cl))
pB <- ggplot(su) +
  geom_ribbon(aes(x, ymin = y_lo, ymax = y_hi, fill = continent), alpha = .13) +
  geom_line(aes(x, y, colour = continent), linewidth = .7) +
  scale_colour_manual(values = COL, name = NULL) + scale_fill_manual(values = COL, guide = "none") +
  facet_wrap(~lab, scales = "free", ncol = 5) +
  labs(x = "Clade richness", y = "Predicted relative richness", title = "B)")

ggsave(file.path(outdir, "FIG_M2_paper.png"),
       (pA / pB) + plot_layout(heights = c(1.2, 1)),
       width = 15, height = 13.0, dpi = 170, bg = "white")
message("DONE -> ", file.path(outdir, "FIG_M2_paper.png"))
