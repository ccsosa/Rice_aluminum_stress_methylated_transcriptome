require(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# ── Load methylation data ─────────────────────────────────────────────────────
load_mC <- function(path, context) {
  dt <- fread(path)[, c(1:15)]
  dt$diff_AZU <- 
    rowMeans(cbind(dt$AC1, dt$AC2, dt$AC3)) -
    rowMeans(cbind(dt$AT1, dt$AT2, dt$AT3))
  dt$diff_BGI <- 
    rowMeans(cbind(dt$BC1, dt$BC2, dt$BC3)) -
    rowMeans(cbind(dt$BT1, dt$BT2, dt$BT3))
  dt <- dt[, c("window_pos", "window_end", "chr", "diff_AZU", "diff_BGI")]
  dt$Context <- context
  return(dt)
}

base_path <- "D:/JJG_OCT25/Data/Epigenome/Number_mCs_per_window100k/"
dmr_path  <- "D:/JJG_OCT25/Data/Epigenome/DMRs_annotated/"

mC_CG  <- load_mC(paste0(base_path, "mCs_windows_sativa_met_CG.csv"),  "CG")
mC_CHG <- load_mC(paste0(base_path, "mCs_windows_sativa_met_CHG.csv"), "CHG")
mC_CHH <- load_mC(paste0(base_path, "mCs_windows_sativa_met_CHH.csv"), "CHH")

mC_all <- rbind(mC_CG, mC_CHG, mC_CHH)

# ── Load DMR data ─────────────────────────────────────────────────────────────
load_dmr <- function(path, source, context) {
  dt <- fread(path)
  dt <- dt[, c("seqnames", "start", "end", "Status")]
  dt$source  <- source
  dt$Context <- context
  dt$chr <- sub("chr0*(\\d+)", "chr\\1", dt$seqnames)
  return(dt)
}

dmr_list <- list(
  load_dmr(paste0(dmr_path, "AZU_CG_AC_AT_StateCalls-filtered.txt_Annotated"),  "AZU", "CG"),
  load_dmr(paste0(dmr_path, "AZU_CHG_AC_AT_StateCalls-filtered.txt_Annotated"), "AZU", "CHG"),
  load_dmr(paste0(dmr_path, "AZU_CHH_AC_AT_StateCalls-filtered.txt_Annotated"), "AZU", "CHH"),
  load_dmr(paste0(dmr_path, "BGI_CG_BC_BT_StateCalls-filtered.txt_Annotated"),  "BGI", "CG"),
  load_dmr(paste0(dmr_path, "BGI_CHG_BC_BT_StateCalls-filtered.txt_Annotated"), "BGI", "CHG"),
  load_dmr(paste0(dmr_path, "BGI_CHH_BC_BT_StateCalls-filtered.txt_Annotated"), "BGI", "CHH")
)
dmr_all <- rbindlist(dmr_list) %>%
  mutate(mid_mb = (start + end) / 2 / 1e6)

# ── Colours ───────────────────────────────────────────────────────────────────
source_colors <- c(AZU = "#D73027", BGI = "#4575B4")
dmr_colors    <- c(hyper = "#E08214", hypo = "#1A9641")

# ── Reshape methylation data ──────────────────────────────────────────────────
df_long <- mC_all %>%
  pivot_longer(cols = c(diff_AZU, diff_BGI),
               names_to = "source", values_to = "diff") %>%
  mutate(source    = sub("diff_", "", source),
         pos_exact = window_pos)

# ── DMR panel: AZU (y = +1 hyper, +2 hypo) / BGI (y = -1 hyper, -2 hypo) ────
make_dmr_panel <- function(d_dmr_azu, d_dmr_bgi, x_range, chrom) {
  
  # y positions:  AZU row at +1.5, BGI row at -1.5
  # arrow direction encodes hyper (up) / hypo (down) within each row
  # We use yend offset of ±0.6 from the row centre
  azu_hyper <- d_dmr_azu %>% filter(Status == "hyper")
  azu_hypo  <- d_dmr_azu %>% filter(Status == "hypo")
  bgi_hyper <- d_dmr_bgi %>% filter(Status == "hyper")
  bgi_hypo  <- d_dmr_bgi %>% filter(Status == "hypo")
  
  y_azu <- 1.5;  y_bgi <- -1.5
  off   <- 0.6   # arrow length offset
  
  ggplot() +
    # ── AZU row ──────────────────────────────────────────────────────────────
    {if (nrow(azu_hyper) > 0)
      geom_segment(data = azu_hyper,
                   aes(x = mid_mb, xend = mid_mb,
                       y = y_azu, yend = y_azu + off),
                   color = dmr_colors["hyper"],
                   arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
                   linewidth = 0.5, alpha = 0.85)} +
    {if (nrow(azu_hypo) > 0)
      geom_segment(data = azu_hypo,
                   aes(x = mid_mb, xend = mid_mb,
                       y = y_azu, yend = y_azu - off),
                   color = dmr_colors["hypo"],
                   arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
                   linewidth = 0.5, alpha = 0.85)} +
    # ── BGI row ──────────────────────────────────────────────────────────────
    {if (nrow(bgi_hyper) > 0)
      geom_segment(data = bgi_hyper,
                   aes(x = mid_mb, xend = mid_mb,
                       y = y_bgi, yend = y_bgi + off),
                   color = dmr_colors["hyper"],
                   arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
                   linewidth = 0.5, alpha = 0.85)} +
    {if (nrow(bgi_hypo) > 0)
      geom_segment(data = bgi_hypo,
                   aes(x = mid_mb, xend = mid_mb,
                       y = y_bgi, yend = y_bgi - off),
                   color = dmr_colors["hypo"],
                   arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
                   linewidth = 0.5, alpha = 0.85)} +
    # ── Row separator & labels ────────────────────────────────────────────────
    geom_hline(yintercept = 0,    linewidth = 0.4, color = "grey60", linetype = "dashed") +
    geom_hline(yintercept = y_azu, linewidth = 0.3, color = "grey80") +
    geom_hline(yintercept = y_bgi, linewidth = 0.3, color = "grey80") +
    # annotate("text", x = x_range[1], y = y_azu, label = "AZU",
    #          hjust = -0.3, vjust = -2, fontface = "bold",
    #          size = 3.2, color = source_colors["AZU"]) +
    # annotate("text", x = x_range[1], y = y_bgi, label = "BGI",
    #          hjust = -0.3, vjust = -2, fontface = "bold",
    #          size = 3.2, color = source_colors["BGI"]) +
    scale_x_continuous(
      limits = x_range,
      labels = function(x) paste0(x, " Mb"),
      expand = expansion(mult = 0.01)
    ) +
    scale_y_continuous(
      limits = c(y_bgi - 1, y_azu + 1),
      breaks = c(y_azu + off/2, y_azu - off/2,
                 y_bgi + off/2, y_bgi - off/2),
      labels = c("hyper (AZU)", "hypo (AZU)", "hyper (BGI)", "hypo (BGI)")
    ) +
    labs(x = paste0("Genomic position · ", chrom), y = "DMR status") +
    theme_bw(base_size = 12) +
    theme(
      panel.grid    = element_blank(),
      axis.text.x   = element_text(angle = 30, hjust = 1, size = 9),
      axis.ticks.y  = element_blank(),
      axis.text.y   = element_text(face = "bold", size = 9,
                                   colour = rep(c(dmr_colors["hyper"],
                                                  dmr_colors["hypo"]), 2)),
      plot.margin   = margin(0, 5, 5, 5)
    )
}

# ── Main plot function ────────────────────────────────────────────────────────
plot_chrom_ctx <- function(chrom, ctx) {
  
  d_line    <- df_long %>% filter(chr == chrom, Context == ctx)
  d_dmr_azu <- dmr_all %>% filter(chr == chrom, Context == ctx, source == "AZU")
  d_dmr_bgi <- dmr_all %>% filter(chr == chrom, Context == ctx, source == "BGI")
  
  if (nrow(d_line) == 0) {
    message("No methylation data for ", chrom, " ", ctx, " — skipping"); return(NULL)
  }
  if (nrow(d_dmr_azu) == 0 && nrow(d_dmr_bgi) == 0) {
    message("No DMR data for ", chrom, " ", ctx, " — skipping");          return(NULL)
  }
  
  x_range <- range(d_line$pos_exact) / 1e6
  
  # ── Panel 1: line plot ──────────────────────────────────────────────────────
  p_line <- ggplot() +
    geom_hline(yintercept = 0, linewidth = 0.4, color = "black") +
    geom_line(
      data = d_line,
      aes(x = pos_exact / 1e6, y = diff, color = source, group = source),
      linewidth = 0.6, alpha = 0.85
    ) +
    geom_point(
      data = d_line,
      aes(x = pos_exact / 1e6, y = diff, color = source),
      size = 0.9, alpha = 0.65
    ) +
    scale_color_manual(values = source_colors, name = "") +
    scale_x_continuous(
      limits = x_range,
      labels = function(x) paste0(x, " Mb"),
      expand = expansion(mult = 0.01)
    ) +
    scale_y_continuous(labels = scales::comma) +
    labs(
      subtitle = "Lines: mean diff per 100 kb window  ·  Arrows: AZU & BGI DMRs",
      x = NULL,
      y = "Mean methylation\ndifference (bp) [Control - Treatment]"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.subtitle    = element_text(size = 9, color = "grey40"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey92"),
      axis.text.x      = element_blank(),
      axis.ticks.x     = element_blank(),
      legend.position  = "top",
      legend.key.width = unit(1.1, "cm"),
      plot.margin      = margin(5, 5, 0, 5)
    )
  
  # ── Panel 2: DMR arrows (AZU + BGI) ────────────────────────────────────────
  p_dmr <- make_dmr_panel(d_dmr_azu, d_dmr_bgi, x_range, chrom)
  
  # ── Combine ─────────────────────────────────────────────────────────────────
  p_line / p_dmr +
    plot_layout(heights = c(4, 1.8), axes = "collect_x")
}

# ── Loop over all chromosomes × contexts ──────────────────────────────────────
chroms   <- sort(unique(df_long$chr))
contexts <- c("CG", "CHG", "CHH")

base_output_dir <- "D:/JJG_OCT25/Figures/diff_meth_plots/"

for (ctx in contexts) {
  ctx_dir <- file.path(base_output_dir, ctx)
  dir.create(ctx_dir, showWarnings = FALSE, recursive = TRUE)
  
  for (chrom in chroms) {
    p <- plot_chrom_ctx(chrom, ctx)
    if (is.null(p)) next
    
    outfile <- file.path(ctx_dir,
                         paste0("diff_meth_AZU_BGI_", chrom, "_", ctx, ".png"))
    ggsave(outfile, plot = p,
           width = 16, height = 7, dpi = 200, bg = "white")
    cat("Saved:", outfile, "\n")
  }
}

cat("\nDone. Total plots:", length(chroms) * length(contexts), "\n")