
CHH <- read.csv("D:/JJG_OCT25/Data/Epigenome/Number_mCs_per_window100k/CHH.csv")
CHG <- read.csv("D:/JJG_OCT25/Data/Epigenome/Number_mCs_per_window100k/CHG.csv")
CG <- read.csv("D:/JJG_OCT25/Data/Epigenome/Number_mCs_per_window100k/CG.csv")


library(ggplot2)
library(ggpubr)
library(ggrepel)
# boxplot(as.numeric(CG$LFC_AZ),as.numeric(CG$LFC_BGI))

uniq <- unique(CHH$chr)
CHR_WILC <- list()

for(i in 1:length(uniq)){
   # i <- 8
  CHH_i <- CHH[which(CHH$chr==uniq[[i]]),]
  CHG_i <- CHG[which(CHG$chr==uniq[[i]]),]
  CG_i <- CG[which(CG$chr==uniq[[i]]),]
  
  # Construcción corregida del data.frame (nombres de columna arreglados)
  x <- data.frame(
    CHR           = uniq[[i]],
    WP            = CHH_i$window_pos,
    WE            = CHH_i$window_end,
    CHH_LFC_AZU   = as.numeric(CHH_i$LFC_AZ),
    CHH_LFC_BGI   = as.numeric(CHH_i$LFC_BGI),
    CHG_LFC_AZU   = as.numeric(CHG_i$LFC_AZ),
    CHG_LFC_BGI   = as.numeric(CHG_i$LFC_BGI),
    CG_LFC_AZU    = as.numeric(CG_i$LFC_AZ),
    CG_LFC_BGI    = as.numeric(CG_i$LFC_BGI)
  )
  x$window_mid <- (x$WP + x$WE) / 2
  
  for(j in 4:10){
  x[,j][which(abs(x[,j])==0)] <- NA
  }
  
  # Asegúrate de que las ventanas estén ordenadas por posición
  x <- x[order(x$window_mid), ]
  
  # Formato largo
  x_long <- data.frame(
    window_mid = rep(x$window_mid, times = 6),
    LFC = c(x$CHH_LFC_AZU, x$CHH_LFC_BGI,
            x$CHG_LFC_AZU, x$CHG_LFC_BGI,
            x$CG_LFC_AZU,  x$CG_LFC_BGI),
    Context = rep(c("CHH", "CHH", "CHG", "CHG", "CG", "CG"), each = nrow(x)),
    Variety = rep(c("AZU", "BGI", "AZU", "BGI", "AZU", "BGI"), each = nrow(x))
  )
  
  # Asegurar tipos correctos
  x_long$window_pos <- x$WP
  x_long$window_end <- x$WE
  x_long$window_mid <- as.numeric(x_long$window_mid)
  x_long$LFC        <- as.numeric(x_long$LFC)
  x_long$Context    <- factor(x_long$Context, levels = c("CG", "CHG", "CHH"))
  x_long$Variety    <- factor(x_long$Variety, levels = c("AZU", "BGI"))
  x_long$LFC[which(is.infinite(x_long$LFC))] <- NA
  
  # Umbral para considerar un punto "outlier" a etiquetar
  outlier_threshold <- 1
  
  x_long$is_outlier <- abs(x_long$LFC) >= outlier_threshold
  x_long$label <- ifelse(x_long$is_outlier,
                         paste0(round(x_long$window_pos/ 1e6,2),"-",
                                round(x_long$window_end/ 1e6,2)," Mb"),
                         NA)
  p <- ggplot(x_long, aes(x = window_mid, y = LFC, fill = Variety)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
    geom_col(position = "identity", alpha = 0.6, width = bar_width, na.rm = TRUE) +
    geom_text_repel(
      data = subset(x_long, is_outlier),
      aes(x = window_mid, y = LFC, label = label, color = Variety),
      size = 2.8,
      show.legend = FALSE,
      max.overlaps = 20,
      segment.size = 0.2,
      min.segment.length = 0
    ) +
    facet_wrap(~ Context, ncol = 3) +
    scale_fill_manual(values = c("AZU" = "#1f78b4", "BGI" = "#e31a1c")) +
    scale_color_manual(values = c("AZU" = "#1f78b4", "BGI" = "#e31a1c")) +
    scale_x_continuous(labels = scales::label_number(scale = 1e-6, suffix = "Mb")) +
    coord_cartesian(ylim = c(-3, 3)) +
    labs(x = paste0("Genomic position (", uniq[[i]], ")"),
         y = "LFC (methylation)",
         fill = "Variety") +
    theme_pubr(legend = "top") +
    theme(strip.background = element_rect(fill = "grey90"),
          axis.text.x = element_text(angle = 0))
  
  # print(p)
  
  ggsave(paste0("D:/JJG_OCT25/RESULTS_MET/LFC_barplot_", uniq[[i]], ".png"),
         plot = p, width = 12, height = 5, dpi = 300)
  
  p_CHH_i <- wilcox.test(as.numeric(CHH_i$LFC_AZ),as.numeric(CHH_i$LFC_BGI),paired = T)
  p_CHG_i <- wilcox.test(as.numeric(CHG_i$LFC_AZ),as.numeric(CHG_i$LFC_BGI),paired = T)
  p_CG_i <- wilcox.test(as.numeric(CG_i$LFC_AZ),as.numeric(CG_i$LFC_BGI),paired = T)
  
  CHR_WILC[[i]] <- data.frame(CHR = uniq[[i]],
             n_AZU_CHH = length(na.omit(as.numeric(CHH_i$LFC_AZ))),
             n_BGI_CHH = length(na.omit(as.numeric(CHH_i$LFC_BGI))),
             n_AZU_CHG = length(na.omit(as.numeric(CHG_i$LFC_AZ))),
             n_BGI_CHG = length(na.omit(as.numeric(CHG_i$LFC_BGI))),
             n_AZU_CG = length(na.omit(as.numeric(CG_i$LFC_AZ))),
             n_BGI_CG = length(na.omit(as.numeric(CG_i$LFC_BGI))),
             #######################
             LFC_AZU_med_CHH = median(as.numeric(CHH_i$LFC_AZ),na.rm = T),
             LFC_BGI_med_CHH = median(as.numeric(CHH_i$LFC_BGI),na.rm = T),
             LFC_AZU_mad_CHH = mad(as.numeric(CHH_i$LFC_AZ),na.rm = T),
             LFC_BGI_mad_CHH = mad(as.numeric(CHH_i$LFC_BGI),na.rm = T),
             LFC_AZU_med_CHG = median(as.numeric(CHG_i$LFC_AZ),na.rm = T),
             LFC_BGI_med_CHG = median(as.numeric(CHG_i$LFC_BGI),na.rm = T),
             LFC_AZU_mad_CHG = mad(as.numeric(CHG_i$LFC_AZ),na.rm = T),
             LFC_BGI_mad_CHG = mad(as.numeric(CHG_i$LFC_BGI),na.rm = T),
             LFC_AZU_med_CG = median(as.numeric(CG_i$LFC_AZ),na.rm = T),
             LFC_BGI_med_CG = median(as.numeric(CG_i$LFC_BGI),na.rm = T),
             LFC_AZU_mad_CG = mad(as.numeric(CG_i$LFC_AZ),na.rm = T),
             LFC_BGI_mad_CG = mad(as.numeric(CG_i$LFC_BGI),na.rm = T),
             #######################
             # LFC_AZU_avg_CHH = mean(as.numeric(CHH_i$LFC_AZ),na.rm = T),
             # LFC_BGI_avg_CHH = mean(as.numeric(CHH_i$LFC_BGI),na.rm = T),
             # LFC_AZU_sd_CHH = sd(as.numeric(CHH_i$LFC_AZ),na.rm = T),
             # LFC_BGI_sd_CHH = sd(as.numeric(CHH_i$LFC_BGI),na.rm = T),
             # LFC_AZU_avg_CHG = mean(as.numeric(CHG_i$LFC_AZ),na.rm = T),
             # LFC_BGI_avg_CHG = mean(as.numeric(CHG_i$LFC_BGI),na.rm = T),
             # LFC_AZU_sd_CHG = sd(as.numeric(CHG_i$LFC_AZ),na.rm = T),
             # LFC_BGI_sd_CHG = sd(as.numeric(CHG_i$LFC_BGI),na.rm = T),
             # LFC_AZU_avg_CG = mean(as.numeric(CG_i$LFC_AZ),na.rm = T),
             # LFC_BGI_avg_CG = mean(as.numeric(CG_i$LFC_BGI),na.rm = T),
             # LFC_AZU_sd_CG = sd(as.numeric(CG_i$LFC_AZ),na.rm = T),
             # LFC_BGI_sd_CG = sd(as.numeric(CG_i$LFC_BGI),na.rm = T),
             #######################
             LFC_med_dir_CHH = median(as.numeric(CHH_i$LFC_AZ),na.rm = T)
                              -
                       median(as.numeric(CHH_i$LFC_BGI),na.rm = T),
             LFC_med_dir_CHG = median(as.numeric(CHG_i$LFC_AZ),na.rm = T)
             -
               median(as.numeric(CHG_i$LFC_BGI),na.rm = T),

             LFC_med_dir_CG = median(as.numeric(CG_i$LFC_AZ),na.rm = T)
             -
               median(as.numeric(CG_i$LFC_BGI),na.rm = T),
             #######################
             # LFC_avg_dir_CHH = mean(as.numeric(CHH_i$LFC_AZ),na.rm = T)
             # - 
             #   mean(as.numeric(CHH_i$LFC_BGI),na.rm = T),
             # LFC_avg_dir_CHG = mean(as.numeric(CHG_i$LFC_AZ),na.rm = T)
             # - 
             #   mean(as.numeric(CHG_i$LFC_BGI),na.rm = T),
             # 
             # LFC_avg_dir_CG = mean(as.numeric(CG_i$LFC_AZ),na.rm = T)
             # - 
             #   mean(as.numeric(CG_i$LFC_BGI),na.rm = T),
             #######################
             CHH = p_CHH_i$p.value,
             CHG = p_CHG_i$p.value,
             CG = p_CG_i$p.value)
}

CHR_WILC <- do.call(rbind,CHR_WILC)

CHR_WILC$CHH_p_adjust <- p.adjust(CHR_WILC$CHH,method = "BH")
CHR_WILC$CHG_p_adjust <- p.adjust(CHR_WILC$CHG,method = "BH")
CHR_WILC$CG_p_adjust <- p.adjust(CHR_WILC$CG,method = "BH")
write.csv(CHR_WILC,"D:/JJG_OCT25/Data/Epigenome/Number_mCs_per_window100k/WILCOXON_mC.CSV",row.names = F)


# wilcox.test(as.numeric(CHH$LFC_AZ),as.numeric(CHH$LFC_BGI),paired = T)
# wilcox.test(as.numeric(CHG$LFC_AZ),as.numeric(CHG$LFC_BGI),paired = T)
# wilcox.test(as.numeric(CG$LFC_AZ),as.numeric(CG$LFC_BGI),paired = T)


