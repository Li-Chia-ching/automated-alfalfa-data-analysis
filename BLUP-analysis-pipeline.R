# =============================================================================
# 项目：苜蓿多年多点数据分析（GGE 双标图 + BLUP + 遗传力）系统
# 数据文件：Alfalfa_5Year_Master_Dataset.csv（需放在 R 项目根目录）
# 功能：自动完成数据清洗、分析、可视化现代化输出、日志记录
# =============================================================================

# ------------------------------ 1. 加载包 ------------------------------------
# 检查并安装必要的包
required_packages <- c("lme4", "lmerTest", "ggplot2", "dplyr", "tidyr", 
                       "tibble", "here", "conflicted", "scales", "ggrepel")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(lme4)
library(lmerTest)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(here)          # 路径管理
library(conflicted)    # 解决函数冲突
library(scales)        # 用于坐标轴格式化
library(ggrepel)       # 用于标签防重叠（GGE图中选用）

# 优先使用 tidyverse 函数
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("mutate", "dplyr")
conflict_prefer("summarise", "dplyr")
conflict_prefer("arrange", "dplyr")
conflict_prefer("lmer", "lmerTest")   # 防止 lmer found in 2 packages

# ------------------------------ 2. 设置输出目录与日志 ------------------------
result_dir <- here("BLUP_Analysis_Results")
if (!dir.exists(result_dir)) dir.create(result_dir)

# 生成时间戳用于日志文件名
timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
log_file <- file.path(result_dir, paste0("Analysis_Log_", timestamp, ".txt"))

# 开启 sink，将后续所有控制台输出同时定向到文件和屏幕
sink(log_file, split = TRUE)

cat("============================================================\n")
cat("苜蓿多年多点数据自动分析系统\n")
cat("分析时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("结果输出目录：", result_dir, "\n")
cat("============================================================\n\n")

# ------------------------------ 3. 读取与预处理数据 --------------------------
file_path <- here("Alfalfa_5Year_Master_Dataset.csv")
if (!file.exists(file_path)) {
  cat("错误：数据文件未找到，请检查路径：", file_path, "\n")
  sink() # 关闭日志
  stop("文件未找到。")
}

df_raw <- read.csv(file_path, stringsAsFactors = FALSE)

cat("--- 数据概况 ---\n")
glimpse(df_raw)

# 将关键列转为因子
df_raw <- df_raw %>%
  mutate(
    Year = as.factor(Year),
    Line = as.factor(Line),
    Rep  = as.factor(Rep)
  )

# 定义性状列表（可根据实际数据调整）
traits <- c("Summer_Height_cm", "Fall_Dormancy_Height_cm", 
            "Fresh_Weight_kg", "Dry_Weight_kg")

cat("\n--- 性状缺失值检查 ---\n")
for (t in traits) {
  if (t %in% names(df_raw)) {
    cat(t, "缺失值个数：", sum(is.na(df_raw[[t]])), "\n")
  } else {
    cat("警告：性状", t, "在数据中未找到。\n")
  }
}

# ------------------------------ 4. GGE 双标图分析（针对产量 Dry_Weight_kg）----
cat("\n========== GGE 双标图分析（性状：Dry_Weight_kg）==========\n")

if ("Dry_Weight_kg" %in% names(df_raw)) {
  # 4.1 聚合产量均值
  yield_met <- df_raw %>%
    filter(!is.na(Dry_Weight_kg)) %>%
    group_by(Line, Year) %>%
    summarise(Yield = mean(Dry_Weight_kg, na.rm = TRUE), .groups = "drop")
  
  # 4.2 构建矩阵并填补缺失
  ge_mat <- yield_met %>%
    pivot_wider(names_from = Year, values_from = Yield) %>%
    column_to_rownames("Line") %>%
    as.matrix()
  
  for(i in 1:ncol(ge_mat)){
    col_mean <- mean(ge_mat[, i], na.rm = TRUE)
    ge_mat[is.na(ge_mat[, i]), i] <- col_mean
  }
  
  # 4.3 中心化与 SVD
  env_means <- colMeans(ge_mat)
  ge_centered <- sweep(ge_mat, 2, env_means, FUN = "-")
  svd_res <- svd(ge_centered)
  pc_var <- svd_res$d^2 / sum(svd_res$d^2) * 100
  PC1_var <- pc_var[1]
  PC2_var <- pc_var[2]
  cat(sprintf("PC1 解释变异: %.1f%%, PC2 解释变异: %.1f%%, 合计: %.1f%%\n",
              PC1_var, PC2_var, PC1_var + PC2_var))
  
  # 4.4 计算得分
  G_scores <- data.frame(
    Line = rownames(ge_centered),
    PC1 = svd_res$u[,1] * sqrt(svd_res$d[1]),
    PC2 = svd_res$u[,2] * sqrt(svd_res$d[2]),
    stringsAsFactors = FALSE
  )
  E_scores <- data.frame(
    Year = colnames(ge_centered),
    PC1 = svd_res$v[,1] * sqrt(svd_res$d[1]),
    PC2 = svd_res$v[,2] * sqrt(svd_res$d[2]),
    stringsAsFactors = FALSE
  )
  
  # 4.5 计算平均产量和 ASV
  ss_ipca1 <- svd_res$d[1]^2 
  ss_ipca2 <- svd_res$d[2]^2
  weight <- ss_ipca1 / ss_ipca2
  G_scores <- G_scores %>%
    mutate(ASV = sqrt(weight * PC1^2 + PC2^2))
  
  mean_yield <- yield_met %>%
    group_by(Line) %>%
    summarise(Mean_Yield = mean(Yield, na.rm = TRUE))
  
  variety_info <- mean_yield %>%
    left_join(G_scores %>% select(Line, ASV, PC1, PC2), by = "Line")
  
  # 4.6 筛选高产稳产 Top 10
  top_stable <- variety_info %>%
    arrange(desc(Mean_Yield)) %>%
    slice_head(n = 20) %>%
    arrange(ASV) %>%
    slice_head(n = 10)
  
  cat("\n===== GGE 推荐的高产稳产品种 Top 10 =====\n")
  print(top_stable %>% select(Line, Mean_Yield, ASV))
  
  # -------------------- GGE 可视化现代化 (Who-Won-Where) --------------------
  # 计算凸包（实现Who-Won-Where）
  chull_indices <- chull(G_scores$PC1, G_scores$PC2)
  G_hull <- G_scores[chull_indices, ]
  G_hull <- rbind(G_hull, G_hull[1,]) # 闭合凸包
  
  # 计算扇区边界线（原点垂直于凸包各边）
  get_lines <- function(hull_df) {
    lines_df <- data.frame(x=numeric(), y=numeric(), xend=numeric(), yend=numeric())
    n <- nrow(hull_df) - 1
    mult <- max(abs(hull_df$PC1), abs(hull_df$PC2)) * 1.5 # 延长线段
    for (i in 1:n) {
      x1 <- hull_df$PC1[i]; y1 <- hull_df$PC2[i]
      x2 <- hull_df$PC1[i+1]; y2 <- hull_df$PC2[i+1]
      slope <- (y2 - y1) / (x2 - x1)
      perp_slope <- -1 / slope
      angle <- atan(perp_slope)
      xend <- mult * cos(angle) * ifelse(x2 > x1, 1, -1)
      yend <- mult * sin(angle) * ifelse(x2 > x1, 1, -1)
      # 处理特殊象限
      if((y2-y1)>0 & xend<0 & slope<0) {xend<--xend; yend<--yend}
      if((y2-y1)<0 & xend>0 & slope<0) {xend<--xend; yend<--yend}
      lines_df <- rbind(lines_df, data.frame(x=0, y=0, xend=xend, yend=yend))
    }
    return(lines_df)
  }
  sector_lines <- get_lines(G_hull)
  
  # 绘制现代化的 Who-Won-Where 双标图
  p_gge <- ggplot() +
    # 绘制扇区参考线
    geom_segment(data = sector_lines, aes(x=x, y=y, xend=xend, yend=yend), 
                 linetype="dotted", color="grey70") +
    # 绘制凸包
    geom_path(data = G_hull, aes(x = PC1, y = PC2), color = "grey50", linetype="dashed") +
    # 品种点与标签
    geom_point(data = G_scores, aes(x = PC1, y = PC2), 
               color = "steelblue", alpha = 0.6, size = 1.8) +
    geom_text_repel(data = G_scores, aes(x = PC1, y = PC2, label = Line), 
                    size = 2.5, color = "steelblue", max.overlaps = 15) +
    # 环境向量
    geom_segment(data = E_scores, aes(x = 0, y = 0, xend = PC1, yend = PC2),
                 arrow = arrow(length = unit(0.2, "cm")), 
                 color = "#D95F02", linewidth = 0.8) + # 使用更清晰的橙红色
    geom_text(data = E_scores, aes(x = PC1, y = PC2, label = Year),
              color = "#D95F02", fontface = "bold", size = 3.5, vjust = 1.2) +
    # 高亮 Top 10 品种
    geom_point(data = top_stable, aes(x = PC1, y = PC2), 
               shape = 23, color = "black", fill = "#E7298A", size = 3, stroke=1) + # 玫红色高亮
    # 参考线与坐标轴
    geom_hline(yintercept = 0, color = "grey40") +
    geom_vline(xintercept = 0, color = "grey40") +
    scale_x_continuous(expand = expansion(mult = 0.1)) +
    scale_y_continuous(expand = expansion(mult = 0.1)) +
    labs(title = "GGE Biplot: Who-Won-Where for Dry Weight Yield",
         subtitle = sprintf("PC1: %.1f%%, PC2: %.1f%% | Pink diamonds: GGE Top 10", PC1_var, PC2_var),
         x = paste0("Primary Component 1 (", round(PC1_var, 1), "%)"),
         y = paste0("Primary Component 2 (", round(PC2_var, 1), "%)")) +
    theme_bw(base_size = 11) + # 统一使用 theme_bw
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5),
          panel.grid.minor = element_blank())
  
  # 保存 GGE 结果
  ggsave(file.path(result_dir, "GGE_Biplot_WhoWonWhere.png"), p_gge, width = 9, height = 7.5, dpi = 300)
  write.csv(variety_info, file.path(result_dir, "Variety_Yield_and_Stability_ASV.csv"), row.names = FALSE)
  write.csv(top_stable, file.path(result_dir, "GGE_Recommended_Top10.csv"), row.names = FALSE)
  
} else {
  cat("\n跳过 GGE 分析：数据中缺少 'Dry_Weight_kg' 列。\n")
}

# ------------------------------ 5. BLUP 与遗传力分析（所有性状）---------------
cat("\n========== BLUP 与广义遗传力分析 ==========\n")

# 准备数据（确保因子类型）
df_met <- df_raw

# 存储结果用于后续绘图
all_blup_list <- list()
h2_results <- data.frame(Trait = character(),
                         V_G = numeric(), V_GE = numeric(), V_e = numeric(),
                         H2 = numeric(), stringsAsFactors = FALSE)

# 循环分析每个性状
for (trait in traits) {
  cat("\n-------------------- 分析性状：", trait, "--------------------\n")
  
  if (!(trait %in% names(df_met)) || all(is.na(df_met[[trait]]))) {
    cat("跳过：性状缺失或无有效数据。\n")
    next
  }
  
  # 拟合混合模型
  formula_str <- paste(trait, "~ Year + (1 | Line) + (1 | Line:Year) + (1 | Rep)")
  # 使用常用优化器控制
  control_params <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))
  
  cat("正在拟合混合模型...\n")
  model <- tryCatch(
    lmer(as.formula(formula_str), data = df_met, control = control_params),
    error = function(e) {
      cat("模型拟合失败：", e$message, "\n")
      return(NULL)
    }
  )
  if (is.null(model)) next
  
  # 查看 Anova
  cat("固定效应 ANOVA 表:\n")
  print(anova(model))
  
  # 收敛性检查
  if (length(model@optinfo$conv$lme4$messages) > 0) {
    cat("警告：模型可能未收敛：", paste(model@optinfo$conv$lme4$messages, collapse = "; "), "\n")
  } else {
    cat("模型收敛成功。\n")
  }
  
  # -------------------- 提取品种 BLUP 与非对称误差（分位数法） --------------------
  blup_raw <- ranef(model, condVar = TRUE)
  
  # 基础 BLUP 值提取
  blup_base <- as.data.frame(blup_raw$Line) %>%
    rownames_to_column("Line") %>%
    rename(BLUP = `(Intercept)`)
  
  # 计算原始数据的偏态分位数偏差（反映真实环境波动的非对称性）
  # 使用 5% 和 95% 分位数来代表 90% 的核心数据波动范围
  quantile_data <- df_met %>%
    filter(!is.na(.data[[trait]])) %>%
    group_by(Line) %>%
    summarise(
      Raw_Mean = mean(.data[[trait]], na.rm = TRUE),
      Q05 = quantile(.data[[trait]], 0.05, na.rm = TRUE),
      Q95 = quantile(.data[[trait]], 0.95, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      # 计算相对于自身均值的真实非对称上下偏差
      dev_lower = Q05 - Raw_Mean,
      dev_upper = Q95 - Raw_Mean
    )
  
  # 将真实偏差映射到 BLUP 上，计算出非对称误差棒边界
  blup_line <- blup_base %>%
    left_join(quantile_data, by = "Line") %>%
    mutate(
      lower_q = BLUP + dev_lower,
      upper_q = BLUP + dev_upper,
      Trait = trait
    ) %>%
    arrange(desc(BLUP))
  
  all_blup_list[[trait]] <- blup_line
  
  # 计算遗传力 H^2
  var_comps <- as.data.frame(VarCorr(model))
  V_G   <- var_comps$vcov[var_comps$grp == "Line"]
  V_GE  <- var_comps$vcov[var_comps$grp == "Line:Year"]
  V_err <- var_comps$vcov[var_comps$grp == "Residual"]
  
  # 稳健计算环境数和重复数
  temp_data <- df_met %>% filter(!is.na(.data[[trait]]))
  e_count <- nlevels(droplevels(temp_data$Year))
  r_count <- temp_data %>%
    group_by(Line, Year) %>%
    summarise(n = n(), .groups = "drop") %>%
    pull(n) %>%
    mean(na.rm = TRUE)
  
  H2 <- V_G / (V_G + (V_GE / e_count) + (V_err / (e_count * r_count)))
  
  h2_results <- rbind(h2_results, data.frame(
    Trait = trait, V_G = V_G, V_GE = V_GE, V_e = V_err, H2 = H2,
    stringsAsFactors = FALSE
  ))
  
  cat(sprintf("方差分量: V_G = %.4f, V_GE = %.4f, V_e = %.4f\n", V_G, V_GE, V_err))
  cat(sprintf("遗传力分析环境数 e = %d, 平均重复数 r = %.2f\n", e_count, r_count))
  cat(sprintf("广义遗传力 H^2 = %.2f%%\n", H2 * 100))
  
  # -------------------- 现代化 BLUP 可视化 (Top 20 非对称误差棒) --------------------
  top20 <- blup_line %>% slice_head(n = 20)
  
  # 创建现代化的条形图
  p_blup <- ggplot(top20, aes(x = reorder(Line, BLUP), y = BLUP, fill = BLUP)) +
    # 使用 linewidth 替代 size 消除警告
    geom_col(width = 0.75, color = "black", linewidth = 0.2, show.legend = FALSE) + 
    # 应用科学的非对称分位数误差棒 (无需截断，反映真实下限)
    geom_errorbar(aes(ymin = lower_q, ymax = upper_q), 
                  width = 0.25, color = "grey30", linewidth = 0.4) +
    scale_fill_viridis_c(option = "mako", begin = 0.2, end = 0.8) + 
    coord_flip() +
    # 稍微扩展 Y 轴，防止两端的偏态误差棒超出画幅
    scale_y_continuous(expand = expansion(mult = c(0.15, 0.15))) +
    labs(title = paste("Top 20 Varieties for", gsub("_", " ", trait)),
         subtitle = "Bars: Genetic BLUP | Error bars: Asymmetric Phenotypic Range (5th-95th Pctl)",
         x = NULL, y = "BLUP Value with Empirical Quantile Variance") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          panel.grid.major.y = element_blank())
  
  # 保存 BLUP 图片和 CSV
  clean_trait_name <- gsub("/", "-", trait) 
  ggsave(filename = file.path(result_dir, paste0("BLUP_Top20_", clean_trait_name, ".png")),
         plot = p_blup, width = 7, height = 6, dpi = 300)
  write.csv(blup_line, file.path(result_dir, paste0("BLUP_Full_Table_", clean_trait_name, ".csv")), row.names = FALSE)
}

# --- 输出遗传力汇总与图片 ---
cat("\n========== 遗传力汇总表 ==========\n")
print(h2_results)
write.csv(h2_results, file.path(result_dir, "Heritability_Summary.csv"), row.names = FALSE)

# 现代化遗传力条形图
h2_plot <- ggplot(h2_results, aes(x = reorder(Trait, H2), y = H2, fill = H2)) +
  geom_col(width = 0.7, color = "black", size=0.3) +
  # 添加参考线（低、中、高遗传力界限）
  geom_hline(yintercept = c(0.3, 0.6), linetype = "dashed", color = "grey60") +
  scale_fill_gradient2(low = "#D7191C", mid = "#FFFFBF", high = "#1A9641", 
                       midpoint = 0.45, limits=c(0, 1), labels=percent,
                       name = expression(H^2)) +
  coord_flip() +
  scale_y_continuous(labels = percent, limits = c(0, 1), expand = c(0, 0)) +
  labs(title = "Broad-sense Heritability (H²) Comparison",
       x = "Agronomic Traits", y = expression(Heritability~(H^2))) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust=0.5),
        legend.position = "right")

ggsave(file.path(result_dir, "Heritability_Comparison_Barplot.png"), h2_plot, width = 7, height = 4.5, dpi = 300)

# -------------------------- 6. 综合分析：BLUP 与 GGE 结果四象限图 -----------
cat("\n========== 绘制高产稳产综合四象限图 ==========\n")

# 直接从内存对象中提取 Dry_Weight_kg 的 BLUP
if (exists("variety_info") && ("Dry_Weight_kg" %in% names(all_blup_list))) {
  
  blup_yield_raw <- all_blup_list[["Dry_Weight_kg"]]
  blup_yield <- blup_yield_raw %>%
    select(Line, BLUP) %>%
    rename(BLUP_Yield = BLUP)
  
  # 合并遗传增量(BLUP)与GGE计算的ASV
  combined_data <- variety_info %>%
    left_join(blup_yield, by = "Line") %>%
    # 标记是否为 Top 10
    mutate(Highlight = ifelse(Line %in% top_stable$Line, "GGE Top10", "Others"))
  
  # 计算中位数作为象限分界线
  med_blup_yield <- median(combined_data$BLUP_Yield, na.rm = TRUE)
  med_asv <- median(combined_data$ASV, na.rm = TRUE)
  
  # 现代化四象限图
  p_quad <- ggplot(combined_data, aes(x = BLUP_Yield, y = ASV)) +
    # 绘制象限参考线
    geom_hline(yintercept = med_asv, linetype = "longdash", color = "grey60") +
    geom_vline(xintercept = med_blup_yield, linetype = "longdash", color = "grey60") +
    # 添加象限标签
    annotate("text", x = Inf, y = -Inf, label = "High Yield\nStable", hjust=1.1, vjust=-0.5, color="grey50", fontface="italic", size=3.5) +
    annotate("text", x = -Inf, y = Inf, label = "Low Yield\nUnstable", hjust=-0.1, vjust=1.5, color="grey50", fontface="italic", size=3.5) +
    # 绘制点，根据Highlight上色
    geom_point(aes(color = Highlight, shape = Highlight, size = Highlight), alpha = 0.8) +
    # 针对Top10添加标签
    geom_text_repel(data = filter(combined_data, Highlight == "GGE Top10"),
                    aes(label = Line), size = 3, color = "black", box.padding = 0.3) +
    # 自定义色彩、形状和大小
    scale_color_manual(values = c("GGE Top10" = "#E7298A", "Others" = "grey70")) +
    scale_shape_manual(values = c("GGE Top10" = 17, "Others" = 16)) +
    scale_size_manual(values = c("GGE Top10" = 3, "Others" = 1.8)) +
    labs(title = "Multi-trait Selection: Breeding Value (BLUP) vs. Stability (ASV)",
         subtitle = "Dry Weight Yield | Dashed lines indicate Medians",
         x = "Genetic Breeding Value (BLUP_Yield)", 
         y = "AMMI Stability Value (ASV, Lower is more stable)") +
    theme_bw() +
    theme(plot.title = element_text(face = "bold", hjust=0.5),
          plot.subtitle = element_text(hjust=0.5),
          legend.position = "bottom", legend.title = element_blank())
  
  # 保存图形
  ggsave(file.path(result_dir, "Comprehensive_Yield_Stability_Quadrant.png"), p_quad, width = 8.5, height = 7, dpi = 300)
  cat("四象限图已保存。\n")
  
} else {
  cat("跳过综合分析：缺少 Dry_Weight_kg 的 GGE 或 BLUP 数据。\n")
}

cat("\n============================================================\n")
cat("全部分析完成！\n")
cat("请在以下目录查阅所有 CSV 表格、PNG 图片和日志文件：\n")
cat(result_dir, "\n")
cat("============================================================\n")

# 关闭 sink 日志记录
sink()
