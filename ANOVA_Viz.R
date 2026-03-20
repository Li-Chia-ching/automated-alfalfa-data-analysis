# =============================================================================
# 5年紫花苜蓿夏季株高联合方差分析与多重比较（最终严谨版）
# 功能：数据清洗、模型拟合、自动处理秩亏、多重比较、高质量图表输出
# 优化：解决由 data.frame() 自动篡改列名导致的 p 值提取 NULL 崩溃问题
# =============================================================================

# 0. 全局设置 ----------------------------------------------------------------
options(stringsAsFactors = FALSE)
options(contrasts = c("contr.sum", "contr.poly"))
verbose <- FALSE

output_dir <- paste0("ANOVA_Results_", format(Sys.Date(), "%Y%m%d"))
if (!dir.exists(output_dir)) dir.create(output_dir)
cat("分析结果将保存至文件夹：", output_dir, "\n")

# 1. 加载必要的包 --------------------------------------------------------------
library(dplyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(multcomp)
library(multcompView)
library(ggplot2)

# 2. 读取数据 ------------------------------------------------------------------
df <- read.csv("Alfalfa_5Year_Master_Dataset.csv")

# 3. 数据清洗与聚合 ------------------------------------------------------------
df_summer <- df %>%
  filter(!is.na(Summer_Height_cm)) %>%
  group_by(Year, Plot, Line, Rep) %>%
  summarise(Summer_Height = mean(Summer_Height_cm, na.rm = TRUE), .groups = "drop") %>%
  mutate(Year = as.factor(Year),
         Line = as.factor(Line),
         Rep  = as.factor(Rep))

# 4. 数据平衡性检查与缺失组合记录 ----------------------------------------------
balance_table <- df_summer %>%
  group_by(Year, Line) %>%
  summarise(n = n(), .groups = "drop")
write.csv(balance_table, file.path(output_dir, "Balance_Check.csv"), row.names = FALSE)

all_combos <- expand.grid(Year = unique(df_summer$Year), Line = unique(df_summer$Line))
missing_combos <- anti_join(all_combos, balance_table, by = c("Year", "Line"))

if (nrow(missing_combos) > 0) {
  cat("\n警告：存在缺失的 Year-Line 组合，共", nrow(missing_combos), "个。\n")
  write.csv(missing_combos, file.path(output_dir, "Missing_Combinations.csv"), row.names = FALSE)
  cat("-> 缺失组合列表已保存为 Missing_Combinations.csv\n")
}

n_ratio <- max(balance_table$n) / min(balance_table$n)
balanced <- n_ratio < 2

# 5. 模型拟合与方差分析 --------------------------------------------------------
if (balanced) {
  cat("\n========== 使用传统方差分析 (aov) ==========\n")
  model_aov <- aov(Summer_Height ~ Year * Line + Error(Rep/Plot), data = df_summer)
  
  aov_summary <- summary(model_aov)
  fixed_table <- as.data.frame(aov_summary$'Error: Rep:Plot'[[1]])
  
  # 【修复】直接提取原 P 值，并创建规范统一的 P_value 列
  p_col_aov <- if ("Pr(>F)" %in% colnames(fixed_table)) "Pr(>F)" else grep("Pr", colnames(fixed_table), value=TRUE)[1]
  fixed_table$P_value <- fixed_table[[p_col_aov]]
  
  # 避免重建 data.frame，直接操作属性
  fixed_table$Effect <- rownames(fixed_table)
  rownames(fixed_table) <- NULL
  
  fixed_table$Signif <- symnum(fixed_table$P_value, 
                               cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                               symbols = c("***", "**", "*", ".", " "))
  
  # 将 Effect 列移到第一位
  fixed_table <- fixed_table[, c("Effect", setdiff(colnames(fixed_table), "Effect"))]
  write.csv(fixed_table, file.path(output_dir, "Anova_Table_aov.csv"), row.names = FALSE)
  
  # 【稳健性提取】使用统一的 P_value 字段
  p_interact <- if("Year:Line" %in% fixed_table$Effect) {
    fixed_table$P_value[fixed_table$Effect == "Year:Line"]
  } else {
    NA_real_
  }
  
  if (length(p_interact) == 1 && !is.na(p_interact) && p_interact < 0.05) {
    emm <- emmeans(model_aov, ~ Line | Year)
  } else {
    emm <- emmeans(model_aov, ~ Line)
  }
  
} else {
  cat("\n========== 使用混合模型 (lmer) ==========\n")
  model_lmer <- lmer(Summer_Height ~ Year * Line + (1|Rep/Plot), data = df_summer)
  
  # 注意：由于数据存在大量缺失组合，采用渐近 Wald 检验 (ddf = "lme4") 
  anova_table <- anova(model_lmer, type = "III", ddf = "lme4")
  anova_df <- as.data.frame(anova_table)
  
  # 【终极修复】极其安全地提取或计算 P 值
  if ("Pr(>F)" %in% colnames(anova_df)) {
    p_vals <- anova_df[["Pr(>F)"]]
  } else if ("Pr(>Chisq)" %in% colnames(anova_df)) {
    p_vals <- anova_df[["Pr(>Chisq)"]]
  } else {
    f_col <- if ("F value" %in% colnames(anova_df)) "F value" else grep("F", colnames(anova_df), value = TRUE)[1]
    df_col <- if ("NumDF" %in% colnames(anova_df)) "NumDF" else if ("npar" %in% colnames(anova_df)) "npar" else "Df"
    p_vals <- pf(as.numeric(anova_df[[f_col]]), as.numeric(anova_df[[df_col]]), Inf, lower.tail = FALSE)
  }
  
  # 显式赋值给 P_value
  anova_df$P_value <- p_vals
  
  # 避免重建 data.frame 篡改列名
  anova_df$Effect <- rownames(anova_df)
  rownames(anova_df) <- NULL
  
  # 生成星号
  anova_df$Signif <- symnum(anova_df$P_value,  
                            cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                            symbols = c("***", "**", "*", ".", " "))
  
  # 调整列顺序，保存
  anova_df <- anova_df[, c("Effect", setdiff(colnames(anova_df), "Effect"))]
  write.csv(anova_df, file.path(output_dir, "Anova_Table_lmer.csv"), row.names = FALSE)
  
  # 【稳健性提取】无论原列名叫什么，统一调用 P_value
  p_interact <- if("Year:Line" %in% anova_df$Effect) {
    anova_df$P_value[anova_df$Effect == "Year:Line"]
  } else {
    NA_real_
  }
  
  if (length(p_interact) == 1 && !is.na(p_interact) && p_interact < 0.05) {
    emm <- emmeans(model_lmer, ~ Line | Year, mode = "asymptotic", 
                   lmerTest.limit = 0, pbkrtest.limit = 0)
  } else {
    emm <- emmeans(model_lmer, ~ Line, mode = "asymptotic", 
                   lmerTest.limit = 0, pbkrtest.limit = 0)
  }
}

# 6. 多重比较与结果保存 --------------------------------------------------------
cat("\n正在计算边际均值和多重比较...\n")

cld_result <- tryCatch({
  multcomp::cld(emm, Letters = letters, adjust = "sidak")
}, error = function(e) {
  warning("cld 计算失败：", e$message)
  NULL
})

if (!is.null(cld_result)) {
  emm_df <- as.data.frame(cld_result)
  if (".group" %in% colnames(emm_df)) emm_df$.group <- trimws(emm_df$.group)
} else {
  emm_df <- as.data.frame(summary(emm))
}

write.csv(emm_df, file.path(output_dir, "Marginal_Means.csv"), row.names = FALSE)

pairs_result <- tryCatch(
  pairs(emm, adjust = "sidak"),
  error = function(e) NULL
)
if (!is.null(pairs_result)) {
  pairs_df <- as.data.frame(summary(pairs_result))
  write.csv(pairs_df, file.path(output_dir, "Pairwise_Comparisons.csv"), row.names = FALSE)
}

# ==================== 可视化部分 ====================
cat("\n正在生成图形...\n")

# 7.1 箱线图
p_box <- ggplot(df_summer, aes(x = Line, y = Summer_Height, fill = Year)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  theme_minimal(base_size = 12) +
  labs(title = "不同品系夏季株高分布（按年份）", x = "品系", y = "株高 (cm)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(output_dir, "Boxplot_Height_by_Line_Year.png"), p_box, width = 12, height = 6, dpi = 300)

# 7.2 交互作用图（Top 10 品系）
means_yr_line <- df_summer %>%
  group_by(Year, Line) %>%
  summarise(Mean_Height = mean(Summer_Height, na.rm = TRUE), .groups = "drop")

top_10_lines <- df_summer %>%
  group_by(Line) %>%
  summarise(Overall_Mean = mean(Summer_Height, na.rm = TRUE)) %>%
  top_n(10, Overall_Mean) %>%
  pull(Line)

means_yr_line_top10 <- means_yr_line %>% filter(Line %in% top_10_lines)

p_interact_top10 <- ggplot(means_yr_line_top10, aes(x = Year, y = Mean_Height, 
                                                    group = Line, color = Line)) +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  theme_minimal(base_size = 12) +
  labs(title = "品系×年份交互作用图 (Top 10 品系)", x = "年份", y = "平均株高 (cm)") +
  theme(legend.position = "right")
ggsave(file.path(output_dir, "Interaction_Plot_Top10.png"), 
       p_interact_top10, width = 8, height = 5, dpi = 300)

# 7.3 带字母标记的柱状图
if (exists("emm_df") && nrow(emm_df) > 0 && ".group" %in% colnames(emm_df)) {
  
  lcl_col <- if("asymp.LCL" %in% colnames(emm_df)) "asymp.LCL" else "lower.CL"
  ucl_col <- if("asymp.UCL" %in% colnames(emm_df)) "asymp.UCL" else "upper.CL"
  
  emm_df_plot <- emm_df[!is.na(emm_df$emmean) & !is.na(emm_df[[lcl_col]]), ]
  
  if ("Year" %in% colnames(emm_df_plot)) {
    p_cld <- ggplot(emm_df_plot, aes(x = Line, y = emmean, fill = Year)) +
      geom_col(position = position_dodge(0.9), width = 0.7) +
      geom_errorbar(aes(ymin = .data[[lcl_col]], ymax = .data[[ucl_col]]),
                    width = 0.2, position = position_dodge(0.9)) +
      geom_text(aes(y = emmean + .data[[ucl_col]]*0.1, label = .group, group = Year),
                position = position_dodge(0.9), vjust = -0.5, size = 3) +
      theme_minimal(base_size = 12) +
      labs(title = "品系间株高多重比较", x = "品系", y = "最小二乘均值 (cm)") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  } else {
    p_cld <- ggplot(emm_df_plot, aes(x = Line, y = emmean)) +
      geom_col(fill = "skyblue", width = 0.7) +
      geom_errorbar(aes(ymin = .data[[lcl_col]], ymax = .data[[ucl_col]]), width = 0.2) +
      geom_text(aes(y = emmean + .data[[ucl_col]]*0.1, label = .group), vjust = -0.5, size = 3) +
      theme_minimal(base_size = 12) +
      labs(title = "品系间株高多重比较", x = "品系", y = "最小二乘均值 (cm)") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  
  ggsave(file.path(output_dir, "CLD_Barplot.png"), p_cld, width = 14, height = 6, dpi = 300)
}

cat("\n所有分析及图表美化完成！结果保存在：", output_dir, "\n")
