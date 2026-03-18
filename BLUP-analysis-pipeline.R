# =============================================================================
# 项目：苜蓿多年多点数据分析（GGE 双标图 + BLUP + 遗传力）
# 数据文件：Alfalfa_5Year_Master_Dataset.csv（需放在 R 项目根目录）
# 说明：本脚本自动完成数据读取、清洗、分析及可视化
# =============================================================================

# ------------------------------ 1. 加载包 ------------------------------------
# 若未安装，请先运行 install.packages(c("lme4", "lmerTest", "ggplot2", "dplyr",
#                                         "tidyr", "tibble", "here", "conflicted"))
library(lme4)
library(lmerTest)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(here)          # 路径管理
library(conflicted)    # 解决函数冲突

# 优先使用 tidyverse 函数
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("mutate", "dplyr")
conflict_prefer("summarise", "dplyr")
conflict_prefer("arrange", "dplyr")
conflict_prefer("lmer", "lmerTest")   # 防止[conflicted] lmer found in 2 packages

# ------------------------------ 2. 读取数据 ------------------------------------
file_path <- here("Alfalfa_5Year_Master_Dataset.csv")
if (!file.exists(file_path)) stop("数据文件未找到，请检查路径：", file_path)

df_raw <- read.csv(file_path, stringsAsFactors = FALSE)

# 查看数据结构
glimpse(df_raw)

# ------------------------------ 3. 数据预处理 ----------------------------------
# 将关键列转为因子
df_raw <- df_raw %>%
  mutate(
    Year = as.factor(Year),
    Line = as.factor(Line),
    Rep  = as.factor(Rep)
  )

# 定义性状列表（根据实际表头，可自由增减）
traits <- c("Summer_Height_cm", "Fall_Dormancy_Height_cm", 
            "Fresh_Weight_kg", "Dry_Weight_kg")

# 检查各性状缺失情况
for (t in traits) {
  cat(t, "缺失值个数：", sum(is.na(df_raw[[t]])), "\n")
}

# 如有需要，可在此处进行缺失值填补（例如用该性状的总体均值），
# 但混合模型本身能处理缺失，通常不填补；GGE 分析需要完整矩阵，后面单独处理。

# ------------------------------ 4. GGE 双标图分析（针对产量 Dry_Weight_kg）----
cat("\n========== GGE 双标图分析（产量性状：Dry_Weight_kg）==========\n")

# 4.1 按品种(Line)和年份(Year)聚合产量均值
yield_met <- df_raw %>%
  filter(!is.na(Dry_Weight_kg)) %>%
  group_by(Line, Year) %>%
  summarise(Yield = mean(Dry_Weight_kg, na.rm = TRUE), .groups = "drop")

# 检查每年观测品种数
print(table(yield_met$Year))

# 4.2 构建品种×年份矩阵
ge_mat <- yield_met %>%
  pivot_wider(names_from = Year, values_from = Yield) %>%
  column_to_rownames("Line") %>%
  as.matrix()

# 4.3 处理缺失值（用该年所有品种的平均值填补）
for(i in 1:ncol(ge_mat)){
  col_mean <- mean(ge_mat[, i], na.rm = TRUE)
  ge_mat[is.na(ge_mat[, i]), i] <- col_mean
}
stopifnot(all(!is.na(ge_mat)))

# 4.4 环境中心化（减去年份均值）
env_means <- colMeans(ge_mat)
ge_centered <- sweep(ge_mat, 2, env_means, FUN = "-")

# 4.5 奇异值分解
svd_res <- svd(ge_centered)
pc_var <- svd_res$d^2 / sum(svd_res$d^2) * 100
PC1_var <- pc_var[1]
PC2_var <- pc_var[2]
cat(sprintf("PC1 解释了 %.1f%%，PC2 解释了 %.1f%%，合计 %.1f%%\n",
            PC1_var, PC2_var, PC1_var + PC2_var))

# 4.6 计算品种得分（对称缩放）
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

# 4.7 计算品种平均产量和 AMMI 稳定性值 (ASV)
mean_yield <- yield_met %>%
  group_by(Line) %>%
  summarise(Mean_Yield = mean(Yield, na.rm = TRUE))

# 计算 ASV
ipca1_score <- svd_res$u[,1] * sqrt(svd_res$d[1])
ipca2_score <- svd_res$u[,2] * sqrt(svd_res$d[2])
ss_ipca1 <- svd_res$d[1]^2          # 特征值平方 = SS_IPCA1
ss_ipca2 <- svd_res$d[2]^2
weight <- ss_ipca1 / ss_ipca2
ASV <- sqrt(weight * ipca1_score^2 + ipca2_score^2)

# 合并到 G_scores
G_scores <- G_scores %>%
  mutate(ASV = ASV)

# 4.8 合并平均产量与稳定性
variety_info <- mean_yield %>%
  left_join(G_scores %>% select(Line, ASV), by = "Line")

# 4.9 筛选高产稳产 Top 10 品种
# 方法：取平均产量前20名，再按 ASV 升序（越小越稳定）取前10
top_by_yield <- variety_info %>%
  arrange(desc(Mean_Yield)) %>%
  slice_head(n = 20)

top_stable <- top_by_yield %>%
  arrange(ASV) %>%
  slice_head(n = 10)

cat("\n===== GGE 推荐的高产稳产品种 Top 10 =====\n")
print(top_stable %>% select(Line, Mean_Yield, ASV))

# 4.10 绘制 GGE 双标图，高亮 Top 10
p_gge <- ggplot() +
  geom_point(data = G_scores, aes(x = PC1, y = PC2), 
             color = "steelblue", alpha = 0.5, size = 2) +
  geom_text(data = G_scores, aes(x = PC1, y = PC2, label = Line), 
            size = 2.5, vjust = -0.8, color = "steelblue") +
  geom_segment(data = E_scores, aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.15, "cm")), 
               color = "red", linewidth = 0.8) +
  geom_text(data = E_scores, aes(x = PC1, y = PC2, label = Year),
            color = "red", fontface = "bold", size = 4, vjust = 1.5) +
  geom_point(data = filter(G_scores, Line %in% top_stable$Line),
             aes(x = PC1, y = PC2), shape = 17, color = "red", size = 4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  labs(title = "GGE Biplot of Dry Weight Yield",
       subtitle = "Red triangles: Top 10 high-yielding & stable varieties",
       x = sprintf("PC1 (%.1f%%)", PC1_var),
       y = sprintf("PC2 (%.1f%%)", PC2_var)) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))

# 保存图形
if (!dir.exists("figures")) dir.create("figures")
ggsave("figures/GGE_Biplot_Top10.png", p_gge, width = 10, height = 8, dpi = 300)

# 保存品种信息表
write.csv(variety_info, "variety_yield_ASV.csv", row.names = FALSE)
write.csv(top_stable, "GGE_Top10_Recommended.csv", row.names = FALSE)

# ------------------------------ 5. BLUP 与遗传力分析（所有性状）---------------
cat("\n========== BLUP 与广义遗传力分析 ==========\n")

# 准备数据（确保因子类型）
df_met <- df_raw %>%
  mutate(
    Year = as.factor(Year),
    Line = as.factor(Line),
    Rep  = as.factor(Rep)
  )

# 存储结果
blup_results <- list()
h2_results <- data.frame(Trait = character(),
                         V_G = numeric(), V_GE = numeric(), V_e = numeric(),
                         H2 = numeric(), stringsAsFactors = FALSE)

# 循环分析每个性状
for (trait in traits) {
  cat("\n-------------------- 性状：", trait, "--------------------\n")
  
  # 检查性状是否存在且有足够数据
  if (!(trait %in% names(df_met))) {
    warning("性状 ", trait, " 不存在于数据框，跳过")
    next
  }
  if (all(is.na(df_met[[trait]]))) {
    warning("性状 ", trait, " 全部缺失，跳过")
    next
  }
  
  # 拟合混合模型：性状 ~ Year + (1|Line) + (1|Line:Year) + (1|Rep)
  formula_str <- paste(trait, "~ Year + (1 | Line) + (1 | Line:Year) + (1 | Rep)")
  # 设置控制参数：使用 bobyqa 优化器，增加最大函数评价次数
  control_params <- lmerControl(optimizer = "bobyqa",
                                optCtrl = list(maxfun = 2e5))
  
  model <- tryCatch(
    lmer(as.formula(formula_str), data = df_met, control = control_params),
    error = function(e) {
      warning("模型拟合失败：", e$message)
      return(NULL)
    }
  )
  if (is.null(model)) next
  
  # 收敛性检查（包含奇异拟合警告）
  if (length(model@optinfo$conv$lme4$messages) > 0) {
    warning("模型可能未收敛或奇异拟合：", 
            paste(model@optinfo$conv$lme4$messages, collapse = "; "))
    # 这里不跳过，继续执行，但提醒用户
  }
  
  # 提取品种 BLUP
  blup_raw <- ranef(model, condVar = TRUE)
  blup_line <- as.data.frame(blup_raw$Line) %>%
    rownames_to_column("Line") %>%
    rename(BLUP = `(Intercept)`) %>%
    mutate(StdErr = sqrt(attr(blup_raw$Line, "postVar")[1, 1, ]),
           lower = BLUP - 1.96 * StdErr,
           upper = BLUP + 1.96 * StdErr) %>%
    arrange(desc(BLUP))
  
  blup_results[[trait]] <- blup_line
  
  # 计算广义遗传力 H^2
  var_comps <- as.data.frame(VarCorr(model))
  V_G   <- var_comps$vcov[var_comps$grp == "Line"]
  V_GE  <- var_comps$vcov[var_comps$grp == "Line:Year"]
  V_err <- var_comps$vcov[var_comps$grp == "Residual"]
  
  e_count <- nlevels(df_met$Year)
  r_count <- df_met %>%
    filter(!is.na(.data[[trait]])) %>%
    group_by(Line, Year) %>%
    summarise(n = n(), .groups = "drop") %>%
    pull(n) %>%
    mean(na.rm = TRUE)
  
  H2 <- V_G / (V_G + (V_GE / e_count) + (V_err / (e_count * r_count)))
  
  h2_results <- rbind(h2_results, data.frame(
    Trait = trait,
    V_G = V_G,
    V_GE = V_GE,
    V_e = V_err,
    H2 = H2,
    stringsAsFactors = FALSE
  ))
  
  cat(sprintf("V_G = %.4f, V_GE = %.4f, V_e = %.4f\n", V_G, V_GE, V_err))
  cat(sprintf("H^2 = %.4f (%.2f%%)\n", H2, H2 * 100))
  
  # 绘制 Top 20 BLUP 图（如果品种数较多）
  if (nrow(blup_line) >= 20) {
    top20 <- blup_line %>% slice_head(n = 20)
  } else {
    top20 <- blup_line
  }
  
  p_blup <- ggplot(top20, aes(x = reorder(Line, BLUP), y = BLUP)) +
    geom_col(fill = "steelblue", width = 0.7) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
    coord_flip() +
    labs(title = paste("Top BLUP -", trait),
         x = "品种", y = "BLUP (95% CI)") +
    theme_minimal(base_size = 11)
  
  ggsave(filename = paste0("figures/BLUP_Top_", trait, ".png"),
         plot = p_blup, width = 8, height = 6, dpi = 300)
  
  # 保存该性状的完整 BLUP 表
  write.csv(blup_line, paste0("BLUP_", trait, ".csv"), row.names = FALSE)
}

# 输出遗传力汇总
cat("\n========== 遗传力汇总 ==========\n")
print(h2_results)
write.csv(h2_results, "Heritability_Summary.csv", row.names = FALSE)

# 绘制遗传力条形图
h2_plot <- ggplot(h2_results, aes(x = reorder(Trait, H2), y = H2)) +
  geom_col(fill = "lightgreen", width = 0.7) +
  geom_hline(yintercept = 0.3, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0.6, linetype = "dashed", color = "orange") +
  coord_flip() +
  labs(title = "广义遗传力 (H²) 比较",
       x = "性状", y = expression(H^2)) +
  theme_minimal(base_size = 12) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1))

ggsave("figures/Heritability_Barplot.png", h2_plot, width = 8, height = 4, dpi = 300)

# -------------------------- 6. 结合 BLUP 与 GGE 结果（直接读取文件）-----------
# 如果产量性状的 BLUP 文件和 ASV 文件都存在，则绘制四象限图
if (file.exists("BLUP_Dry_Weight_kg.csv") && file.exists("variety_yield_ASV.csv")) {
  
  cat("\n========== 绘制高产稳产四象限图 ==========\n")
  
  # 读取 BLUP 数据（产量性状）
  blup_yield <- read.csv("BLUP_Dry_Weight_kg.csv") %>%
    dplyr::select(Line, BLUP) %>%
    dplyr::rename(BLUP_Yield = BLUP)
  
  # 读取 ASV 数据（包含平均产量和 ASV）
  asv_data <- read.csv("variety_yield_ASV.csv")
  
  # 合并数据
  combined <- asv_data %>%
    left_join(blup_yield, by = "Line") %>%
    # 标记 Top 10 品种（需从 GGE 结果中获取，若无则用分位数法）
    mutate(Highlight = ifelse(Line %in% top_stable$Line, "Top10", "Others"))
  
  # 若 top_stable 不存在（例如未运行 GGE 部分），改用中位数划分
  if (!exists("top_stable")) {
    med_yield <- median(combined$Mean_Yield, na.rm = TRUE)
    med_asv <- median(combined$ASV, na.rm = TRUE)
    combined <- combined %>%
      mutate(Highlight = case_when(
        Mean_Yield > med_yield & ASV < med_asv ~ "Ideal",
        TRUE ~ "Others"
      ))
  }
  
  # 绘制四象限图
  p_quad <- ggplot(combined, aes(x = Mean_Yield, y = ASV, color = Highlight)) +
    geom_point(size = 2.5, alpha = 0.8) +
    geom_hline(yintercept = median(combined$ASV, na.rm = TRUE), 
               linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = median(combined$Mean_Yield, na.rm = TRUE), 
               linetype = "dashed", color = "grey40") +
    scale_color_manual(values = c("Top10" = "red", "Ideal" = "red", "Others" = "grey60")) +
    labs(title = "高产稳产四象限图（平均产量 vs ASV）",
         x = "平均产量 (kg)", y = "ASV (越小越稳定)") +
    theme_minimal() +
    theme(legend.position = "none")
  
  # 保存图形
  ggsave("figures/Yield_Stability_Quadrant.png", p_quad, width = 8, height = 6, dpi = 300)
  cat("高产稳产四象限图已保存至 figures/ 文件夹\n")
  
} else {
  cat("\n跳过四象限图：缺少 BLUP_Dry_Weight_kg.csv 或 variety_yield_ASV.csv 文件。\n")
}

cat("\n================= 全部分析完成！结果保存在当前工作目录 =================\n")
