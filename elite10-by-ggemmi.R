# =============================================================================
# 基于GGE双标图与AMMI思路的高产稳产品种筛选脚本（MET数据：5年×多个品种）
# 目标：计算每个品种的平均产量和稳定性指标，综合排序，输出Top 10推荐品种
# 数据要求：包含列 Line (品种), Year (年份), 以及目标性状 (如 Dry_Weight_kg)
# 注：本脚本采用环境中心化（GGE）进行SVD，并用PC2绝对值作为稳定性度量，
#     也可以替换为其他稳定性指标（如ASV）。
# =============================================================================

# 1. 加载必要包 ----------------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)

# 2. 数据准备：从原始数据框 df 中提取产量性状并按品种×年份聚合 --------------------
#   假设原始数据框 df 包含：Line, Year, Rep, Plot, Dry_Weight_kg 等
#   如果使用其他性状（如鲜重），请替换下面代码中的列名

# 确保年份为因子（后续需用于环境中心化）
df <- df %>% mutate(Year = as.factor(Year))

# 按品种(Line)和年份(Year)聚合产量（取各小区均值，若有多重复）
yield_met <- df %>%
  filter(!is.na(Dry_Weight_kg)) %>%          # 剔除产量缺失的记录
  group_by(Line, Year) %>%
  summarise(Yield = mean(Dry_Weight_kg, na.rm = TRUE), .groups = "drop")

# 检查数据量
cat("聚合后的数据维度：", dim(yield_met), "\n")
table(yield_met$Year)  # 查看每年有多少品种

# 3. 构建基因型 × 环境（年份）矩阵 ------------------------------------------------
ge_mat <- yield_met %>%
  pivot_wider(names_from = Year, values_from = Yield) %>%
  column_to_rownames("Line") %>%
  as.matrix()

# 查看矩阵概况
print(dim(ge_mat))
head(ge_mat)

# 4. 处理缺失值 ----------------------------------------------------------------
#   如果某年某品种缺失，用该年所有品种的平均值填补（也可用其他方法）
for(i in 1:ncol(ge_mat)){
  col_mean <- mean(ge_mat[, i], na.rm = TRUE)
  ge_mat[is.na(ge_mat[, i]), i] <- col_mean
}
# 检查是否还有缺失
stopifnot(all(!is.na(ge_mat)))

# 5. 环境中心化（GGE模型：减去每个环境（年份）的均值） ----------------------------
env_means <- colMeans(ge_mat)                # 各年份平均产量
ge_centered <- sweep(ge_mat, 2, env_means, FUN = "-")  # 中心化

# 6. 奇异值分解 SVD ------------------------------------------------------------
svd_res <- svd(ge_centered)

# 计算各主成分解释的方差比例
pc_var <- svd_res$d^2 / sum(svd_res$d^2) * 100
PC1_var <- pc_var[1]
PC2_var <- pc_var[2]
cat(sprintf("PC1解释了%.1f%%的变异，PC2解释了%.1f%%的变异，合计%.1f%%\n",
            PC1_var, PC2_var, PC1_var + PC2_var))

# 7. 计算品种得分（G-score）和环境得分（E-score）---------------------------------
#   采用对称缩放（使双标图可同时展示品种和环境）
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

# 8. 计算每个品种的平均产量（原始值）和稳定性指标 ---------------------------------
#   平均产量直接从原始聚合数据计算，确保准确
mean_yield <- yield_met %>%
  group_by(Line) %>%
  summarise(Mean_Yield = mean(Yield, na.rm = TRUE))

# 稳定性指标：这里采用PC2得分的绝对值（|PC2|越小，表示品种在各环境下的表现越一致）
#   也可使用类似AMMI稳定性值（ASV）的计算方法，但PC2绝对值简单直观
stability <- G_scores %>%
  dplyr::select(Line, PC1, PC2) %>%
  mutate(Stability = abs(PC2))

# 合并数据
variety_info <- mean_yield %>%
  left_join(stability, by = "Line")

# 9. 综合排名筛选Top 10高产稳产品种 --------------------------------------------
#   方法：先按产量降序排列，然后从高产群体中挑选稳定性较好的品种。
#   这里采用两个步骤：取产量前20的品种，再按稳定性（Stability从小到大）排序取前10。
#   您也可以构建一个综合指数，例如：Index = Mean_Yield - 2*Stability，然后排序。

top_by_yield <- variety_info %>%
  arrange(desc(Mean_Yield)) %>%
  slice_head(n = 20)               # 先选出产量最高的20个品种

top_stable <- top_by_yield %>%
  arrange(Stability) %>%            # 稳定性越好（Stability越小）越靠前
  slice_head(n = 10)                # 最终Top 10

cat("\n===== 推荐的高产稳产品种 Top 10 =====\n")
print(top_stable %>% dplyr::select(Line, Mean_Yield, Stability))

# 如果您想看到所有品种的排名，可输出全部排序
# variety_info %>% arrange(desc(Mean_Yield), Stability) %>% print(n=20)

# 10. 绘制GGE双标图，并高亮Top 10品种 ------------------------------------------
#     图中红色箭头代表年份，蓝色点为品种，Top 10品种用红色三角形突出显示

p <- ggplot() +
  # 品种点（全部）
  geom_point(data = G_scores, aes(x = PC1, y = PC2), color = "steelblue", alpha = 0.5, size = 2) +
  # 品种标签（全部）
  geom_text(data = G_scores, aes(x = PC1, y = PC2, label = Line), 
            size = 2.5, vjust = -0.8, color = "steelblue") +
  # 环境向量（从原点出发的箭头）
  geom_segment(data = E_scores, aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.15, "cm")), color = "red", linewidth = 0.8) +
  geom_text(data = E_scores, aes(x = PC1, y = PC2, label = Year),
            color = "red", fontface = "bold", size = 4, vjust = 1.5) +
  # 高亮Top 10品种（红色三角形）
  geom_point(data = filter(G_scores, Line %in% top_stable$Line),
             aes(x = PC1, y = PC2), shape = 17, color = "red", size = 4) +
  # 参考线
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  labs(title = "GGE Biplot of Yield (Dry Weight)",
       subtitle = "Red triangles indicate Top 10 high-yielding and stable varieties",
       x = sprintf("PC1 (%.1f%%)", PC1_var),
       y = sprintf("PC2 (%.1f%%)", PC2_var)) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))

print(p)

# 保存图形
ggsave("GGE_Biplot_Top10.png", p, width = 10, height = 8, dpi = 300)

# 11. 输出品种排名表格（含产量和稳定性）------------------------------------------
write.csv(variety_info, "variety_yield_stability.csv", row.names = FALSE)
cat("\n完整品种信息已保存至 variety_yield_stability.csv\n")
