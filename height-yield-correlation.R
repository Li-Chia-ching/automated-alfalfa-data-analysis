# 确保 showtext 已加载并启用
library(showtext)
showtext_auto()

# 生成热图（假设 cor_df 已准备好）
p <- ggplot(cor_df, aes(x = var1, y = var2, fill = cor)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label), size = 6, family = "Arial") + 
  scale_fill_gradient2(
    low = "#2166ac", mid = "#f7f7f7", high = "#b2182b",
    midpoint = 0, limit = c(-1, 1), space = "Lab",
    name = "Pearson\nCorrelation"
  ) +
  theme_minimal() +
  theme(
    text = element_text(family = "Arial"),
    axis.text.x = element_text(angle = 14, hjust = 1, vjust = 1, 
                               size = 16, family = "Arial"),           
    axis.text.y = element_text(size = 16, family = "Arial"),          
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.text = element_text(size = 10),                              # 图例数值略小一点
    legend.title = element_text(size = 12)                              # 图例标题与坐标轴一致
  ) +
  coord_fixed()
p

# 保存为PNG，增大图片尺寸至8x7英寸，为文字留足空间
ggsave("correlation_heatmap.png", plot = p, width = 8, height = 7, dpi = 300)
