# ============================================================================
# Spring_Phenotype.R
# 春季返青期苜蓿表型数据分析脚本（最终优化版 v4）
# 改进：
#   - 死亡率条形图改为水平条形图，文字横向
#   - 节间数与分枝数条形图使用 Duncan's method (P<0.05) 进行多重比较，标小写字母
#   - 基于帕累托前沿综合四个性状，全局筛选200个极端单株
#   - 使用PCA双标图展示选中个体在四维空间中的分布
#   - 所有图表均输出对应的原始绘图数据 (CSV格式)
# ============================================================================

# 清理环境
rm(list = ls())

# 加载必要的包
library(dplyr)
library(ggplot2)
library(tidyr)
library(agricolae)      # 用于Duncan多重比较
library(showtext)       # 支持字体
library(rPref)          # 帕累托前沿计算
library(FactoMineR)     # PCA分析
library(factoextra)     # PCA可视化

# 设置 ggplot2 默认字体
font_family <- "Arial"
tryCatch({
  showtext_auto()
  font_add("Arial", "arial.ttf")
}, error = function(e) {
  message("Arial 字体不可用，将使用系统默认 sans 字体")
  font_family <- "sans"
})

# 核心分析函数------------------------------------------------------------------

spring_phenotype_analysis <- function(data,
                                      group_col = "Group",
                                      matrix_col = "Matrix",
                                      traits = c("Plant_Height", "Internode", 
                                                 "Branch_Number", "Multifoliate_Score"),
                                      death_col = "Death_Code",
                                      extreme_percent = 0.03) {
  
  required_cols <- c(group_col, matrix_col, traits, death_col)
  if (!all(required_cols %in% colnames(data))) {
    stop("数据框中缺少必要的列：", 
         paste(setdiff(required_cols, colnames(data)), collapse = ", "))
  }
  
  data[[group_col]] <- as.factor(data[[group_col]])
  
  extreme_list <- list()
  groups <- unique(data[[group_col]])
  
  # 单性状极端个体提取
  for (g in groups) {
    sub_data <- data[data[[group_col]] == g, ]
    n <- nrow(sub_data)
    k <- ceiling(extreme_percent * n)
    
    if (k == 0) {
      warning(paste("Group", g, "样本量过小，无法提取极端个体"))
      next
    }
    
    for (trait in traits) {
      if (!is.numeric(sub_data[[trait]])) next
      
      valid_idx <- !is.na(sub_data[[trait]])
      if (sum(valid_idx) < 2 * k) next
      
      order_idx <- order(sub_data[[trait]][valid_idx])
      low_indices <- which(valid_idx)[order_idx[1:k]]
      high_indices <- which(valid_idx)[order_idx[(n - k + 1):n]]
      
      low_matrix <- sub_data[[matrix_col]][low_indices]
      high_matrix <- sub_data[[matrix_col]][high_indices]
      
      if (length(low_matrix) > 0) {
        extreme_list <- append(extreme_list, list(data.frame(Group = g, Trait = trait, Side = "Low", Matrix = low_matrix, stringsAsFactors = FALSE)))
      }
      if (length(high_matrix) > 0) {
        extreme_list <- append(extreme_list, list(data.frame(Group = g, Trait = trait, Side = "High", Matrix = high_matrix, stringsAsFactors = FALSE)))
      }
    }
  }
  
  candidate_extreme <- if (length(extreme_list) > 0) do.call(rbind, extreme_list) else data.frame()
  
  # 死亡率统计
  mortality_df <- data %>%
    group_by(.data[[group_col]]) %>%
    summarise(Total = n(),
              Deaths = sum(.data[[death_col]] %in% c(1, 2), na.rm = TRUE),
              Mortality_Rate = Deaths / Total, .groups = "drop")
  
  total_mortality <- data %>%
    summarise(Group = "Overall", Total = n(),
              Deaths = sum(.data[[death_col]] %in% c(1, 2), na.rm = TRUE),
              Mortality_Rate = Deaths / Total)
  mortality_df <- bind_rows(mortality_df, total_mortality)
  
  return(list(candidate_extreme = candidate_extreme, mortality = mortality_df))
}

# 描述性统计--------------------------------------------------------------------

compute_descriptive_stats <- function(data, traits, group_col = "Group") {
  stats_list <- list()
  for (trait in traits) {
    df <- data %>%
      group_by(.data[[group_col]]) %>%
      summarise(Trait = trait, N = sum(!is.na(.data[[trait]])),
                Mean = mean(.data[[trait]], na.rm = TRUE),
                SD = sd(.data[[trait]], na.rm = TRUE),
                SE = SD / sqrt(N), .groups = "drop")
    stats_list[[trait]] <- df
  }
  return(bind_rows(stats_list))
}

# 多重比较（Duncan's method, P < 0.05）与小写字母标记---------------------------

add_duncan_letters <- function(data, trait, group_col = "Group") {
  df <- data[!is.na(data[[trait]]), c(group_col, trait)]
  if (nrow(df) == 0) return(NULL)
  
  aov_res <- aov(as.formula(paste(trait, "~", group_col)), data = df)
  duncan_res <- duncan.test(aov_res, group_col, group = TRUE, console = FALSE)
  
  letters_df <- data.frame(Group = rownames(duncan_res$groups), 
                           Letters = tolower(as.character(duncan_res$groups$groups)),
                           stringsAsFactors = FALSE)
  return(letters_df)
}

# 可视化函数（含原始数据导出）--------------------------------------------------

## 1. 死亡率水平条形图-----------------------------------------------------------
plot_mortality <- function(mortality_df, output_path) {
  plot_data <- mortality_df %>% filter(Group != "Overall")
  
  p <- ggplot(plot_data, aes(x = Group, y = Mortality_Rate, fill = Group)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(aes(label = paste0(round(Mortality_Rate * 100, 1), "%")), 
              hjust = -0.2, size = 3.5, family = font_family) +
    labs(title = "Mortality rate by group", x = "Group", y = "Mortality rate") +
    scale_y_continuous(labels = scales::percent_format(), expand = expansion(mult = c(0, 0.2))) +
    scale_fill_viridis_d(option = "D", begin = 0.2, end = 0.8) +
    coord_flip() +
    theme_classic() +
    theme(text = element_text(family = font_family), legend.position = "none")
  
  ggsave(file.path(output_path, "1_Mortality_barplot.png"), p, width = 8, height = 5, dpi = 300)
  write.csv(plot_data, file.path(output_path, "1_Mortality_plot_data.csv"), row.names = FALSE)
  return(p)
}

## 2. 株高箱线图-----------------------------------------------------------------
plot_plant_height <- function(data, output_path) {
  plot_data <- data %>% filter(!is.na(Plant_Height))
  
  p <- ggplot(plot_data, aes(x = Group, y = Plant_Height, fill = Group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 0.8, alpha = 0.5, color = "gray40") +
    labs(title = "Plant height distribution", x = "Group", y = "Plant height (cm)") +
    scale_fill_viridis_d(option = "D", begin = 0.2, end = 0.8) +
    theme_classic() +
    theme(text = element_text(family = font_family),
          axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  
  ggsave(file.path(output_path, "2_Plant_Height_boxplot.png"), p, width = 8, height = 5, dpi = 300)
  write.csv(plot_data %>% select(Matrix, Group, Plant_Height), 
            file.path(output_path, "2_Plant_Height_plot_data.csv"), row.names = FALSE)
  return(p)
}

## 3. 节间数与分枝数的均值±SE条形图 + Duncan小写字母标记-------------------------
plot_count_traits <- function(data, output_path) {
  stats <- data %>%
    group_by(Group) %>%
    summarise(
      Internode_mean = mean(Internode, na.rm = TRUE),
      Internode_se = sd(Internode, na.rm = TRUE) / sqrt(sum(!is.na(Internode))),
      Branch_Number_mean = mean(Branch_Number, na.rm = TRUE),
      Branch_Number_se = sd(Branch_Number, na.rm = TRUE) / sqrt(sum(!is.na(Branch_Number))),
      .groups = "drop"
    ) %>%
    pivot_longer(cols = c(Internode_mean, Branch_Number_mean), names_to = "Trait", values_to = "Mean") %>%
    mutate(SE = ifelse(Trait == "Internode_mean", Internode_se, Branch_Number_se),
           Trait = ifelse(Trait == "Internode_mean", "Internode", "Branch number"))
  
  letters_internode <- add_duncan_letters(data, "Internode")
  letters_branch <- add_duncan_letters(data, "Branch_Number")
  
  letters_internode$Trait <- "Internode"
  letters_branch$Trait <- "Branch number"
  letters_all <- bind_rows(letters_internode, letters_branch)
  
  plot_data <- left_join(stats, letters_all, by = c("Group", "Trait"))
  
  p <- ggplot(plot_data, aes(x = Group, y = Mean, fill = Group)) +
    geom_bar(stat = "identity", width = 0.7, alpha = 0.7) +
    geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2) +
    geom_text(aes(y = Mean + SE, label = Letters), vjust = -0.8, size = 4, family = font_family) +
    facet_wrap(~ Trait, scales = "free_y", ncol = 1) +
    labs(title = "Internode and branch number (Duncan's test, P<0.05)", x = "Group", y = "Value") +
    scale_fill_viridis_d(option = "D", begin = 0.2, end = 0.8) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
    theme_classic() +
    theme(text = element_text(family = font_family),
          axis.text.x = element_text(angle = 45, hjust = 1), 
          legend.position = "none",
          strip.text = element_text(face = "bold", size = 11))
  
  ggsave(file.path(output_path, "3_Count_traits_barchart.png"), p, width = 10, height = 8, dpi = 300)
  write.csv(plot_data, file.path(output_path, "3_Count_traits_plot_data.csv"), row.names = FALSE)
  return(p)
}

## 4. 多叶评分堆叠图-------------------------------------------------------------
plot_multifoliate_score <- function(data, output_path) {
  score_data <- data %>%
    filter(!is.na(Multifoliate_Score)) %>%
    group_by(Group, Multifoliate_Score) %>%
    summarise(Freq = n(), .groups = "drop") %>%
    group_by(Group) %>%
    mutate(Percent = Freq / sum(Freq) * 100)
  
  p <- ggplot(score_data, aes(x = Group, y = Percent, fill = as.factor(Multifoliate_Score))) +
    geom_bar(stat = "identity", position = "fill", width = 0.7) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(title = "Multifoliate score composition", x = "Group", y = "Proportion", fill = "Score") +
    scale_fill_viridis_d(option = "D", begin = 0.2, end = 0.8) +
    theme_classic() +
    theme(text = element_text(family = font_family),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(output_path, "4_Multifoliate_Score_stackedbar.png"), p, width = 9, height = 6, dpi = 300)
  write.csv(score_data, file.path(output_path, "4_Multifoliate_Score_plot_data.csv"), row.names = FALSE)
  return(p)
}

## 5. 多指标帕累托综合筛选（200个单株）+ PCA双标图可视化-------------------------

select_pareto_top200 <- function(data, traits, output_path, n_select = 200, seed = 123) {
  set.seed(seed)
  
  # 1. 过滤缺失值，确保帕累托计算有效
  df_valid <- data %>% 
    filter(if_all(all_of(traits), ~ !is.na(.)))
  
  if(nrow(df_valid) < n_select) {
    warning("有效完整数据量少于预期的选择数量！将返回所有有效数据。")
    n_select <- nrow(df_valid)
  }
  
  # 2. 定义帕累托偏好 (四个性状均为越大越好)
  pref <- high(Plant_Height) * high(Internode) * high(Branch_Number) * high(Multifoliate_Score)
  
  # 3. 计算所有个体的帕累托前沿层级
  pareto_res <- psel(df_valid, pref, top = nrow(df_valid), show_level = TRUE)
  
  # 4. 逐层提取，直到刚好满足 n_select 个
  selected_df <- data.frame()
  current_level <- 1
  needed <- n_select
  
  while(needed > 0 && current_level <= max(pareto_res$.level)) {
    front_data <- pareto_res %>% filter(.level == current_level)
    
    if (nrow(front_data) <= needed) {
      selected_df <- bind_rows(selected_df, front_data)
      needed <- needed - nrow(front_data)
    } else {
      # 按 Group 比例在当前层中分层抽样
      sampled_front <- front_data %>%
        group_by(Group) %>%
        sample_frac(size = needed / nrow(front_data)) %>% 
        ungroup()
      
      # 修正四舍五入导致的数量偏差
      diff <- needed - nrow(sampled_front)
      if (diff > 0) {
        remain <- anti_join(front_data, sampled_front, by = "Matrix")
        sampled_front <- bind_rows(sampled_front, remain %>% sample_n(diff))
      } else if (diff < 0) {
        sampled_front <- sampled_front %>% sample_n(needed)
      }
      
      selected_df <- bind_rows(selected_df, sampled_front)
      needed <- 0
    }
    current_level <- current_level + 1
  }
  
  # 去除 rPref 自动生成的 .level 列
  selected_df <- selected_df %>% select(-.level)
  
  # 5. PCA 双标图可视化（使用所有四个性状）
  pca_data <- df_valid %>% select(all_of(traits))
  pca_res <- PCA(pca_data, scale.unit = TRUE, graph = FALSE)
  
  # 提取个体坐标
  indiv_coord <- as.data.frame(pca_res$ind$coord[, 1:2])
  colnames(indiv_coord) <- c("Dim1", "Dim2")
  indiv_coord$Matrix <- df_valid$Matrix
  indiv_coord$Group <- df_valid$Group
  indiv_coord$Status <- ifelse(indiv_coord$Matrix %in% selected_df$Matrix, "Selected", "Unselected")
  
  # 提取变量载荷（前两个主成分）
  var_loadings <- as.data.frame(pca_res$var$coord[, 1:2])
  colnames(var_loadings) <- c("Dim1", "Dim2")
  var_loadings$Variable <- rownames(var_loadings)
  
  # 缩放载荷以便在图上显示（避免箭头过长）
  scale_factor <- max(abs(indiv_coord$Dim1), abs(indiv_coord$Dim2)) / max(abs(var_loadings[,1:2])) * 0.8
  
  # 绘制 PCA 双标图
  p_pca <- ggplot() +
    # 未选中个体（灰色）
    geom_point(data = filter(indiv_coord, Status == "Unselected"),
               aes(x = Dim1, y = Dim2), color = "gray85", alpha = 0.6, size = 1.5) +
    # 选中个体（按 Group 着色）
    geom_point(data = filter(indiv_coord, Status == "Selected"),
               aes(x = Dim1, y = Dim2, color = Group), size = 3, alpha = 0.9) +
    # 变量箭头（载荷）
    geom_segment(data = var_loadings,
                 aes(x = 0, y = 0, xend = Dim1 * scale_factor, yend = Dim2 * scale_factor),
                 arrow = arrow(length = unit(0.2, "cm")), color = "red", linewidth = 0.8) +
    geom_text(data = var_loadings,
              aes(x = Dim1 * scale_factor * 1.1, y = Dim2 * scale_factor * 1.1, label = Variable),
              color = "red", size = 4, family = font_family, hjust = 0.5, vjust = 0.5) +
    labs(title = paste("PCA biplot of selected individuals (Top", n_select, ")"),
         subtitle = paste("Based on:", paste(traits, collapse = ", ")),
         x = paste0("PC1 (", round(pca_res$eig[1,2], 1), "%)"),
         y = paste0("PC2 (", round(pca_res$eig[2,2], 1), "%)")) +
    scale_color_viridis_d(option = "D", begin = 0.2, end = 0.8) +
    theme_classic() +
    theme(text = element_text(family = font_family),
          legend.position = "right")
  
  ggsave(file.path(output_path, "5_Pareto_Selection_PCA.png"), p_pca, width = 8, height = 6, dpi = 300)
  write.csv(indiv_coord, file.path(output_path, "5_Pareto_Selection_PCA_data.csv"), row.names = FALSE)
  
  # 可选：同时保留原来的二维散点图（株高 vs 节间数）作为补充
  plot_data_scatter <- df_valid %>%
    mutate(Status = ifelse(Matrix %in% selected_df$Matrix, "Selected (Pareto Top)", "Unselected"))
  p_scatter <- ggplot(plot_data_scatter, aes(x = Plant_Height, y = Internode)) +
    geom_point(data = filter(plot_data_scatter, Status == "Unselected"), 
               color = "gray85", alpha = 0.6, size = 1.5) +
    geom_point(data = filter(plot_data_scatter, Status != "Unselected"), 
               aes(color = Group), size = 3, alpha = 0.9) +
    labs(title = paste("Pareto Front Selection (Top", n_select, ") - 2D view"),
         subtitle = "Plant Height vs Internode",
         x = "Plant Height (cm)", y = "Internode Count") +
    scale_color_viridis_d(option = "D", begin = 0.2, end = 0.8) +
    theme_classic() +
    theme(text = element_text(family = font_family),
          legend.position = "right")
  ggsave(file.path(output_path, "5_Pareto_Selection_Scatter.png"), p_scatter, width = 8, height = 6, dpi = 300)
  write.csv(plot_data_scatter, file.path(output_path, "5_Pareto_Selection_Scatter_data.csv"), row.names = FALSE)
  
  return(list(pca_plot = p_pca, scatter_plot = p_scatter, selected_df = selected_df))
}

# 主程序------------------------------------------------------------------------

analysis_date <- Sys.Date()
folder_name <- paste0("Phenotype_Viz_Upgrade_", analysis_date)
if (!dir.exists(folder_name)) dir.create(folder_name)
cat("创建文件夹:", folder_name, "\n")

data_file <- "Raw_SprGreenUp.csv"
if (!file.exists(data_file)) {
  stop("数据文件不存在：", data_file)
}
data <- read.csv(data_file, stringsAsFactors = FALSE)

cat("数据读取成功，共", nrow(data), "行记录。\n")

# 执行核心分析
results <- spring_phenotype_analysis(data)

# 执行多指标帕累托综合筛选（200个单株）并生成PCA双标图
traits_for_pareto <- c("Plant_Height", "Internode", "Branch_Number", "Multifoliate_Score")
selection_results <- select_pareto_top200(data, traits_for_pareto, folder_name, n_select = 200)

selected_extreme_200 <- selection_results$selected_df

# 保存常规统计结果
write.csv(results$candidate_extreme, file.path(folder_name, "candidate_extreme_individuals.csv"), row.names = FALSE)
write.csv(selected_extreme_200, file.path(folder_name, "selected_extreme_200.csv"), row.names = FALSE)
write.csv(results$mortality, file.path(folder_name, "mortality.csv"), row.names = FALSE)

traits <- c("Plant_Height", "Internode", "Branch_Number", "Multifoliate_Score")
descriptive_stats <- compute_descriptive_stats(data, traits)
write.csv(descriptive_stats, file.path(folder_name, "descriptive_statistics.csv"), row.names = FALSE)

cat("\n正在生成其余图形...\n")
plot_mortality(results$mortality, folder_name)
plot_plant_height(data, folder_name)
plot_count_traits(data, folder_name)
plot_multifoliate_score(data, folder_name)

cat("\n===== 分析摘要 =====\n")
cat("图表及对应原始数据 CSV 均已生成在文件夹:", folder_name, "内。\n")
cat("最终选定覆盖各组的综合极端个体数:", nrow(selected_extreme_200), "\n")
cat("帕累托筛选 PCA 双标图已保存为 5_Pareto_Selection_PCA.png\n")
cat("分析完成！\n")
