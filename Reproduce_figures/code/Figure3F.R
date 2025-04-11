###画图（条形图）
data<-read.csv("C:\\Figure3F.csv", header = TRUE, stringsAsFactors = FALSE)
library(ggplot2)
library(dplyr)
## 1. 定义你想要的顺序
method_order <- c( "DeepCCDS", "DeepCCDS265", "DeepDR", 
                   "DeepTTA","Precily", "BANDRP","DeepCDR", "DrugCell", 
                   "lasso", "ridge","elastic",  "SVM")

## 2. 将method列转换为因子并指定顺序
data$method <- factor(data$method, levels = method_order)
# 绘制图形
###PCC
ggplot(data, aes(x = method, y = PCC_mean, fill = class)) +
  geom_bar(stat = "identity", 
           position = position_dodge(width = 0.7), 
           width = 0.7) +
  geom_errorbar(aes(ymin = PCC_mean - PCC_standard, 
                    ymax = PCC_mean + PCC_standard),
                position = position_dodge(width = 0.7), 
                width = 0.25,
                linewidth = 0.4) +
  labs(x = "Method",
       y = "PCC in GDSC",
       fill = "") +  # 清空图例标题
  scale_fill_manual(values = c("cell" = "#ecb761", "drug" = "#5abed0"),
                    labels = c("Cell", "Drug")) +  # 自定义图例标签
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
    axis.ticks.y = element_line(),  # 添加y轴刻度线
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.position = "top",  # 图例放在上方
    legend.justification = "center",  # 图例居中
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(-10, 0, -5, 0),  # 调整图例与图之间的距离
    plot.margin = unit(c(1, 1, 1, 1), "cm")  # 调整图形边距
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_x_discrete(expand = expansion(add = 0.5))  # 增加方法间的间隔

