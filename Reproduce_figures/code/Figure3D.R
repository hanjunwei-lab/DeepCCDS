
library(ggplot2)
library(dplyr)
data<-read.csv("C:\\Figure3D.csv" ,header = TRUE, stringsAsFactors = FALSE)

avg_data <- data %>%
  group_by(method) %>%
  summarise(mean_PCC = mean(PCC),
            mean_RMSE = mean(RMSE)) %>%
  ungroup()

deep_avg <- avg_data %>%
  filter(method %in% c("DeepCCDS", "DeepCCDS265")) %>%
  summarise(method = "DeepCCDS &\nDeepCCDS265",
            mean_PCC = mean(mean_PCC),
            mean_RMSE = mean(mean_RMSE))

other_avg <- avg_data %>%
  filter(!method %in% c("DeepCCDS", "DeepCCDS265")) %>%
  summarise(method = "Other\nMethods",
            mean_PCC = mean(mean_PCC),
            mean_RMSE = mean(mean_RMSE))


#PCC
sorted_main <- avg_data %>%
  arrange(desc(mean_PCC))


final_data <- bind_rows(
  sorted_main,
  deep_avg,
  other_avg
)


final_data$method <- factor(final_data$method, 
                            levels = c(as.character(sorted_main$method), 
                                       "DeepCCDS &\nDeepCCDS265", 
                                       "Other\nMethods"))


ggplot(final_data, aes(x = method, y = mean_PCC)) +
  geom_bar(stat = "identity", fill = c(rep("steelblue", nrow(sorted_main)), "orange", "gray")) +
  geom_text(aes(label = round(mean_PCC, 3)), vjust = -0.3, size = 3) +
  labs(title = "Average PCC by Method", 
       x = "Method", 
       y = "Average PCC") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



#RMSE
sorted_main_RMSE <- avg_data %>%
  arrange(desc(mean_RMSE))


final_data_RMSE <- bind_rows(
  sorted_main_RMSE,
  deep_avg,
  other_avg
)


final_data_RMSE$method <- 
  factor(final_data_RMSE$method, 
         levels = c(as.character(sorted_main_RMSE$method), 
                    "DeepCCDS &\nDeepCCDS265", 
                    "Other\nMethods"))

ggplot(final_data_RMSE, aes(x = method, y = mean_RMSE)) +
  geom_bar(stat = "identity", fill = c(rep("steelblue", nrow(sorted_main_RMSE)), "orange", "gray")) +
  geom_text(aes(label = round(mean_RMSE, 3)), vjust = -0.3, size = 3) +
  labs(title = "Average RMSE by Method", 
       x = "Method", 
       y = "Average RMSE") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
