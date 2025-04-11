
library(tidyverse)
library(cowplot)
library(export)

data<- read.csv("C:\\Figure3E.csv")

emptybar<- 1 #每组之间的间隔为1（小间隙）
empty<- data.frame(matrix(NA, emptybar*nlevels(factor(data$group))+3, 
                          ncol(data)))

colnames(empty)<- colnames(data)
data<-data[order(data$auc,decreasing = T),]
empty$group<- c(rep(levels(factor(data$group)), each = emptybar),
                rep("P", 3))

data<- rbind(data, empty) %>% arrange(group)
data$id<- c(c(14:26),c(1:13),c(27:42))
data<-data[order(data$id),] 


p1<-ggplot(data, aes(x = as.factor(id)))+ 
  geom_bar(aes(y = auc), stat="identity", alpha = 1);p1



group_data1<-data[c(1,12,14,25,27,38),]
group_data2<- data.frame(matrix(group_data1$id, byrow = T,
                                ncol = 2, nrow = 3),
                         Group = paste0("group", 1:3))
group_data2$group<-c("GDSC","CCLE","NCI60")
colnames(group_data2)<-c("x", "xend", "Group")
group_data1[c(1,3,5),5]<-group_data1[c(1,3,5),5]-0.8
group_data1[c(2,4,6),5]<-group_data1[c(2,4,6),5]+0.8

p2<- p1 +
  geom_segment(data = group_data1, aes(x = id, y = auc, xend = id, yend = 1), 
               color = "#1C7875", alpha = 0.8, size = 0.7, lty = 2)+
  geom_segment(data = group_data2, aes(x = x-0.8, y = 1, xend = xend+0.8, yend = 1), 
               color = "#1C7875", alpha = 0.8, size = 0.7, lty = 2)+
  ylim(-1,1)+

  theme_minimal()+
  theme(legend.position = "none",
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank())+
  coord_polar();p2#将平面直角坐标系转化为极坐标系



