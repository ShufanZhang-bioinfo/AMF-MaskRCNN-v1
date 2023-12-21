library(ggplot2)
library(viridis)
library(gridExtra)

#visualize overall colonization using barplots
factor(col.overall$region, levels = c("TOP", "MID", "BOT")) -> col.overall$region
p1 <- ggplot(col.overall, aes(x=region, y=count_density^0.5, fill = region)) + 
  geom_boxplot()+
  geom_jitter(width = 0.1)+
  theme_classic()+
  theme(text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(margin = margin(r = 20)),
        legend.position = "none")+
  xlab("")+
  ylab("count density ^ 1/2")


p2 <- ggplot(col.overall, aes(x=region, y=percent_colonization^(1/3), fill = region)) + 
  geom_boxplot()+
  geom_jitter(width = 0.1)+
  theme_classic()+
  theme(text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(margin = margin(r = 20)),
        legend.position = "none")+
  xlab("")+
  ylab("percent colonization ^ 1/3")

grid.arrange(p2, p1, nrow = 1)

#visualize class size
factor(seg.amf.sum.class$region, levels = c("TOP", "MID", "BOT")) -> seg.amf.sum.class$region
ggplot(seg.amf.sum.class, aes(x=region, y=log(avg.size), fill = class)) +
  scale_fill_viridis_d()+
  geom_boxplot()+
  facet_grid(~class)+
  theme_classic()+
  theme(text = element_text(size = 16),
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(margin = margin(r = 20), vjust = -0.25),
        legend.position = "none")+
  xlab("")+
  ylab("log(class size)")

factor(col.class.count2$region, levels = c("TOP", "MID", "BOT")) -> col.class.count2$region
ggplot(col.class.count2, aes(x=region, y=count_density^(1/4), fill = class2)) +
  scale_fill_viridis_d()+
  geom_boxplot()+
  facet_grid(~class2)+
  theme_classic()+
  theme(text = element_text(size = 16),
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(margin = margin(r = 20), vjust = -0.25),
        legend.position = "none")+
  xlab("")+
  ylab("count density^1/4")
