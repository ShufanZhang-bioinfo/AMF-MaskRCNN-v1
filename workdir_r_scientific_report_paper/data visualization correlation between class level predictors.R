
library(GGally)
ggpairs(col.overall.class.predictors.scaled, columns = 7:16, aes(colour = region))
full_opc_class <- lmer(percent_colonization^(1/3) ~  region + (1 | accession) + (1 | accession:region), data=col_overall, REML = FALSE)  #insensitive to the order of strain and region in the equation 
summary(full_opc)

cbind(colscale(col_overall_sample_predictors[,7:16], center = TRUE, scale = FALSE))

for (i in 2:6){
  hist(col_class_count_spread[,i])
}

### correlation analysis of class-level predictors
library(reshape)
library(ggplot2)
#col.overall.class.predictors.scaled[is.na(col.overall.class.predictors.scaled)] <- 0
col.overall.class.predictors.scaled[,c(7:16)] -> col_corr
colnames(col_corr) <- c("arb_count", "exH_count", "inH_count",
                        "sp_count",  "ves_count", "arb_size", 
                        "exH_size",  "inH_size",  "sp_size",  
                        "ves_size" )
#col_class_corr <- col_class_corr[,c(2,1,3)]
cormat <- cor(col_corr,use ="complete.obs")

# Create a ggheatmap
ggheatmap <- ggplot(as.data.frame.table(cormat), aes(Var2, Var1, fill = Freq))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(text = element_text(size = 16))+
  coord_fixed()

# Print the heatmap
print(ggheatmap)
ggheatmap + 
  #geom_text(aes(Var2, Var1, label = cortest), color = "black", size = 8) +
  theme(
    axis.text = element_text(size = 20),
    axis.text.x = element_text(angle = 90,hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    #legend.position = c(0.6, 0.7),
    legend.position = "none",
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))

colnames(col_corr) <- c("arb_count", "exH_count", "inH_count",
                        "sp_count",  "ves_count", "arb_size", 
                        "exH_size",  "inH_size",  "sp_size",  
                        "ves_size" )
ggpairs(col_corr) +
  theme_minimal()

library(leaps)
col_overall_sample_predictors[,c(5,7:16)] -> tmp
modelselect <- regsubsets(percent_colonization~., nbest = 2, data = tmp)
summary(modelselect)


plot(pc.arb.slope, type = c("p", "smooth"))
library(lattice)
qqmath(pc.arb.slope, id = 0.05)
