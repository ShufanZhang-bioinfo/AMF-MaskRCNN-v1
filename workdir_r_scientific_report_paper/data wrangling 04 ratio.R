logitTransform <- function(p) { log(p/(1-p)) }

# creating %mutualistic structure
#col.overall.class.predictors.scaled$ratio <- (col.overall.class.predictors.scaled$arb_count_scaled + col.overall.class.predictors.scaled$exH_count_scaled) / (col.overall.class.predictors.scaled$arb_count_scaled + col.overall.class.predictors.scaled$exH_count_scaled + col.overall.class.predictors.scaled$inH_count_scaled + col.overall.class.predictors.scaled$sp_count_scaled + col.overall.class.predictors.scaled$ves_count_scaled)
factor(col.overall.class.predictors$region, levels = c("MID", "TOP", "BOT")) -> col.overall.class.predictors$region
col.overall.class.predictors$ratio <- (col.overall.class.predictors$arb_count + col.overall.class.predictors$exH_count) / (col.overall.class.predictors$arb_count + col.overall.class.predictors$exH_count + col.overall.class.predictors$inH_count + col.overall.class.predictors$sp_count + col.overall.class.predictors$ves_count)

# transform proportion data using logit transformation
t.ratio <- logitTransform(col.overall.class.predictors$ratio)
hist(t.ratio,30)
col.overall.class.predictors$logit_ratio <- t.ratio


library(emmeans) 
emmeans(pc.rin.arbcs.rg.slope.arbc, pairwise ~ region)
emmeans(cd.rin.5fixcd.slope.exh, pairwise ~ region)

hist(col.overall.class.predictors.scaled$ratio,30)
range(col.overall.class.predictors.scaled$ratio)

hist(col.overall.class.predictors$ratio,20)
range(col.overall.class.predictors$ratio)

# building model 
ratio.empty <- lmer(logit_ratio ~ (1|accession), data = col.overall.class.predictors, REML = TRUE)
ratio.rg <- lmer(logit_ratio ~ region + (1|accession), data = col.overall.class.predictors, REML = TRUE)
#ratio.rg.slope <- lmer(logit_ratio ~ region + (region|accession), data = col.overall.class.predictors, REML = TRUE)
anova(ratio.empty,ratio.rg, test = "LRT")
#anova(ratio.rg, ratio.rg.slope, test = "LRT")
summary(ratio.rg)
performance::icc(ratio.rg)
r.squaredGLMM(ratio.rg)

# tuckey test for between region differences 
emmeans(ratio.rg, pairwise ~ region)

# calculate accession order for ratio using transformed data
#seg.amf %>% group_by(accession, class) %>% summarise(total.count = n()) -> tmp
#col.class.count2 %>% group_by(accession, class2) %>% summarise(total.count = sum(count)) -> tmp
#spread(tmp, key="class2", value="total.count") -> tmp2 
#tmp2$ratio <- (tmp2$arb + tmp2$exH) / (tmp2$arb+tmp2$exH+tmp2$inH+tmp2$sp+tmp2$ves)
#tmp2[order(tmp2$ratio, decreasing = TRUE),] -> col.accession.ratio
#col.accession.ratio$accession -> accession_order_ratio

library(dplyr)
col.overall.class.predictors %>% group_by(accession) %>% summarise(avg_logit_ratio = mean(logit_ratio)) -> tmp
tmp[order(tmp$avg_logit_ratio,decreasing = TRUE),] -> tmp
tmp$accession -> accession_order_ratio

merge(col.accession.avg, tmp, by = "accession") -> col.accession.avg2
col.accession.avg2$avg_trans_cd2 <- col.accession.avg2$avg_trans_cd*10^3
gather(col.accession.avg2[,c(1,2,6,5)], key = "phenotype", value = "value",-accession) -> col.accession.avg.long

# visualize ranking of lines 
ggplot(col.accession.avg.long, aes(x=factor(accession, levels = accession_order_pc), y=value,fill=phenotype))+   #arbuscule presents either in forest or isolated bushes, which give rise to a large variation in pixel density
  geom_col(position = position_dodge(width = 0.8))+
  #facet_grid(~region)+
  #scale_fill_viridis_d()+
  theme_classic()+
  #my_theme2+
  xlab("")

ggpairs(col.accession.avg2, columns = c(2,6,5), upper = list(continuous = wrap("cor", method = "spearman")))

col.overall.class.predictors.scaled[,c(1:6,17)] -> col.overall2
col.overall2$percent_colonization^(1/3) -> col.overall2$trans_pc
col.overall2$count_density^(1/2) -> col.overall2$trans_cd
t.ratio -> col.overall2$trans_ratio

ggpairs(col.overall.class.predictors, columns = 18:20, 
        upper = list(continuous = wrap("cor", method = "pearson")),
        columnLabels = c("Percent colonization", "Count density", "Proportion NE"))+
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank()  #remove y axis ticks
        )
cor.test(col.overall.class.predictors$trans_percent_colonization, col.overall.class.predictors$trans_count_density,method = "pearson")
cor.test(col.overall.class.predictors$trans_percent_colonization, col.overall.class.predictors$logit_ratio, method = "pearson")
cor.test(col.overall.class.predictors$trans_count_density, col.overall.class.predictors$logit_ratio, method = "pearson")


cor.test(col.accession2$total.pc, col.accession2$total.cd,method = "spearman")
cor.test(col.accession2$total.pc, col.accession2$ratio, method = "spearman")
cor.test(col.accession2$total.cd, col.accession2$ratio, method = "spearman")

cor.test(col.accession.avg2$avg_trans_pc, col.accession.avg2$avg_trans_cd,method = "spearman")
cor.test(col.accession.avg2$avg_trans_pc, col.accession.avg2$avg_logit_ratio, method = "spearman")
cor.test(col.accession.avg2$avg_trans_cd, col.accession.avg2$avg_logit_ratio, method = "spearman")

#=======================================================================
#       multiple comparison table
#=======================================================================
as.data.frame(emmeans(pc.rin.arbcs.rg.slope.arbc, pairwise ~ region)$contrast)[,c(1,2,6)] -> tmp1
as.data.frame(emmeans(cd.rin.5fixcd.slope.exh, pairwise ~ region)$contrast)[,c(1,2,6)] -> tmp2
as.data.frame(emmeans(ratio.rg, pairwise ~ region)$contrast)[,c(1,2,6)] -> tmp3

merge(tmp1,tmp2,by = "contrast") -> tmp1
merge(tmp1,tmp3,by = "contrast") -> modeltable.region
#modeltable.class.traits <- as.data.frame(matrix(rep(1,9), nrow = 3, ncol = 3))
#modeltable.class.traits[,1] <- c("BOT - MID", "BOT - TOP", "MID - TOP")
for (i in class_trait_list2){
  tmp <- col.overall.class.predictors.scaled.long[col.overall.class.predictors.scaled.long$traits == i,]
  tmp.model <- lmer(value ~ region + (1|accession), data = tmp, REML = TRUE)
  tmp2 <- as.data.frame(emmeans(tmp.model, pairwise ~ region)$contrast)[,c(2,6)]
  modeltable.class.traits <- cbind(modeltable.class.traits, tmp2)
}
modeltable.class.traits[,c(1:11)] -> modeltable.class.count.traits
modeltable.class.traits[,c(1,12:21)] -> modeltable.class.size.traits
