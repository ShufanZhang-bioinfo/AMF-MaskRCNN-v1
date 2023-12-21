###########################################################################################
# function for summarizing random intercept models
###########################################################################################
rintercept.model.stat <- function(lmermodel, nullmodel){
  fxef <- as.data.frame(fixef(lmermodel))
  rownames(fxef) -> fxef$stat.names
  fxef <- fxef[,c(2,1)]
  
  #intercept <- fixef(lmermodel)[1]
  #arb.count <- fixef(lmermodel)[2]
  #arb.size <- fixef(lmermodel)[3]
  #reg.top <- fixef(lmermodel)[4]
  #reg.bot <- fixef(lmermodel)[5]
  
  var.accession.intercept <- c()
  var.sample <- c()
  var.fixed <- c()
  
  pcv.accession <- c()
  pcv.sample <- c()
  
  intracc <- c()
  
  rsqr.mrg <- c()
  rsqr.cnd <- c()
  
  aic <- c()
  bic <- c()
  dev <- c()
  
  #random effects
  var.accession.intercept <- VarCorr(lmermodel)$accession[1]
  var.sample <- sigma(lmermodel)^2
  var.fixed <- var(as.vector(fixef(lmermodel) %*% t(model.matrix(lmermodel))))
  
  pcv.accession <- ((VarCorr(nullmodel)$accession[1] - VarCorr(lmermodel)$accession[1]) / VarCorr(nullmodel)$accession[1])
  pcv.sample <- ((sigma(nullmodel)^2 - sigma(lmermodel)^2)/sigma(nullmodel)^2) #pcv.sample is the same as pcv.residual
  
  intracc <- performance::icc(lmermodel)[[1]]
  # While the adjusted ICC only relates to the random effects, the unadjusted ICC also takes the fixed effects variances into account, more precisely, the fixed effects variance is added to the denominator of the formula to calculate the ICC (see Nakagawa et al. 2017).
  
  a <- anova(nullmodel, lmermodel, test = "LRT")
  aic <- a$AIC[2]
  bic <- a$BIC[2]
  dev <- a$deviance[2]
  
  rsqr.mrg <- r.squaredGLMM(lmermodel)[1]
  rsqr.cnd <- r.squaredGLMM(lmermodel)[2]
  
  stat.names <- c("Variance of fixed effects", "Variance between accessions", "Variance between root samples",   
                  "Between accessions", "Between root samples", 
                  " ","Marginal R2", "Conditional R2","AIC", "BIC", "Deviance")
  stat.values <- c(var.fixed, var.accession.intercept, var.sample,  
                   pcv.accession, pcv.sample,
                   intracc, rsqr.mrg, rsqr.cnd, aic, bic, dev)
  output <- data.frame(stat.names, stat.values)
  colnames(fxef) <- colnames(output)
  output <- rbind(fxef, output)
  #colnames(output)[2] <- deparse(substitute(lmermodel))
  output[output == 0] <- NA
  print(output)
}

###########################################################################################
# function for summarizing random slope models
###########################################################################################

rslope.model.stat <- function(lmermodel, nullmodel){
  ##fix effects
  fxef <- as.data.frame(fixef(lmermodel))
  rownames(fxef) -> fxef$stat.names
  fxef <- fxef[,c(2,1)]
  
  #random effects
  v <- as.data.frame(VarCorr(lmermodel))
  vn <- as.data.frame(VarCorr(nullmodel))
  
  pcv.accession <- c()
  pcv.sample <- c()
  vpc.accession <- c()
  rsqr.mrg <- c()
  rsqr.cnd <- c()
  aic <- c()
  bic <- c()
  dev <- c()
  
  median.x <- median(model.matrix(lmermodel)[,2]) # calling the value of arbc
  var.accession <- v$vcov[1] + v$vcov[2]*(median.x^2) + 2*v$vcov[3]*median.x
  var.accession.null <- vn$vcov[1]
  var.sample <- v$vcov[4]
  var.fixed <- var(as.vector(fixef(lmermodel) %*% t(model.matrix(lmermodel))))
  
  #pcv
  pcv.accession <- ((var.accession.null - var.accession) / var.accession.null)
  pcv.sample <- ((vn$vcov[2] - v$vcov[4])/vn$vcov[2]) #pcv.sample is the same as pcv.residual
  
  # variance partition coef
  vpc.accession <- var.accession / (var.accession + sigma(lmermodel)^2)
  
  a <- anova(nullmodel, lmermodel, test = "LRT")
  aic <- a$AIC[2]
  bic <- a$BIC[2]
  dev <- a$deviance[2]
  
  rsqr.mrg <- r.squaredGLMM(lmermodel)[1]
  rsqr.cnd <- r.squaredGLMM(lmermodel)[2]
  
  #output statistics
  stat.names <- c("Variance of fixed effects", "Variance between accessions", "Variance between root samples", 
                  "Between accessions", "Between root samples", 
                  " ","Marginal R2", "Conditional R2","AIC", "BIC", "Deviance")
  stat.values <- c(var.fixed, var.accession, var.sample, pcv.accession, pcv.sample,
                   vpc.accession, rsqr.mrg, rsqr.cnd, aic, bic, dev)
  
  output <- data.frame(stat.names, stat.values)
  colnames(fxef) <- colnames(output)
  output <- rbind(fxef, output)
  #colnames(output)[2] <- deparse(substitute(lmermodel))
  output[output == 0] <- NA
  print(output)
}

###########################################################################################
# build models for overall percent colonization
###########################################################################################
library(lme4)
library(MuMIn)
# model fit percent colonization random intercept
pc.empty <- lmer(percent_colonization^(1/3) ~ (1|accession), data = col.overall.class.predictors.scaled, REML = TRUE)
pc.rin.rg <- lmer(percent_colonization^(1/3) ~ region + (1|accession), data = col.overall.class.predictors.scaled, REML = TRUE)
pc.rin.arbc.rg <- lmer(percent_colonization^(1/3) ~ arb_count_scaled + region + (1|accession), data = col.overall.class.predictors.scaled, REML = TRUE)
pc.rin.arbs.rg <- lmer(percent_colonization^(1/3) ~ arb_size_scaled + region + (1|accession), data = col.overall.class.predictors.scaled, REML = TRUE)
pc.rin.arbcs.rg <- lmer(percent_colonization^(1/3) ~ arb_count_scaled + arb_size_scaled + region + (1|accession), data = col.overall.class.predictors.scaled, REML = TRUE)

# model fit percent colonization random intercept and slope
pc.rin.arbcs.rg.slope.arbc <- lmer(percent_colonization^(1/3) ~ arb_count_scaled + arb_size_scaled + region + (arb_count_scaled|accession), data = col.overall.class.predictors.scaled, REML = TRUE)
# model with unscaled arb count as fix effect failed to converge
#pc.rin.arbcs.rg.slope.rg <- lmer(percent_colonization^(1/3) ~ arb_count_scaled + arb_size_scaled + region + (region|accession), data = col.overall.class.predictors.scaled, REML = TRUE)
#anova(pc.rin.arbcs.rg,pc.rin.arbcs.rg.slope.rg,test = "LRT")
# root region effect is not different between sibling lines

###########################################################################################
# build models summary tables
###########################################################################################

# build model statistic table for random intercept models
modeltable <- rintercept.model.stat(pc.empty,pc.empty)
#tmp <- rintercept.model.stat(pc.rin.rg,pc.empty)

#tmp <- rintercept.model.stat(pc.rin.arbs.rg,pc.empty)
#tmp$id <- seq(1:nrow(tmp))
#modeltable <- merge(modeltable, tmp, by = "stat.names", all.y = TRUE)
#modeltable <- modeltable[order(modeltable$id),]
n=0
for (i in list(pc.rin.rg, pc.rin.arbs.rg, pc.rin.arbc.rg, pc.rin.arbcs.rg)){
  n=n+1
  print(n)
  tmp <- rintercept.model.stat(i, pc.empty)
  if (n == 4){
    tmp$id <- seq(1:nrow(tmp))
    modeltable <- merge(modeltable, tmp, by = "stat.names", sort = FALSE, all = TRUE)
    modeltable <- modeltable[order(modeltable$id),]
    modeltable <- modeltable[,-7]
  }else{
    modeltable <- merge(modeltable, tmp, by = "stat.names", sort = FALSE, all = TRUE)
  }
}

colnames(modeltable) <- c("stat.names","Null model", "Model with one fixed effect", "Model A with two fixed effects",
                          "Model B with two fixed effects", "Model A with three fixed effects")

# build model statistic table for random intercept + slope models
#modeltable2 <- rintercept.model.stat(pc.empty,pc.empty)
#modeltable2[7,1] <- "Variance Partition Coefficient (VPC)"

modeltable2 <- rslope.model.stat(pc.rin.arbcs.rg.slope.arbc, pc.empty)
modeltable3 <- merge(modeltable, modeltable2, by = "stat.names", sort = FALSE)

colnames(modeltable3) <- c("","Null model", "Model with one fixed effect", "Model A with two fixed effects",
                           "Model B with two fixed effects", "Model A with three fixed effects", "Model B with three fixed effects")
modeltable3[1,1] <- "Intercept"

###########################################################################################
# get coefficients to plot fitted models
###########################################################################################

# get model coefficients for three var random intercept models
coef.pc.rin.arbcs.rg <- as.data.frame(coef(pc.rin.arbcs.rg)$accession)
coef.pc.rin.arbcs.rg$TOP <- coef.pc.rin.arbcs.rg$`(Intercept)` + coef.pc.rin.arbcs.rg$regionTOP
#coef.pc.rin.arbcs.rg$BOT <- coef.pc.rin.arbcs.rg$`(Intercept)` + coef.pc.rin.arbcs.rg$regionBOT
coef.pc.rin.arbcs.rg$accession <- row.names(coef.pc.rin.arbcs.rg)

gather(coef.pc.rin.arbcs.rg[,-c(4,5)], key = "int.type", value = "intercept", -accession, -arb_count_scaled, -arb_size_scaled) -> tmp
gather(tmp, key = "slope.type", value = "slope",-accession, -int.type, -intercept) -> coef.pc.rin.arbcs.rg.long
factor(coef.pc.rin.arbcs.rg.long$int.type, levels = c("TOP", "(Intercept)")) -> coef.pc.rin.arbcs.rg.long$int.type
levels(coef.pc.rin.arbcs.rg.long$int.type) <- c("TOP","Intercept")  

# get model coefficients for random intercept + slope models
coef.pc.rin.arbcs.rg.slope.arbc <- as.data.frame(coef(pc.rin.arbcs.rg.slope.arbc)$accession)
#coef.pc.rin.arbcs.rg.slope.arbc$BOT <- coef.pc.rin.arbcs.rg.slope.arbc$`(Intercept)` + coef.pc.rin.arbcs.rg.slope.arbc$regionBOT
coef.pc.rin.arbcs.rg.slope.arbc$TOP <- coef.pc.rin.arbcs.rg.slope.arbc$`(Intercept)` + coef.pc.rin.arbcs.rg.slope.arbc$regionTOP
coef.pc.rin.arbcs.rg.slope.arbc$accession <- row.names(coef.pc.rin.arbcs.rg.slope.arbc)

gather(coef.pc.rin.arbcs.rg.slope.arbc[,-c(4,5)], key = "int.type", value = "intercept", -accession, -arb_count_scaled, -arb_size_scaled) -> tmp
gather(tmp, key = "slope.type", value = "slope",-accession, -int.type, -intercept) -> coef.pc.rin.arbcs.rg.slope.arbc.long
#model.predicted.pc.rin.arb.rg.long$yhat = model.predicted.pc.rin.arb.rg.long$intercept + model.predicted.pc.rin.arb.rg.long$x * model.predicted.pc.rin.arb.rg.long$arb_count_scaled 
factor(coef.pc.rin.arbcs.rg.slope.arbc.long$int.type, levels = c("TOP", "(Intercept)")) -> coef.pc.rin.arbcs.rg.slope.arbc.long$int.type
levels(coef.pc.rin.arbcs.rg.slope.arbc.long$int.type) <- c("TOP", "Intercept") 

# order accession
coef.pc.rin.arbcs.rg.long$accession <- factor(coef.pc.rin.arbcs.rg.long$accession, levels = accession_order_pc)
coef.pc.rin.arbcs.rg.slope.arbc.long$accession <- factor(coef.pc.rin.arbcs.rg.slope.arbc.long$accession, levels = accession_order_pc)

# plot the accession variance as a function of x variable
v.pc.slope <- as.data.frame(VarCorr(pc.rin.arbcs.rg.slope.arbc))
range.pc.slope <- range(model.matrix(pc.rin.arbcs.rg.slope.arbc)[,2])
x.arbc <- seq(range.pc.slope[1], range.pc.slope[2]+0.1, 0.1)
var.accession.pc <- v.pc.slope$vcov[1] + v.pc.slope$vcov[2]*(x.arbc^2) + 2*v.pc.slope$vcov[3]*x.arbc
#plot(x.arbc, var.accession.pc)


save.image("workspace_modeltable.rds")


