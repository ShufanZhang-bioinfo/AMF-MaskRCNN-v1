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
  
  median.x <- median(model.matrix(lmermodel)[,2])
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
# build models for overall count density
###########################################################################################

# random intercept models
cd.empty <- lmer(count_density^0.5 ~ (1|accession), data = col.overall.class.predictors.scaled)
cd.rin.rg <- lmer(count_density^0.5 ~ region + (1|accession), data = col.overall.class.predictors.scaled)
cd.rin.5fix <- lmer(count_density^0.5 ~ region + arb_count_scaled + exH_count_scaled +
                      sp_count_scaled + ves_count_scaled + (1|accession), data = col.overall.class.predictors.scaled)

# random intercept models testing interaction effects
cd.rin.int.arb <- lmer(count_density^0.5 ~ region*arb_count_scaled + exH_count_scaled +
                         sp_count_scaled + ves_count_scaled + (1|accession), data = col.overall.class.predictors.scaled)
cd.rin.int.exH <- lmer(count_density^0.5 ~ region*exH_count_scaled + arb_count_scaled + 
                         sp_count_scaled + ves_count_scaled + (1|accession), data = col.overall.class.predictors.scaled)
cd.rin.int.sp <- lmer(count_density^0.5 ~ region*sp_count_scaled + arb_count_scaled + 
                        exH_count_scaled + ves_count_scaled + (1|accession), data = col.overall.class.predictors.scaled)
cd.rin.int.ves <- lmer(count_density^0.5 ~ region*ves_count_scaled + arb_count_scaled + 
                         exH_count_scaled + sp_count_scaled + (1|accession), data = col.overall.class.predictors.scaled)
#none of the interactions were significant

# random intercept and slope models 
#cd.slope.rg <- lmer(count_density^0.5 ~ region + arb_count_scaled + exH_count_scaled +
#                       sp_count_scaled + ves_count_scaled + (region|accession), data = col.overall.class.predictors.scaled) #is.singular
cd.slope.arb <- lmer(count_density^0.5 ~ region + arb_count_scaled + exH_count_scaled +
                       sp_count_scaled + ves_count_scaled + (arb_count_scaled|accession), data = col.overall.class.predictors.scaled)
cd.slope.exh <- lmer(count_density^0.5 ~ region + arb_count_scaled + exH_count_scaled +
                       sp_count_scaled + ves_count_scaled + (exH_count_scaled|accession), data = col.overall.class.predictors.scaled)
cd.slope.sp <- lmer(count_density^0.5 ~ region + arb_count_scaled + exH_count_scaled +
                      sp_count_scaled + ves_count_scaled + (sp_count_scaled|accession), data = col.overall.class.predictors.scaled)
cd.slope.ves <- lmer(count_density^0.5 ~ region + arb_count_scaled + exH_count_scaled +
                       sp_count_scaled + ves_count_scaled + (ves_count_scaled|accession), data = col.overall.class.predictors.scaled)
# model with random slope for exH is the only one significantly different 

###########################################################################################
# build models summary tables
###########################################################################################

# build model statistic table for random intercept models
modeltable.cd <- rintercept.model.stat(cd.empty,cd.empty)
#tmp <- rintercept.model.stat(cd.rin.rg,cd.empty)

#tmp <- rintercept.model.stat(cd.rin.arbs.rg,cd.empty)
#tmp$id <- seq(1:nrow(tmp))
#modeltable <- merge(modeltable, tmp, by = "stat.names", all.y = TRUE)
#modeltable <- modeltable[order(modeltable$id),]
n=0
for (i in list(cd.rin.rg, cd.rin.5fix)){
  n=n+1
  #print(n)
  tmp <- rintercept.model.stat(i, cd.empty)
  if (n == 2){
    tmp$id <- seq(1:nrow(tmp))
    modeltable.cd <- merge(modeltable.cd, tmp, by = "stat.names", sort = FALSE, all = TRUE)
    modeltable.cd <- modeltable.cd[order(modeltable.cd$id),]
    modeltable.cd <- modeltable.cd[,-5]
  }else{
    modeltable.cd <- merge(modeltable.cd, tmp, by = "stat.names", sort = FALSE, all = TRUE)
  }
}

colnames(modeltable.cd) <- c("stat.names","Null model", "Model with one fixed effect", "Model with five fixed effects")

# build model statistic table for random intercept + slope models
#modeltable2 <- rintercept.model.stat(cd.empty,cd.empty)
#modeltable2[7,1] <- "Variance Partition Coefficient (Vcd)"

modeltable.cd2 <- rslope.model.stat(cd.slope.exh, cd.empty)
modeltable.cd3 <- merge(modeltable.cd, modeltable.cd2, by = "stat.names", sort = FALSE)

# change the scale of estimated variables for display
#power9 <- function(x) {x*(10^9)}
#apply(modeltable.cd3[1:10,2:5], 2, "*")
modeltable.cd3[1:10,] %>% 
  mutate_at(vars(-stat.names),
            .funs = funs(. * 10^9)) -> modeltable.cd3[1:10,]

colnames(modeltable.cd3) <- c("","Null model", "Model with one fixed effect", "Model A with five fixed effects","Model B with five fixed effects")

modeltable.cd3[1,1] <- "Intercept"
###########################################################################################
# get coefficients to plot fitted models
###########################################################################################

# get model coefficients for three var random intercept models
coef.cd.rin.5fix <- as.data.frame(coef(cd.rin.5fix)$accession)
coef.cd.rin.5fix$TOP <- coef.cd.rin.5fix$`(Intercept)` + coef.cd.rin.5fix$regionTOP
coef.cd.rin.5fix$accession <- row.names(coef.cd.rin.5fix)

gather(coef.cd.rin.5fix[,-c(2,3)], key = "int.type", value = "intercept", -accession, -arb_count_scaled, -exH_count_scaled,
       -sp_count_scaled,-ves_count_scaled) -> tmp
gather(tmp, key = "slope.type", value = "slope",-accession, -int.type, -intercept) -> coef.cd.rin.5fix.long
factor(coef.cd.rin.5fix.long$int.type, levels = c("TOP", "(Intercept)")) -> coef.cd.rin.5fix.long$int.type
levels(coef.cd.rin.5fix.long$int.type) <- c("TOP","Intercept")  

# get model coefficients for random intercept + slope models
coef.cd.slope.exh <- as.data.frame(coef(cd.slope.exh)$accession)
coef.cd.slope.exh$TOP <- coef.cd.slope.exh$`(Intercept)` + coef.cd.slope.exh$regionTOP
coef.cd.slope.exh$accession <- row.names(coef.cd.slope.exh)

gather(coef.cd.slope.exh[,-c(2,3)], key = "int.type", value = "intercept", -accession, -arb_count_scaled, -exH_count_scaled,
       -sp_count_scaled,-ves_count_scaled) -> tmp
gather(tmp, key = "slope.type", value = "slope",-accession, -int.type, -intercept) -> coef.cd.slope.exh.long
#model.predicted.cd.rin.arb.rg.long$yhat = model.predicted.cd.rin.arb.rg.long$intercept + model.predicted.cd.rin.arb.rg.long$x * model.predicted.cd.rin.arb.rg.long$arb_count_scaled 
factor(coef.cd.slope.exh.long$int.type, levels = c("TOP", "(Intercept)")) -> coef.cd.slope.exh.long$int.type
levels(coef.cd.slope.exh.long$int.type) <- c("TOP", "Intercept") 

# order accession
coef.cd.rin.5fix.long$accession <- factor(coef.cd.rin.5fix.long$accession, levels = accession_order_pc)
coef.cd.slope.exh.long$accession <- factor(coef.cd.slope.exh.long$accession, levels = accession_order_pc)

# plot the accession variance as a function of x variable
v.cd.slope <- as.data.frame(VarCorr(cd.slope.exh))
range.cd.slope <- range(model.matrix(cd.slope.exh)[,5])
x.exh <- seq(range.cd.slope[1], range.cd.slope[2]+0.1, 0.1)
var.accession.cd <- v.cd.slope$vcov[1] + v.cd.slope$vcov[2]*(x.exh^2) + 2*v.cd.slope$vcov[3]*x.exh

#plot(x.exh, var.accession.cd)

save.image("workspace_modeltable.rds")


