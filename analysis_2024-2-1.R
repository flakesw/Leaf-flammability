#*******************************************************************************
#* Analysis of flammability
#* Sam Flake, 2 Feb 2024
#*******************************************************************************

# Packages and options ---------------------------------------------------
library("effects")
library("MuMIn")
library("lme4")
library("tidyverse")
library("lmerTest")

options(warnPartialMatchDollar = FALSE)
options(warn = 1)

# Import data ------------------------------------------------------------------

leaf_data <- read.csv("./clean data/leaf_level_data.csv")
leaf_sp_level <- read.csv("./clean data/leaf_sp_level_data.csv")
burn_level <- read.csv("./clean data/burn_level_data.csv")
burn_sp <- read.csv("./clean data/burn_sp_level_data.csv")


#*******************************************************************************
# Exploration of fire measurements ---------------------------------------------
#*******************************************************************************

summary(lm(burn_level$Fireline.intensity.kWm.1 ~ burn_level$litter_spec_vol + burn_level$SLAcm2g.1 + burn_level$FG))

fire_correlations <- as.matrix(cor(burn_level[, c("litter_spec_vol", "Average.Flame.Height", "PercentCombustion", 
                                                "HeatFlux_SBG01_Max", "Heatflux_flat_Max", "Temp_C_15_10_Max",        
                                                "Temp_C_30_10_Max", "Temp_C_45_10_Max", "Temp_C_30_5_Max",         
                                                "Temp_C_30_20_Max", "Temp_C_IR_Max", "Mean.ROS.ms.1",
                                                "Fireline.intensity.kWm.1")]))

write.table(fire_correlations, "./Outputs/fire_correlations.csv", sep = ",", quote = FALSE, row.names = TRUE, col.names=NA)

#PCA
species_burn <- burn_level %>%
  group_by(Species) %>%
  summarise(across(where(is.numeric), mean)) %>%
  select(c("litter_spec_vol", "Average.Flame.Height", "PercentCombustion", 
           "HeatFlux_SBG01_Max", "Heatflux_flat_Max", "Temp_C_15_10_Max",        
           "Temp_C_30_10_Max", "Temp_C_45_10_Max", "Temp_C_30_5_Max",         
           "Temp_C_30_20_Max", "Temp_C_IR_Max", "Mean.ROS.ms.1",
           "Fireline.intensity.kWm.1")) %>%
  scale()

species_burn

fire_intensity_pca <- princomp(species_burn)
summary(fire_intensity_pca)
biplot(fire_intensity_pca)


#run the PCA
pca_results <- vegan::rda(species_burn, scale = TRUE)
trait_names <- gsub("_", " ", rownames(pca_results$CA$v))
sp_scores <- vegan::scores(pca_results, choices = c(1:10), display = c("sites"))
trait_scores <- vegan::scores(pca_results, choices = c(1:10), display = c("species"))

summary(pca_results)

biplot(pca_results, choices = c(1,2))

set.seed(45750766)

#----------------------
#Create Figure S2 (PCA of flammability traits)

svg(file ="./Outputs/Figure S2 pca.svg", 
          width = 7, 
          height=7,
          pointsize = 9)
plot(sp_scores[, 2] ~ sp_scores[, 1],
     xlab = "PC1", ylab = "PC2",
     xlim = c( (min(c(sp_scores[, 1], trait_scores[, 1]) - .5)), 
               (max(c(sp_scores[, 1], trait_scores[, 1])))),
     ylim = c( (min(c(sp_scores[, 2], trait_scores[, 2]))), 
               (max(c(sp_scores[, 2], trait_scores[, 2])) + 0.4)))

abline(h = 0)
abline(v = 0)

segments(x0 = 0, y0 = 0, x1 = trait_scores[, 1], y1 = trait_scores[, 2])
yoffsets <- c(-0.1, 0, 0, -0.03, 0.03, 0, 0.03, 0, 0, 0, 0, 0.1, 0)
xoffsets <- c(-0.05, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.05, 0)

text(x = trait_scores[, 1]+xoffsets, y = trait_scores[, 2]+yoffsets, labels = trait_names, cex = 0.8)

dev.off()

#*******************************************************************************
#* Scatterplots to demonstrate overall relationships
#*******************************************************************************
#------------------------------------------------
# Fig xx: relationship between litter specific volume and packing ratio
mod <- lm(log(litter_spec_vol) ~ log(packing_ratio), data = burn_sp)
pred <- predict(mod, 
                newdata = data.frame(packing_ratio = seq(0.01, max(burn_sp$packing_ratio), length.out = 100)),
                se.fit = TRUE)
pred.df <- data.frame(litter_spec_vol = exp(pred$fit), #this name has to match the name in burn_sp, no idea why
                      hi = exp(pred$fit + pred$se.fit),
                      low = exp(pred$fit - pred$se.fit),
                      packing_ratio = data.frame(packing_ratio = seq(0.01, max(burn_sp$packing_ratio), length.out = 100)))

g <- ggplot(data = burn_sp, aes(y = litter_spec_vol, x = packing_ratio)) +
  geom_point(color = "steelblue") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  ylab(expression("Litter specific volume (V"["s, litter"]*") (m"^3*" kg"^-1*")")) +
  xlab(expression("Packing ratio (unitless)")) +
  geom_line(data = pred.df, aes(y = litter_spec_vol, x = packing_ratio)) + 
  geom_line(data = pred.df, aes(y = hi, x = packing_ratio), color = "grey", alpha = 0.4) + 
  geom_line(data = pred.df, aes(y = low, x = packing_ratio), color = "grey", alpha = 0.4) +
  geom_ribbon(data = pred.df, aes(ymin = low, ymax = hi, x = packing_ratio), color = "grey", alpha = 0.3)
plot(g)


# Fig xx: Relationship between leaf specific volume and litter specific volume -----------------------------------------------------------
plot(burn_sp$litter_spec_vol ~ burn_sp$hull_specific_vol)
abline(0,1)
abline(0, coef(lm(litter_spec_vol ~ hull_specific_vol + 0, data = burn_sp)))
summary(lm(litter_spec_vol ~ hull_specific_vol + 0, data = burn_sp, offset = hull_specific_vol))

ggplot(data = burn_sp, aes(x = hull_specific_vol, y = litter_spec_vol)) +
  geom_point(color = "steelblue") + 
  geom_smooth(method='lm', color = "black", formula= (y ~ x + 0)) + 
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  ylab(expression("Litter specific volume (V"["s, litter"]*") (m"^3*" kg"^-1*")")) +
  xlab(expression("Convex hull specific volume (V"["s, leaf"]*") (m"^3*" kg"^-1*")"))



#*******************************************************************************
# Hypothesis 1
#********************************************************************************

# Exploration of leaf size measurements

# Estimating volume/

# For each version, how well does SLA*S*C predict bulk specific volume (i.e. 1/bulk density)? 

plot(leaf_data$S1 ~ leaf_data$C1)
plot(leaf_data$S1 ~ leaf_data$C1H)
plot(leaf_data$C1 ~ leaf_data$C1H) #two measures of leaf curl maybe exponential/sigmoid relationship?
plot(leaf_data$S1 ~ leaf_data$C1H_new)
plot(leaf_data$C1H ~ leaf_data$C1H_new)
plot(leaf_data$C1 ~ leaf_data$C1H_new)

plot(leaf_data$dryHull ~ log(leaf_data$C1H))
plot(log(leaf_data$dryHull) ~ log(leaf_data$C1H_new))


plot(leaf_data$S2 ~ leaf_data$C2)
plot(leaf_data$S2 ~ leaf_data$C2H)
plot(leaf_data$C2 ~ I(leaf_data$C2H/10)) #two measures maybe sigmoid relationship?

plot(leaf_data$S3 ~ leaf_data$C3)
plot(leaf_data$S3 ~ leaf_data$C3H)
plot(log(leaf_data$C3) ~ log(leaf_data$C3H/10)) #two measures mostly linearly related
summary(lm(log(leaf_data$C3) ~ log(leaf_data$C3H/10))) #two measures mostly linearly related



plot(leaf_data$S1 ~ leaf_data$S2) #leaf area is pretty well correlated with L*W
plot(leaf_data$dryHull ~ I(leaf_data$DryL * leaf_data$DryW * leaf_data$DryH)) #convex hull is pretty well correlated with L*W*H
plot(leaf_data$dryHull ~ leaf_data$particle_volume)

#-------------------
## original models
## These have balanced units
mod1 <- lm(log(litter_spec_vol) ~ log(Areacm2^0.5) + log(DryH / (Areacm2^0.5)) + log(SLAcm2g.1), data = burn_sp)
mod1 <- glm(litter_spec_vol ~ log(Areacm2^0.5) + log(DryH / (Areacm2^0.5)) + log(SLAcm2g.1), 
            data = burn_sp,
            family = gaussian(link = "log"))

summary(mod1)
car::vif(mod1)
AIC(mod1)
plot(residuals(mod1) ~ fitted(mod1))

mod2 <- lm(log(litter_spec_vol) ~ log(Areacm2^0.5) + log(((dryHull/1000)^(1/3))/(Areacm2^0.5)) + log(SLAcm2g.1), data = burn_sp)
mod2 <- glm(litter_spec_vol ~ log(Areacm2^0.5) + log(((dryHull/1000)^(1/3))/(Areacm2^0.5)) + log(SLAcm2g.1), 
           data = burn_sp,
           family = gaussian(link = "log"))


summary(mod2)
car::vif(mod2)
AIC(mod2)
plot(residuals(mod2) ~ fitted(mod2))


mod3 <- lm(log(litter_spec_vol) ~ log((DryL * DryW)^0.5) + I(log(DryH / ((DryL * DryW)^0.5))) + log(SLAcm2g.1), data = burn_sp)
summary(mod3)
car::vif(mod3)
AIC(mod3)
plot(residuals(mod3) ~ fitted(mod3))
 

mod4 <- lm(log(litter_spec_vol) ~ log((DryL * DryW)^0.5) + I(log(((dryHull/1000)^(1/3)) / ((DryL * DryW)^0.5))) + log(SLAcm2g.1), data = burn_sp)
summary(mod4)
car::vif(mod4)
AIC(mod4)
plot(residuals(mod4) ~ fitted(mod4))


mod5 <- lm(log(litter_spec_vol) ~ log(DryL) + log(I(DryH/DryL)) + log(SLAcm2g.1) + 0, data = burn_sp)
summary(mod5)
car::vif(mod5)
AIC(mod5)
plot(residuals(mod5) ~ fitted(mod5))
 
mod6 <- lm(log(litter_spec_vol) ~ log(DryL) + log(((dryHull/1000)^(1/3))/DryL) + log(SLAcm2g.1), data = burn_sp)
summary(mod6)
car::vif(mod6)
AIC(mod6)
plot(residuals(mod6) ~ fitted(mod6))


#-----------------------------------------------------
# Figure 5: Make effects plot for volume
#---------------------------------------------------------

temp <- burn_sp %>%
  mutate(C = ((dryHull/1000)^(1/3))/DryL)
mod6 <- lm(log(litter_spec_vol) ~ log(DryL) + log(C) + log(SLAcm2g.1), data = temp)
mod6 <- glm(litter_spec_vol ~ scale(log(DryL)) + scale(log(C)) + scale(log(SLAcm2g.1)), data = temp,
            family = gaussian(link = "log"))
summary(mod6)
car::vif(mod6)
AIC(mod6)
plot(residuals(mod6) ~ fitted(mod6))

eff <- Effect(focal.predictors = c("DryL"), 
              mod = mod6, xlevels = 100,
              partial.residuals = TRUE,
              se = TRUE)
resids <- data.frame(x = eff$data$DryL,
                     resid = exp(eff$residuals + predict(mod6)))
dat <- data.frame(y = exp(eff$fit),
                  lower = exp(eff$lower),
                  upper = exp(eff$upper),
                  var = eff[["x"]][["DryL"]])
poly <- data.frame(x = c(dat$var, rev(dat$var)),
                   y = c(dat$upper, rev(dat$lower)))

S_eff <- ggplot() +
  geom_line(data = dat, aes(x = var, y = y), color = "darkblue") +
  geom_polygon(data = poly, aes(x = x, y = y), fill = "darkblue", alpha = 0.3) +
  geom_point(data = resids, aes(x = x, y = resid), color = "darkblue") +
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=9)) +
  xlab(label = "Leaf size (Dry leaf length, cm)") +
  ylab(label = expression("Litter specific volume (m"^3*" kg"^-1*")"))
plot(S_eff)

eff <- Effect(focal.predictors = c("C"), 
              mod = mod6, xlevels = 100,
              partial.residuals = TRUE,
              se = TRUE
              )
resids <- data.frame(x = eff$data$C,
                     resid = exp(eff$residuals + predict(mod6)))
dat <- data.frame(y = exp(eff$fit),
                  lower = exp(eff$lower),
                  upper = exp(eff$upper),
                  var = eff[["x"]][["C"]])
poly <- data.frame(x = c(dat$var, rev(dat$var)),
                   y = c(dat$upper, rev(dat$lower)))

C_eff <- ggplot() +
  geom_line(data = dat, aes(x = var, y = y), color = "darkblue") +
  geom_polygon(data = poly, aes(x = x, y = y), fill = "darkblue", alpha = 0.3) +
  geom_point(data = resids, aes(x = x, y = resid), color = "darkblue") +
  xlab(expression("Leaf curl (V"^(1/3)*" L"^-1*", unitless)"))  +
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=9))
plot(C_eff)



eff <- Effect(focal.predictors = c("SLAcm2g.1"), 
              mod = mod6, xlevels = 100,
              partial.residuals = TRUE,
              se = TRUE)
resids <- data.frame(x = eff$data$SLAcm2g.1,
                     resid = exp(eff$residuals + predict(mod6)))
dat <- data.frame(y = exp(eff$fit),
                  lower = exp(eff$lower),
                  upper = exp(eff$upper),
                  var = eff[["x"]][["SLAcm2g.1"]])
poly <- data.frame(x = c(dat$var, rev(dat$var)),
                   y = c(dat$upper, rev(dat$lower)))

SLA_eff <- ggplot() +
  geom_line(data = dat, aes(x = var, y = y), color = "darkblue") +
  geom_polygon(data = poly, aes(x = x, y = y), fill = "darkblue", alpha = 0.3) +
  geom_point(data = resids, aes(x = x, y = resid), color = "darkblue")+
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=9)) +
  xlab(expression("Specific leaf area (cm"^2*" g"^-1*")"))
plot(SLA_eff)

theme_set(theme_bw())
theme_update(panel.grid.major = element_blank(),
             panel.grid.minor = element_blank())
effect_combined <- cowplot::plot_grid(S_eff, 
                   C_eff + theme(axis.text.y = element_blank(),
                              axis.title.y = element_blank()), 
                   SLA_eff + theme(axis.text.y = element_blank(),
                             axis.title.y = element_blank()),
                   nrow = 1, ncol = 3,
                   rel_widths = c(1.3, 1, 1),
                   label_x = c(0.3,0.1,0.1),
                   labels = "auto",
                   align = "h")
plot(effect_combined)

ggsave(effect_combined,
       filename = "./Outputs/Fig 5 effects_traits.png",
       device = "png",
       units = "in",
       width = 6,
       height = 2.3)

#-------------------------------------------------------------------------------
# Hypothesis 2: How are fresh leaf traits related to dry traits?
#-----------------------------------------------------------------------------

# correlations
burn_level$change_in_hull <- burn_level$dryHull/burn_level$freshHull
burn_level$aspect_ratio <- burn_level$FreshL/burn_level$FreshW

leaf_sp_level$change_in_hull <- leaf_sp_level$dryHull / leaf_sp_level$freshHull
m <- mean(leaf_sp_level$change_in_hull)
hist((leaf_sp_level$change_in_hull - m)/m)
range((leaf_sp_level$change_in_hull - m)/m)

plot(burn_level$change_in_hull ~ burn_level$SLAcm2g.1)
    summary(lm(burn_level$change_in_hull ~ burn_level$SLAcm2g.1))
plot(burn_level$change_in_hull ~ burn_level$FreshTH) #strong relationship
     summary(lm(burn_level$change_in_hull ~ burn_level$FreshTH)) #significant
plot(burn_level$change_in_hull ~ burn_level$FreshL)
     summary(lm(burn_level$change_in_hull ~ burn_level$FreshL)) #significant
plot(burn_level$change_in_hull ~ burn_level$FreshW)
     summary(lm(burn_level$change_in_hull ~ burn_level$FreshW)) #significant
plot(burn_level$change_in_hull ~ burn_level$aspect_ratio)
     summary(lm(burn_level$change_in_hull ~ burn_level$aspect_ratio)) #not significant
plot(burn_level$change_in_hull ~ burn_level$C1Hfresh) # not very strong correlation
    summary(lm(burn_level$change_in_hull ~ burn_level$C1Hfresh))
plot(burn_level$change_in_hull ~ burn_level$C2Hfresh) #not very strong correlation
    summary(lm(burn_level$change_in_hull ~ burn_level$C2Hfresh))
plot(burn_level$change_in_hull ~ burn_level$C3Hfresh) #stronger
  summary(lm(burn_level$change_in_hull ~ burn_level$C3Hfresh))

change_model <- lm(change_in_hull ~ FreshTH + FreshL + FreshW + C1Hfresh*SLAcm2g.1, data = burn_level)

summary(change_model)

lm.without<-update(change_model, ~. - FreshTH)

rsq::rsq.partial(change_model)


plot(burn_level$change_in_hull ~ burn_level$SLAcm2g.1)
  m1 <- summary(lm(burn_level$change_in_hull ~ burn_level$SLAcm2g.1))
plot(burn_level$change_in_hull ~ burn_level$FreshTH) #strong relationship
  m2 <- summary(lm(burn_level$change_in_hull ~ burn_level$FreshTH)) #significant
plot(burn_level$change_in_hull ~ burn_level$FreshL)
  m3 <- summary(lm(burn_level$change_in_hull ~ burn_level$FreshL)) #significant
plot(burn_level$change_in_hull ~ burn_level$FreshW)
  m4 <- summary(lm(burn_level$change_in_hull ~ burn_level$FreshW)) #significant
plot(burn_level$change_in_hull ~ burn_level$aspect_ratio)
  m5 <- summary(lm(burn_level$change_in_hull ~ burn_level$aspect_ratio)) #not significant
plot(burn_level$change_in_hull ~ burn_level$C3Hfresh) #stronger
  m6 <- summary(lm(burn_level$change_in_hull ~ burn_level$C3Hfresh))

#------------------------------------------
## Figure 6: change in traits effects plots
#--------------------------------------------
  
p1 <- ggplot(data = burn_level, aes(x = SLAcm2g.1, y = change_in_hull)) +
  geom_point(color = "steelblue") + 
  #geom_smooth(method='lm', color = "black", formula= (y ~ x)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  xlab(expression("Specific leaf area (cm"^2*" g"^-1*")")) +
  ylab(expression("Change in volume (proportion)"))+
  # annotate("text", x = -Inf, y = Inf, label = paste("(a)"), vjust = 2, hjust = -1) +
  annotate("text", x = Inf, y = Inf, label = paste("p =", round(m1$coefficients[2,4], 2)), vjust = 2, hjust = 2) + 
  geom_hline(yintercept = 1, linetype = 2)
p2 <- ggplot(data = burn_level, aes(x = FreshTH, y = change_in_hull)) +
  geom_point(color = "steelblue") + 
  geom_smooth(method='lm', color = "black", formula= (y ~ x)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  xlab(expression("Leaf thickness (mm)")) +
  ylab(expression("Change in volume (proportion)"))+
  # annotate("text", x = -Inf, y = Inf, label = paste("(b)"), vjust = 2, hjust = -1) +
  annotate("text", x = Inf, y = Inf, label = paste("p <0.001"), vjust = 2, hjust = 2)+ 
  geom_hline(yintercept = 1, linetype = 2)
p3 <- ggplot(data = burn_level, aes(x = FreshL, y = change_in_hull)) +
  geom_point(color = "steelblue") + 
  geom_smooth(method='lm', color = "black", formula= (y ~ x)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  xlab(expression("Leaf length (cm)")) +
  ylab(expression("Change in volume (proportion)"))+
  # annotate("text", x = -Inf, y = Inf, label = paste("(c)"), vjust = 2, hjust = -1) +
  annotate("text", x = Inf, y = Inf, label = paste("p <0.001"), vjust = 2, hjust = 2)+ 
  geom_hline(yintercept = 1, linetype = 2)
p4 <- ggplot(data = burn_level, aes(x = FreshW, y = change_in_hull)) +
  geom_point(color = "steelblue") + 
  geom_smooth(method='lm', color = "black", formula= (y ~ x)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  xlab(expression("Leaf width (cm)"))+
  ylab(expression("Change in volume (proportion)"))+
  # annotate("text", x = -Inf, y = Inf, label = paste("(d)"), vjust = 2, hjust = -1) +
  annotate("text", x = Inf, y = Inf, label = paste("p <0.001"), vjust = 2, hjust = 2)+ 
  geom_hline(yintercept = 1, linetype = 2)
p5 <- ggplot(data = burn_level, aes(x = aspect_ratio, y = change_in_hull)) +
  geom_point(color = "steelblue") + 
  geom_smooth(method='lm', color = "black", formula= (y ~ x)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  xlab(expression("Aspect ratio (unitless)")) +
  ylab(expression("Change in volume (proportion)"))+
  # annotate("text", x = -Inf, y = Inf, label = paste("(e)"), vjust = 2, hjust = -1) +
  annotate("text", x = Inf, y = Inf, label = paste("p = 0.002"), vjust = 2, hjust = 2)+ 
  geom_hline(yintercept = 1, linetype = 2)
p6 <- ggplot(data = burn_level, aes(x = C3Hfresh, y = change_in_hull)) +
  geom_point(color = "steelblue") + 
  geom_smooth(method='lm', color = "black", formula= (y ~ x)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  xlab(expression("Leaf curl (unitless)")) +
  ylab(expression("Change in volume (proportion)"))+
  # annotate("text", x = -Inf, y = Inf, label = paste("(f)"), vjust = 2, hjust = -1) +
  annotate("text", x = Inf, y = Inf, label = paste("p < 0.001"), vjust = 2, hjust = 2)+ 
  geom_hline(yintercept = 1, linetype = 2)

cowplot::plot_grid(p1, 
                   p2 + theme(axis.text.y = element_blank(),
                                  axis.title.y = element_blank() ), 
                   p3+ theme(axis.text.y = element_blank(),
                             axis.title.y = element_blank() ),
                   p4, 
                   p5+ theme(axis.text.y = element_blank(),
                             axis.title.y = element_blank() ), 
                   p6+ theme(axis.text.y = element_blank(),
                             axis.title.y = element_blank() ),
                   nrow = 2, ncol = 3,
                   labels = "auto",
                   align = "v")


#-------------------------------------------------------------------------------
# Hypothesis 3:

# SEM of leaf traits
#########################################
library("lavaan")
library("semPlot")
library("sem")
# remotes::install_github("maksimrudnev/LittleHelpers")
library("LittleHelpers")



sem_data <- burn_level[burn_level$success, ] %>%
  group_by(Species) %>%
  summarise(across(.cols = where(is.numeric),.fns = mean))

#curl selected from models above
sem_data$curl <- sem_data$dryHull^(1/3)/sem_data$DryL
sem_data$logCurl <- log(sem_data$curl)
sem_data$logL <- log(sem_data$`DryL`)
sem_data$other_ash <- sem_data$All_ash - sem_data$Al.Conc.
sem_data$log_ash <- log(sem_data$All_ash)
sem_data$log_other_ash <- log(sem_data$other_ash)
sem_data$log_al <- log(sem_data$Al.Conc.)
sem_data$S3 <- log(sem_data$S3)
sem_data$C3H <- log(sem_data$C3H)
sem_data$logSLA <- log(sem_data$SLAcm2g.1)

plot(sem_data$SLAcm2g.1 ~ sem_data$S3)
plot(sem_data$SLAcm2g.1 ~ sem_data$C3H)
plot(sem_data$S3 ~ sem_data$C3H)
plot(burn_level$S3 ~ burn_level$C3H)
plot(sem_data$logL ~ sem_data$S3)

sem_data <- sem_data %>%
  mutate(across(where(is.numeric), scale))

# sem_data$FG <- as.numeric(as.factor(sem_data$FG)) - 1

leaf_model <- ' 

            log_litter_spec_vol ~  a*logSLA + d*S3 + f*C3H
            Fireline.intensity.kWm.1 ~ b*log_litter_spec_vol + c*logSLA + e*S3 + g*C3H
            
            #SLA total effect
            SLA_total := c + (a*b)
            SLA_indirect := a*b
            
            #S3 total effect
            S3_total := e + (d*b)
            S3_indirect := d*b
            
            #C3H total effect
            C3H_total := g + (f*b)
            C3H_indirect := f*b
            
            '

fit <- lavaan::sem(leaf_model, data = sem_data, estimator = "ML", fixed.x=FALSE)
semPaths(object = fit, what = "std", layout = "tree", residuals = FALSE, nCharNodes = 10, 
         sizeMan = 12, sizeLat = 7, label.cex = 1.2, edge.label.cex = 1.2)
summary(fit, fit.measures = TRUE)
fitMeasures(fit)

resid(fit)
LittleHelpers::lav_to_graph(fit, file = "Fig 7 leaf_sem.svg")


#-------------------------------------------------------------------------------
#Hypothesis 4 -- Aluminum


######
# aluminum and probability of ignition

binom_mod <- glm(prop_success ~ scale(log(Al.Conc.)) + scale(log(litter_spec_vol)),
                 data = burn_sp,
                 family = binomial(link = "logit"),
                 weights = rep(3,nrow(burn_sp)))
summary(binom_mod)

#ANCOVA
binom_mod <- glm(prop_success ~ accumulator + scale(log(litter_spec_vol)),
                 data = burn_sp,
                 family = binomial(link = "logit"),
                 weights = rep(3,nrow(burn_sp)))
summary(binom_mod)

newdat <- expand.grid(accumulator = c(TRUE, FALSE),
                      litter_spec_vol = seq(0, 0.2, 0.001))
preds_binom <- cbind(preds = predict(binom_mod, newdat), newdat)
plot(boot::inv.logit(preds_binom$preds) ~ newdat$litter_spec_vol)
plot(burn_sp$prop_success ~ burn_sp$litter_spec_vol, col = ifelse(burn_sp$accumulator, "red", "blue"))

##
#intensity relationship
mod2 <- glm(Fireline.intensity.kWm.1 ~ accumulator + scale(log(litter_spec_vol)), 
            data = burn_sp[!is.na(burn_sp$Al.Conc.) & burn_sp$Fireline.intensity.kWm.1 > 0, ],
            family = gaussian(link = "log"))
summary(mod2)

newdat <- expand.grid(Al.Conc. = c(100, 5000),
                      litter_spec_vol = seq(0, 0.2, by = 0.001))
preds <- predict(mod, newdata = newdat)
plot(preds ~ newdat$litter_spec_vol)

#more or less equivalent using binary accumulator variable
newdat <- expand.grid(accumulator = c(TRUE, FALSE),
                      litter_spec_vol = seq(0, 0.2, by = 0.001))
preds_int <- cbind(preds = predict(mod2, newdata = newdat, type = "response"), newdat)
plot((preds_int$preds) ~ newdat$litter_spec_vol)

#------------------
#Fig 8: Aluminum relationships

zero_mod <-  ggplot(data = preds_binom, aes(x = litter_spec_vol, 
                                            y = boot::inv.logit(preds), 
                                            col = accumulator,
                                            linetype = accumulator,
                                            shape = accumulator)) +
  geom_line(linewidth = 1, alpha = 0.5) +
  geom_point(data = burn_sp[!is.na(burn_sp$accumulator), ], 
             aes(x = litter_spec_vol, y = prop_success, col = accumulator),
             alpha = 0.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.title = element_blank()) +
  scale_color_discrete(labels = c("Non-accumulator", "Al accumulator")) +
  scale_linetype_discrete(labels = c("Non-accumulator", "Al accumulator")) +
  scale_shape_discrete(labels = c("Non-accumulator", "Al accumulator")) +
  ylab(expression("Proportion of ignitions succcessful")) +
  xlab("")
plot(zero_mod)

intensity_mod <- ggplot(data = preds_int, aes(x = litter_spec_vol, 
                                              y = preds, 
                                              col = accumulator,
                                              linetype = accumulator,
                                              shape = accumulator)) +
  geom_line(linewidth = 1, alpha = 0.5) +
  geom_point(data = burn_level[complete.cases(burn_level[, c("Fireline.intensity.kWm.1", "litter_spec_vol", "accumulator")]) &
                                 burn_level$Fireline.intensity.kWm.1 > 0, ],
             aes(y = Fireline.intensity.kWm.1,
                 x = litter_spec_vol, 
                 col = accumulator),
             alpha = 0.5)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none") + 
  ylab(expression("Fireline intensity (kW m"^-1*")")) +
  # xlab(expression("Litter specific volume ("*"kg m"^-3*")"))
  xlab("")
plot(intensity_mod)

al_plots <- cowplot::plot_grid(zero_mod + 
                                 theme(legend.position = c(0.7, 0.2))+
                                 theme(axis.text=element_text(size=7),
                                       axis.title=element_text(size=9)), 
                               intensity_mod + 
                                 theme(legend.position="none") +
                                 theme(axis.text=element_text(size=7),
                                       axis.title=element_text(size=9)),
                               align='v', vjust=1, scale = 1, 
                               labels = "auto", nrow = 2, axis = "l")

al_plots <- cowplot::ggdraw(cowplot::add_sub(al_plots, expression("Litter specific volume ("*"kg m"^-3*")"),
                size = 9.5, vpadding = grid::unit(0.5, "lines"), vjust = -0.85, hjust = 0.33))
plot(al_plots)
ggsave(plot = al_plots, 
       filename = "./outputs/al_plots.png",
       device = "png",
       units = "in",
       scale = 1,
       width = 3.2,
       height = 6)

# SEM for Al---------------------------------------------------
al_model <- ' 
            log_litter_spec_vol ~  k*logSLA + l*S3 + m*C3H + b*log_al
            Fireline.intensity.kWm.1 ~ c*log_litter_spec_vol + n*logSLA + a*log_al
            
            #Correlations
            log_al ~~ h*logSLA + i*S3 + j*C3H
            
            #Al effects
            Al_direct := a
            Al_indirect := b*c
            Al_shared := h*n + h*k*c + i*l*c + j*m*c
            Al_total := a + b*c
            '
fit_al <- lavaan::sem(al_model, data = sem_data)
semPaths(object = fit_al, what = "std", layout = "spring", residuals = FALSE, nCharNodes = 10, 
         sizeMan = 9, sizeLat = 7, label.cex = 0.9, edge.label.cex = 0.9)
summary(fit_al, fit.measures= TRUE)



LittleHelpers::lav_to_graph(fit_al, file = "al_sem.svg")

summary(lm(sem_data$Fireline.intensity.kWm.1 ~ sem_data$log_al))


#-------------------------------------------------------------------------------
# Figure 5 differences in functional type
library("multcompView")

#significant -- pseudoreplicated
boxplot(burn_level$Fireline.intensity.kWm.1 ~ burn_level$FG)
summary(aov(burn_level$Fireline.intensity.kWm.1 ~ burn_level$FG))

#significant in mixed-effects model with three types of species
mod <- lmer(Fireline.intensity.kWm.1 ~ FG + (1|Species), data = burn_level[burn_level$success, ])
summary(mod)
anova(mod)
m.emm <- emmeans::emmeans(mod, "FG", )


boxplot(burn_sp$Fireline.intensity.kWm.1 ~ burn_sp$FG)
model <- aov(Fireline.intensity.kWm.1 ~ FG, data = burn_sp)
summary(model)
fire_int_tukey <- TukeyHSD(model, "FG", conf.level=0.95)

boxplot(burn_sp$SLAcm2g.1 ~ burn_sp$FG)
model <- aov(SLAcm2g.1 ~ FG, data = burn_sp)
summary(model)
sla_tukey <- TukeyHSD(model, "FG", conf.level=0.95)

boxplot(burn_sp$dryHull ~ burn_sp$FG)
model <- aov(dryHull ~ FG, data = burn_sp)
summary(model)
hull_tukey <- TukeyHSD(model, conf.level=0.95)

boxplot(burn_sp$S3 ~ burn_sp$FG)
model <- aov(S3 ~ FG, data = burn_sp)
summary(model)
s3_tukey <- TukeyHSD(model, "FG", conf.level=0.95)

boxplot(burn_sp$C3H ~ burn_sp$FG)
model <- aov(C3H ~ FG, data = burn_sp)
summary(model)
c3H_tukey <- TukeyHSD(model, "FG", conf.level=0.95)

boxplot(log(burn_sp$Al.Conc.) ~ burn_sp$FG)
model <- aov(log(Al.Conc.) ~ FG, data = burn_sp)
summary(model)
al_tukey <- TukeyHSD(model, conf.level=0.95)

boxplot(log(litter_spec_vol) ~ FG, data = burn_sp)
model <- aov(log(litter_spec_vol) ~ FG, data = burn_sp)
summary(model)
spec_vol_tukey <- TukeyHSD(model, conf.level=0.95)

generate_label_df <- function(TUKEY, variable){
  # Extract labels and factor levels from Tukey post-hoc
  Tukey.levels <- TUKEY[[variable]][,4]
  Tukey.labels <- data.frame(multcompView::multcompLetters(Tukey.levels, reversed = TRUE)['Letters'])
  
  #I need to put the labels in the same order as in the boxplot :
  Tukey.labels$treatment=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$treatment) , ]
  return(Tukey.labels)
}


#----------------
#Make multipanel figure

theme_update(plot.margin = unit(c(0.5, 0, 0,0.25), "cm"))

# Tukey test to study each pair of treatment :
TUKEY <- fire_int_tukey
labels<-generate_label_df(TUKEY , "FG")#generate labels using function
names(labels)<-c('Letters','FG')#rename columns for merging
yvalue<-aggregate(Fireline.intensity.kWm.1 ~ FG, data = burn_sp, mean)# obtain letter position for y axis using means
final<-merge(labels,yvalue) #merge dataframes

fire_int_fg <- ggplot(data = burn_sp[!is.na(burn_sp$FG), ], 
                      aes(x = FG, y = Fireline.intensity.kWm.1)) +
  geom_boxplot(color = "steelblue") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  ylab(expression("Fireline intensity (kW m"^-1*")")) +
  xlab(expression("Functional group")) + 
  geom_text(data = final, aes(label = Letters),vjust=-5,hjust=-.5)
plot(fire_int_fg)

## panel for specivif volume
TUKEY <- spec_vol_tukey
labels<-generate_label_df(TUKEY , "FG")#generate labels using function
names(labels)<-c('Letters','FG')#rename columns for merging
yvalue<-aggregate(litter_spec_vol ~ FG, data = burn_sp, mean)# obtain letter position for y axis using means
final<-merge(labels,yvalue) #merge dataframes

spec_vol_fg <- ggplot(data = burn_sp[!is.na(burn_sp$FG), ], aes(x = FG, y = litter_spec_vol)) +
  geom_boxplot(color = "steelblue") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  ylab(expression("Litter specific volume (m"^3*"g"^-1*")")) +
  xlab(expression("Functional group")) + 
  geom_text(data = final, aes(label = Letters),vjust=-5,hjust=-.5)
plot(spec_vol_fg)

#panel for SLA
TUKEY <- sla_tukey
labels<-generate_label_df(TUKEY , "FG")#generate labels using function
names(labels)<-c('Letters','FG')#rename columns for merging
yvalue<-aggregate(SLAcm2g.1 ~ FG, data = burn_sp, mean)# obtain letter position for y axis using means
final<-merge(labels,yvalue) #merge dataframes

sla_fg <- ggplot(data = burn_sp[!is.na(burn_sp$FG), ], aes(x = FG, y = SLAcm2g.1)) +
  geom_boxplot(color = "steelblue") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  ylab(expression("Specific leaf area (cm"^2*" g"^-1*")")) +
  xlab(expression("Functional group"))+ 
  geom_text(data = final, aes(label = Letters),vjust=-5,hjust=-.5)
plot(sla_fg)


# Tukey test to study each pair of treatment :
TUKEY <- s3_tukey
labels<-generate_label_df(TUKEY , "FG")#generate labels using function
names(labels)<-c('Letters','FG')#rename columns for merging
yvalue<-aggregate(S3 ~ FG, data = burn_sp, mean)# obtain letter position for y axis using means
final<-merge(labels,yvalue) #merge dataframes


s3_fg <- ggplot(data = burn_sp[!is.na(burn_sp$FG), ], aes(x = FG, y = S3)) +
  geom_boxplot(color = "steelblue") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  ylab(expression("Leaf size (cm)")) +
  xlab(expression("Functional group"))+ 
  geom_text(data = final, aes(label = Letters),vjust=-4.5,hjust=-.5)
plot(s3_fg)


# Tukey test to study each pair of treatment :
TUKEY <- c3H_tukey
labels<-generate_label_df(TUKEY , "FG")#generate labels using function
names(labels)<-c('Letters','FG')#rename columns for merging
yvalue<-aggregate(C3H ~ FG, data = burn_sp, mean)# obtain letter position for y axis using means
final<-merge(labels,yvalue) #merge dataframes

C3H_fg <- ggplot(data = burn_sp[!is.na(burn_sp$FG), ], aes(x = FG, y = C3H)) +
  geom_boxplot(color = "steelblue") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  ylab(expression("Leaf curl (unitless)")) +
  xlab(expression("Functional group"))+ 
  geom_text(data = final, aes(label = Letters),vjust=-4.5,hjust=-.5)
plot(C3H_fg)

# gridExtra::grid.arrange(fire_int_fg, spec_vol_fg, sla_fg, s3_fg, C3H_fg, nrow = 2)

fg_plots <- cowplot::plot_grid(fire_int_fg +
                                 theme(axis.text=element_text(size=7),
                                       axis.title=element_text(size=8)), 
                               spec_vol_fg+
                                 theme(axis.text=element_text(size=7),
                                       axis.title=element_text(size=8)),
                               sla_fg+
                                 theme(axis.text=element_text(size=7),
                                       axis.title=element_text(size=8)), 
                               s3_fg+
                                 theme(axis.text=element_text(size=7),
                                       axis.title=element_text(size=8)), 
                               C3H_fg+
                                 theme(axis.text=element_text(size=7),
                                       axis.title=element_text(size=8)),
                               align='v', vjust=1, hjust = -0.5, scale = 1, 
                               labels = "auto", 
                               nrow = 2, ncol = 3, axis = "l")
plot(fg_plots)
ggsave(plot = fg_plots, 
       filename = "./outputs/fg_plots.png",
       device = "png",
       units = "in",
       scale = 1,
       width = 6,
       height = 4)

#-------------------------------------------------------------------------------
# Figure 1
fig_intensity_specvol <- ggplot(data = burn_level[burn_level$success, ], aes(x = litter_spec_vol, y = Fireline.intensity.kWm.1)) +
  geom_point(color = "steelblue") + 
  geom_smooth(method='lm', color = "black", formula= (y ~ x)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  ylab(expression("Fireline intensity (kW m"^-1*")")) +
  xlab(expression("Litter specific volume (m"^3*" kg"^-1*")"))
ggsave(fig_intensity_specvol,
       device = "png",
       units = "in",
       filename = "Outputs/Figure4-intensity vs specific volume.png",
       height= 3,
       width = 3)

summary(lm(Fireline.intensity.kWm.1 ~ litter_spec_vol, data = burn_level[burn_level$success, ]))

ggplot(data = burn_level[burn_level$success, ], aes(x = litter_bulk_density, y = Fireline.intensity.kWm.1)) +
  geom_point(color = "steelblue") + 
  geom_smooth(method='lm', color = "black", formula= (y ~ x)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  ylab(expression("Fireline intensity (kW m"^-1*")")) +
  xlab(expression("Litter bulk density (kg"^-1*"m"^3*")"))



#Figure S1
#correlations among volumes

library(GGally)
options(scipen = 999)

leaf_data <- leaf_data %>% filter(leaf_data$dryHull > 0.002)

lowerFn <- function(data, mapping, method = "lm", ...) {
  p <- ggplot(data = data, mapping = mapping) +
    geom_point(colour = "steelblue") +
    geom_smooth(method = method, color = "black", ...) +
    scale_x_continuous(trans = "log10", 
                       breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000, 100000),
                       labels = function(x) as.character(x)) +
    scale_y_continuous(trans = "log10", 
                       breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000, 100000),
                       labels = function(x) as.character(x)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2)
  p
}

diagFn <- function(data, mapping, method = "lm", ...) {
  p <- ggplot(data = data, mapping = mapping) +
    geom_density() +
    scale_x_continuous(trans = "log10",
                       breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000, 100000),
                       labels = function(x) as.character(x)) +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          )
  p
}

ggpairs(data    = leaf_data, 
        lower   = list(continuous = wrap(lowerFn, method = "lm")),
        diag    = list(continuous = wrap(diagFn)),
        columns = c("freshHull", "particle_volume", "rectangular_volume"),
        columnLabels = c("Convex hull", "Particle volume", "Rectangular volume"))

