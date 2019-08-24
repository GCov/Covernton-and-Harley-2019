
snail1 <- read.csv("sal.snail.csv", header = TRUE)

snail1 <- subset(snail1, Tank != 'NA')

head(snail1)#view data

library(survival)
library(coxme)
library(ggplot2)
library(ggfortify)
library(survminer)
library(cowplot)
library(MuMIn)

npalette <- c('#D95A43', '#246378', '#FFC524', '#8C7C79')

## Fix data

snail1$death <- as.numeric(snail1$X)
snail1$status <- as.numeric(snail1$X.1)
snail1$length <- as.numeric(snail1$Snail.Length)
snail1$treat <- as.character(snail1$Treatment)
snail1$treat <- as.factor(snail1$treat)
snail1$tank <- as.factor(snail1$Tank)
snail1$tub <- as.character(snail1$Tub)
snail1$tub <- as.factor(snail1$tub)

names(snail1)

snail1 <- snail1[c('treat', 'tank', 'tub', 'length', 'death', 'status')]
head(snail1)
snail1$pop <- factor(snail1$tub, levels = c('WRJ', 'WRA', 'TS'))

snail1$death <- snail1$death-1 ## Make 1st day t=0
head(snail1)
summary(snail1)

## Plot

fit1 <- survfit(Surv(death, status) ~ 1, 
                data = subset(snail1, tub == 'WRJ' & treat == 'High'))
fit2 <- survfit(Surv(death, status) ~ 1, 
                data = subset(snail1, tub == 'WRA' & treat == 'High'))
fit3 <- survfit(Surv(death, status) ~ 1, 
                data = subset(snail1, tub == 'TS' & treat == 'High'))

fit4 <- survfit(Surv(death, status) ~ 1, 
                data =  subset(snail1, tub == 'WRJ' & treat == 'Low'))
fit5 <- survfit(Surv(death, status) ~ 1, 
                data = subset(snail1, tub == 'WRA' & treat == 'Low'))
fit6 <- survfit(Surv(death, status) ~ 1, 
                data = subset(snail1, tub == 'TS' & treat == 'Low'))

fit7 <- survfit(Surv(death, status) ~ 1, 
                data = subset(snail1, tub == 'WRJ' & treat == 'Rescue'))
fit8 <- survfit(Surv(death, status) ~ 1, 
                data = subset(snail1, tub == 'WRA' & treat == 'Rescue'))
fit9 <- survfit(Surv(death, status) ~ 1, 
                data =  subset(snail1, tub == 'TS' & treat == 'Rescue'))

A <- ggsurvplot(fit1, color = 'black', linetype = 'solid', 
                size = 0.5, censor.size = 4,
                conf.int = FALSE, conf.int.fill = npalette[1],
                xlim = c(0,22))$plot +
  xlab('') +
  ylab('') +
  scale_x_continuous(limits = c(0,22), breaks = c(0, 7, 14, 21)) +
  theme_bw() +
  theme(
    text = element_text(size = 10),
    axis.text = element_text(size = 6),
    strip.background = element_blank(),
    strip.text =  element_text(size = 8),
    panel.grid = element_blank(),
    plot.margin = margin(c(0.5,0,0,0))
  )

B <- ggsurvplot(fit2, color = 'black', linetype = 'solid', 
                size = 0.5, censor.size = 4,
                conf.int = TRUE, conf.int.fill = npalette[1],
                xlim = c(0,22))$plot +
  xlab('') +
  ylab('') +
  theme_bw() +
  scale_x_continuous(limits = c(0,22), breaks = c(0, 7, 14, 21)) +
  theme(
    text = element_text(size = 10),
    axis.text = element_text(size = 6),
    strip.background = element_blank(),
    strip.text =  element_text(size = 8),
    panel.grid = element_blank(),
    plot.margin = margin(c(0.5,0,0,0))
  )

C <- ggsurvplot(fit3, color = 'black', linetype = 'solid', 
                size = 0.5, censor.size = 4,
                conf.int = FALSE, conf.int.fill = npalette[1],
                xlim = c(0,22))$plot +
  xlab('') +
  ylab('') +
  theme_bw() +
  scale_x_continuous(limits = c(0,22), breaks = c(0, 7, 14, 21)) +
  theme(
    text = element_text(size = 10),
    axis.text = element_text(size = 6),
    strip.background = element_blank(),
    strip.text =  element_text(size = 8),
    panel.grid = element_blank(),
    plot.margin = margin(c(0.5,0,0,0))
    )

D <- ggsurvplot(fit4, color = 'black', linetype = 'solid', 
                size = 0.5, censor.size = 4,
                conf.int = TRUE, conf.int.fill = npalette[1],
                xlim = c(0, 22))$plot +
  xlab('') +
  ylab('Proportion Alive') +
  theme_bw() +
  scale_x_continuous(limits = c(0,22), breaks = c(0, 7, 14, 21)) +
  theme(
    text = element_text(size = 8),
    axis.text = element_text(size = 6),
    strip.background = element_blank(),
    strip.text =  element_text(size = 7),
    panel.grid = element_blank(),
    plot.margin = margin(c(0,0,0,0))
    )

E <- ggsurvplot(fit5, color = 'black', linetype = 'solid', 
                size = 0.5, censor.size = 4,
                conf.int = TRUE, conf.int.fill = npalette[1],
                xlim = c(0,22))$plot +
  xlab('') +
  ylab('') +
  theme_bw() +
  scale_x_continuous(limits = c(0,22), breaks = c(0, 7, 14, 21)) +
  theme(
    text = element_text(size = 10),
    axis.text = element_text(size = 6),
    strip.background = element_blank(),
    strip.text =  element_text(size = 8),
    panel.grid = element_blank(),
    plot.margin = margin(c(0,0,0,0))
    )

G <- ggsurvplot(fit6, color = 'black', linetype = 'solid', 
                size = 0.5, censor.size = 4,
                conf.int = TRUE, conf.int.fill = npalette[1],
                xlim = c(0,22))$plot +
  xlab('') +
  ylab('') +
  theme_bw() +
  scale_x_continuous(limits = c(0,22), breaks = c(0, 7, 14, 21)) +
  theme(
    text = element_text(size = 10),
    axis.text = element_text(size = 6),
    strip.background = element_blank(),
    strip.text =  element_text(size = 8),
    panel.grid = element_blank(),
    plot.margin = margin(c(0,0,0,0))
    )

H <- ggsurvplot(fit7, color = 'black', linetype = 'solid', 
                size = 0.5, censor.size = 4,
                conf.int = TRUE, conf.int.fill = npalette[1],
                xlim = c(0,22))$plot +
  xlab('') +
  ylab('') +
  theme_bw() +
  scale_x_continuous(limits = c(0,22), breaks = c(0, 7, 14, 21)) +
  theme(
    text = element_text(size = 10),
    axis.text = element_text(size = 6),
    strip.background = element_blank(),
    strip.text =  element_text(size = 8),
    panel.grid = element_blank(),
    plot.margin = margin(c(0,0,0,0))
    )

I <- ggsurvplot(fit8, color = 'black', linetype = 'solid', 
                size = 0.5, censor.size = 4,
                conf.int = TRUE, conf.int.fill = npalette[1],
                xlim = c(0,22))$plot +
  xlab('Day') +
  ylab('') +
  theme_bw() +
  scale_x_continuous(limits = c(0,22), breaks = c(0, 7, 14, 21)) +
  theme(
    text = element_text(size = 8),
    axis.text = element_text(size = 6),
    strip.background = element_blank(),
    strip.text =  element_text(size = 7),
    panel.grid = element_blank(),
    plot.margin = margin(c(0,0,0,0))
    )

J <- ggsurvplot(fit9, color = 'black', linetype = 'solid', 
                size = 0.5, censor.size = 4,
                conf.int = TRUE, conf.int.fill = npalette[1],
                xlim = c(0,22))$plot +
  xlab('') +
  ylab('') +
  theme_bw() +
  scale_x_continuous(limits = c(0,22), breaks = c(0, 7, 14, 21)) +
  theme(
    text = element_text(size = 10),
    axis.text = element_text(size = 6),
    strip.background = element_blank(),
    strip.text =  element_text(size = 8),
    panel.grid = element_blank(),
    plot.margin = margin(c(0,0,0,0))
    )

xoff <- 0.1 # relative x position of label, within one plot
yoff <- 1 # relative y position of label, within one plot
labels <- c('WRJ High', 'WRA High', 'TS High',
            'WRJ Low', 'WRA Low', 'TS Low',
            'WRJ Variable', 'WRA Variable', 'TS Variable')

surv_plot <-
  plot_grid(
    A,
    B,
    C,
    D,
    E,
    G,
    H,
    I,
    J,
      align = 'hv',
    ncol = 3,
    nrow = 3,
    labels = labels,
    label_size = 7,
    hjust = 0,
    label_x = 0.22,
    label_y = 0.52,
    label_fontface = 'plain'
  ) + 
  theme(plot.margin = margin(c(0.5,0,0,0)))

surv_plot

## Stats

## First reduce data set to be only rescue and low treatments

snail2 <- subset(snail1, treat != 'High')
snail2$treat <- as.character(snail2$treat)
snail2$treat <- as.factor(snail2$treat)

## Coxme model with interactions between population, treatment, and length
coxmod1 <- coxme(Surv(death, status) ~ pop*treat*length + (1|tank), 
                 data = snail2)
summary(coxmod1)

## Without length interactions
coxmod2 <- coxme(Surv(death, status) ~ pop*treat + length + (1|tank), 
                 data = snail2)

## Without population interactions
coxmod3 <- coxme(Surv(death, status) ~ pop + treat*length + (1|tank), 
                 data = snail2)

## Without treatment interactions
coxmod4 <- coxme(Surv(death, status) ~ pop*length + treat + (1|tank), 
                 data = snail2)

## Without any interactions 
coxmod5 <- coxme(Surv(death, status) ~ pop + length + treat + (1|tank), 
                 data = snail2)

AICc(coxmod1, coxmod2, coxmod3, coxmod4, coxmod5)
anova(coxmod1, coxmod2, coxmod3, coxmod4, coxmod5)

## model with no interactions (coxmod 5) seems to be best according to AICc

summary(coxmod5)

## Test the overall significance of population and treatment

coxmod6 <- coxme(Surv(death, status) ~ length + treat + (1|tank), 
                 data = snail2)
coxmod7 <- coxme(Surv(death, status) ~ pop + length + (1|tank), 
                 data = snail2)
coxmod8 <- coxme(Surv(death, status) ~ pop + treat + (1|tank), 
                 data = snail2)

anova(coxmod5, coxmod6)  # population significant (p = 0.009)
anova(coxmod5, coxmod7)  # treatment significant (p < 0.001)
anova(coxmod5, coxmod8)  # length significant (p < 0.001)

## Model validation

coxmod9 <- coxph(Surv(death, status) ~ pop + length + treat + strata(tank),
                 data = snail2)
fit.test <- cox.zph(coxmod9)
fit.test
plot(fit.test, var = 1)
plot(fit.test, var = 2)
plot(fit.test, var = 3)
plot(fit.test, var = 4)

## Plot model predictions 

hazard_plot <- ggplot(snail2) +
  geom_line(aes(
    x = length,
    y = predict(coxmod5, type = 'risk'),
    linetype = treat
  )) +
  facet_wrap( ~ pop) +
  scale_linetype_discrete(c('solid', 'dashed')) +
  labs(x = 'Shell Length (mm)', y = 'Hazard Ratio') +
  theme_bw() +
  theme(
    text = element_text(size = 8),
    axis.text = element_text(size = 7),
    strip.background = element_blank(),
    strip.text = element_text(size = 8),
    panel.grid = element_blank(),
    legend.title = element_blank(),
    legend.margin = margin(c(0,0,0,0))
  )

png(
  filename = "Cox Model Predictions.png",
  width = 15.24,
  height = 10,
  units = "cm",
  pointsize = 10,
  res = 600
)

plot_grid(surv_plot, hazard_plot, 
          ncol = 1,
          labels = c('A', 'B'),
          label_size = 8,
          rel_heights = c(2,1))

dev.off()

## Try GLM to determine the effect of length vs. population

## First subset just for low treatment

lowtreat <- subset(snail2, treat == 'Low')

hist(lowtreat$death)
hist(log(lowtreat$death))

## separate only by population

library(plyr)
lowtreat$pop <- mapvalues(lowtreat$pop, 
                          from = c('WRA','WRJ'),
                          to = c('WR', 'WR'))

summary(lowtreat)

## GLM time
library(glmmTMB)
library(MuMIn)
library(nlme)

## Scale and centre length

lowtreat$cen.length <- scale(lowtreat$length, center = TRUE, scale = TRUE)

popmod1 <- glmer(death ~ pop*cen.length + (1 | tank), 
                   family = poisson(link = 'log'),
                   weights = length,
                   data = lowtreat)
summary(popmod1)

plot(resid(popmod1, type = 'pearson') ~ 
       predict(popmod1, type = 'link'))
abline(0,0)  ## Odd pattern for those that survive long

plot(resid(popmod1, type = 'pearson') ~ 
       lowtreat$cen.length)
abline(0,0)  ## lower variance for larger snails

plot(resid(popmod1, type = 'pearson') ~ 
       lowtreat$pop)
abline(0,0)  ## lower variance for TS snails

## Try NB distribution

popmod2 <- glmmTMB(death ~ pop*cen.length + (1 | tank), 
                   family = nbinom1(link = 'logit'),
                   data = lowtreat)
summary(popmod2)

## all p values are 1, something messed up going on

popmod3 <- lmer(death ~ pop*cen.length + (1 | tank), 
                weights = length,
                data = lowtreat)
summary(popmod3)

plot(resid(popmod3, type = 'pearson') ~ 
       predict(popmod3, type = 'link'))
abline(0,0)  ## Odd pattern for those that survive long

plot(resid(popmod3, type = 'pearson') ~ 
       lowtreat$cen.length)
abline(0,0)  ## lower variance for larger snails

plot(resid(popmod3, type = 'pearson') ~ 
       lowtreat$pop)
abline(0,0)  ## lower variance for TS snails

AICc(popmod1, popmod3)  # regular ANCOVA is better fit

## Try log-transforming day of death

popmod4 <- lmer(log(death) ~ pop*cen.length + (1 | tank),
                weights = length,
                data = lowtreat)
summary(popmod4)

plot(resid(popmod4) ~ 
       predict(popmod4))
abline(0,0)  ## Odd pattern for those that survive long

plot(resid(popmod4) ~ 
       lowtreat$cen.length)
abline(0,0)  ## lower variance for larger snails

plot(resid(popmod4) ~ 
       lowtreat$pop)
abline(0,0)  # lower variance for TS snails

AICc(popmod1, popmod3, popmod4)  # now popmod4 is the best fit

plot(log(death) ~ cen.length, data = lowtreat, ylim = c(1,3))
points(predict(popmod1) ~ lowtreat$cen.length, col = 'red')

## Try with variance structure

popmod5 <- lme(death ~ cen.length*pop, 
               random = (~1 | tank), 
               weights = varIdent(form=~1 | pop),
               data = lowtreat,
               method = 'ML')
summary(popmod5)

plot(resid(popmod5, type = 'pearson') ~ 
       predict(popmod5, type = 'link'))
abline(0,0) 

plot(resid(popmod5, type = 'pearson') ~ 
       lowtreat$cen.length)
abline(0,0)

plot(resid(popmod5, type = 'pearson') ~ 
       lowtreat$pop)
abline(0,0)  # best looking in terms of assumptions

AICc(popmod1, popmod3, popmod4, popmod5)  # 2nd lowest AICc for popmod5

plot(death ~ cen.length, data = lowtreat)
points(predict(popmod5) ~ lowtreat$cen.length, col = 'red')

summary(popmod5)
test1<- lme(death ~ cen.length + pop, 
            random = (~1 | tank), 
            weights = varIdent(form=~1 | pop),
            data = lowtreat,
            method = 'ML')
test2 <- lme(death ~ cen.length, 
             random = (~1 | tank), 
             weights = varIdent(form=~1 | pop),
             data = lowtreat,
             method = 'ML')
test3 <- lme(death ~ pop, 
             random = (~1 | tank), 
             weights = varIdent(form=~1 | pop),
             data = lowtreat,
             method = 'ML')

anova(popmod5, test1)  # interaction sig at p=0.01
anova(test1, test2)  # population sig at p<0.01
anova(test1, test3)  # length sig at p<0.01
