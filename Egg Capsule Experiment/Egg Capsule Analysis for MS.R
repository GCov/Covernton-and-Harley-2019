library(plyr)
library(survival)
library(coxme)
library(ggplot2)
library(ggfortify)
library(coxme)
library(MuMIn)
library(cowplot)

## Load dataset

eggcaps <-
  read.csv("eggcaps.csv", header = TRUE)
head(eggcaps)

npalette <- c('#D95A43', '#246378', '#FFC524', '#8C7C79')

## Explore data

summary(eggcaps)
eggcaps$block <- as.factor(eggcaps$block)
eggcaps$sal <- as.factor(eggcaps$sal)
eggcaps$stat <- as.numeric(eggcaps$stat)


hist(eggcaps$purpdays)
plot(purpdays ~ sal, data = eggcaps)

eggcaps$sal <- mapvalues(
  eggcaps$sal,
  from = c('9', '12', '15'),
  to = c('9 psu', '12 psu', '15 psu')
)

with(eggcaps, tapply(purpdays, sal, mean))
with(eggcaps, tapply(purpdays, sal, sd))

## Survival analysis

## Create dat survival object

eggfit <- survfit(Surv(purpdays, stat) ~ sal, data = eggcaps)


survplot <-
  autoplot(
    eggfit,
    surv.geom = 'line',
    surv.colour = 'black',
    surv.size = 0.3,
    conf.int = TRUE,
    conf.int.fill = npalette[1],
    censor = TRUE,
    censor.shape = '+',
    censor.size = 0.5,
    facets = TRUE,
    ncol = 3,
    xlab = 'Day',
    ylab = 'Proportion of Non-Purple \nEgg Capsules',
    ylim = c(0, 1)
  ) +
  theme_bw() +
  scale_x_continuous(limits = c(0, 56),
                     breaks = seq(from = 0, to = 56, by = 14)) +
  theme(
    text = element_text(size = 8),
    axis.text = element_text(size = 6),
    strip.background = element_blank(),
    strip.text =  element_text(size =),
    panel.grid = element_blank(),
    plot.margin = margin(c(0,0.2,0,0), unit = 'cm')
  )


## Mixed effects coxph model withblock as random effect

fullmodel <- coxme(eggsurv ~ sal + (1 | block / beaker), data = eggcaps)
summary(fullmodel)

## Block has low variance, try removing

mem2 <- coxme(eggsurv ~ sal + (1 | beaker), data = eggcaps)
summary(mem2)
anova(fullmodel, mem2) 
AICc(fullmodel, mem2)
eggcaps$sal <- relevel(eggcaps$sal, '12')
eggsurv <- Surv(eggcaps$purpdays, eggcaps$stat)

## Try with no random effects

m1 <- coxph(eggsurv ~ sal, data = eggcaps)
summary(m1)
anova(mem2, m1)
AIC(mem2, m1) ## CANNOT take out the effect of beaker, mem2 is best model
cox.zph(m1)## looks like assumption of proportional hazard not violated in simplified model
plot(cox.zph(m1))
sfit1 <- survfit(m1)
plot(sfit1)

fixed.effects(mem2)
random.effects(mem2)

## Use full model

summary(fullmodel)
nullmodel <- coxme(eggsurv ~ (1 | block / beaker), data = eggcaps) 
anova(fullmodel, nullmodel)

## Look at original survival object

eggfit <- survfit(eggsurv ~ sal, data = eggcaps)
summary(eggfit)
m1fitted <- predict(m1, type = "risk", se.fit = TRUE)
m1fitted <- as.data.frame(m1fitted)
m1fit <- m1fitted$fit
plot(m1fit ~ eggcaps$purpdays)
plot(resid(m1) ~ m1fit)


## Plot fullmodel predictions

predictplot <-
  ggplot(eggcaps) +
  geom_boxplot(
    aes(x = sal,
        y = predict(fullmodel, type = 'risk')),
    fill = npalette[1],
    alpha = 0.5,
    size = 0.5,
    outlier.size = 0.5
  ) +
  labs(x = 'Treatment', y = 'Hazard Ratio') +
  theme_bw() +
  theme(
    text = element_text(size = 8),
    axis.text = element_text(size = 7),
    strip.background = element_blank(),
    strip.text = element_text(size = 8),
    panel.grid = element_blank(),
    legend.title = element_blank(),
    legend.margin = margin(c(0, 0, 0, 0))
  )

png(
  filename = "Egg Capsule Stress.png",
  width = 7.62,
  height = 7.62,
  units = "cm",
  pointsize = 10,
  res = 600
)

plot_grid(survplot, predictplot, labels = c('A', 'B'),
          label_size = 8, ncol = 1, rel_heights = c(1.5,1))

dev.off()
