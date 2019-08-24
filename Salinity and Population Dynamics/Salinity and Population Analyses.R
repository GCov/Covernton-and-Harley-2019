library(plyr)
library(nlme)
library(lme4)
library(MuMIn)
library(ggplot2)
library(zoo)
library(cowplot)
library(glmmTMB)


npalette <- c('#D95A43', '#246378', '#FFC524', '#8C7C79')

## Salinity analysis

salinity <-
  read.csv("salinity.csv", header = TRUE)

summary(salinity)
head(salinity$Date)

salinity$Site <- as.factor(salinity$Site)
salinity$sal_reg <- as.factor(salinity$salinity.regime)
salinity$sal.dat <- as.Date(salinity$Date)
head(salinity$sal.dat)

salinity$Site <-
  mapvalues(
    salinity$Site,
    from = c(levels(salinity$Site)),
    to = c(
      'Dunbar',
      'Figurehead Point',
      'Kitsilano Point',
      'Point Atikinson',
      'Navvy Jack Point',
      'Tower Beach North',
      'Tower Beach South',
      'Waterloo',
      'Weston Park',
      'White Rock',
      'Yellow Point'
    )
  )

salinity$Site <-
  factor(
    salinity$Site,
    levels = c(
      'Tower Beach South',
      'Tower Beach North',
      'Dunbar',
      'Waterloo',
      'Point Atikinson',
      'Kitsilano Point',
      'Weston Park',
      'Navvy Jack Point',
      'Figurehead Point',
      'White Rock',
      'Yellow Point'
    )
  )

## What is the best predictor of site salinity?

## Figure out what measure of FR outflow best predicts salinities at our sites

## First pull out different FR outflow measures and pair them with salinity data
FRoutflow <- 
  read.csv("FRoutflow.csv", header = TRUE)

summary(salinity)
summary(FRoutflow)
head(FRoutflow)

FRoutflow$Date <- as.Date(FRoutflow$Date)

FRoutflow2 <-
  subset(FRoutflow, PARAM == 1)  # only want outflow data
salinity$Date <- salinity$sal.dat

## Now calculate total outflow over preceding n days

salinity_grab <- list()

fourdaymean <- rollmeanr(FRoutflow2$Value, k = 4, align = 'right')
fourdayframe <- data.frame(FRoutflow2[4:15545, ], fourdaymean)

twoweekmean <- rollmeanr(FRoutflow2$Value, k = 14, align = 'right')
twoweekframe <- data.frame(FRoutflow2[14:15545, ], twoweekmean)

for (i in 1:length(salinity$Date)) {
  dateminus1 <-
    FRoutflow2$Value[FRoutflow2$Date == (salinity$Date[i] - 1)]
  dateminus4 <-
    FRoutflow2$Value[FRoutflow2$Date == (salinity$Date[i] - 4)]
  fourdayave <- fourdayframe$fourdaymean[fourdayframe$Date ==
                                           salinity$Date[i]]
  twoweekave <- twoweekframe$twoweekmean[twoweekframe$Date ==
                                           salinity$Date[i]]
  data <-
    data.frame(salinity[i, ], dateminus1, dateminus4, fourdayave,
               twoweekave)
  salinity_grab[[i]] <- data
}

salinity_lag <- data.frame(rbind.fill(salinity_grab))

head(salinity_lag)

plot(salinity ~ log(dateminus1), data = salinity_lag)
plot(salinity ~ log(fourdayave), data = salinity_lag)
plot(salinity ~ log(twoweekave), data = salinity_lag)

## Model salinity as a function of Fraser River outflow

fit0 <- gls(
  salinity ~ log(dateminus1) * Site,
  correlation = corAR1(form = ~ Date | Site),
  data = salinity_lag
)

fit0.1 <- gls(
  salinity ~ log(fourdayave) * Site,
  correlation = corAR1(form = ~ Date | Site),
  data = salinity_lag
)

fit0.2 <- gls(
  salinity ~ log(twoweekave) * Site,
  correlation = corAR1(form = ~ Date | Site),
  data = salinity_lag
)

fit0.3 <- gls(
  salinity ~ log(dateminus4) * Site,
  correlation = corAR1(form = ~ Date | Site),
  data = salinity_lag
)

fit0.4 <- gls(
  salinity ~ log(fourdayave) * Site,
  correlation = corCAR1(form = ~ Date | Site),
  data = salinity_lag
)

AICc(fit0, fit0.1, fit0.2, fit0.3, fit0.4)

## corCAR1 variance structure and log(twoweekave) best fit
## However oceanographic models say that 4 days best predicts salinity
## in Burrard Inlet, so we'll use 'dateminus4' (fit0.3)

coef(fit0.3)
summary(fit0.3)

plot(salinity ~ dateminus4, data = salinity_lag)
points(fitted(fit0.3) ~ salinity_lag$dateminus4, col = 'red')

plot(resid(fit0.3) ~ fitted(fit0.3))
abline(0, 0)
plot(resid(fit0.3) ~ salinity_lag$Site)

## Get predicted values and se and change order of sites

salinity_lag$predict <- predict(fit0.3)

salinity_lag$Site <-
  factor(
    salinity$Site,
    levels = c(
      'Tower Beach South',
      'Tower Beach North',
      'Dunbar',
      'Waterloo',
      'Point Atikinson',
      'Kitsilano Point',
      'Weston Park',
      'Navvy Jack Point',
      'Figurehead Point',
      'White Rock',
      'Yellow Point'
    )
  )

## Plot salinity by site against FR outflow 4 days prior

Plot1 <-
  ggplot(
    subset(
      salinity_lag,
      Site != 'Yellow Point' &
        Site != 'Navvy Jack Point' &
        Site != 'Figurehead Point' &
        Site != 'White Rock'
    ),
    aes(x = dateminus4, y = salinity)
  ) +
  geom_point(size = 0.3) +
  geom_line(aes(y = predict), linetype = 'dashed', size = 0.5) +
  geom_abline(slope = 0,
              intercept = 15,
              linetype = 'dotted') +
  facet_wrap(~ Site, scales = "fixed", ncol = 2) +
  xlab('Fraser River Discharge (-4 Day Lag)') +
  ylab('Salinity') +
  coord_cartesian(ylim = c(0, 30), xlim = c(0, 11500)) +
  ggtitle('Low Salinity Sites') +
  theme_bw() +
  theme(
    legend.position = 'none',
    plot.title = element_text(face = 'bold', hjust = 0.5, size = 10),
    text = element_text(size = 8),
    axis.text = element_text(size = 6),
    strip.background = element_blank(),
    strip.text = element_text(size = 8),
    panel.grid = element_blank(),
    plot.margin = margin(c(0, 0.5, 0.2, 0.2), unit = 'cm')
  )

Plot2 <-
  ggplot(
    subset(
      salinity_lag,
      Site == 'Yellow Point' |
        Site == 'Navvy Jack Point' |
        Site == 'Figurehead Point' |
        Site == 'White Rock'
    ),
    aes(x = dateminus4, y = salinity)
  ) +
  geom_point(size = 0.3) +
  geom_line(aes(y = predict), linetype = 'dashed', size = 0.5) +
  geom_abline(slope = 0,
              intercept = 15,
              linetype = 'dotted') +
  facet_wrap(~ Site, scales = "fixed", ncol = 1) +
  xlab('Fraser River Discharge (-4 Day Lag)') +
  ylab('Salinity') +
  coord_cartesian(ylim = c(0, 30), xlim = c(0, 11500)) +
  ggtitle('High Salinity Sites') +
  theme_bw() +
  theme(
    legend.position = 'none',
    plot.title = element_text(face = 'bold', hjust = 0.5, size = 10),
    text = element_text(size = 8),
    axis.title.y = element_blank(),
    axis.text = element_text(size = 6),
    strip.background = element_blank(),
    strip.text = element_text(size = 8),
    panel.grid = element_blank(),
    plot.margin = margin(c(0, 0.5, 0.2, 0.5), unit = 'cm'),
    axis.title.x = element_blank()
  )

salplot <- plot_grid(Plot1,
                     Plot2,
                     align = 'h',
                     rel_widths = c(1.8, 1))

FRplot <-
  ggplot(
    subset(
      FRoutflow2,
      Date >= '2010-01-01' & Date < '2017-01-01'
    ),
    aes(x = Date, y = Value)
  ) +
  geom_line(size = 0.2) +
  xlab('') +
  ylab(expression(atop(NA, 
                       atop('Fraser River', paste('Discharge (',m^3~s^-1*')'))))) +
  scale_x_date(date_breaks = '3 months', 
               date_labels = '%b %Y',
               expand = c(0,0)) +
  theme_bw() +
  theme(
    legend.position = 'none',
    plot.title = element_text(face = 'bold', hjust = 0.5, size = 10),
    text = element_text(size = 10),
    axis.text = element_text(size = 6),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(),
    strip.text = element_text(size = 8),
    panel.grid = element_blank(),
    plot.margin = margin(c(0.2, 0.2, 0, 0), unit = 'cm')
  )

png(
  filename = "Salinity Measurements By Site.png",
  width = 15.24,
  height = 12,
  units = "cm",
  pointsize = 10,
  res = 600
)

plot_grid(FRplot,
          salplot,
          align = 'h',
          axis = 'left',
          nrow = 2,
          rel_heights = c(1, 2.5))

dev.off()

## Plot Fraser outflow

png(
  filename = "Fraser Outflow.png",
  width = 11.43,
  height = 9,
  units = "cm",
  pointsize = 10,
  res = 600
)



dev.off()

## Now feed the model historical Fraser River outflow data

head(FRoutflow)

FRoutflow$Date <- as.Date(FRoutflow$Date)

FRoutflow86to16 <- subset(FRoutflow,
                          PARAM == 1 &
                            Date >= "1985-12-28" &
                            Date < "2019-01-01")
summary(FRoutflow86to16)

list1 <- list()

for (i in 5:length(FRoutflow86to16$Date)) {
  fourdayflow <- FRoutflow86to16$Value[i - 4]
  data <- data.frame(FRoutflow86to16[i, ], fourdayflow)
  list1[[i]] <- data
}

temp_data <- data.frame(rbind.fill(list1))

temp_data <- subset(temp_data, Date >= '1986-01-01')

head(temp_data)

sim_data <- list()

for (i in 1:length(levels(salinity_lag$Site))) {
  Site <- rep(levels(salinity_lag$Site)[i],
              length(temp_data$fourdayflow))
  sal.dat <- temp_data$Date
  dateminus4 <- temp_data$fourdayflow
  dat <- data.frame(Site, sal.dat, dateminus4)
  sim_data[[i]] <- dat
}

sim_dataframe <- data.frame(rbind.fill(sim_data))
summary(sim_dataframe)

sim_dataframe$sal.pred <- predict(fit0.3, newdata = sim_dataframe)
head(sim_dataframe)
summary(sim_dataframe)

sim_dataframe$Year <- format(as.Date(sim_dataframe$sal.dat), '%Y')

## Plot salinity predictions

Plot3 <-
  ggplot(
    subset(
      subset(
        sim_dataframe,
        Site != 'Yellow Point' &
          Site != 'Navvy Jack Point' &
          Site != 'Figurehead Point' &
          Site != 'White Rock'
      ),
      sal.dat >= '2010-01-01' & sal.dat <= '2015-01-01'
    ),
    aes(x = sal.dat, y = sal.pred)
  ) +
  geom_line() +
  geom_abline(slope = 0,
              intercept = 15,
              linetype = 'dotted') +
  facet_wrap( ~ Site, scales = "fixed", ncol = 2) +
  xlab('Date') +
  ylab('Salinity') +
  coord_cartesian(ylim = c(0, 32)) +
  scale_x_date(
    expand = c(0, 0),
    date_breaks = '6 months',
    date_labels = c('%b-%Y')
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  ggtitle('Low Salinity Sites') +
  theme(
    legend.position = 'none',
    plot.title = element_text(face = 'bold', hjust = 0.5, size = 10),
    text = element_text(size = 8),
    axis.text = element_text(size = 6),
    strip.background = element_blank(),
    strip.text = element_text(size = 8),
    panel.grid = element_blank(),
    plot.margin = margin(c(0.2, 0.7, 0.2, 0.2), unit = 'cm'),
    panel.spacing = unit(1, "lines"),
    axis.text.x = element_text(
      angle = -90,
      vjust = 0,
      hjust = 0
    )
  )

Plot4 <-
  ggplot(
    subset(
      subset(
        sim_dataframe,
        Site == 'Yellow Point' |
          Site == 'Navvy Jack Point' |
          Site == 'Figurehead Point' |
          Site == 'White Rock'
      ),
      sal.dat >= '2010-01-01' & sal.dat <= '2015-01-01'
    ),
    aes(x = sal.dat, y = sal.pred)
  ) +
  geom_line() +
  geom_abline(slope = 0,
              intercept = 15,
              linetype = 'dotted') +
  facet_wrap( ~ Site, scales = "fixed", ncol = 1) +
  xlab('Date') +
  ylab('Salinity') +
  coord_cartesian(ylim = c(0, 32)) +
  scale_x_date(
    expand = c(0, 0),
    date_breaks = '6 months',
    date_labels = c('%b-%Y')
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  ggtitle('High Salinity Sites') +
  theme(
    legend.position = 'none',
    plot.title = element_text(face = 'bold', hjust = 0.5, size = 10),
    text = element_text(size = 8),
    axis.text = element_text(size = 6),
    strip.background = element_blank(),
    strip.text = element_text(size = 8),
    panel.grid = element_blank(),
    plot.margin = margin(c(0.2, 0.7, 0.2, 0.2), unit = 'cm'),
    panel.spacing = unit(1, "lines"),
    axis.text.x = element_text(
      angle = -90,
      vjust = 0,
      hjust = 0
    ),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )


png(
  filename = "Salinity Predictions By Site.png",
  width = 15.24,
  height = 9,
  units = "cm",
  pointsize = 10,
  res = 600
)

plot_grid(Plot3,
          Plot4,
          align = 'h',
          rel_widths = c(1.9, 1))

dev.off()


## Create a measure of 'salinity stress'

library(plyr)

sal_summary <- ddply(
  sim_dataframe,
  c('Site', 'Year'),
  summarize,
  under20 = length(sal.pred[sal.pred <= 20]),
  under15 = length(sal.pred[sal.pred <= 15]),
  under9 = length(sal.pred[sal.pred <= 9]),
  min.sal = min(sal.pred),
  mean.sal = mean(sal.pred),
  min.date = format((sal.dat[sal.pred == min(sal.pred)])[1], '%j'),
  under15.date = format((sal.dat[sal.pred <= 15])[1], '%j'),
  stress = log(sum(1 / sal.pred)),
  variance = sd(sal.pred) ^ 2
)

head(sal_summary)
summary(sal_summary)

sal_summary$Year <- as.numeric(as.character(sal_summary$Year))
sal_summary$min.date <-
  as.numeric(as.character(sal_summary$min.date))
sal_summary$under15.date <-
  as.numeric(as.character(sal_summary$under15.date))

## Subset for the years of interest
sal_summary2 <- subset(sal_summary,  Year >= 2010 & Year < 2015)

## scale and center stress index

sal_summary2$stressnew <- scale(sal_summary2$stress, 
                                scale = TRUE, center = TRUE)
summary(sal_summary2$stressnew)

png(
  filename = "Low Salinity Days By Site.png",
  width = 16.9,
  height = 10,
  units = "cm",
  pointsize = 10,
  res = 600
)

ggplot(
  subset(
    sal_summary2,
    Site != 'Figurehead Point' &
      Site != 'Navvy Jack Point' &
      Site != 'White Rock' & Site != 'Yellow Point'
  ),
  aes(x = Year, y = under20)
) +
  geom_line(aes(x = Year, y = under15, colour = "15 psu or less"), size = 0.5) +
  geom_line(aes(x = Year, y = under9, colour = "9 psu or less"), size = 0.5) +
  facet_wrap(~ Site, scales = "fixed") +
  xlab('Year') +
  ylab('Number of Days') +
  theme_bw() +
  guides(colour = guide_legend(override.aes = list(size = 6))) +
  labs(size = "n", colour = "Salinities of:") +
  scale_colour_brewer(type = 'qual') +
  theme(
    text = element_text(size = 10),
    panel.spacing = unit(2, "lines"),
    axis.text = element_text(size = 6),
    strip.background = element_blank(),
    strip.text.x = element_text(size = 10),
    panel.grid = element_blank(),
    legend.position = "top",
    legend.justification = "centre",
    legend.margin = margin(0, 0, 0, 0)
  )

dev.off()

## Plot stress for the same years

Plot5 <-
  ggplot(
    subset(
      sal_summary2,
      Site != 'Yellow Point' &
        Site != 'Navvy Jack Point' &
        Site != 'Figurehead Point' &
        Site != 'White Rock'
    ),
    aes(x = Year, y = stressnew)
  ) +
  geom_abline(
    slope = 0,
    intercept = 0,
    linetype = 'dotted',
    size = 0.3
  ) +
  geom_abline(
    slope = 0,
    intercept = 1,
    linetype = 'dotted',
    size = 0.3
  ) +
  geom_abline(
    slope = 0,
    intercept = -1,
    linetype = 'dotted',
    size = 0.3
  ) +
  geom_line(size = 0.5) +
  geom_point(size = 1) +
  facet_wrap(~ Site, scales = "fixed", ncol = 2) +
  xlab('Year') +
  ylab('Stress Index \n (Scaled and Centered)') +
  scale_y_continuous(limits = c(-1.5, 2.5),
                     breaks = seq(from = -2, to = 2.5, by = 1)) +
  ggtitle('Low Salinity Sites') +
  theme_bw() +
  theme(
    legend.position = 'none',
    plot.title = element_text(face = 'bold', hjust = 0.5, size = 10),
    text = element_text(size = 8),
    axis.text = element_text(size = 6),
    strip.background = element_blank(),
    strip.text = element_text(size = 8),
    panel.grid = element_blank(),
    plot.margin = margin(c(0.2, 0.7, 0.2, 0.2), unit = 'cm'),
    panel.spacing = unit(1, "lines")
  )

Plot6 <-
  ggplot(subset(
    sal_summary2,
    Site == 'Yellow Point' |
      Site == 'Navvy Jack Point' |
      Site == 'Figurehead Point' |
      Site == 'White Rock'
  ),
  aes(x = Year, y = stressnew)
  ) +
  geom_abline(
    slope = 0,
    intercept = 0,
    linetype = 'dotted',
    size = 0.3
  ) +
  geom_abline(
    slope = 0,
    intercept = 1,
    linetype = 'dotted',
    size = 0.3
  ) +
  geom_abline(
    slope = 0,
    intercept = -1,
    linetype = 'dotted',
    size = 0.3
  ) +
  geom_line(size = 0.5) +
  geom_point(size = 1) +
  facet_wrap( ~ Site, scales = "fixed", ncol = 1) +
  xlab('Year') +
  ylab('Stress Index \n (Scaled and Centered)') +
  scale_y_continuous(limits = c(-1.5, 2.5),
                     breaks = seq(from = -2, to = 2.5, by = 1)) +
  ggtitle('High Salinity Sites') +
  theme_bw() +
  theme(
    legend.position = 'none',
    plot.title = element_text(face = 'bold', hjust = 0.5, size = 10),
    text = element_text(size = 8),
    axis.text = element_text(size = 6),
    strip.background = element_blank(),
    strip.text = element_text(size = 8),
    panel.grid = element_blank(),
    plot.margin = margin(c(0.2, 0.7, 0.2, 0.2), unit = 'cm'),
    axis.title.y = element_blank(),
    panel.spacing = unit(1, "lines")
  )


png(
  filename = "Annual Salinity Stress By Site.png",
  width = 15.24,
  height = 9,
  units = "cm",
  pointsize = 10,
  res = 600
)

plot_grid(Plot5,
          Plot6,
          align = 'h',
          rel_widths = c(2, 1))

dev.off()


## Do some modelling of factors that could have been historically affecting
## N. lamellosa populations

highstress <- subset(sal_summary, stress > 2.71)  ## Excludes high salinity sites
highstress$Year <- as.numeric(highstress$Year)

sal.mod <- lmer(stress ~ Year + (1 | Site), data = highstress)
summary(sal.mod)  ## Each year stress goes up by ~ 0.001
plot(resid(sal.mod) ~ fitted(sal.mod))  ## Fit looks good
plot(resid(sal.mod) ~ highstress$Year)  ## Not too bad
plot(resid(sal.mod) ~ highstress$Site)  ## Slightly decreasing variance with decreasing stress at a site

null.sal.mod <- lmer(stress ~ 1 + (1 | Site), data = highstress)

anova(sal.mod, null.sal.mod)  ## p = 0.003

highstress$predict <- predict(sal.mod)

variancemod <- lmer(variance ~ Year + (1 | Site), data = highstress)
summary(variancemod)
plot(resid(variancemod) ~ fitted(variancemod))
variancenull <- lmer(variance ~ (1 | Site), data = highstress)
anova(variancemod, variancenull)  # p = 0.0.4

highstress$predictvariance <- predict(variancemod)

FRtheme <-   
  theme_bw() +
  theme(
    text = element_text(size = 8),
    axis.text = element_text(size = 6),
    panel.grid = element_blank(),
    plot.margin = margin(c(0.2, 0.7, 0.2, 0.5), unit = 'cm')
  )

TSsal <- subset(sal_summary, 
                Year <= 2016 & 
                  Site == 'Tower Beach South')

mod1 <- lm(under15 ~ Year, data = TSsal)
plot(mod1)  # seems fine
summary(mod1)  # p = 0.73, R2 = -0.03

under15plot <-
  ggplot(TSsal) +
  geom_point(aes(x = Year,
                 y = under15),
             size = 0.5) +
  geom_ribbon(
    aes(
      x = Year,
      ymin = predict(mod1) - predict(mod1, se.fit = TRUE)$se.fit,
      ymax = predict(mod1) + predict(mod1, se.fit = TRUE)$se.fit
    ),
    fill = npalette[1],
    size = 0.5,
    alpha = 0.5
  ) +
  geom_line(aes(x = Year,
                y = predict(mod1)),
            linetype = 'dashed', size = 0.5) +
  xlab('Year') +
  ylab('Days Under 15 psu') +
  FRtheme

mod2 <- glm(under9 ~ Year, data = TSsal, family = 'poisson')
plot(mod2)  # seems decent
summary(mod2)  # p < 0.001

under9plot <- 
  ggplot(TSsal) +
  geom_point(aes(x = Year,
                 y = under9),
             size = 0.5) +
  geom_ribbon(
    aes(
      x = Year,
      ymin = exp(predict(mod2) - predict(mod2, se.fit = TRUE)$se.fit),
      ymax = exp(predict(mod2) + predict(mod2, se.fit = TRUE)$se.fit)
    ),
    fill = npalette[1],
    size = 0.5,
    alpha = 0.5
  ) +
  geom_line(aes(x = Year,
                y = exp(predict(mod2))),
            linetype = 'dashed', size = 0.5) +
  xlab('Year') +
  ylab('Days Under 9 psu') +
  theme_classic() + 
  FRtheme +
  annotate('text', x = 1987, y = 65, label = '*', size = 10)

mod3 <- lm(log(min.sal) ~ Year, data = TSsal)
plot(mod3)  # seems decent
summary(mod3)  # p = 0.56, R2 = -0.02

minsalplot <- 
  ggplot(TSsal) +
  geom_point(aes(x = Year,
                 y = min.sal),
             size = 0.5) +
  geom_ribbon(
    aes(
      x = Year,
      ymin = exp(predict(mod3) - predict(mod3, se.fit = TRUE)$se.fit),
      ymax = exp(predict(mod3) + predict(mod3, se.fit = TRUE)$se.fit)
    ),
    fill = npalette[1],
    size = 0.5,
    alpha = 0.5
  ) +
  geom_line(aes(x = Year,
                y = exp(predict(mod3))),
            linetype = 'dashed', size = 0.5) +
  xlab('Year') +
  ylab('Annual Salinity\n Minimum') +
  theme_classic() + 
  FRtheme

mod4 <- lm(mean.sal ~ Year, data = TSsal)
plot(mod4)  # looks pretty good
summary(mod4)  # p = 0.12, R2 = 0.05

meansalplot <- 
  ggplot(TSsal) +
  geom_point(aes(x = Year,
                 y = mean.sal),
             size = 0.5) +
  geom_ribbon(
    aes(
      x = Year,
      ymin = predict(mod4) - predict(mod4, se.fit = TRUE)$se.fit,
      ymax = predict(mod4) + predict(mod4, se.fit = TRUE)$se.fit
    ),
    fill = npalette[1],
    size = 0.5,
    alpha = 0.5
  ) +
  geom_line(aes(x = Year,
                y = predict(mod4)),
            linetype = 'dashed', size = 0.5) +
  xlab('Year') +
  ylab('Annual Salinity Mean') +
  theme_classic() + 
  FRtheme

mod5 <- lm(min.date ~ Year, data = TSsal)
plot(mod5)  # looks pretty good
summary(mod5)  # p = 0.59, R2 = -0.02

mindateplot <- 
  ggplot(TSsal) +
  geom_point(aes(x = Year,
                 y = min.date),
             size = 0.5) +
  geom_ribbon(
    aes(
      x = Year,
      ymin = predict(mod5) - predict(mod5, se.fit = TRUE)$se.fit,
      ymax = predict(mod5) + predict(mod5, se.fit = TRUE)$se.fit
    ),
    fill = npalette[1],
    size = 0.5,
    alpha = 0.5
  ) +
  geom_line(aes(x = Year,
                y = predict(mod5)),
            linetype = 'dashed', size = 0.5) +
  xlab('Year') +
  ylab('Date of Annual\n Salinity Minimum') +
  theme_classic() + 
  FRtheme

mod6 <- lm(stress ~ Year, data = TSsal)
plot(mod6)  # looks pretty good
summary(mod6)  # p = 0.2, R2 = 0.01

stressplot <- 
  ggplot(TSsal) +
  geom_point(aes(x = Year,
                 y = stress),
             size = 0.5) +
  geom_ribbon(
    aes(
      x = Year,
      ymin = predict(mod6) - predict(mod6, se.fit = TRUE)$se.fit,
      ymax = predict(mod6) + predict(mod6, se.fit = TRUE)$se.fit
    ),
    fill = npalette[1],
    size = 0.5,
    alpha = 0.5
  ) +
  geom_line(aes(x = Year,
                y = predict(mod6)),
            linetype = 'dashed', size = 0.5) +
  xlab('Year') +
  ylab('Annual Salinity Stress Metric') +
  theme_classic() + 
  FRtheme

mod7 <- lm(variance ~ Year, data = TSsal)
plot(mod7)  # looks pretty good
summary(mod7)  # p = 0.32, R2 = 0.00

varianceplot <- 
  ggplot(TSsal) +
  geom_point(aes(x = Year,
                 y = variance),
             size = 0.5) +
  geom_ribbon(
    aes(
      x = Year,
      ymin = predict(mod7) - predict(mod7, se.fit = TRUE)$se.fit,
      ymax = predict(mod7) + predict(mod7, se.fit = TRUE)$se.fit
    ),
    fill = npalette[1],
    size = 0.5,
    alpha = 0.5
  ) +
  geom_line(aes(x = Year,
                y = predict(mod7)),
            linetype = 'dashed', size = 0.5) +
  xlab('Year') +
  xlab('Year') + 
  ylab('Annual Salinity Variance') +
  theme_classic() + 
  FRtheme

png(
  filename = "Patterns Since 1986.png",
  width = 15.24,
  height =11,
  units = "cm",
  pointsize = 10,
  res = 600
)


plot_grid(
  under15plot,
  under9plot,
  minsalplot,
  meansalplot,
  mindateplot,
  varianceplot,
  stressplot,
  labels = c('A', 'B', 'C', 'D', 'E', 'F', 'G'),
  nrows = 3,
  ncol = 3,
  label_size = 12
)

dev.off()


## Population analyses

# Load required packages.

library(ggplot2)
library(plyr)
library(lme4)
library(cowplot)
library(glmmTMB)
library(MuMIn)
library(glmmTMB)
library(colorspace)

# Set up data.

nuc.abund <-
  read.csv("~nuc.abund.csv", header = TRUE)

summary(nuc.abund)
nuc.abund$date <- as.Date(nuc.abund$date, "%d-%b-%y")
nuc.abund$live.tot <- as.numeric(nuc.abund$live.tot)

## Convert to be in terms of m^2 (density)

names(nuc.abund)
nuc.abund[,6:11] <- nuc.abund[,6:11]/10
head(nuc.abund[,6:11])

nuc.abund$site <-
  mapvalues(
    nuc.abund$site,
    from = c(levels(nuc.abund$site)),
    to = c(
      'Dunbar',
      'Figurehead Point',
      'Kitsilano Point',
      'Point Atikinson',
      'Tower Beach North',
      'Tower Beach South',
      'Waterloo',
      'Weston Park',
      'Navvy Jack Point',
      'White Rock',
      'Yellow Point'
    )
  )

nuc.abund$site <-
  factor(
    nuc.abund$site,
    levels = c(
      'Tower Beach South',
      'Tower Beach North',
      'Dunbar',
      'Waterloo',
      'Point Atikinson',
      'Kitsilano Point',
      'Weston Park',
      'Navvy Jack Point',
      'Figurehead Point',
      'White Rock',
      'Yellow Point'
    )
  )

# Average abundance data over the transects that were surveyed.

nuc.abund2 <- subset(nuc.abund, live.tot != 'NA')

ave.abund <-
  ddply(
    nuc.abund2,
    c('year', 'date', 'salinity.regime', 'site'),
    summarize,
    live.tot.mean = mean(live.tot),
    live.tot.se = sd(live.tot)/sqrt(length(live.tot)),
    live.adults.mean = mean(live.adults),
    live.adults.se = sd(live.adults)/sqrt(length(live.adults)),
    live.juv.mean = mean(live.juveniles),
    live.juv.se = sd(live.juveniles)/sqrt(length(live.juveniles)),
    empty.tot.mean = mean(empty.tot, na.rm = TRUE),
    empty.tot.se = sd(empty.tot, na.rm = TRUE)/sqrt(length(empty.tot)),
    egg.mass.mean = mean(egg.mass, na.rm = TRUE),
    egg.mass.se = sd(egg.mass, na.rm = TRUE)/sqrt(length(egg.mass))
  )
head(ave.abund)
summary(ave.abund)

# Make some plots of abundances.

## Plot abundances over time at different sites

A <-
  ggplot(subset(ave.abund, salinity.regime == 'fresh'),
         aes(x = date, y = live.tot.mean)) +
  geom_ribbon(
    aes(ymin = live.tot.mean - live.tot.se,
        ymax = live.tot.mean + live.tot.se),
    alpha = 0.5,
    colour = 'black',
    fill = npalette[2],
    size = 0.3,
    show.legend = FALSE
  ) +
  geom_line(size = 0.3) +
  geom_point(size = 1) +
  facet_wrap( ~ site, scales = 'fixed', ncol = 2) +
  xlab('Date') +
  ylab(expression(paste('Abundance (ind '~m^-2*')'))) +
  scale_x_date(
    expand = c(0, 0),
    date_breaks = '6 months',
    date_labels = '%b. %Y',
    limits = as.Date(c('2011-04-01',
                       '2014-10-01'))
  ) +
  scale_y_continuous(limits = c(0,25)) +
  ggtitle('Low Salinity Sites') +
  theme_bw() +
  theme(
    plot.title = element_text(face = 'bold', hjust = 0.5, size = 10),
    text = element_text(size = 8),
    axis.text = element_text(size = 6),
    axis.text.x = element_text(
      angle = -45,
      vjust = -0.5,
      hjust = 0.3
    ),
    strip.background = element_blank(),
    strip.text = element_text(size = 8),
    panel.grid = element_blank(),
    plot.margin = margin(c(0.2, 0.7, 0.2, 0.2), unit = 'cm')
  )

B <-
  ggplot(subset(ave.abund, salinity.regime == 'salty'),
         aes(x = date, y = live.tot.mean)) +
  geom_ribbon(
    aes(ymin = live.tot.mean - live.tot.se,
        ymax = live.tot.mean + live.tot.se),
    alpha = 0.5,
    colour = 'black',
    fill = npalette[2],
    size = 0.3,
    show.legend = FALSE
  ) +
  geom_line(size = 0.3) +
  geom_point(size = 1) +
  facet_wrap( ~ site, scales = 'fixed', ncol = 1) +
  xlab('Date') +
  ylab(expression(paste('Abundance (ind '~m^-2))) +
  scale_x_date(
    expand = c(0, 0),
    date_breaks = '6 months',
    date_labels = '%b. %Y',
    limits = as.Date(c('2011-04-01',
                       '2014-10-01'))
  ) +
  scale_y_continuous(limits = c(0,25)) +
  ggtitle('High Salinity Sites') +
  theme_bw() +
  theme(
    plot.title = element_text(face = 'bold', hjust = 0.5, size = 10),
    text = element_text(size = 8),
    axis.text = element_text(size = 6),
    axis.text.x = element_text(
      angle = -45,
      vjust = -0.5,
      hjust = 0.3
    ),
    strip.background = element_blank(),
    strip.text =  element_text(size = 8),
    axis.title.y = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(c(0.2, 0.7, 0.2, 0.2), unit = 'cm')
  )

png(
  filename = "Nucella abundances by site.png",
  width = 15.24,
  height = 13,
  units = "cm",
  pointsize = 10,
  res = 600
)

plot_grid(A,
          B,
          align = 'h',
          rel_widths = c(1.8, 1))

dev.off()

## Plot number dead on top of live abundances

nuc.abund2.1 <- subset(nuc.abund2, empty.tot != 'NA')

summary(nuc.abund2.1)

nuc.abund2.1$prop.dead <-
  with(nuc.abund2.1, empty.tot / (live.tot + empty.tot))

prop_dead <- ddply(
  nuc.abund2.1,
  c('year', 'date', 'salinity.regime', 'site'),
  summarize,
  mean.prop.dead = mean(prop.dead, na.rm = TRUE),
  se.prop.dead = sd(prop.dead, na.rm = TRUE)/sqrt(length(prop.dead))
)

summary(prop_dead)

prop_dead$mean.prop.dead[is.na(prop_dead$mean.prop.dead)] <- 0

C <-
  ggplot(subset(prop_dead, salinity.regime == 'fresh'),
         aes(x = date, y = mean.prop.dead)) +
  geom_ribbon(
    aes(ymin = mean.prop.dead - se.prop.dead,
        ymax = mean.prop.dead + se.prop.dead),
    alpha = 0.1,
    colour = 'black',
    fill = npalette[4],
    show.legend = FALSE
  ) +
  geom_line(aes(x = date, y = mean.prop.dead), size = 0.3) +
  geom_point(aes(x = date, y = mean.prop.dead), size = 1) +
  facet_wrap( ~ site, scales = 'fixed', ncol = 2) +
  xlab('Date') +
  ylab('Proportion of Shells that are Dead') +
  scale_x_date(
    expand = c(0, 0),
    date_breaks = '6 months',
    date_labels = '%b. %Y',
    limits = as.Date(c('2011-04-01',
                       '2014-10-01'))
  ) +
  scale_y_continuous(limits = c(0,1.1), breaks = c(0,0.25,0.5,0.75,1)) +
  ggtitle('Low Salinity Sites') +
  theme_bw() +
  theme(
    plot.title = element_text(face = 'bold', hjust = 0.5, size = 10),
    text = element_text(size = 8),
    axis.text = element_text(size = 6),
    axis.text.x = element_text(
      angle = -45,
      vjust = -0.5,
      hjust = 0.3
    ),
    strip.background = element_blank(),
    strip.text = element_text(size = 8),
    panel.grid = element_blank(),
    plot.margin = margin(c(0.2, 0.7, 0.2, 0.2), unit = 'cm')
  )

D<-
  ggplot(subset(prop_dead, salinity.regime == 'salty'),
         aes(x = date, y = mean.prop.dead)) +
  geom_ribbon(
    aes(ymin = mean.prop.dead - se.prop.dead,
        ymax = mean.prop.dead + se.prop.dead),
    alpha = 0.1,
    colour = 'black',
    fill = npalette[4],
    show.legend = FALSE
  ) +
  geom_line(aes(x = date, y = mean.prop.dead), size = 0.3) +
  geom_point(aes(x = date, y = mean.prop.dead), size = 1) +
  facet_wrap( ~ site, scales = 'fixed', ncol = 1) +
  xlab('Date') +
  ylab('Proportion of Shells that are Dead') +
  scale_x_date(
    expand = c(0, 0),
    date_breaks = '6 months',
    date_labels = '%b. %Y',
    limits = as.Date(c('2011-04-01',
                       '2014-10-01'))
  ) +
  scale_y_continuous(limits = c(0,1.1), breaks = c(0,0.25,0.5,0.75,1)) +
  ggtitle('High Salinity Sites') +
  theme_bw() +
  theme(
    plot.title = element_text(face = 'bold', hjust = 0.5, size = 10),
    text = element_text(size = 8),
    axis.text = element_text(size = 6),
    axis.text.x = element_text(
      angle = -45,
      vjust = -0.5,
      hjust = 0.3
    ),
    strip.background = element_blank(),
    strip.text =  element_text(size = 8),
    axis.title.y = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(c(0.2, 0.7, 0.2, 0.2), unit = 'cm')
  )

png(
  filename = "Nucella dead by site.png",
  width = 15.24,
  height = 10,
  units = "cm",
  pointsize = 10,
  res = 600
)

plot_grid(C,
          D,
          align = 'h',
          rel_widths = c(1.8, 1))

dev.off()

## Analysis of proportion juvenile

## Pull out the measurements we have, by transect, for 2012-2014

nuc.abund3 <- subset(nuc.abund2, year != '2011' &
                       site != 'Tower South' & live.juveniles != 'NA')


nuc.abund3$prop.live.juvies <-
  with(nuc.abund3, live.juveniles / live.tot)

prop.ave <- ddply(
  nuc.abund3,
  c('year', 'date', 'salinity.regime', 'site'),
  summarize,
  mean.prop.live.juvies = mean(prop.live.juvies, na.rm = TRUE),
  se.prop.live.juvies = 
    sd(prop.live.juvies, na.rm = TRUE)/sqrt(length(prop.live.juvies))
)

summary(prop.ave)

# Plot.

## Proportion live, adults vs. juvies

prop.ave[is.na(prop.ave)] <- 0

E <-
  ggplot(
    subset(prop.ave, salinity.regime == 'fresh'),
    aes(x = date, y = mean.prop.live.juvies)
  ) +
  geom_ribbon(
    aes(
      ymin = mean.prop.live.juvies - se.prop.live.juvies,
      ymax = mean.prop.live.juvies + se.prop.live.juvies
    ),
    alpha = 0.5,
    colour = 'black',
    fill = npalette[3],
    show.legend = FALSE
  ) +
  geom_line(size = 0.3) +
  geom_point(size = 1) +
  facet_wrap( ~ site, scales = 'fixed', ncol = 2) +
  xlab('Date') +
  ylab('Proportion Juveniles') +
  ggtitle('Low Salinity Sites') +
  scale_y_continuous(
    limits = c(-0.2, 1.1),
    breaks = c(0, 0.25, 0.50, 0.75, 1),
    expand = c(0, 0)
  ) +
  scale_x_date(
    expand = c(0, 0),
    date_breaks = '6 months',
    date_labels = '%b. %Y',
    limits = as.Date(c('2012-04-01',
                       '2014-10-01'))
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = 'bold', hjust = 0.5, size = 10),
    text = element_text(size = 8),
    axis.text = element_text(size = 6),
    strip.background = element_blank(),
    strip.text = element_text(size = 8),
    panel.grid = element_blank(),
    plot.margin = margin(c(0.2, 0.7, 0.2, 0.2), unit = 'cm'),
    axis.text.x = element_text(
      angle = -45,
      vjust = -0.5,
      hjust = 0.3
    )
  )

G <-
  ggplot(
    subset(prop.ave, salinity.regime == 'salty'),
    aes(x = date, y = mean.prop.live.juvies)
  ) +
  geom_ribbon(
    aes(
      ymin = mean.prop.live.juvies - se.prop.live.juvies,
      ymax = mean.prop.live.juvies + se.prop.live.juvies
    ),
    alpha = 0.5,
    colour = 'black',
    fill = npalette[3],
    show.legend = FALSE
  ) +
  geom_line(size = 0.3) +
  geom_point(size = 1) +
  facet_wrap( ~ site, scales = 'fixed', ncol = 1) +
  xlab('Date') +
  ylab('Proportion Juveniles') +
  ggtitle('High Salinity Sites') +
  scale_y_continuous(
    limits = c(-0.2, 1.1),
    breaks = c(0, 0.25, 0.50, 0.75, 1),
    expand = c(0, 0)
  ) +
  scale_x_date(
    expand = c(0, 0),
    date_breaks = '6 months',
    date_labels = '%b. %Y',
    limits = as.Date(c('2012-04-01',
                       '2014-10-01'))
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = 'bold', hjust = 0.5, size = 10),
    text = element_text(size = 8),
    axis.text = element_text(size = 6),
    strip.background = element_blank(),
    strip.text =  element_text(size = 8),
    axis.title.y = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(c(0.2, 0.7, 0.2, 0.2), unit = 'cm'),
    axis.text.x = element_text(
      angle = -45,
      vjust = -0.5,
      hjust = 0.3
    )
  )

png(
  filename = "Proportion live juveniles by site.png",
  width = 15.24,
  height = 10,
  units = "cm",
  pointsize = 10,
  res = 600
)

plot_grid(E,
          G,
          align = 'h',
          rel_widths = c(1.9, 1))

dev.off()


## Bring in salinity data
## Model which aspects of salinity most affect mortality in N. lamellosa total
## populations and proportion of adults vs. juveniles

## Use sal_summary, but subset for 2010-2013, for only low salinity sites

sal_char <-
  subset(
    sal_summary,
    Year > 2009 & Year < 2014
  )

head(sal_char)
summary(sal_char)

sal_char$Year <- as.character(sal_char$Year)
sal_char$Year <- as.numeric(sal_char$Year)
sal_char$Year <- sal_char$Year + 1
sal_char$Year <- as.character(sal_char$Year)
sal_char$Year <- format(as.Date(sal_char$Year, format = '%Y'), '%Y')

# Add next year abundances
summary(nuc.abund2)
nuc.abund2$year <- format(as.Date(as.character(nuc.abund2$year),
                                  format = '%Y'), '%Y')

join.abund <- nuc.abund2[c('year',
                           'site',
                           'transect',
                           'num.tran',
                           'live.tot',
                           'live.adults',
                           'live.juveniles',
                           'empty.tot',
                           'empty.adult',
                           'empty.juv',
                           'egg.mass')]
head(join.abund)

join.abund$Year <- join.abund$year
join.abund$Site <- join.abund$site
join.abund <- join.abund[c('Site',
                           'Year',
                           'transect',
                           'num.tran',
                           'live.tot',
                           'live.adults',
                           'live.juveniles',
                           'empty.tot',
                           'empty.adult',
                           'empty.juv',
                           'egg.mass')]

library(dplyr)


pop_mod <- left_join(sal_char, join.abund, by = c('Site', 'Year'))

head(pop_mod)
summary(pop_mod)

pop_mod2 <- subset(pop_mod, live.tot != 'NA')
summary(pop_mod2)
pop_mod2$Site <- as.factor(pop_mod2$Site)
pop_mod2$Year <- as.character(pop_mod2$Year)
pop_mod2$Year <- format(as.Date(pop_mod2$Year, format = '%Y'), '%Y')

## First centre and scale predictors

head(pop_mod2)

pop_mod2$under15new <-
  scale(pop_mod2$under15, center = TRUE, scale = TRUE)
pop_mod2$min.salnew <-
  scale(pop_mod2$min.sal, center = TRUE, scale = TRUE)
pop_mod2$mean.salnew <-
  scale(pop_mod2$mean.sal, center = TRUE, scale = TRUE)
pop_mod2$min.datenew <-
  scale(pop_mod2$min.date, center = TRUE, scale = TRUE)
pop_mod2$under15.datenew <-
  scale(pop_mod2$under15.date, center = TRUE, scale = TRUE)
pop_mod2$stressnew <-
  scale(pop_mod2$stress, center = TRUE, scale = TRUE)
pop_mod2$variancenew <-
  scale(pop_mod2$variance, center = TRUE, scale = TRUE)

pop_mod2$ID <- with(pop_mod2, interaction(Site,  transect))  # Add unique ID for transect

## Specify overdispersion function

dispfun <- function(m) {
  r <- residuals(m,type="pearson")
  n <- df.residual(m)
  dsq <- sum(r^2)
  c(dsq=dsq,n=n,disp=dsq/n)
}

## ZI Poisson model
m1 <-
  glmmTMB(
    live.tot ~  stressnew * Site +
      (1 + 1 | Year),
    family = poisson(link = log),
    ziformula = ~ 1,
    data = pop_mod2
  )
summary(m1)  # convergence problem

## Poisson model
m2 <-
  glmmTMB(live.tot ~  stressnew * Site +
            (1 + 1 | Year),
          family = 'poisson',
          data = pop_mod2)
summary(m2)
dispfun(m2)  # overdispersed

## ZI Negative Binomial model 1
m3 <-
  glmmTMB(
    live.tot ~  stressnew * Site +
      (1 + 1 | Year),
    family = nbinom1(),
    ziformula = ~ 1,
    data = pop_mod2
  )
summary(m3)  # dispersion is better

## ZI Negative Binomial model 2
m3.1 <-
  glmmTMB(
    live.tot ~  stressnew * Site +
      (1 + 1 | Year),
    family = nbinom2(),
    ziformula = ~ 1,
    data = pop_mod2
  )
summary(m3.1) # overdispersed

## Negative Binomial model
m4 <-
  glmmTMB(live.tot ~  stressnew * Site +
            (1 + 1 | Year),
          family = nbinom2(),
          data = pop_mod2)
summary(m4)
dispfun(m4)  # good

m5 <-   glmmTMB(live.tot ~  stressnew * Site +
                  (1 + 1 | Year),
                family = nbinom1(),
                data = pop_mod2)
dispfun(m5)

AICc(m1, m2, m3.1, m3, m4, m5)

## m4, the NB model gives the best fit according to AICc
## And the dispersion looks good

plot(resid(m4, type = 'pearson') ~ fitted(m4))
abline(0,0)

## Try without random effects

m6 <-
  glmmTMB(live.tot ~  stressnew * Site,
          family = 'nbinom2',
          data = pop_mod2)

m7 <-
  glmmTMB(live.tot ~  stressnew * Site,
          family = 'nbinom1',
          data = pop_mod2)

AICc(m4, m6, m7) # better fit with no random effects and log link

summary(m7)
dispfun(m7)  # dispersion parameter = 1.26

plot(resid(m7, type = 'pearson') ~ fitted(m7)) # pattern in variance
abline(0, 0)

plot(resid(m7) ~ pop_mod2$Site)
plot(resid(m7) ~ pop_mod2$stressnew)

m8<- update(m7, . ~ . - stressnew:Site)
m9 <- update(m8, . ~ . - stressnew)
m10 <- update(m8, . ~ . - Site)

AICc(m7, m8, m9, m10)  # doesn't improve fit to take any terms out

anova(m7, m8)  # interaction is significant at p < 0.001
anova(m8, m9) # stress significant at p < 0.001
anova(m7, m10) # site significant at p < 0.001

library(DHARMa)

simulateResiduals(
  fittedModel = m7,
  n = 1000,
  plot = T,
  seed = 1
)

# fit looks awesome

## Plot model predictions

pop_mod2$min.pred <- (predict(m7, type = 'response') - 
                        predict(m7, type = 'response',
                                se.fit = TRUE)$se.fit)

pop_mod2$max.pred <- (predict(m7, type = 'response') + 
                        predict(m7, type = 'response',
                                se.fit = TRUE)$se.fit)

pop_mod2$pred <- predict(m7, type = 'response')

###
pred1 <-
  ggplot(subset(pop_mod2, 
                Site != 'Yellow Point' &
                  Site != 'Navvy Jack Point' &
                  Site != 'Figurehead Point' &
                  Site != 'White Rock'), 
         aes(x = stressnew, y = live.tot)) +
  geom_jitter(size = 0.3) +
  geom_ribbon(
    aes(
      ymin = min.pred,
      ymax = max.pred
    ),
    alpha = 0.5,
    fill = npalette[1]
  ) +
  geom_line(aes(x = stressnew,
                y = pred),
            linetype = 'dashed',
            size = 0.5) +
  facet_wrap( ~ Site, scales = 'fixed', ncol = 2) +
  xlab('Stress Index \n (Scaled and Centered)') +
  ylab(expression(paste('Abundance (ind '~m^-2*')'))) +
  scale_y_continuous(trans = 'log1p', breaks = c(0,1,5,10,20)) +
  theme_bw() +
  theme(
    legend.position = 'none',
    text = element_text(size = 8),
    axis.text = element_text(size = 7),
    strip.background = element_blank(),
    strip.text = element_text(size = 8),
    panel.grid = element_blank()
  )

pred2 <-
  ggplot(subset(pop_mod2, 
                Site == 'Yellow Point' |
                  Site == 'Navvy Jack Point' |
                  Site == 'Figurehead Point' |
                  Site == 'White Rock'), 
         aes(x = stressnew, y = live.tot)) +
  geom_jitter(size = 0.3) +
  geom_ribbon(
    aes(
      ymin = min.pred,
      ymax = max.pred
    ),
    alpha = 0.5,
    fill = npalette[1]
  ) +
  geom_line(aes(x = stressnew,
                y = pred),
            linetype = 'dashed',
            size = 0.5) +
  facet_wrap( ~ Site, scales = 'fixed', ncol = 1) +
  xlab('Stress Index \n (Scaled and Centered)') +
  ylab('') +
  scale_y_continuous(trans = 'log1p', breaks = c(0,1,5,10,20)) +
  theme_bw() +
  theme(
    legend.position = 'none',
    text = element_text(size = 8),
    axis.text = element_text(size = 7),
    strip.background = element_blank(),
    strip.text = element_text(size = 8),
    panel.grid = element_blank()
  )

png(
  filename = "Population Model Predictions.png",
  width = 15.24,
  height = 10,
  units = "cm",
  pointsize = 10,
  res = 600
)

plot_grid(pred1,
          pred2,
          align = 'h',
          rel_widths = c(2, 1))

dev.off()

## Hind-cast population estimates

nucella_pop <-
  subset(
    sal_summary,
    Site != 'Figurehead Point' &
      Site != 'Navvy Jack Point' &
      Site != 'White Rock' & Site != 'Yellow Point'
  )

nucella_pop$stressnew <-
  scale(nucella_pop$stress, center = TRUE, scale = TRUE)

nucella_pop$snailpop <- predict(m7, newdata = nucella_pop,
                                type = 'response')
nucella_pop$snailpopse <- predict(m7,
                                  newdata = nucella_pop,
                                  type = 'response',
                                  se.fit = TRUE)$se.fit

## Plot prediction since 1986

library(MASS)

ggplot(subset(nucella_pop, Site != 'Tower Beach North'), 
       aes(x = Year, y = floor(snailpop))) +
  geom_line(colour = 'skyblue') +
  facet_wrap( ~ Site, scales = 'fixed') +
  xlab('Year') +
  ylab('Predicted Abundance at Year + 1') +
  scale_y_continuous(limits = c(0,100)) +
  scale_x_continuous(
    limits = c(1986, 2016),
    breaks = c(seq(
      from = 1985, to = 2015, by = 5
    ))
  ) +
  theme_bw() +
  theme(
    legend.position = 'none',
    text = element_text(size = 10),
    panel.spacing = unit(2, "lines"),
    axis.text = element_text(size = 6),
    strip.background = element_blank(),
    strip.text.x = element_text(size = 10),
    panel.grid = element_blank(),
    plot.margin = margin(0, 1, 0, 1, 'cm')
  )

## Now model effects on ratio juveniles to adults

pop_mod3 <- subset(pop_mod, live.adults != 'NA' & Year > 2011)
summary(pop_mod3)
pop_mod3$Year <- as.factor(pop_mod3$Year)

pop_mod3$Site <- as.factor(pop_mod3$Site)
pop_mod3$under15new <-
  scale(pop_mod3$under15, center = TRUE, scale = TRUE)
pop_mod3$min.salnew <-
  scale(pop_mod3$min.sal, center = TRUE, scale = TRUE)
pop_mod3$mean.salnew <-
  scale(pop_mod3$mean.sal, center = TRUE, scale = TRUE)
pop_mod3$min.datenew <-
  scale(pop_mod3$min.date, center = TRUE, scale = TRUE)
pop_mod3$under15.datenew <-
  scale(pop_mod3$under15.date, center = TRUE, scale = TRUE)
pop_mod3$stressnew <-
  scale(pop_mod3$stress, center = TRUE, scale = TRUE)
pop_mod3$ID <- 
  with(pop_mod3, interaction(Site,  transect))  # Add unique ID for transect

head(pop_mod3)

## Now model

pop_mod3$percent.juveniles <- with(pop_mod3, live.juveniles / live.tot) * 100

pop_mod3$percent.juveniles[is.na(pop_mod3$percent.juveniles)] <- 0

pop_mod3$prop.juveniles <- with(pop_mod3, live.juveniles / live.tot) 

pop_mod3$prop.juveniles[is.na(pop_mod3$prop.juveniles)] <- 0

summary(pop_mod3)

## cbind model 1
model1 <-
  glmmTMB(
    cbind(live.juveniles, live.adults) ~  stressnew +
      (1 | Year/Site),
    family = binomial(),
    data = pop_mod3,
    REML = TRUE
  )
summary(model1)
plot(resid(model1, type = 'pearson') ~ fitted(model1))
abline(0,0)  # doesn't look too bad

## Try with different random effect structure

model2 <- 
  glmmTMB(
    cbind(live.juveniles, live.adults) ~  stressnew +
      (1 | Year),
    family = binomial(),
    data = pop_mod3,
    REML = TRUE
  )

model3 <- 
  glmmTMB(
    cbind(live.juveniles, live.adults) ~  stressnew +
      (1 | Site),
    family = binomial(),
    data = pop_mod3,
    REML = TRUE
  )

AICc(model1, model2, model3)  ## Better fit with all random effects

## Compare with no random effects

model4 <- 
  glmmTMB(
    cbind(live.juveniles, live.adults) ~  stressnew +
      (1 | Year/Site),
    family = binomial(),
    data = pop_mod3,
    REML = FALSE
  )

model5 <- 
  glm(
    cbind(live.juveniles, live.adults) ~  stressnew,
    family = binomial(),
    data = pop_mod3
  )

AICc(model4, model5)  # best fit with no random effects

plot(resid(model5, type = 'pearson') ~ pop_mod3$stressnew)
abline(0,0)

plot(resid(model5, type = 'pearson') ~ pop_mod3$Site)

## Plot predictions

juvmodplot <-
  ggplot(pop_mod3,
         aes(x = stressnew, y = prop.juveniles)) +
  geom_ribbon(
    aes(
      x = stressnew,
      ymin = predict(model5, type = 'response') -
        predict(model5, type = 'response', se.fit = TRUE)$se.fit,
      ymax = predict(model5, type = 'response') +
        predict(model5, type = 'response', se.fit = TRUE)$se.fit
    ),
    size = 0.5,
    alpha = 0.5,
    fill = npalette[1]
  ) +
  geom_line(aes(x = stressnew, y = predict(model5, type = 'response')),
            size = 0.3,
            linetype = 'dashed') +
  geom_jitter(size = 0.3) +
  xlab('Stress Index \n (Scaled & Centered)') +
  ylab('Proportion Juveniles') +
  coord_cartesian(ylim = c(0, 1),
                  xlim = c(-1.5, 2.5)) +
  theme_bw() +
  theme(
    plot.title = element_text(face = 'bold', hjust = 0.5, size = 10),
    text = element_text(size = 8),
    axis.text = element_text(size = 7),
    strip.background = element_blank(),
    strip.text = element_text(size = 8),
    panel.grid = element_blank()
  )

## Now model proportion live vs. dead

pop_mod4 <- subset(pop_mod, empty.tot != 'NA')

summary(pop_mod4)
pop_mod4$Year <- as.character(pop_mod4$Year)
pop_mod4$Year <- as.factor(pop_mod4$Year)
pop_mod4$Site <- as.factor(pop_mod4$Site)

pop_mod4$under15new <-
  scale(pop_mod4$under15, center = TRUE, scale = TRUE)
pop_mod4$min.salnew <-
  scale(pop_mod4$min.sal, center = TRUE, scale = TRUE)
pop_mod4$mean.salnew <-
  scale(pop_mod4$mean.sal, center = TRUE, scale = TRUE)
pop_mod4$min.datenew <-
  scale(pop_mod4$min.date, center = TRUE, scale = TRUE)
pop_mod4$under15.datenew <-
  scale(pop_mod4$under15.date, center = TRUE, scale = TRUE)
pop_mod4$stressnew <-
  scale(pop_mod4$stress, center = TRUE, scale = TRUE)

M1 <-
  glmmTMB(
    cbind(empty.tot, live.tot) ~ stressnew +
      (1 | Year / Site),
    family = binomial(link = 'cloglog'),
    data = pop_mod4
  )
summary(M1)
plot(resid(M1, type = 'pearson') ~ fitted(M1))
abline(0,0)

M2 <-
  glmmTMB(
    cbind(empty.tot, live.tot) ~  stressnew +
      (1 + 1 | Year / Site),
    family = binomial(link = 'logit'),
    data = pop_mod4
  )

summary(M2)
AICc(M1, M2) # better fit with cloglog

## Try removing individual random effects

M3 <-
  glmmTMB(
    cbind(empty.tot, live.tot) ~  stressnew +
      (1 + 1 | Year),
    family = binomial(link = 'cloglog'),
    data = pop_mod4,
    REML = TRUE
  )

M4 <-
  glmmTMB(
    cbind(empty.tot, live.tot) ~  stressnew +
      (1 + 1 | Site),
    family = binomial(link = 'cloglog'),
    data = pop_mod4,
    REML = TRUE
  )

## refit M2 with REML

M1.1 <- update(M1, REML = TRUE)

AICc(M1.1, M3, M4)  # better fit to keep in all random effects

## Try with no random effects

M5 <-
  glm(
    cbind(empty.tot, live.tot) ~  stressnew,
    family = binomial(link = cloglog),
    data = pop_mod4
  )

AICc(M1, M5)  # M1 is best fit

summary(M1)

plot(resid(M1, type = 'pearson') ~ fitted(M1))  # does look a bit overdispersed

## Try betabinomial distribution

M6 <- 
  glmmTMB(
    cbind(empty.tot, live.tot) ~  stressnew +
      (1 | Year / Site),
    family = betabinomial(link = 'cloglog'),
    data = pop_mod4
  )

AICc(M1, M6)  # better as regular binomial model

summary(M1)

## Plot predictions

pop_mod4$prop.dead <- 
  with(pop_mod4, empty.tot / (live.tot + empty.tot))

pop_mod4$prop.dead[is.na(pop_mod4$prop.dead)] <- 0

summary(pop_mod4$prop.dead)

## Generate model predictions

simdead <- list()
sedead <- list()

stressnew <- seq(from = min(pop_mod4$stressnew), 
                 to = max(pop_mod4$stressnew),
                 length.out = length(pop_mod4$prop.dead))
Year <- as.factor(rep('2013', length(pop_mod4$prop.dead)))
fake <- data.frame(stressnew, Year)

for (i in 1:length(unique(pop_mod4$Site))) {
  fake$Site <- 
    as.factor(rep(unique(pop_mod4$Site)[i], length(pop_mod4$prop.dead)))
  predict <- predict(M1, newdata = fake, type = 'response',
                     allow.new.levels = FALSE)
  se <- as.numeric(predict(M1, newdata = fake, type = 'response',
                           allow.new.levels = FALSE, se.fit = TRUE)$se.fit)
  simdead[[i]] <- predict
  sedead[[i]] <- se
}

simdead2 <- as.data.frame(do.call(cbind, simdead))
head(simdead2)

sedead2 <- as.data.frame(do.call(cbind, sedead))
head(sedead2)

deadmin <- simdead2 - sedead2
deadmax <- simdead2 + sedead2

deadmodplot <- ggplot(pop_mod4,
                      aes(x = stressnew, y = prop.dead)) +
  geom_line(
    aes(x = fake$stressnew, y = simdead2$V1),
    size = 0.3,
    linetype = 'dashed',
    colour = npalette[1]
  ) +
  geom_line(
    aes(x = fake$stressnew, y = simdead2$V2),
    size = 0.3,
    linetype = 'dashed',
    colour = npalette[1]
  ) +
  geom_line(
    aes(x = fake$stressnew, y = simdead2$V3),
    size = 0.3,
    linetype = 'dashed',
    colour = npalette[1]
  ) +
  geom_line(
    aes(x = fake$stressnew, y = simdead2$V4),
    size = 0.3,
    linetype = 'dashed',
    colour = npalette[1]
  ) +
  geom_line(
    aes(x = fake$stressnew, y = simdead2$V5),
    size = 0.3,
    linetype = 'dashed',
    colour = npalette[1]
  ) +
  geom_line(
    aes(x = fake$stressnew, y = simdead2$V6),
    size = 0.3,
    linetype = 'dashed',
    colour = npalette[1]
  ) +
  geom_line(
    aes(x = fake$stressnew, y = simdead2$V7),
    size = 0.3,
    linetype = 'dashed',
    colour = npalette[1]
  ) +
  geom_line(
    aes(x = fake$stressnew, y = simdead2$V8),
    size = 0.3,
    linetype = 'dashed',
    colour = npalette[1]
  ) +
  geom_line(
    aes(x = fake$stressnew, y = simdead2$V9),
    size = 0.3,
    linetype = 'dashed',
    colour = npalette[1]
  ) +
  geom_line(
    aes(x = fake$stressnew, y = simdead2$V10),
    size = 0.3,
    linetype = 'dashed',
    colour = npalette[1]
  ) +
  geom_line(
    aes(x = fake$stressnew, y = simdead2$V11),
    size = 0.3,
    linetype = 'dashed',
    colour = npalette[1]
  ) +
  geom_jitter(size = 0.3) +
  xlab('Stress Index \n (Scaled & Centered)') +
  ylab('Proportion Dead') +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw() +
  theme(
    plot.title = element_text(face = 'bold', hjust = 0.5, size = 10),
    text = element_text(size = 8),
    axis.text = element_text(size = 7),
    strip.background = element_blank(),
    strip.text = element_text(size = 8),
    panel.grid = element_blank()
  )

png(
  filename = "Juv and Death Model Predictions.png",
  width = 7.62,
  height = 4,
  units = "cm",
  pointsize = 10,
  res = 600
)

plot_grid(juvmodplot, deadmodplot)

dev.off()


