#final plots for thesis draft
#4/24/22
#load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load("tidyverse", "lubridate", "MASS", "mosaic", "ReIns", "patchwork", 
               "igraph", "caret", "scales", "goftest", "gofgamma",
               "CITAN", "cmstatr", "eva", "mclust", "emdist", "transport", "stats",
               "gpdAd")

#made in section 11 of final1ksample.R
data <- read.csv("HortonRatios/formatted_allSources.csv")
source('paretoFunctions.R')

#color palette
ms_color <- "#ffbe0b"
mc_color <- "#8338ec"
kz_color <- "#ff006e"
################################
#1. histograms at each spatial scale
################################

#all distributions, stacked
big <- ggplot(data, aes(width_m)) +
  geom_histogram(binwidth = 0.1,
                 color = "black",
                 aes(fill = source),
                 position = "identity") +
  #facet_grid(rows = vars(source))+
  theme_classic() +
  labs(y = "Frequency",
       x = "Width (m)")+
  # geom_line(data = d, aes(x = lineSeq, y = value, color = fit))+
  scale_y_log10(expand = expansion(mult = c(0, .1)),
                breaks = c(1, 10, 100),
                labels = c("1", "10", "100"))+
  scale_x_log10(expand = expansion(mult = c(0, .1)),
                breaks = c(0.1, 10, 1000),
                labels = c("0.1", "10", "1,000"))+
  scale_fill_manual(
    name = "",
    values = c(
      "Konza: Headwaters" = kz_color,
      "McDowell: Local" = mc_color,
      "Mississippi: Continental" = ms_color
    ),
    labels = c("Konza: Headwaters",
               "McDowell: Local",
               "Mississippi: Continental"))+
  theme(legend.position = "top")

#ms = mississippi, continental
ms <- filter(data, source == "Mississippi: Continental")

ms_widths <- ggplot(ms, aes(width_m)) +
  geom_histogram(binwidth = 0.1,
                 color = "black",
                 aes(fill = source)) +
  theme_classic() +
  labs(y = "Frequency",
       x = "Width (m)") +
  # geom_line(data = d, aes(x = lineSeq, y = value, color = fit))+
  scale_y_log10(expand = expansion(mult = c(0, .1)),
                breaks = c(1, 10, 100),
                labels = c("1", "10", "100")) +
  scale_x_log10(expand = expansion(mult = c(0, .1)),
                breaks = c(1, 10, 100, 1000),
                labels = c("1", "10", "100", "1,000")) +
  scale_fill_manual(
    name = "",
    values = c("Mississippi: Continental" = "#ffbe0b"),
    labels = c("Mississippi: Continental")
  )+
  theme(legend.position = "")

ms_widths
#mc = mcdowell, regional
mc <- filter(data, source == "McDowell: Local")

mc_widths <- ggplot(mc, aes(width_m)) +
  geom_histogram(binwidth = 0.1,
                 color = "black",
                 aes(fill = source)) +
  theme_classic() +
  labs(y = "Frequency",
       x = "Width (m)") +
  # geom_line(data = d, aes(x = lineSeq, y = value, color = fit))+
  scale_y_log10(expand = expansion(mult = c(0, .1)),
                breaks = c(1, 10, 100),
                labels = c("1", "10", "100")) +
  scale_x_log10(expand = expansion(mult = c(0, .1)),
                breaks = c(0.03, 1, 3, 10, 30),
                labels = c("0.03", "1", "3", "10", "30")) +
  scale_fill_manual(
    name = "",
    values = c("McDowell: Local" = "#8338ec"),
    labels = c("McDowell: Local")
  )+
  theme(legend.position = "")

#kz = konza, headwaters
kz <- filter(data, source == "Konza: Headwaters")

kz_widths <- ggplot(kz, aes(width_m)) +
  geom_histogram(binwidth = 0.1,
                 color = "black",
                 aes(fill = source)) +
  theme_classic() +
  labs(y = "Frequency",
       x = "Width (m)") +
  # geom_line(data = d, aes(x = lineSeq, y = value, color = fit))+
  scale_y_log10(expand = expansion(mult = c(0, .1)),
                breaks = c(1, 10, 100),
                labels = c("1", "10", "100")) +
  scale_x_log10(expand = expansion(mult = c(0, .1)),
                breaks = c(0.01, 0.1, 1, 10),
                labels = c("0.01", "0.1", "1", "10")) +
  scale_fill_manual(
    name = "",
    values = c("Konza: Headwaters" = kz_color),
    labels = c("Konza: Headwaters")
  )+
  theme(legend.position = "")


# scale_x_log10(expand = expansion(mult = c(0, .1)),
#               breaks = c(1, 10, 100, 1000),
#               labels = c("1", "10", "100", "1,000"),
#               limits = c(1, 2000)
# )+
#   scale_y_log10(expand = expansion(mult = c(0, .1)),
#                 breaks = c(0.1, 1, 10, 100),
#                 labels = c("0.1" ,"1", "10", "100"),
#                 limits = c(1, 100)

#all <- 
  
  big / (kz_widths | mc_widths | ms_widths) + plot_annotation(tag_levels = 'A') + plot_layout(heights = c(2, 1))

# ggsave(
#   filename = "allfits.png",
#   plot = last_plot(),
#   width = 10,
#   height = 12,
#   units = "in"
# )
################################################
#2. fitting stat distributions to each spatial scale
#################################################
#code from section 2 of final1ksample.R

#trying to make a for loop what will produce all of the fits
sources <- unique(data$source)
for(i in 1:length(sources)){
  s = sources[i]
  t = dplyr::filter(data, source == s)
  
  w = t$width_m
  
  minW = min(w)
  int = 10
  
  dlnInt = 0.1
  lineSeq = seq(0, ceiling(max(w)/dlnInt)*dlnInt, dlnInt)
  breaks = seq(0, (max(w)+int), int)
  
  
  parfit = pareto.mle(w)#consider running using mode W, rn using min W
  dpar = dpareto(lineSeq, parfit[[1]], parfit[[2]])
  ppar = ppareto(lineSeq, parfit[[1]], parfit[[2]])
  
  gfit = suppressWarnings(fitdistr(w, "gamma")$estimate)
  dg = dgamma(lineSeq, gfit[1], gfit[2])
  pg = pgamma(lineSeq, gfit[1], gfit[2])
  # logNormal
  lnfit = fitdistr(w, "log-normal")$estimate
  dln = dlnorm(lineSeq, lnfit[1], lnfit[2])
  pln = plnorm(lineSeq, lnfit[1], lnfit[2])
  # weibull
  wfit = suppressWarnings(fitdistr(w, "weibull", lower = c(0,0))$estimate)
  dw = dweibull(lineSeq, wfit[1], wfit[2])
  pw = pweibull(lineSeq, wfit[1], wfit[2])
  
  #Anderson Darling test, many packages
  #lognormal, goftest package
  adl = ad.test(w, "plnorm", estimated = TRUE)
  adl.p = adl$p.value
  adl.t = unname(adl$statistic)
  #gamma, gofgamma package
  adg = test.AD(w)
  #pvalue not produced from this package
  adg.t = adg$T.value
  #weibull, cmstatr
  adw = anderson_darling_weibull(x = w)
  adw.p = adw$osl
  adw.t = adw$A
  #pareto, agop
  adp = pareto2_test_ad(w, s = parfit$xm)
  adp.p = adp$p.value
  adp.t = unname(adp$statistic)
  
  #Kolmogorov Smirnov test
  jw = jitter(w) # to remove ties
  # gamma
  gamks = ks.test(jw, "pgamma", gfit[1], gfit[2], alternative="two.sided")
  # logNormal
  lnks = ks.test(jw, "plnorm", lnfit[1], lnfit[2], alternative="two.sided")
  # weibull
  weibks = ks.test(jw, "pweibull", wfit[1], wfit[2], alternative="two.sided")
  # pareto
  parks = pareto.test(jw, 2e2)
  
  #Earth Movers Distance
  rpar = rpareto(lineSeq, parfit[[1]], parfit[[2]])
  rpar = rpar[rpar < max(w)] #filter out crazy big random values
  
  rg = rgamma(lineSeq, gfit[1], gfit[2])
  rln = rlnorm(lineSeq, lnfit[1], lnfit[2])
  rw = rweibull(lineSeq, wfit[1], wfit[2])
  
  #comparing to weibull
  emd.w = wasserstein1d(w, rw, p=1)
  #comparing to ln
  emd.ln = wasserstein1d(w, rln, p=1)
  #gamma
  emd.g = wasserstein1d(w, rg, p=1)
  #pareto
  emd.p = wasserstein1d(w, rpar, p=1)
  
  p1 = c((parfit[1]), (gfit[1]), (lnfit[1]), (wfit[1]))
  results = data.frame("scale" = s,
                       "fit" = c("Pareto", "Gamma", "log-normal", "Weibull"),
                       "p1" = unlist(unname(p1)),
                       "p2" = unlist(unname(c(parfit[2], gfit[2], lnfit[2], wfit[2]))),
                       "adp" = c(adp.p, NA, adl.p, adw.p),
                       "adt" = c(adp.t, adg.t, adl.t, adw.t),
                       "ksp" = c(parks$p, gamks$p.value, lnks$p.value, weibks$p.value),
                       "ksd" = c(parks$D, gamks$statistic, lnks$statistic, weibks$statistic),
                       "emd" = c(emd.p, emd.g, emd.ln, emd.w)
                       )
  if(i == 1) allresults = results
  if(i > 1) allresults = bind_rows(allresults, results) 
  #return(p1)
}
allresults2 <- signif(allresults, 2)
write.csv(allresults, "stattestresults_6_23_22.csv")

names = c("bill", "pete", "sally")
ages = c(18, 26, 25)
data.frame("people" = names, "experience" = ages, "fits" = p1)

#testing fits for other data
w <- widths
w = filter(data, source == sources[3]) 

w = w$width_m
w <- w[w > 21]

minW <- min(w)
int = 0.1
dlnInt = 0.1
lineSeq = seq(0, ceiling(max(w)/dlnInt)*dlnInt, dlnInt)
breaks = seq(0, (max(w)+int), int)
h = hist(w, breaks, plot=F)


parfit = pareto.mle(w)#consider running using mode W, rn using min W
dpar = dpareto(lineSeq, parfit[[1]], parfit[[2]])
ppar = ppareto(lineSeq, parfit[[1]], parfit[[2]])

gfit = suppressWarnings(fitdistr(w, "gamma")$estimate)
dg = dgamma(lineSeq, gfit[1], gfit[2])
pg = pgamma(lineSeq, gfit[1], gfit[2])
# logNormal
lnfit = fitdistr(w, "log-normal")$estimate
dln = dlnorm(lineSeq, lnfit[1], lnfit[2])
pln = plnorm(lineSeq, lnfit[1], lnfit[2])
# weibull
wfit = suppressWarnings(fitdistr(w, "weibull", lower = c(0,0))$estimate)
dw = dweibull(lineSeq, wfit[1], wfit[2])
pw = pweibull(lineSeq, wfit[1], wfit[2])

#Anderson Darling test, many packages
#lognormal, goftest package
ad.test(w, "plnorm", estimated = TRUE)
w <- jitter(w)
ks.test(w, "plnorm")
#gamma, gofgamma package
test.AD(w, alpha = 0.1)
test.KS(w) #same answer as my code
#weibull, cmstatr
anderson_darling_weibull(x = w, alpha = 0.1)
#pareto
pareto2_test_ad(w, s = parfit$xm) #agop
gpdAd(w)

#Earth Movers Distance
rpar = rpareto(lineSeq, parfit[[1]], parfit[[2]])
rg = rgamma(lineSeq, gfit[1], gfit[2])

rln = rlnorm(lineSeq, lnfit[1], lnfit[2])
pln = plnorm(lineSeq, lnfit[1], lnfit[2])

rw = rweibull(lineSeq, wfit[1], wfit[2])
breaks = seq(0, (max(widths)+20), 20)
#comparing to weibull
hist(widths, breaks,  freq = F)
hist(rw, breaks, freq = F, )
plot(density(widths), log = "x")
lines(density(rw), col = "blue")
wasserstein1d(w, rw, p=1)
#comparing to ln
breaks = seq(0, (max(widths)+20), 20)
hist(widths, breaks,  freq = T)
breaks <- seq(0, (max(rln)+20), 20)
hist(pln, breaks, freq = F, )
plot(density(widths), log = "x")
lines(density(rln), col = "blue")
wasserstein1d(w, rln, p=1)
#gamma
plot(density(widths), log = "x")
lines(density(rg), col = "blue")
wasserstein1d(w, rg, p=1)
#pareto
plot(density(widths), log = "x")
rpar <- rpar[rpar < max(w)]
lines(density(rpar), col = "blue")
wasserstein1d(w, rpar, p=1)

#################################
#Splitting up distributions, but with AD test
######################################



#plot of min filterd width versus p-value
#giga function
Pvals <- function(w){
  n = unique(c(seq(1, 100, 1), seq(100, 1000, 10))) 
  
  #setting up progress bar
  n_iter = length(n)
  init = numeric(n_iter)
  end = numeric(n_iter)
  
  #initiating progress bar
  pb = txtProgressBar(min = 0,
                      max = n_iter,
                      style = 3,
                      width = n_iter, # Needed to avoid multiple printings
                      char = "=") 
  
  for(i in 1:length(n)){
    init[i] = Sys.time()
    #filtering out to certain size
    
    widths.gre = w[which(w > i)]
    widths.les = w[which(w < i)]
    
    
    #widths greater than
    #KS test for pareto
    #t.gre = pareto.test(widths.gre, B = 1e3)
    #KS test for log-normal
    # jg = jitter(w)
    # lnfit.gre = fitdistr(widths.gre, "log-normal")$estimate
    # lnks.gre = ks.test(jg, "plnorm", lnfit.gre[1], lnfit.gre[2], alternative="two.sided")
    # #KS test for gamma
    # gfit.gre = suppressWarnings(fitdistr(widths.gre, "gamma", lower = c(1,1))$estimate)
    # gamks.gre = ks.test(jg, "pgamma", gfit.gre[1], gfit.gre[2], alternative="two.sided")
    # #KS test for weibull
    # wfit.gre = suppressWarnings(fitdistr(widths.gre, "weibull", lower = c(1,1))$estimate)
    # weibks.gre = ks.test(jg, "pweibull", wfit.gre[1], wfit.gre[2], alternative="two.sided")
    # 
    parfit.gre = pareto.mle(widths.gre)
    adp.gre = pareto2_test_ad(widths.gre, s = parfit.gre$xm)
    adp.p.gre = adp.gre$p.value
    adp.t.gre = unname(adp.gre$statistic)
    
    #lognormal, goftest package
    adl.gre = ad.test(widths.gre, "plnorm", estimated = TRUE)
    adl.p.gre = adl.gre$p.value
    adl.t.gre = unname(adl.gre$statistic)
    #gamma, gofgamma package
    # adg.gre = test.AD(widths.gre)
    # #pvalue not produced from this package
    # adg.t.gre = adg.gre$T.value
    #weibull, cmstatr
    # adw.gre = anderson_darling_weibull(x = widths.gre)
    # adw.p.gre = adw.gre$osl
    # adw.t.gre = adw.gre$A
    
    #widths less than
    #KS test for pareto
    #t.les = pareto.test(widths.les, B = 1e3)
    #KS test for log-normal
    # jl = jitter(w)
    # lnfit.les = fitdistr(widths.les, "log-normal")$estimate
    # lnks.les = ks.test(jl, "plnorm", lnfit.les[1], lnfit.les[2], alternative="two.sided")
    # #KS test for gamma
    # gfit.les = suppressWarnings(fitdistr(widths.les, "gamma", lower = c(1,1))$estimate)
    # gamks.les = ks.test(jl, "pgamma", gfit.les[1], gfit.les[2], alternative="two.sided")
    # #KS test for weibull
    # wfit.les = suppressWarnings(fitdistr(widths.les, "weibull", lower = c(1,1))$estimate)
    # weibks.les = ks.test(jl, "pweibull", wfit.les[1], wfit.les[2], alternative="two.sided")
    # 
    parfit.les = pareto.mle(widths.les)
    adp.les = pareto2_test_ad(widths.les, s = parfit.les$xm)
    adp.p.les = adp.les$p.value
    adp.t.les = unname(adp.les$statistic)
    
    #lognormal, goftest package
    adl.les = ad.test(widths.les, "plnorm", estimated = TRUE)
    adl.p.les = adl.les$p.value
    adl.t.les = unname(adl.les$statistic)
    #gamma, gofgamma package
    # adg.les = test.AD(widths.les)
    # #pvalue not produced from this package
    # adg.t.les = adg.les$T.value
    #weibull, cmstatr
    # adw.les = anderson_darling_weibull(x = widths.les)
    # adw.p.les = adw.les$osl
    # adw.t.les = adw.les$A
    
    raw = data.frame("width" = n[i],
                     #widths greater than
                     # "Ppval.gre" = t.gre$p, "PD.gre" = t.gre$D, 
                     # "LNpval.gre" = lnks.gre$p.value, "LND.gre" = lnks.gre$statistic,
                     # "Gpval.gre" = gamks.gre$p.value, "GD.gre" = gamks.gre$statistic,
                     # "Wpval.gre" = weibks.gre$p.value, "WD.gre" = weibks.gre$p.value,
                     # #Widths less than
                     # "Ppval.les" = t.les$p, "PD.les" = t.les$D,
                     # "LNpval.les" = lnks.les$p.value, "LND.les" = lnks.les$statistic,
                     # "Gpval.les" = gamks.les$p.value, "GD.les" = gamks.les$statistic,
                     # "Wpval.les" = weibks.les$p.value, "WD.les" = weibks.les$statistic,
                     #number of widths
                     "numberof>" = length(widths.gre),
                     "numberof<" = length(widths.les),
                     #AD results
                     #pareto
                     "adp.p.gre" = adp.p.gre,
                     "adp.t.gre" = adp.t.gre,
                     "adp.p.les" = adp.p.les,
                     "adp.t.les" = adp.t.les,
                     #lognormal
                     "adl.p.gre" = adl.p.gre,
                     "adl.t.gre" = adl.t.gre,
                     "adl.p.les" = adl.p.les,
                     "adl.t.les" = adl.t.les
                     #gamma
                     # "adg.gre" = adg.gre$Decision,
                     # "adg.t.gre" = adg.t.gre,
                     # "adg.les" = adg.les$Decision,
                     # "adg.t.les" = adg.t.les,
                     #weibull
                     # "adw.p.gre" = adw.p.gre,
                     # "adw.t.gre" = adw.t.gre,
                     # "adw.p.les" = adw.p.les,
                     # "adw.t.les" = adw.t.les
                     )
    
    if(i == 1) alldat = raw
    if(i > 1) alldat = bind_rows(alldat, raw)
    
    #Progress bar stuff
    end[i] = Sys.time()
    setTxtProgressBar(pb, i)
    time = round(seconds_to_period(sum(end - init)), 0)
    
    # Estimated remaining time based on the
    # mean time that took to run the previous iterations
    est = n_iter * (mean(end[end != 0] - init[init != 0])) - time
    remainining = round(seconds_to_period(est), 0)
    
    cat(paste(" // Execution time:", time,
              " // Estimated time remaining:", remainining), "")
  }
  
  close(pb)
  
  return(alldat)
}

pv <- Pvals(widths)
write.csv(pv, "fittingDistr_ADresults.csv", row.names = FALSE)
pv <- read.csv("fittingDistr_ADresults.csv")


test = w[which(w > i)]
test = w[which(w < 100)]
w3 <- sample(w, 100, FALSE)
test <- rpar(100, )

test = w[w > 22]
parfit = pareto.mle(w)#consider running using mode W, rn using min W
pareto.test(jitter(test), 2e2) #modeW

pareto.mle(w)
pareto2_test_ad(w, s = 0.332)
adp.p.gre = adp.gre$p.value
adp.t.gre = unname(adp.gre$statistic)

#lognormal, goftest package
lnfit = fitdistr(w3, "log-normal")$estimate
dln = dlnorm(seq(1,100, 1), lnfit[1], lnfit[2])


ad.test(w3, "plnorm", estimated = TRUE)
rln = rlnorm(1000, 2, 1)
ad.test(rln, "plnorm", estimated = TRUE)
ks.test(w3, "plnorm",lnfit[1], lnfit[2], alternative="two.sided")

adl.p.gre = adl.gre$p.value
adl.t.gre = unname(adl.gre$statistic)

#log-normal
dfLN.gre <- pv %>%
  dplyr::select(width, LNpval.gre, LND.gre, numberof.)%>%
  mutate(Value = derivedFactor(
    "p < 0.05" = LNpval.gre < 0.05,
    "p > 0.05" = LNpval.gre >= 0.05))%>%
  filter(numberof. > 53)

ln.gre <- ggplot(dfLN.gre)+
  geom_point(aes(x = width, y = LND.gre, color = Value), size = 2)+
  theme_classic()+
  labs(title = "Log-Normal fit K-S test results (Widths greater than n)",
       x = "Narrowest Width n (m)",
       y = "D")+
  lims(x = c(0, 100))+
  scale_color_manual(#values = c("p < 0.05", "p < 0.05"), 
    values = c("grey", "#7570b3"),
    name = "P-value")
#############################
#3. Plotting fits (if it is possible)
#3.1 continental scale
###########################

e4 <- read.csv("every4_measured.csv")
e41 <- na.omit(e4)

rivers2 <- filter(e41, lake == 0)

rivers <- read.csv("GE_widths.csv")
widths <- c(rivers$width_m, rivers2$width_m)

#widths <- rivers$width_m
w = widths
w[c(123, 323)] <- 1
w[c(888, 897, 1015, 1075, 1164, 1185)] <- 1

int = 2
breaks = seq(1, (max(w)+int), int)
counts = hist(w, breaks = breaks, plot = F)$counts

lnfit = fitdistr(w, "log-normal")$estimate
dln = dlnorm(breaks, lnfit[1], lnfit[2])
parfit = pareto.mle(w)
dpar = dpareto(breaks, parfit[[1]], parfit[[2]])
gfit = suppressWarnings(fitdistr(w, "gamma")$estimate)
dg = dgamma(breaks, gfit[1], gfit[2])
wfit = suppressWarnings(fitdistr(w, "weibull")$estimate)
dw = dweibull(breaks, wfit[1], wfit[2])


listx <- c()
listy <- c()
for(i in 1: length(breaks)){
  
  h = counts[i]
  coordsx = c(breaks[i], breaks[i], breaks[i+1], breaks[i+1])
  coordsy = c(0, h, h, 0)
  
  
  listx = append(listx, coordsx)
  listy = append(listy, coordsy)
}

coords <- data.frame(x = listx,
                     y = listy)

fits <- c("Pareto, Log-normal", "Weibull", "Gamma")
pdfs <- c(dpar, dln, dw, dg)

for(i in 1:length(fits)){
  f = pdfs[i]
  temp = data.frame(x = breaks,
                     y = dens <- (length(f)/sum(f)) * f)%>% 
  mutate(x2 = lag(x),
         y2 = lag(y), 
         fit = fits[i])
  if(i == 1) alldat = temp
  if(i > 1) alldat = bind_rows(alldat, temp)
}

points <- data.frame(x = breaks,
                     y = dens <- (length(dpar)/sum(dpar)) * dpar
)%>% 
  mutate(x2 = lag(x),
         y2 = lag(y),
         fit = "Pareto")
points.dl <- data.frame(x = breaks,
                     y = dens <- (length(dln)/sum(dln)) * dln
)%>% 
  mutate(x2 = lag(x),
         y2 = lag(y),
         fit = "Log-normal")
points.g <- data.frame(x = breaks,
                        y = dens <- (length(dg)/sum(dg)) * dg
)%>% 
  mutate(x2 = lag(x),
         y2 = lag(y),
         fit = "Gamma")
points.w <- data.frame(x = breaks,
                       y = dens <- (length(dw)/sum(dw)) * dw
)%>% 
  mutate(x2 = lag(x),
         y2 = lag(y),
         fit = "Weibull")

points3 <- bind_rows(points, points.dl, points.g, points.w)

ggplot(coords, aes(x = x, y = y)) +
  geom_polygon(fill = "lightgrey")+
  scale_x_log10(expand = expansion(mult = c(0, .1)),
                breaks = c(1, 10, 100, 1000),
                labels = c("1", "10", "100", "1,000"),
                limits = c(1, 2000)
                )+
  scale_y_log10(expand = expansion(mult = c(0, .1)),
                breaks = c(0.1, 1, 10, 100),
                labels = c("0.1" ,"1", "10", "100"),
                limits = c(1, 100)
  )+
  geom_segment(data = points3, aes(x = x, y = y,
                                   xend = x2, yend = y2,
                                   color = fit))+
  labs(x = "River Width (m)",
       y = "Frequency")+
  theme_classic()+
  scale_color_manual(name = "",
                     values = c("#1b9e77",
                                "#7570b3",
                                "#e7298a",
                                "#d95f02"),
                     labels = c("Gamma (EMD = 36.41)",
                                "Log-normal (EMD = 20.73)",
                                "Pareto (EMD = 25.29)",
                                "Weibull (EMD = 28.83)"))+
  theme(legend.position = "right")

##above is code to generate proper plot, below is code I am tweaking to fix first bin
cont <- ggplot(coords, aes(x = x, y = y)) +
  geom_polygon(fill = "lightgrey")+
  scale_x_log10(expand = expansion(mult = c(0, .1)),
                breaks = c(1, 10, 100, 1000),
                labels = c("1", "10", "100", "1,000"),
                limits = c(1, 1500)
  )+
  scale_y_log10(expand = expansion(mult = c(0, .1)),
                breaks = c(0.1, 1, 10, 100),
                labels = c("0.1" ,"1", "10", "100"),
                limits = c(0.1, 150)
  )+
  geom_segment(data = points3, aes(x = x, y = y,
                                   xend = x2, yend = y2,
                                   color = fit))+
  labs(x = "River Width (m)",
       y = "Frequency")+
  theme_classic()+
  scale_color_manual(name = "",
                     values = c("#1b9e77",
                                "#7570b3",
                                "#e7298a",
                                "#d95f02"),
                     labels = morewidthslabel)

morewidthslabel = c("Gamma (EMD = 30.86)",
                    "Log-normal (EMD = 16.22)",
                    "Pareto (EMD = 37.62)",
                    "Weibull (EMD = 23.89)")
oldlabels = c("Gamma",
              "Log-normal",
              "Pareto",
              "Weibull")
#trying to fix y axis limits
cont <- ggplot(coords, aes(x = x, y = y)) +
  geom_polygon(fill = "lightgrey")+
  scale_x_log10(expand = expansion(mult = c(0, .1)),
                breaks = c(1, 10, 100, 1000),
                labels = c("1", "10", "100", "1,000"),
                limits = c(1, 1500)
  )+
  scale_y_log10(expand = expansion(mult = c(0, .1)),
                breaks = c(0.5, 1, 10, 100),
                labels = c("0.5", "1", "10", "100"),
                limits = c(0.5, 150)
  )+
  geom_segment(data = points3, aes(x = x, y = y,
                                   xend = x2, yend = y2,
                                   color = fit))+
  labs(x = "River Width (m)",
       y = "Frequency",
       title = "Continental")+
  theme_classic()+
  scale_color_manual(name = "",
                     values = c("#1b9e77",
                                "#7570b3",
                                "#e7298a",
                                "#d95f02"),
                     labels = oldlabels)

################
#implementing for narrow widths
#3.2 headwater scale
###############
w2 = filter(data, source == "Konza: Headwaters")
#write.csv(w2, "headwaterwidths.csv", row.names = FALSE)

w = filter(data, source == "Konza: Headwaters")
w = w$width_m

int = 0.1
breaks = seq(0.01, (max(w)+ int), int)
counts = hist(w, breaks = breaks, plot = F)$counts
lnfit = fitdistr(w, "log-normal")$estimate
dln = dlnorm(breaks, lnfit[1], lnfit[2])

parfit = pareto.mle(w)

dpar = dpareto(breaks, parfit[[1]], parfit[[2]])

gfit = suppressWarnings(fitdistr(w, "gamma")$estimate)
dg = dgamma(breaks, gfit[1], gfit[2])
wfit = suppressWarnings(fitdistr(w, "weibull")$estimate)
dw = dweibull(breaks, wfit[1], wfit[2])


# points <- data.frame(x = breaks,
#                      dln = dln,
#                      dln_l = length(dln),
#                      dln_sum = sum(dln),
#                      counts = c(counts, 0),
#                      test = dens <- length(dln)/sum(dln) * dln,
#                      y = dens <- length(dln) * dln * int
# )

# points2 <- points %>% 
#   mutate(x2 = lag(x),
#          y2 = lag(y))
listx <- c()
listy <- c()
for(i in 1: length(breaks)){
  
  h = counts[i]
  coordsx = c(breaks[i], breaks[i], breaks[i+1], breaks[i+1])
  coordsy = c(0, h, h, 0)
  
  
  listx = append(listx, coordsx)
  listy = append(listy, coordsy)
}

coords <- data.frame(x = listx,
                     y = listy)

fits <- c("Pareto, Log-normal", "Weibull", "Gamma")
pdfs <- c(dpar, dln, dw, dg)

for(i in 1:length(fits)){
  f = pdfs[i]
  temp = data.frame(x = breaks,
                    y = dens <- length(f)/sum(f) * f)%>% 
    mutate(x2 = lag(x),
           y2 = lag(y), 
           fit = fits[i])
  if(i == 1) alldat = temp
  if(i > 1) alldat = bind_rows(alldat, temp)
}

points <- data.frame(x = breaks,
                     y = dens <- (length(dpar)/sum(dpar)) * dpar *20 
)%>% 
  mutate(x2 = lag(x),
         y2 = lag(y),
         fit = "Pareto")
points.dl <- data.frame(x = breaks,
                        y = dens <- (length(dln)/sum(dln)) * dln *20
)%>% 
  mutate(x2 = lag(x),
         y2 = lag(y),
         fit = "Log-normal")
points.g <- data.frame(x = breaks,
                       y = dens <- length(dg)/sum(dg) * dg*20
)%>% 
  mutate(x2 = lag(x),
         y2 = lag(y),
         fit = "Gamma")
points.w <- data.frame(x = breaks,
                       y = dens <- length(dw)/sum(dw) * dw*20
)%>% 
  mutate(x2 = lag(x),
         y2 = lag(y),
         fit = "Weibull")

coords2 <- coords[-1:-4,]

points3 <- bind_rows(points, points.dl, points.g, points.w)
head <- ggplot(coords2, aes(x = x, y = y)) +
  geom_polygon(fill = "lightgrey")+
  # scale_x_log10(expand = expansion(mult = c(0, .1)),
  #               breaks = c(0.01, 0.1, 1, 10),
  #               labels = c("0.01", "0.1", "1", "10"),
  #               limits = c(0.01, 10)
  # )+
  # scale_y_log10(expand = expansion(mult = c(0, .1)),
  #               breaks = c(0.01, 1, 100),
  #               ,
  #               limits = c(0.01, 100)
  # )+
  scale_x_log10(expand = expansion(mult = c(0, .1)),
                breaks = c(0.1, 1, 10),
                labels = c("0.1", "1", "10"),
                limits = c(0.1, 10))+
  scale_y_log10(expand = expansion(mult = c(0, .1)),
                breaks = c(1, 10, 100),
                labels = c("1", "10", "100"),
                limits = c(1, 200))+
  geom_segment(data = points3, aes(x = x, y = y,
                                   xend = x2, yend = y2,
                                   color = fit))+
  scale_color_manual(name = "",
                     values = c("#1b9e77",
                                "#7570b3",
                                "#e7298a",
                                "#d95f02"),
                     labels = oldlabels)+
  theme_classic()+
  labs(x = "River Width (m)",
       y = "Frequency",
       title = "Headwater")


################
#implementing for medium widths
#3.3 local scale
###############
w = filter(data, source == "McDowell: Local")

#write.csv(w, "localwidths.csv", row.names = FALSE)
w = w$width_m


int = 0.5
breaks = seq(0.01, (max(w)+int), int)
counts = hist(w, breaks = breaks, plot = F)$counts

lnfit = fitdistr(w, "log-normal")$estimate
dln = dlnorm(breaks, lnfit[1], lnfit[2])
parfit = pareto.mle(w)
dpar = dpareto(breaks, parfit[[1]], parfit[[2]])
gfit = suppressWarnings(fitdistr(w, "gamma")$estimate)
dg = dgamma(breaks, gfit[1], gfit[2])
wfit = suppressWarnings(fitdistr(w, "weibull")$estimate)
dw = dweibull(breaks, wfit[1], wfit[2])


listx <- c()
listy <- c()
for(i in 1: length(breaks)){
  
  h = counts[i]
  coordsx = c(breaks[i], breaks[i], breaks[i+1], breaks[i+1])
  coordsy = c(0, h, h, 0)
  
  
  listx = append(listx, coordsx)
  listy = append(listy, coordsy)
}

coords <- data.frame(x = listx,
                     y = listy)

fits <- c("Pareto, Log-normal", "Weibull", "Gamma")
pdfs <- c(dpar, dln, dw, dg)

for(i in 1:length(fits)){
  f = pdfs[i]
  temp = data.frame(x = breaks,
                    y = dens <- length(f)/sum(f) * f)%>% 
    mutate(x2 = lag(x),
           y2 = lag(y), 
           fit = fits[i])
  if(i == 1) alldat = temp
  if(i > 1) alldat = bind_rows(alldat, temp)
}

points <- data.frame(x = breaks,
                     y = dens <- length(dpar)/sum(dpar) * dpar *10
)%>% 
  mutate(x2 = lag(x),
         y2 = lag(y),
         fit = "Pareto")
points.dl <- data.frame(x = breaks,
                        y = dens <- length(dln)/sum(dln) * dln *10
)%>% 
  mutate(x2 = lag(x),
         y2 = lag(y),
         fit = "Log-normal")
points.g <- data.frame(x = breaks,
                       y = dens <- length(dg)/sum(dg) * dg*10
)%>% 
  mutate(x2 = lag(x),
         y2 = lag(y),
         fit = "Gamma")
points.w <- data.frame(x = breaks,
                       y = dens <- length(dw)/sum(dw) * dw*10
)%>% 
  mutate(x2 = lag(x),
         y2 = lag(y),
         fit = "Weibull") %>% filter()

points3 <- bind_rows(points, points.dl, points.g, points.w)

ggplot(coords, aes(x = x, y = y)) +
  geom_polygon(fill = "lightgrey")+
  scale_x_log10(expand = expansion(mult = c(0, .1))
                ,
                breaks = c(0.01, 0.1, 1, 10),
                labels = c("0.01", "0.1", "1", "10"),
                limits = c(0.01, 200)
  )+
  scale_y_log10(expand = expansion(mult = c(0, .1))
                ,
                breaks = c(0.1, 1, 10, 100),
                labels = c("0.1" ,"1", "10", "100"),
                limits = c(1, 100)
  )+
  geom_segment(data = points3, aes(x = x, y = y,
                                   xend = x2, yend = y2,
                                   color = fit))+
  labs(x = "River Width (m)",
       y = "Frequency")+
  theme_classic()+
  scale_color_manual(name = "",
                     values = c("#1b9e77",
                                "#7570b3",
                                "#e7298a",
                                "#d95f02"),
                     labels = c("Gamma (EMD = 0.55)",
                                "Log-normal (EMD = 1.2)",
                                "Pareto (EMD = 2.8)",
                                "Weibull (EMD = 0.63)"))+
  theme(legend.position = "right")

##modified to remove first bin
coords2 <- coords[-1:-4,]
points4 <- filter(points3, x > 1, x2 > 1)

ggplot(coords2, aes(x = x, y = y)) +
  geom_polygon(fill = "lightgrey")+
  scale_x_log10(limits = c(0.51, 100))+
  scale_y_log10()

loc <- ggplot(coords2, aes(x = x, y = y)) +
  geom_polygon(fill = "lightgrey")+
   scale_x_log10(expand = expansion(mult = c(0, .1))
                ,
                breaks = c(0.1, 1, 10, 100),
                labels = c("0.1","1", "10", "100"),
                limits = c(0.51, 100)
  )+
  scale_y_log10(expand = expansion(mult = c(0, .1))
                ,
                breaks = c(0.1, 1, 10, 100),
                labels = c("0.1" ,"1", "10", "100"),
                limits = c(1, 100)
  )+
  geom_segment(data = points3, aes(x = x, y = y,
                                   xend = x2, yend = y2,
                                   color = fit))+
  labs(x = "River Width (m)",
       y = "Frequency",
       title = "Local")+
  theme_classic()+
  scale_color_manual(name = "",
                     values = c("#1b9e77",
                                "#7570b3",
                                "#e7298a",
                                "#d95f02"),
                     labels = oldlabels)+
  theme(legend.position = "right")




# logNormal
w <- widths
minW <- min(w)
int = 10
dlnInt = 0.1
lineSeq = seq(0, ceiling(max(w)/dlnInt)*dlnInt, dlnInt)
breaks = seq(0, (max(w)+int), int)
h = hist(w, breaks, plot=F)

lnfit = fitdistr(w, "log-normal")$estimate
dln = dlnorm(breaks, lnfit[1], lnfit[2])
pln = plnorm(lineSeq, lnfit[1], lnfit[2])
plot(breaks, dln, log = 'x')

k=1000
set.seed(104)
val = rlnorm(k)
dens = density(val, n=512)

# Convert to counts
dens$y = k/sum(dens$y) * dens$y
dens$y = k/sum(dens$y) * dens$y
dens <- length(dln)/sum(dln) * dln
plot(dens)

ggplot()

fitd <- data.frame(lineSeq, dpar, dw, dln, dg)
t <- gather(fitd, key = "fit", value = "value", -lineSeq) %>% 
  filter(value > 0)

t$value <- t$value * 10 #* length(lineSeq[dpar<maxY&dpar!=0])
d <- filter(t, value > 0)
d <- na.omit(d)

lineSeq = seq(0, ceiling(max(w)/dlnInt)*dlnInt, dlnInt)

rpar = rpareto(lineSeq, parfit[[1]], parfit[[2]])
dt <- data.frame(lineSeq, dpar, value = v)
v <- dpar * seq(0.1, ) * length(lineSeq[dpar<maxY&dpar!=0])
breaks = c(1,1^(1:8))

ms_widths + geom_line(data = d, aes(x = lineSeq, y = value, color = fit))
                

ggplot() +
  geom_density(ms, mapping = aes(width_m)) +
  theme_classic()+
  labs(y = "N measurements",
       x = "Width (m)") +
  # geom_line(data = d, aes(x = lineSeq, y = value, color = fit))+
  scale_x_log10(expand = expansion(mult = c(0, .1)),
                limits = c(1, 1500))+ 
  geom_line(data = t, aes(x = lineSeq, y = value, color = fit))
##############################
# Make combined plot, showing all distributions and fits
#############################

patch <- cont / loc / head 

patch + plot_annotation(tag_levels = 'A')

ggsave(
  filename = "allfits.png",
  plot = last_plot(),
  width = 10,
  height = 12,
  units = "in"
)

##########################
#Doing for pareto fit, widths greater than 21
##########################
#pareto greater than 21 m wide
#weibull less than 49 m wide
e4 <- read.csv("every4_measured.csv")
e41 <- na.omit(e4)

rivers2 <- filter(e41, lake == 0)

rivers <- read.csv("GE_widths.csv")
widths <- c(rivers$width_m, rivers2$width_m)

w <- widths

int = 10
dlnInt = 0.1


#pareto
pw <- jitter(w[w > 22])

lineSeq.p = seq(0, ceiling(max(pw)/dlnInt)*dlnInt, dlnInt)
breaks.p = seq(0, (max(pw)+int), 5)
h.p = hist(pw, breaks.p, plot=F)

maxX = max(lineSeq.p)
maxY =  max(c(h$density))#, dg))#removed gamma dist max
breaks = seq(0, (max(w)+int), 2)

parfit = pareto.mle(pw)#consider running using mode W, rn using min W
dpar = dpareto(lineSeq.p, parfit[[1]], parfit[[2]])
ppar = ppareto(lineSeq.p, parfit[[1]], parfit[[2]])
parks <- pareto.test(pw, 2e2) #modeW

#pareto plot
h = hist(pw, breaks.p, freq=F, plot=T,
         ylim=c(0, 0.04), 
         xlim=c(0,200), xlab="Width (m)", 
         main='Pareto fit (W > 21 m)', yaxt='n',
         xaxt = 'n',
         col="light gray", border=0, ylab="Density", las = 1)
axis(2, seq(0, 0.04, 0.01), las = 2)
axis(1, c(0, 50, 100, 150, 200))
lines(lineSeq.p[dpar<maxY&dpar!=0], dpar[dpar<maxY&dpar!=0], col="#e7298a", lw = 2)

w <- jitter(widths[widths > 22])

int = 5
breaks = seq(1, (max(w)+int), int)
counts = hist(w, breaks = breaks, plot = F)$counts

parfit = pareto.mle(pw)#consider running using mode W, rn using min W
dpar = dpareto(breaks, parfit[[1]], parfit[[2]])
ppar = ppareto(breaks, parfit[[1]], parfit[[2]])


listx <- c()
listy <- c()
for(i in 1: length(breaks)){
  
  h = counts[i]
  coordsx = c(breaks[i], breaks[i], breaks[i+1], breaks[i+1])
  coordsy = c(0, h, h, 0)
  
  
  listx = append(listx, coordsx)
  listy = append(listy, coordsy)
}

coords <- data.frame(x = listx,
                     y = listy)

fits <- c("Pareto")
pdfs <- c(dpar)

for(i in 1:length(fits)){
  f = pdfs[i]
  temp = data.frame(x = breaks,
                    y = dens <- length(f)/sum(f) * f)%>% 
    mutate(x2 = lag(x),
           y2 = lag(y), 
           fit = fits[i])
  if(i == 1) alldat = temp
  if(i > 1) alldat = bind_rows(alldat, temp)
}

points <- data.frame(x = breaks,
                     y = dens <- length(dpar)/sum(dpar) * dpar * 1.9
)%>% 
  mutate(x2 = lag(x),
         y2 = lag(y),
         fit = "Pareto")


points3 <- (points[-1:-6,])

ggplot(coords, aes(x = x, y = y)) +
  geom_polygon(fill = "lightgrey")+
  scale_x_log10(expand = expansion(mult = c(0, .1)),
                #breaks = c(1, 10, 100, 1000),
                #labels = c("1", "10", "100", "1,000"),
                limits = c(10, 2000)
  )+
  scale_y_log10(expand = expansion(mult = c(0, .1)),
                #breaks = c(0.1, 1, 10, 100),
                #labels = c("0.1" ,"1", "10", "100"),
                limits = c(0.5, 100)
  )+
  geom_segment(data = points3, aes(x = x, y = y,
                                   xend = x2, yend = y2,
                                   color = fit),
               size = 2)+
  labs(x = "River Width (m)",
       y = "Frequency")+
  theme_classic()+
  scale_color_manual(name = "",
                     values = c("#e7298a"),
                     labels = c("Pareto fit"))+
  theme(legend.position = "right")



#Everything above works, doing with expanded dataset
rivers <- read.csv("GE_widths.csv")
widths <- rivers$width_m
w <- widths

int = 10
dlnInt = 0.1


#pareto
pw <- jitter(w[w > 22])

lineSeq.p = seq(0, ceiling(max(pw)/dlnInt)*dlnInt, dlnInt)
breaks.p = seq(0, (max(pw)+int), 5)
h.p = hist(pw, breaks.p, plot=F)

maxX = max(lineSeq.p)
maxY =  max(c(h.p$density))#, dg))#removed gamma dist max
breaks = seq(0, (max(w)+int), 2)

parfit = pareto.mle(pw)#consider running using mode W, rn using min W
dpar = dpareto(lineSeq.p, parfit[[1]], parfit[[2]])
ppar = ppareto(lineSeq.p, parfit[[1]], parfit[[2]])
parks <- pareto.test(pw, 2e2) #modeW

#pareto plot
h = hist(pw, breaks.p, freq=F, plot=T,
         ylim=c(0, 0.04), 
         xlim=c(0,200), xlab="Width (m)", 
         main='Pareto fit (W > 21 m)', yaxt='n',
         xaxt = 'n',
         col="light gray", border=0, ylab="Density", las = 1)
axis(2, seq(0, 0.04, 0.01), las = 2)
axis(1, c(0, 50, 100, 150, 200))
lines(lineSeq.p[dpar<maxY&dpar!=0], dpar[dpar<maxY&dpar!=0], col="#e7298a", lw = 2)

w <- jitter(widths[widths > 22])

int = 5
breaks = seq(1, (max(w)+int), int)
counts = hist(w, breaks = breaks, plot = F)$counts

parfit = pareto.mle(pw)#consider running using mode W, rn using min W
dpar = dpareto(breaks, parfit[[1]], 0.96)
ppar = ppareto(breaks, parfit[[1]], 0.96)


listx <- c()
listy <- c()
for(i in 1: length(breaks)){
  
  h = counts[i]
  coordsx = c(breaks[i], breaks[i], breaks[i+1], breaks[i+1])
  coordsy = c(0.1, h, h, 0.1)
  
  
  listx = append(listx, coordsx)
  listy = append(listy, coordsy)
}

coords <- data.frame(x = listx,
                     y = listy)

fits <- c("Pareto")
pdfs <- c(dpar)

for(i in 1:length(fits)){
  f = pdfs[i]
  temp = data.frame(x = breaks,
                    y = dens <- (length(f)/sum(f)) * f)%>% 
    mutate(x2 = lag(x),
           y2 = lag(y), 
           fit = fits[i])
  if(i == 1) alldat = temp
  if(i > 1) alldat = bind_rows(alldat, temp)
}

points <- data.frame(x = breaks,
                     y = dens <- (length(dpar)/sum(dpar)) * dpar 
)%>% 
  mutate(x2 = lag(x),
         y2 = lag(y),
         fit = "Pareto")


points3 <- (points[-1:-6,])

ggplot(coords, aes(x = x, y = y)) +
  geom_polygon(fill = "lightgrey")+
  scale_x_log10(expand = expansion(mult = c(0, .1)),
                breaks = c(10, 22, 100, 1000),
                labels = c("10", "22", "100", "1,000"),
                limits = c(10, 2000)
  )+
  scale_y_log10(expand = expansion(mult = c(0, .1)),
                breaks = c(0.1, 1, 10, 100),
                labels = c("0.1" ,"1", "10", "100"),
                limits = c(0.1, 100)
  )+
  geom_segment(data = points3, aes(x = x, y = y,
                                   xend = x2, yend = y2,
                                   color = fit),
               size = 1)+
  
  labs(x = "River Width (m)",
       y = "Frequency")+
  theme_classic()+
  scale_color_manual(name = "",
                     values = c("#e7298a"),
                     labels = c("This Study"))+
  theme(legend.position = "right")

##################################################
#4. Validation of measurements
##################################################
#using metrics from Ryan's RODEO paper, and George's NARWidth paper
#read in data (from section 12 of final1ksample.R)
jk <- read.csv("edited_GRWL_validation.csv")
jk <- filter(jk, width_m >= 90)
gs <- read.csv("USGS_insitu_validation.csv")
gs <- na.omit(gs)
fi <- read.csv("field_validation.csv")

sd(jk$diff)
sd(gs$width_diff_m)
sd(fi$diff)

SE <- function(x) sd(x) / sqrt(length(x))

SE(jk$diff)
SE(gs$width_diff_m)
SE(fi$diff)
#calculate Bias
bias <- function(xih, xi){
  if(length(xih) != length(xi)) {
    print("lengths differ")
    }
  else {
    sum(xih - xi)/length(xih)
    }
}

bias(fi$GE_width, fi$field_width)
bias(gs$GE_width, gs$USGS_width)
bias(jk$GE_width, jk$width_m)

#calculate RMSE
rootmse <- function(xih, xi){
  if(length(xih) != length(xi)) {
    print("lengths differ")
  }
  else {
    sqrt(sum( ((xi - xih)^2) /length(xih)))
  }
}

rootmse(fi$GE_width, fi$field_width)
rootmse(gs$GE_width, gs$USGS_width)
rootmse(jk$GE_width, jk$width_m)

#calculate Mean Absolute Error (MAE)
mae <- function(xih, xi){
  if(length(xih) != length(xi)) {
    print("lengths differ")
  }
  else {
    sum(abs(xi - xih))/length(xih)
  }
}

mae(fi$GE_width, fi$field_width)
mae(gs$GE_width, gs$USGS_width)
mae(jk$GE_width, jk$width_m)


g <- data.frame(m = measw, r = refw)
g <- na.omit(g)

bias(g$m, g$r)
rootmse(g$m, g$r)
mae(g$m, g$r)

########################
#5. redoing all analysis with bigger sample of widths at continental scale
######################
e4 <- read.csv("every4_measured.csv")
e41 <- na.omit(e4)

rivers2 <- filter(e41, lake == 0)
rivers <- read.csv("GE_widths.csv")

widths <- c(rivers$width_m, rivers2$width_m)

#simple histogram
hist(widths, breaks = seq(0, max(widths)+100, 100))
hist(log10(widths), breaks = 100)

# row.names(UAA) <- 1:nrow(UAA)
##########################
#6. residuals of USGS validation dataset dates
###########################
gs <- read.csv("USGS_insitu_validation.csv")

d <- as.numeric(difftime(gs$GE_date, gs$closest_USGS_date, units = "days"))

hist(d, xlim = c(-100, 400), breaks = 100, xlab = "Difference in Days", main = "")

###############
#calculating grwl estimate

#######################
# Revision of boxplots dataframe
#######################
#fixing name of McDowell data
data <- read.csv("HortonRatios/formatted_allSources.csv")
data2 <- filter(data, source == "McDowell: Regional") %>% 
  mutate(source = "McDowell: Local")
data3 <- filter(data, source != "McDowell: Regional")
data4 <- rbind(data2, data3)

#write.csv(data4, "HortonRatios/formatted_allSources.csv", row.names = FALSE)

#fixing stream order mismatch
data <- read.csv("HortonRatios/formatted_allSources.csv")
data2 <- filter(data, order == 13) %>% 
  mutate(order = 12)
data3 <- filter(data, order != 13)
data4 <- rbind(data2, data3)
# write.csv(data4, "HortonRatios/formatted_allSources.csv", row.names = FALSE)


#######################
# AD and KS test for revisions
#######################
#Adapted from George's small stream code
#testing for 

w <- widths
medFoW = median(w, na.rm=T)
modeFoW = density(w, bw=10, kernel="gaussian", na.rm=T)
modeFoW = modeFoW$x[which.max(modeFoW$y)]
minW = modeFoW

w <- widths
minW <- min(w)
int = 10
dlnInt = 0.1
lineSeq = seq(0, ceiling(max(w)/dlnInt)*dlnInt, dlnInt)
breaks = seq(0, (max(w)+int), int)
h = hist(w, breaks, plot=F)


parfit = pareto.mle(w)#consider running using mode W, rn using min W
dpar = dpareto(lineSeq, parfit[[1]], parfit[[2]])
ppar = ppareto(lineSeq, parfit[[1]], parfit[[2]])

gfit = suppressWarnings(fitdistr(w, "gamma")$estimate)
dg = dgamma(lineSeq, gfit[1], gfit[2])
pg = pgamma(lineSeq, gfit[1], gfit[2])
# logNormal
lnfit = fitdistr(w, "log-normal")$estimate
dln = dlnorm(lineSeq, lnfit[1], lnfit[2])
pln = plnorm(lineSeq, lnfit[1], lnfit[2])
# weibull
wfit = suppressWarnings(fitdistr(w, "weibull", lower = c(0,0))$estimate)
dw = dweibull(lineSeq, wfit[1], wfit[2])
pw = pweibull(lineSeq, wfit[1], wfit[2])

# Two sided One sample KS GOF test:
jw = jitter(widths) # to remove ties
# gamma
gamks = ks.test(jw, "pgamma", gfit[1], gfit[2], alternative="two.sided")
# logNormal
lnks = ks.test(jw, "plnorm", lnfit[1], lnfit[2], alternative="two.sided")
# weibull
weibks = ks.test(jw, "pweibull", wfit[1], wfit[2], alternative="two.sided")
# pareto
parks = pareto.test(jw, 2e2) #modeW

#Anderson Darling test, many packages
#lognormal, goftest package
ad.test(w, "plnorm", estimated = TRUE)

x <- rnorm(1000, mean=2, sd=1)
plot(density(x))
ad.test(x, "pnorm", mean=2, sd=1)

rln = rlnorm(200, lnfit[1], lnfit[2])
plot(density(w), log = 'x')
plot(rln, log = 'xy')
points(dln)

# Line-shaded polygons
plot(c(1, 9), c(1,9), type = "n")
polygon(x = c(1,1,4,4), y = c(1,4,4,1))
polygon(1:9, c(2,1,2,1,NA,2,1,2,1),
        density = c(10, 20), angle = c(-45, 45))


#gamma, gofgamma package
test.AD(w)
test.KS(w) #same answer as my code
#weibull, cmstatr
anderson_darling_weibull(x = w)
#pareto
pareto2_test_ad(w, s = 0.8906805) #agop
gpdAd(w)

# add histogram: 
maxX = max(lineSeq)
maxY =  max(c(h$density, dln))#, dg))#removed gamma dist max
breaks = seq(0, (max(w)+int), 2)

#log width
l <- log(w)
lb <- log(breaks)

###Final Plot
# png(file="differentFits2.png",
#     width=8.75, height=5, units = "in", res = 500)
h = hist(w, breaks = breaks, freq=F, plot=T,
         #ylim=c(0, maxY), 
         xlim=c(1,200), 
         xlab="River Width (m)", 
         main='', yaxt='n', 
         col="light gray", border=0, ylab="Density", las = 1)
axis(2, seq(0, 0.03, 0.01), las = 2)
lines(lineSeq[dpar<maxY&dpar!=0], dpar[dpar<maxY&dpar!=0], col="#e7298a", lw = 2)
lines(lineSeq, dw, col="#d95f02", lw = 2)
lines(lineSeq, dln, col="#7570b3", lw = 2)
lines(lineSeq, dg, col="#1b9e77", lw = 2)
legend(150, 0.02, legend=c("gamma",
                           "weibull",
                           "log-normal",
                           "pareto"),
       col = c("#1b9e77",
               "#d95f02",
               "#7570b3",
               "#e7298a"), 
       lty=1, lwd = 2, cex=0.8)

#dev.off()

