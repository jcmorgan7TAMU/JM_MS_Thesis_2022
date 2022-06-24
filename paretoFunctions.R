###This script just contains support Pareto functions
#12/20/2021
#Adapted from George's small stream width analysis script, and my previous analysis

#example of loading a script that contains a function
#source(paste0(workingDir, '/EDfig2_EDtab3_GOF.R'))


#host of pareto functions
#ensure that there are no NAs in width dataset
# hist(as.numeric(ten5$width))
# wid <- as.numeric(ten5$width)
# wid2 <- wid[!is.na(wid)]

# calculate mode first order stream width for Pareto fit:
# find first order mean stream width:
# Don't actually use this I think
# fOw_raw = w_raw[tab$stream_order == 1]
# notZero = (fOw_raw!=0 & !is.na(fOw_raw))
# fOw = fOw_raw[notZero]#*2.54# inches to cm convert
# medFoW = median(fOw, na.rm=T)
# modeFoW = density(fOw, bw=10, kernel="gaussian", na.rm=T)
# modeFoW = modeFoW$x[which.max(modeFoW$y)]
# minW = modeFoW

# distribution, cdf, quantile and random functions for Pareto distributions
#pdf
dpareto <- function(x, xm, alpha) ifelse(x > xm , alpha*xm**alpha/(x**(alpha+1)), 0)
#cdf
ppareto <- function(q, xm, alpha) ifelse(q > xm , 1 - (xm/q)**alpha, 0 )
qpareto <- function(p, xm, alpha) ifelse(p < 0 | p > 1, NaN, xm*(1-p)**(-1/alpha))
rpareto <- function(n, xm, alpha) qpareto(runif(n), xm, alpha)


pareto.mle <- function(x)
{
  xm <- min(x)
  alpha <- length(x)/(sum(log(x))-length(x)*log(xm))
  return( list(xm = xm, alpha = alpha))
}


pareto.test <- function(x, B = 1e3)
{
  a <- pareto.mle(x)
  
  # KS statistic
  D <- ks.test(x, function(q) ppareto(q, a$xm, a$alpha))$statistic
  
  # estimating p value with parametric bootstrap
  #B <- 1e5
  n <- length(x)
  emp.D <- numeric(B)
  for(b in 1:B)
  {
    xx <- rpareto(n, a$xm, a$alpha);
    aa <- pareto.mle(xx)
    emp.D[b] <- ks.test(xx, function(q) ppareto(q, aa$xm, aa$alpha))$statistic
  }
  
  return(list(xm = a$xm, alpha = a$alpha, D = D, p = sum(emp.D > D)/B))
}


#plot
# curve(dpareto(x, 27.443,1.03), xlim = c(10, 1000), ylim = c(1e-04, 1e+01),
#       col = "black", xlab = c("Width (m)"), ylab = c("Density"), log = "xy", lwd = 2)#all, mode as min
# curve(dpareto(x, 90, 0.83), add = TRUE, col = "blue", lwd = 2, lty = 3)#GRWL, global
# curve(dpareto(x, 90, 0.96), add = TRUE, col = "lightblue", lwd = 2, lty = 3)#GRWL Mississippi
# curve(dpareto(x, 92.19, 1.18), add = TRUE, col = "salmon1", lwd = 2)#bigger than 90 m
# curve(dpareto(x, 0.61, 0.29), add = TRUE, col = "salmon3", lwd = 2)#all, true min
# legend(150, 10, legend=c("This study",
#                          "GRWL, global",
#                          "GRWL, Mississippi Basin",
#                          "This study, > 90m",
#                          "This study, < 90m"),
#        col = c("black",
#                "blue",
#                "lightblue",
#                "salmon1",
#                "salmon3"), 
#        lty=c(1, 3,3,1,1), lwd = 2, cex=0.8)


