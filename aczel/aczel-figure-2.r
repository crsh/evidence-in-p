library("dplyr")
library("baymedr")
library("BayesFactor")

n_jab <- function(n1, n2 = NULL) {
  if(is.null(n2)) n2 <- 0
  os <- n2 == 0
  n1 * (n2 + os) / ((n1 + n2) - (n1 + n2 - 1) * os)
}

### t.test
B01_paired = function(t.stat, N, log.prior.dens, lo = -Inf, up = Inf, ...){
  dt(t.stat, N - 1)/p.y.alt_p(t,N,log.prior.dens,lo,up,...)
}

B01_independent = function(t.stat, n1,n2, log.prior.dens, lo = -Inf, up = Inf, ...){
  dt(t.stat, n1+n2 - 2)/p.y.alt_i(t,n1,n2,log.prior.dens,lo,up,...)
}

## Bayes factor
### Source: https://gist.github.com/richarddmorey/3dbd911466389f7f263e
### for paired samples
p.y.alt_p = function(t.stat, N, log.prior.dens,lo=-Inf,up=Inf,...){
  normalize = integrate(function(delta,...){
    exp(log.prior.dens(delta, ...))
  }, lower = lo, upper = up, ...)[[1]]
  py = integrate(function(delta,t.stat,N,...){
    exp(
      dt(t, N - 1, ncp = delta*sqrt(N), log = TRUE) +
        log.prior.dens(delta, ...)  
    )
  },lower = lo, upper = up, t.stat = t.stat, N = N,stop.on.error = FALSE, ...)[[1]]
  py/normalize
}

### for independent samples / ess = effective smaple size
p.y.alt_i = function(t.stat, n1, n2, log.prior.dens,lo=-Inf,up=Inf,...){
  normalize = integrate(function(delta,...){
    exp(log.prior.dens(delta, ...))
  }, lower = lo, upper = up, ...)[[1]]
  py = integrate(function(delta,t.stat,n1,n2,...){
    ess = n1*n2 / (n1+n2)
    exp(
      dt(t, n1+n2 - 2, ncp = delta*sqrt(ess), log = TRUE) +
        log.prior.dens(delta, ...)  
    )
  },lower = lo, upper = up, t.stat = t.stat, n1=n1, n2=n2,stop.on.error = FALSE, ...)[[1]]
  py/normalize
}

## Prior of the Bayes factor
### Normal and Cauchy prior functions 
# normal.prior = function(delta, mu=0, sd=1){
#   dnorm(delta, mu, sd, log = TRUE)
# }

# cauchy.prior = function(delta, rscale){
#   dcauchy(delta, scale = rscale, log = TRUE)
# }

sample_nonsig <- read.csv("./aczel/data.csv", row.names = 1)

sample_nonsig <- 
  sample_nonsig |> mutate(
    Reference=as.character(Reference),
    Link=as.character(Link),
    t=as.numeric(as.character(t)),
    Associated.statistics=as.character(Associated.statistics),
    p=as.numeric(as.character(p)),
    Negative.statement=as.character(Negative.statement))

nonsigB<- subset(sample_nonsig, is.na(sample_nonsig$t)!=TRUE&is.na(sample_nonsig$N1)!=TRUE)
#### Clearing target variables
nonsigB$BF01default <- NULL
# nonsigB$BF01informed <- NULL
# nonsigB$BF01normal <- NULL

### Transforming NAs to 0
nonsigB$N2 <- ifelse(is.na(nonsigB$N2)==TRUE,0,nonsigB$N2)
#### Sample size
# nonsigB$N <- nonsigB$N1 + nonsigB$N2
#### Design: 0 = paired, 1 = independent
# nonsigB$Design <- ifelse(nonsigB$N2==0, 0,1)

nonsigB$p.calculated <- 2*pt(abs(nonsigB$t), ifelse(nonsigB$N2==0, nonsigB$N1-1, (nonsigB$N1+nonsigB$N2 -2)), lower = FALSE)
nonsigB$OneTailed <- ifelse(0.01+nonsigB$p.calculated/2<nonsigB$p, 0,1)

### Calculating Bayes factor
for (i in 1:nrow(nonsigB)){
  nonsigB[i,"BF01default"] <-1/exp(ttest.tstat(t = nonsigB[i, "t"],
                                      n1 = nonsigB[i, "N1"],
                                      n2 = nonsigB[i, "N2"],
                                      nullInterval = NULL,
                                      rscale = "medium",
                                      complement = FALSE,
                                      simple = FALSE)[['bf']])
    
### Bayes factor with normal prior
# t = nonsigB[i, "t"]
# N1 = nonsigB[i, "N1"]
# N2 = nonsigB[i, "N2"]
# N = nonsigB[i, "N1"]
  
# ifelse(test = nonsigB[i, "Design"]==0,
#         nonsigB[i, "BF01normal"] <- B01_paired(t, N, normal.prior,-Inf,Inf,mu=0,sd=0.5),
#         nonsigB[i, "BF01normal"] <- B01_independent(t, N1,N2, normal.prior,-Inf,Inf,mu=0,sd=0.5))
  }

### Bayes factor with informed prior
#### Informed prior
# prior.location1 <- - 0.34999
# prior.location2 <- 0.34999
# prior.scale <- 0.1021
# prior.df <- 3

# for (i in seq_len(nrow(nonsigB))) {
  
#   print(i)
  
#   if (!is.na(nonsigB$t[i]) && !is.na(nonsigB$N1[i]) &&
#       nonsigB$Type.of.statistics[i] %in% c("One Sample T-Test", "Paired Samples T-Test")) {
    
#     BF10_1 <- baymedr:::bf10_t(nonsigB$t[i], n1 = nonsigB$N1[i],
#                      prior_loc = prior.location1, prior_scale = prior.scale,
#                      prior_df = prior.df)$bf_10
#     BF10_2 <- baymedr:::bf10_t(nonsigB$t[i], n1 = nonsigB$N1[i],
#                      prior_loc = prior.location2, prior_scale = prior.scale,
#                      prior_df = prior.df)$bf_10
#     nonsigB$BF01informed[i] <- 1/mean(c(BF10_1, BF10_2))
    
#   } else if (!is.na(nonsigB$t[i]) && !is.na(nonsigB$N1[i]) && !is.na(nonsigB$N2[i]) &&
#              nonsigB$Type.of.statistics[i] == "Independent Samples T-Test") {
    
#     BF10_1 <- baymedr:::bf10_t(nonsigB$t[i], n1 = nonsigB$N1[i], n2 = nonsigB$N2[i], ind_samples = TRUE,
#                      prior_loc = prior.location1, prior_scale = prior.scale,
#                      prior_df = prior.df)$bf_10
#     BF10_2 <- baymedr:::bf10_t(nonsigB$t[i], n1 = nonsigB$N1[i], n2 = nonsigB$N2[i], ind_samples = TRUE,
#                      prior_loc = prior.location2, prior_scale = prior.scale,
#                      prior_df = prior.df)$bf_10
#     nonsigB$BF01informed[i] <- 1/mean(c(BF10_1, BF10_2))
    
#   }}

### Producing log Bayes factors with reciprocal values
#### ln Bayes factors
nonsigB$logBF01default <- log(nonsigB$BF01default)
# nonsigB$logBF01normal <- log(nonsigB$BF01normal)
# nonsigB$logBF01informed <- log(nonsigB$BF01informed)

# nonsigB$logN <- log(nonsigB$N)

par(cex.main = 1.3, mar = c(4.5, 6, 4, 7) + 0.1, mgp = c(3, 1, 0), cex.lab = 1.3, 
    font.lab = 2, cex.axis = 1.3, las = 1)

plot(
  x= log(jab::jab_p(nonsigB[,'p'] * (2 - !nonsigB[,'OneTailed']), n_jab(nonsigB$N1, nonsigB$N2))), y=nonsigB[,'logBF01default']
  # x= nonsigB[,'p'], y=nonsigB[,'logBF01default'] # Original plot
  # , xlim = c(0, 1)
  # , xlim = c(0, 1)
  , ylim = c(-1 * log(4), log(30)), 
     xlab = "", ylab = "", cex.lab = 1.3, cex.axis = 1.3, las = 1, yaxt = "n",xaxt="n", cex = 1, 
     bty = "n", type = "p", pch = 21, bg = "grey")


labelsUpper = log(c(30, 10, 3, 1))#bebebe
lablesUppercorx = log(c(10, 3, 1))
labelsLower = -1 * log(c(3, 1))
criticalP = c(labelsLower, 0, labelsUpper)
criticalPcorx = c(lablesUppercorx,0,labelsLower)
abline(h = 0)
abline(a = 0, b = 1)
axis(side = 4, at = labelsUpper + 0.602, tick = FALSE, cex.axis = 1, labels = c ( "", "Strong H0", "Moderate H0", "Anecdotal H0"))
axis(side = 4, at = labelsLower - 0.602, tick = FALSE, cex.axis = 1, labels = c( "", "Anecdotal H1"))
axis(side = 2, at = c(criticalP), tick = TRUE, las = 2, cex.axis = 1, labels = c( "1/3", "1", "",  "30", "10", "3", ""))
axis(side=1, cex.axis= 0.7)
# axis(side=1, at= c(0,0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,1), cex.axis= 0.7) # Original plot
mtext("Bayes Factors in Favor of the Null", side = 2, line = 2.5, las = 0, cex = 1)
mtext("log JAB", side = 1, line = 2.5, las = 1, cex = 1)
# mtext("P-values", side = 1, line = 2.5, las = 1, cex = 1)
grid::grid.text("", 0.97, 0.5, rot = 270, gp = grid::gpar(cex = 1.3))
for (idx in 1:length(criticalPcorx)) {
  abline(h = criticalPcorx[idx], col = "darkgrey", lwd = 1, lty = 2)
}
