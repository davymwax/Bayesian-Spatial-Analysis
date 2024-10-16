## ----"Setting global chunk options", echo=FALSE----
rm(list = ls())
knitr::opts_chunk$set(echo = FALSE, 
                      message = FALSE, 
                      warning = FALSE, 
                      fig.align = "center", 
                      out.extra = "")


## ----"Loading the packages", include=FALSE--------
library("MASS")
library("maps")
library("fields")
library("sf")
library("akima")
library("spBayes")
library("RColorBrewer")
library("classInt")
library("gstat")
library("sp")
library("tidyverse")
library("tidyr")
library("ggplot2")
library("geoR")
library("kableExtra")
library("nlme")
library("lme4")
library("CARBayes")
library("maptools")
library("spatstat")
library("SpatialEpi")
library("spatialreg")
library("spdep")
library("foreign")
library("splancs")


## ----"Defining the grid", cache=TRUE, include=FALSE----
x.grid <- seq(0, 10, length=30)
y.grid <- seq(0, 10, length=30)
n.grid <- length(x.grid)
coord.pts <- matrix(0, nrow = n.grid*n.grid,
                    ncol = 2)
coord.pts[,1] <- rep(x.grid, each=n.grid)
coord.pts[,2] <- rep(y.grid, n.grid)
locations <- plot(coord.pts[,1], coord.pts[,2], 
     xlab="X", ylab="Y",
     main="Locations where the spatial process \n is simulated at",
     type="p",
     pch=21,
     col="gray")


## ----"Distance Matrix"----------------------------
N <- dim(coord.pts)[1]
dist.mat <- rdist(coord.pts)


## ----"defining the variance-covariance matrices", include=FALSE----
exp.cov <- function(sigma2, phi, tau2, dist.pts){
  cov.value <- tau2 + (sigma2*(exp(-(dist.pts/phi))))
  return(cov.value)
}

Sigma.exp <- matrix(0, N, N)
sigma2.exp <- 3
phi.exp <- 2
tau2.exp <- 1

for (i in 1:N){
  for (j in 1:N){
    Sigma.exp[i, j] <- exp.cov(sigma2.exp, phi.exp, tau2.exp, dist.mat[i, j])
  }
}


## ----"simulating the processes"-------------------
set.seed(183123)

## design matrix
X <- matrix(0, nrow=N, ncol=2)
for (i in 1:N){
  for (j in 1:2){
    X[i, j] <- runif(1, 0, 1)
  }
}

## mu(s)
mu <- rep(0, N)
for (i in 1:N){
  mu[i] <- 4.5 + 3.8*X[i, 1] + 8*X[i, 2]
}

## simulation 2 of Gaussian process (exponential covariance)
field.exp <- mvrnorm(1, mu, Sigma.exp)

image.plot(x.grid,
           y.grid,
           matrix(field.exp, n.grid, n.grid),
           col=terrain.colors(100),
           xlab="X",
           ylab="Y",
           main="Simulated Spatial Process with \n Exponential Covariance Function")


## ----"setting1", include=FALSE, cache=TRUE--------
set.seed(21233)
y <- field.exp
X_1 <- X[,1]
X_2 <- X[,2]
n <- dim(X)[1]

beta.starting <- rep(0, 2)
sigma2.starting <- 0.2
phi.starting <- 0.4
tau2.starting <- 0.16

setting1.fit <- spLM(y ~ X_1 + X_2,
                         coords = coord.pts,
                         knots = c(6, 6, 0.1),
                         starting=list("beta"= beta.starting,
                                       "phi"= phi.starting,
                                       "sigma.sq"= sigma2.starting,
                                       "tau.sq"=tau2.starting),
                         tuning=list("phi"= 0.006,
                                     "sigma.sq"= 0.0035,
                                     "tau.sq"= 0.001),
                        priors=list("beta.Flat",
                                    "phi.Unif"= c(0.15, 3),
                                    "sigma.sq.IG"=c(2, 1),
                                    "tau.sq.IG"=c(2, 0.5)),
                        cov.model="exponential",
                        n.samples=30000,
                        verbose=TRUE, n.report = 5000)


## ----"fig1", fig.cap="Trace Plots of Mean Model Parameters Setting 1",cache=TRUE----
n.samples <- 30000
burn.in <- 0.4*n.samples
setting1fit.other.params <- spRecover(setting1.fit,
                                          start=burn.in,
                                          verbose = FALSE)
par(mai=rep(0.25, 4), mfrow = c(2, 2))
plot(setting1fit.other.params$p.beta.recover.samples[, 1:3])


## ----"fig2", fig.cap="Trace Plots of the Covariance Parameters Setting 1"----
par(mai=rep(0.25, 4), mfrow = c(2, 2))
plot(setting1.fit$p.theta.samples)


## ----"tab1"---------------------------------------
beta_coeff.q <- round(apply(setting1fit.other.params$p.beta.recover.samples,
                          2, quantile,c(0.025, 0.5, 0.975)), 3)
beta_coeff.mean <- round(apply(setting1fit.other.params$p.beta.recover.samples,2,mean),3)
beta_coeff <- rbind(beta_coeff.mean, beta_coeff.q)
colnames(beta_coeff) <-c("Intercept", "X_1",
                         "X_2")
rownames(beta_coeff) <-c("Mean", "2.5%", "Median", "97.5%")
beta_coeff %>%
  kable(format = "latex",
        align = "c",
        digits = 4,
        caption = "Bayesian Estimates of Mean Model Parameters Under Knot Setting 1",
        booktabs = TRUE) %>%
  kable_styling(latex_options = c("striped", "scale_down"))


## ----"tab2"---------------------------------------
theta.q <- round(apply(setting1fit.other.params$p.theta.samples,
                          2, quantile,c(0.025, 0.5, 0.975)), 3)
theta.mean <- round(apply(setting1fit.other.params$p.theta.samples, 2, mean),3)
theta.estimates <- rbind(theta.mean, theta.q)
rownames(theta.estimates) <-c("Mean", "2.5%", "Median", "97.5%")
theta.estimates %>%
  kable(format = "latex",
        align = "c",
        digits = 4,
        caption = "Bayesian Estimates of Covariance Function Parameters Under Knot Setting 1",
        booktabs = TRUE) %>%
  kable_styling(latex_options = c("striped", "scale_down"))


## ----"setting2", include=FALSE, cache=TRUE--------
set.seed(212)
setting2.fit <- spLM(y ~ X_1 + X_2,
                         coords = coord.pts,
                         knots = c(12, 12, 0.1),
                         starting=list("beta"= beta.starting,
                                       "phi"= phi.starting,
                                       "sigma.sq"= sigma2.starting,
                                       "tau.sq"=tau2.starting),
                         tuning=list("phi"= 0.006,
                                     "sigma.sq"= 0.0035,
                                     "tau.sq"= 0.001),
                        priors=list("beta.Flat",
                                    "phi.Unif"= c(0.15, 3),
                                    "sigma.sq.IG"=c(2, 1),
                                    "tau.sq.IG"=c(2, 0.5)),
                        cov.model="exponential",
                        n.samples=30000,
                        verbose=TRUE, n.report = 5000)


## ----"fig3", fig.cap="Trace Plots of Mean Model Parameters Setting 2",cache=TRUE----
n.samples <- 30000
burn.in <- 0.4*n.samples
setting2fit.other.params <- spRecover(setting2.fit,
                                          start=burn.in,
                                          verbose = FALSE)
par(mai=rep(0.25, 4), mfrow = c(2, 2))
plot(setting2fit.other.params$p.beta.recover.samples[, 1:3])


## ----"fig4", fig.cap="Trace Plots of the Covariance Parameters Setting 2"----
par(mai=rep(0.25, 4), mfrow = c(2, 2))
plot(setting2.fit$p.theta.samples)


## ----"tab3"---------------------------------------
beta_coeff.q.2 <- round(apply(setting2fit.other.params$p.beta.recover.samples,
                          2, quantile,c(0.025, 0.5, 0.975)), 3)
beta_coeff.mean.2 <- round(apply(setting2fit.other.params$p.beta.recover.samples,2,mean),3)
beta_coeff.2 <- rbind(beta_coeff.mean.2, beta_coeff.q.2)
colnames(beta_coeff.2) <-c("Intercept", "X_1",
                         "X_2")
rownames(beta_coeff.2) <-c("Mean", "2.5%", "Median", "97.5%")
beta_coeff.2 %>%
  kable(format = "latex",
        align = "c",
        digits = 4,
        caption = "Bayesian Estimates of Mean Model Parameters Under Knot Setting 2",
        booktabs = TRUE) %>%
  kable_styling(latex_options = c("striped", "scale_down"))


## ----"tab4"---------------------------------------
theta.q.2 <- round(apply(setting2fit.other.params$p.theta.samples,
                          2, quantile,c(0.025, 0.5, 0.975)), 3)
theta.mean.2 <- round(apply(setting2fit.other.params$p.theta.samples, 2, mean),3)
theta.estimates.2 <- rbind(theta.mean.2, theta.q.2)
rownames(theta.estimates.2) <-c("Mean", "2.5%", "Median", "97.5%")
theta.estimates.2 %>%
  kable(format = "latex",
        align = "c",
        digits = 4,
        caption = "Bayesian Estimates of Covariance Function Parameters Under Knot Setting 2",
        booktabs = TRUE) %>%
  kable_styling(latex_options = c("striped", "scale_down"))


## ----"324 knots", include=FALSE, cache=TRUE-------
knots.x.grid <- runif(18, 0, 10)
knots.y.grid <- runif(18, 0, 10)
knots.n.grid <- length(knots.x.grid)
knots.pts <- matrix(0, nrow = knots.n.grid*knots.n.grid,
                    ncol = 2)
knots.pts[,1] <- rep(knots.x.grid, each=knots.n.grid)
knots.pts[,2] <- rep(knots.y.grid, knots.n.grid)


## ----"setting3", include=FALSE, cache=TRUE--------
set.seed(1223)
setting3.fit <- spLM(y ~ X_1 + X_2,
                         coords = coord.pts,
                         knots = knots.pts,
                         starting=list("beta"= beta.starting,
                                       "phi"= phi.starting,
                                       "sigma.sq"= sigma2.starting,
                                       "tau.sq"= tau2.starting),
                         tuning=list("phi"= 0.006,
                                     "sigma.sq"= 0.0035,
                                     "tau.sq"= 0.001),
                        priors=list("beta.Flat",
                                    "phi.Unif"= c(0.15, 3),
                                    "sigma.sq.IG"=c(2, 1),
                                    "tau.sq.IG"=c(2, 0.5)),
                        cov.model="exponential",
                        n.samples=30000,
                        verbose=TRUE, n.report = 5000)


## ----"fig5", fig.cap="Trace Plots of Mean Model Parameters Setting 3",cache=TRUE----
n.samples <- 30000
burn.in <- 0.4*n.samples
setting3.fit.other.params <- spRecover(setting3.fit,
                                          start=burn.in,
                                          verbose = FALSE)
par(mai=rep(0.25, 4), mfrow = c(2, 2))
plot(setting3.fit.other.params$p.beta.recover.samples[, 1:3])


## ----"fig6", fig.cap="Trace Plots of the Covariance Parameters Setting 3"----
par(mai=rep(0.25, 4), mfrow = c(2, 2))
plot(setting3.fit$p.theta.samples)


## ----"tab5"---------------------------------------
beta_coeff.q.3 <- round(apply(setting3.fit.other.params$p.beta.recover.samples,
                          2, quantile,c(0.025, 0.5, 0.975)), 3)
beta_coeff.mean.3 <- round(apply(setting3.fit.other.params$p.beta.recover.samples,2,mean),3)
beta_coeff.3 <- rbind(beta_coeff.mean.3, beta_coeff.q.3)
colnames(beta_coeff.3) <-c("Intercept", "X_1",
                         "X_2")
rownames(beta_coeff.3) <-c("Mean", "2.5%", "Median", "97.5%")
beta_coeff.3 %>%
  kable(format = "latex",
        align = "c",
        digits = 4,
        caption = "Bayesian Estimates of Mean Model Parameters Under Knot Setting 3",
        booktabs = TRUE) %>%
  kable_styling(latex_options = c("striped", "scale_down"))


## ----"tab6"---------------------------------------
theta.q.3 <- round(apply(setting3.fit.other.params$p.theta.samples,
                          2, quantile,c(0.025, 0.5, 0.975)), 3)
theta.mean.3 <- round(apply(setting3.fit.other.params$p.theta.samples, 2, mean),3)
theta.estimates.3 <- rbind(theta.mean.3, theta.q.3)
rownames(theta.estimates.3) <-c("Mean", "2.5%", "Median", "97.5%")
theta.estimates.3 %>%
  kable(format = "latex",
        align = "c",
        digits = 4,
        caption = "Bayesian Estimates of Covariance Function Parameters Under Knot Setting 3",
        booktabs = TRUE) %>%
  kable_styling(latex_options = c("striped", "scale_down"))

