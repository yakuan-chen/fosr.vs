## ----, echo=FALSE, message=FALSE, warning=FALSE--------------------------

library(fosr.vs)
library(refund)
library(BayesFoSR)
library(ggplot2)
library(MASS)

###############################################################
## set simulation design elements
###############################################################

set.seed(100)

I = 100
p = 20
D = 50
grid = seq(0, 1, length = D)

## coefficient functions
beta.true = matrix(0, p, D)
beta.true[1,] = sin(2*grid*pi)
beta.true[2,] = cos(2*grid*pi)
beta.true[3,] = 2

## basis functions for correlated errors
psi.true = matrix(NA, 2, D)
psi.true[1,] = sin(4*grid*pi)
psi.true[2,] = cos(4*grid*pi)
lambda = c(3,1)

## predictors and PC scores
X = matrix(rnorm(I*p), I, p); X[,1] = 1
C = cbind(rnorm(I, mean = 0, sd = lambda[1]), rnorm(I, mean = 0, sd = lambda[2]))

## fixed effects, PC effects, uncorrelated errors
fixef = X%*%beta.true
pcaef = C %*% psi.true
error = matrix(rnorm(I*D), I, D)

Yi.true = fixef
Yi.pca = fixef + pcaef
Yi.obs = fixef + pcaef + error

data = as.data.frame(X)
data$Y = Yi.obs

## response
Yi.obs = fixef + pcaef + error

## ----, echo = FALSE, fig.align='center', fig.width=6---------------------
## observed data
obs.data = cbind(as.vector(t(Yi.obs)), rep(1:I, each = D), rep(grid, I))
obs.data = as.data.frame(obs.data)
colnames(obs.data) = c("y", "curve", "grid")

ggplot(obs.data, aes(x = grid, y = y, group = curve)) + geom_path(alpha = .3) +
  theme_bw() + labs(x = "Grid", y = "Observed Data")

## ----, echo = FALSE, results='hide', fig.align='center', fig.width=6-----
coefs.mcp = fit.fosr.vs$coefficients
coefs.vb = fit.vb$beta.pm
coefs.pffr = matrix(NA, nrow = p, ncol = D)
for(i in 1:p){
  coefs.pffr[i,] = coef(fit.pffr, n1 = D)$smterms[[paste0("x", i, "(yindex)")]]$value
}

plot.dat = data.frame(est = c(as.vector(t(beta.true[1:3,])), as.vector(t(coefs.mcp[1:3,])), as.vector(t(coefs.pffr[1:3,])), as.vector(t(coefs.vb[1:3,]))),
                      grid = rep(grid, 3*4),
                      curve = rep(rep(1:3, each = D), 4),
                      Method = rep(c("Truth", "VarSelect", "PFFR", "VB"), each = 3*D))

ggplot(plot.dat, aes(x = grid, y = est, group = Method, color = Method)) + geom_path(alpha = .5) + 
  theme_bw() + labs(x = "Grid", y = "Beta") + facet_grid(.~curve)

## ----, echo = FALSE, fig.align='center', fig.width=4---------------------
plot.dat = data.frame(est = c(as.vector(t(beta.true[4:20,])), as.vector(t(coefs.mcp[4:20,])), as.vector(t(coefs.pffr[4:20,])), as.vector(t(coefs.vb[4:20,]))),
                      grid = rep(grid, 17*4),
                      curve = rep(1:(17*4), each = D),
                      Method = rep(c("Truth", "VarSelect", "PFFR", "VB"), each = 17*D))

ggplot(plot.dat, aes(x = grid, y = est, group = curve, color = Method)) + geom_path(alpha = .5) + 
  theme_bw() + labs(x = "Grid", y = "Beta")

