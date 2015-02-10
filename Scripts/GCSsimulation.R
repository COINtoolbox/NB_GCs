# Simulation study for count data ##############################################
#
# Simulate data from 3 distributions: gaussian, poisson and neg.binomial.
# Fit models to each of these 3 assuming again gaussian, 
# poisson and neg. binomial.
# Assess model fit for all 9 combinations. 
#
# Bart Buelens, 11 Dec 2014.

library(ggplot2)
library(MASS)
library(ggthemes)

### Generate 3 data sets #######################################################
n = 200
x = runif(n, min=2, max=5)
yLi = x + rnorm(n, mean=0, sd=0.5)
summary(yLi)
yLi[yLi<0] = 0
summary(exp(yLi))
plot(x, yLi)
plot(x, exp(yLi))

yPo = rpois(n, lambda=20*x)
summary(yPo)
plot(x, yPo)
plot(x, log(yPo))

yNb = rnbinom(n, size=3, mu=20*x)
summary(yNb)
plot(x, yNb)
plot(x, log(yNb))

D = data.frame(x = rep(x,3), 
               count = c(exp(yLi), yPo, yNb), 
               logcount = c(yLi, log(yPo), log(yNb)),
               key = rep(1:n,3))
D$dgm = factor(rep(c("true.linm","true.pois","true.negb"), each=n),
               levels=c("true.linm","true.pois","true.negb"),
               ordered=TRUE)
D$logcount[D$logcount==-Inf] = 0

ggplot(D, aes(x=x, y=count)) + 
   geom_point(aes(color=dgm)) + 
   theme(legend.position="none") +
    facet_wrap(~ dgm)

ggplot(D, aes(x=x, y=logcount)) + 
   geom_point(aes(color=dgm)) + 
   theme(legend.position="none") +
  facet_wrap(~ dgm)

### Code for Prediction interval estimation ####################################
predInterval = function(Data, bN = 100) {
   N = nrow(Data)
   xdf = data.frame(x = Data$x)
   predN = length(Data$x)
   Blinm = Bpois = Bnegb = matrix(data = NA, nrow = bN, ncol = predN)
   BlinmU = BlinmL = BpoisU = BpoisL = BnegbU = BnegbL = Blinm
   qlwr = 0.025
   qupr = 0.975
   for (i in 1:bN) {
      G = Data[sample(1:N, N, replace = TRUE),]
      linmFit = glm(logcount ~ x, G, family="gaussian")
      poisFit = glm(count ~ x, G, family="poisson")
      negbFit = glm.nb(count ~ x, G)
      linmP = predict(linmFit, xdf)
      poisP = predict(poisFit, xdf)
      negbP = predict(negbFit, xdf)
      Blinm[i,] = linmP
      BlinmL[i,] = qnorm(qlwr, mean = linmP, sd = sd(residuals(linmFit)))
      BlinmU[i,] = qnorm(qupr, mean = linmP, sd = sd(residuals(linmFit)))
      Bpois[i,] = poisP
      BpoisL[i,] = log(qpois(qlwr, lambda = exp(poisP)))
      BpoisU[i,] = log(qpois(qupr, lambda = exp(poisP)))
      Bnegb[i,] = negbP
      BnegbL[i,] = log(qnbinom(qlwr, size = negbFit$theta, mu = exp(negbP)))
      BnegbU[i,] = log(qnbinom(qupr, size = negbFit$theta, mu = exp(negbP)))
   } 
   returnDataF = data.frame(pred = c(colMeans(Blinm), colMeans(Bpois), 
                                     colMeans(Bnegb)),
                         lwr = c(colMeans(BlinmL), colMeans(BpoisL), 
                                 colMeans(BnegbL)),
                         upr = c(colMeans(BlinmU), colMeans(BpoisU), 
                                 colMeans(BnegbU)),
                         key = rep(Data$key, 3))
   returnDataF$fit = factor(rep(c("fit.linm","fit.pois","fit.negb"), 
                                each=nrow(Data)),
                            levels=c("fit.linm", "fit.pois", "fit.negb"),
                            ordered=TRUE)  
   return(returnDataF)
}

### Fit each of the three with all 3 ###########################################
D = rbind(D, D, D)
D$fit = factor(rep(c("fit.linm", "fit.pois", "fit.negb"), each = 3 * n),
               levels=c("fit.linm", "fit.pois", "fit.negb"),
               ordered=TRUE)
D$yhat = NA
trueModels = c("true.linm", "true.pois", "true.negb")
for (thisDgm in trueModels) {
   D[D$fit=="fit.linm" & D$dgm==thisDgm, "yhat"] = 
      lm(logcount ~ x, subset(D, dgm==thisDgm & fit=="fit.linm"))$fitted.values
   
   D[D$fit=="fit.pois" & D$dgm==thisDgm, "yhat"] = 
      log(glm(count ~ x, subset(D, dgm==thisDgm & fit=="fit.pois"), 
              family="poisson")$fitted.values)
   
   D[D$fit=="fit.negb" & D$dgm==thisDgm, "yhat"] = 
      log(glm.nb(count ~ x, 
                 subset(D, dgm==thisDgm & fit=="fit.negb"))$fitted.values)
   
}

ggplot(D, aes(x=x, y=logcount)) + 
   geom_point(alpha=0.3) + 
   geom_line(aes(x=x, y=yhat), colour="blue") + 
   facet_grid(dgm ~ fit)

### Compute prediction intervals ###############################################

bootstrapN = 50

piLinm = predInterval(D[D$dgm=="true.linm" & D$fit=="fit.linm",], bN=bootstrapN)
piLinm$dgm = "true.linm"
piPois = predInterval(D[D$dgm=="true.pois" & D$fit=="fit.linm",], bN=bootstrapN)
piPois$dgm = "true.pois"
piNegb = predInterval(D[D$dgm=="true.negb" & D$fit=="fit.linm",], bN=bootstrapN)
piNegb$dgm = "true.negb"
piAll = rbind(piLinm, piPois, piNegb)

piAll$dgm = factor(piAll$dgm,levels=c("true.linm","true.pois","true.negb"),
                   ordered=TRUE)

D = merge(D, piAll, by=c("key", "dgm", "fit"))
dim(D)
head(D)
summary(D$lwr)
D$lwr[D$lwr == -Inf] = 0

# plot with prediction intervals, log count scale
ggplot(D, aes(x=x, y=logcount)) + 
   geom_point(alpha=0.5) + 
  geom_line(aes(x=x, y=yhat), colour="blue") + 
  geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=0.3, fill="blue") +
   facet_grid(dgm ~ fit)

# same on the count scale 
ggplot(D, aes(x=x, y=count)) + 
  geom_point(alpha=0.5) + 
  geom_line(aes(x=x, y=exp(yhat)), colour="blue") + 
  geom_ribbon(aes(ymin=exp(lwr), ymax=exp(upr)), alpha=0.3, fill="blue") +
  facet_grid(dgm ~ fit)+
  theme_economist_white(gray_bg = F, base_size = 11, base_family = "sans")

### Bias diagnostic ############################################################
A = aggregate(D$count, by = list(D$dgm, D$fit), FUN = sum)
names(A) = c("true","fitted","count")
B = aggregate(exp(D$yhat), by = list(D$dgm, D$fit), FUN = sum)
A[,"estimate"] = B$x
A$relBias = (A$estimate - A$count) / A$count
A 
# Conclusion: linear model is biased, others OK

### Coverage diagnostic ########################################################
D$cvrg = D$lwr <= D$logcount & D$logcount <= D$upr
AA = aggregate(D$cvrg, by = list(D$dgm, D$fit), FUN = function(x) sum(x)/n)
names(AA) = c("true","fitted","cvrg")
AA 
# Conclusion: coverage fitting poisson is only good when data are true poisson

# Overall conclusion: fitting negbin to data from all three distributions gives 
# reasonable results. Linear is biased and Poisson underestimates variance.

