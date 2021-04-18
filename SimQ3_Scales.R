#clean

rm(list = ls())

#-----------------------

source("FastTau/FastTau.R")

n = 100

#-----------


library(MASS)
library(latex2exp)

makeData = function(p, eps, m){
  
  # '.g' for good, '.b' for bad
  
  n.g = (1-eps)*n
  n.b = n - n.g
  
  mu.g = rep(0, p-1)
  mu.b = c(100, rep(0,p-2))
  I = diag(rep(1, p-1))
  
  #good data
  X.g = mvrnorm(n = n.g, mu = mu.g, Sigma = I)
  X.g.p = as.matrix(rep(1, n.g))          #pth covariate : x_ip = 1
  X.g = cbind(X.g, X.g.p) 
  
  err = rnorm(n = n.g, mean = 0, sd = 1)
  B_0 = as.matrix(rep(0, p))              #coefficient vector
  y.g =  X.g %*% (B_0) + err
  
  #bad data : high leverage points
  X.b = mvrnorm(n = n.b, mu = mu.b, Sigma = 0.1^2 * I)
  X.b = cbind(X.b, as.matrix(rep(1, n.b)))  #pth covariate : x_ip = 1
  
  y.b = as.matrix(rnorm(n = n.b, mean = m, sd = 0.1^2))
  
  
  #full dataset
  
  X = rbind(X.g, X.b)         # n X (p-1) matrix
  y = rbind(y.g, y.b)
  
  return(list("X" = X, "y" = y))
}  

# p = {2, 5, 10}
# eps = {0.1, 0.2}
# m = {100, 120, ... 200}



#############################################
#               ESTIMATION                  #
#############################################
# Res. Scale Estimate in different settings #  
#############################################
ptm <- proc.time()
iter = 20
M = c(100,120,140, 160, 180,200)
P = c(2,5,10,20)

#-------------------------------------------#
# Case 1 - k = 2, t = 5, N = 500            #
#-------------------------------------------#

N2 = 500

#matrices for all (p,m) combinations
ARR.tau = matrix(0.0, nrow = length(P), ncol = length(M))         #Final
ARR.M = matrix(0.0, nrow = length(P), ncol = length(M))
ARR.LS = matrix(0.0, nrow = length(P), ncol = length(M))
ARR.MM = matrix(0.0, nrow = length(P), ncol = length(M))

for (m in 1:length(M)){
  for (p in 1:length(P)){
    
    sig.tau = rep(0, iter)
    sig.M = rep(0, iter)
    sig.LS = rep(0, iter)
    sig.MM = rep(0, iter)
    #sig.LTS = rep(0, iter)
    #sig.S = rep(0, iter)
    
    
    for (i in 1:iter){
      if ((i%%50==0)){print(i)}  
      model = makeData(P[p], 0.2, M[m])
      y = model$y
      X1 = model$X
      
      #bringing the column of 1's to the front
      X1 = X1[,c(P[p], 1:(P[p]-1))]
      
      #X without the 1's column
      X1. = as.matrix(X1[,-1])
      
      TauModel = FastTau(x=X1, y=y, N=N2, kk=2, tt=5, approximate=0)
      sig.tau[i] = TauModel$scale 
      
      LSModel <- lm(y ~ X1.)
      sig.LS[i] = summary(LSModel)$sigma
      
      #LTSModel <- ltsReg(y ~ X1.) #intercept automatically added
      #sig.LTS[i] = LTSModel$scale
      
      #m1 <- lmrob..M..fit(X1, y, beta.initial = coef(LTSModel), scale = LTSModel$scale,
                          #control = lmrob.control(copute.rd=T)) #dataset should have an intercept column
      m1 = rlm(X1, y, method = "M", maxit = 10000)
      sig.M[i] = summary(m1)$sigma
      
      #Ops <- lmrob.control(copute.rd=T)
      #Ops$max.it <- N2
      #MMModel <- lmrob(y ~ X1., control = Ops)
      MMModel = rlm(X1, y, method = "MM", maxit = 10000)
      sig.MM[i] = summary(MMModel)$sigma
      

    }
    
    ARR.tau[p,m] = mean(sig.tau)
    ARR.M[p,m] = mean(sig.M)
    ARR.MM[p,m] = mean(sig.MM)
    ARR.LS[p,m] = mean(sig.LS)

  }
  print(N2)
}

#-------------------------------------------#
# Case 2 - k = 2, t = 5, N = 200            #
#-------------------------------------------#

N1 = 200

#matrices for all (p,m) combinations
ARR.tau.2 = matrix(0.0, nrow = length(P), ncol = length(M))         #Final
ARR.M.2 = matrix(0.0, nrow = length(P), ncol = length(M))
ARR.LS.2 = matrix(0.0, nrow = length(P), ncol = length(M))
ARR.MM.2 = matrix(0.0, nrow = length(P), ncol = length(M))


for (m in 1:length(M)){
  for (p in 1:length(P)){
    
    sig.tau = rep(0, iter)
    sig.M = rep(0, iter)
    sig.LS = rep(0, iter)
    sig.MM = rep(0, iter)
    sig.LTS = rep(0, iter)
    sig.S = rep(0, iter)
    
    
    for (i in 1:iter){
      if ((i%%50==0)){print(i)}  
      
      model = makeData(P[p], 0.2, M[m])
      y = model$y
      X1 = model$X
      
      #bringing the column of 1's to the front
      X1 = X1[,c(P[p], 1:P[p]-1)]
      
      #X without the 1's column
      X1. = as.matrix(X1[,-1])
      
      TauModel = FastTau(x=X1, y=y, N=N1, kk=2, tt=5, approximate=0)
      sig.tau[i] = TauModel$scale 
      
      LSModel <- lm(y ~ X1.)
      sig.LS[i] = summary(LSModel)$sigma
      
      #LTSModel <- ltsReg(y ~ X1.) #intercept automatically added
      #sig.LTS[i] = LTSModel$scale
      
      #m1 <- lmrob..M..fit(X1, y, beta.initial = coef(LTSModel), scale = LTSModel$scale,
      #control = lmrob.control(copute.rd=T)) #dataset should have an intercept column
      m1 = rlm(X1, y, method = "M", maxit = 10000)
      sig.M[i] = summary(m1)$sigma
      
      #Ops <- lmrob.control(copute.rd=T)
      #Ops$max.it <- N2
      #MMModel <- lmrob(y ~ X1., control = Ops)
      MMModel = rlm(X1, y, method = "MM", maxit = 10000)
      sig.MM[i] = summary(MMModel)$sigma
      
      
    }
    
    ARR.tau.2[p,m] = mean(sig.tau)
    ARR.M.2[p,m] = mean(sig.M)
    ARR.MM.2[p,m] = mean(sig.MM)
    ARR.LS.2[p,m] = mean(sig.LS)
    
  }
  print(N1)
}

#-------------------------------------------#
# Case 3 - k = 2, t = 5, N = 1000           #
#-------------------------------------------#

N3 = 1500

#matrices for all (p,m) combinations
ARR.tau.10 = matrix(0.0, nrow = length(P), ncol = length(M))         #Final
ARR.M.10 = matrix(0.0, nrow = length(P), ncol = length(M))
ARR.LS.10 = matrix(0.0, nrow = length(P), ncol = length(M))
ARR.MM.10 = matrix(0.0, nrow = length(P), ncol = length(M))


for (m in 1:length(M)){
  for (p in 1:length(P)){
    
    sig.tau = rep(0, iter)
    sig.M = rep(0, iter)
    sig.LS = rep(0, iter)
    sig.MM = rep(0, iter)
    #sig.LTS = rep(0, iter)
    #sig.S = rep(0, iter)
    
    
    for (i in 1:iter){
      if ((i%%50==0)){print(i)}  
      
      model = makeData(P[p], 0.2, M[m])
      y = model$y
      X1 = model$X
      
      #bringing the column of 1's to the front
      X1 = X1[,c(P[p], 1:(P[p]-1))]
      
      #X without the 1's column
      X1. = as.matrix(X1[,-1])
      
      TauModel = FastTau(x=X1, y=y, N=N3, kk=2, tt=5, approximate=0)
      sig.tau[i] = TauModel$scale 
      
      LSModel <- lm(y ~ X1.)
      sig.LS[i] = summary(LSModel)$sigma
      
      #LTSModel <- ltsReg(y ~ X1.) #intercept automatically added
      #sig.LTS[i] = LTSModel$scale
      
      #m1 <- lmrob..M..fit(X1, y, beta.initial = coef(LTSModel), scale = LTSModel$scale,
      #control = lmrob.control(copute.rd=T)) #dataset should have an intercept column
      m1 = rlm(X1, y, method = "M", maxit = 10000)
      sig.M[i] = summary(m1)$sigma
      
      #Ops <- lmrob.control(copute.rd=T)
      #Ops$max.it <- N2
      #MMModel <- lmrob(y ~ X1., control = Ops)
      MMModel = rlm(X1, y, method = "MM", maxit = 10000)
      sig.MM[i] = summary(MMModel)$sigma
    }
    
    ARR.tau.10[p,m] = mean(sig.tau)
    ARR.M.10[p,m] = mean(sig.M)
    ARR.MM.10[p,m] = mean(sig.MM)
    ARR.LS.10[p,m] = mean(sig.LS)
    
  }
  print(N3)
}
TotTime  = proc.time() - ptm


#-----------------------------------------------------------------------


colnames(ARR.tau) = c("m100","m120", "m140", "m160", "m180", "m200")
rownames(ARR.tau) = c("p2", "p5", "p10", "p20")

par(mfrow = c(1,2))
corrplot(ARR.tau, is.corr = F, method = "color", addgrid.col = "grey", 
         tl.col = "black", cl.lim = c(min(ARR.tau, ARR.MM), max(ARR.tau, ARR.MM)))
corrplot(ARR.MM, is.corr = F, method = "color", addgrid.col = "grey",
         tl.col = "black", cl.lim = c(min(ARR.tau, ARR.MM), max(ARR.tau, ARR.MM)))

#########################################
#--------eps = 0.1----------------------#
#########################################


#averaging over the 6 sample types:eps=0.1

arr.tau = rowMeans(ARR.tau)
arr.LS = rowMeans(ARR.LS)
arr.M = rowMeans(ARR.M)
arr.MM = rowMeans(ARR.MM)

arr.tau.2 = rowMeans(ARR.tau.2)
arr.LS.2 = rowMeans(ARR.LS.2)
arr.M.2 = rowMeans(ARR.M.2)
arr.MM.2 = rowMeans(ARR.MM.2)

arr.tau.10 = rowMeans(ARR.tau.10)
arr.LS.10 = rowMeans(ARR.LS.10)
arr.M.10 = rowMeans(ARR.M.10)
arr.MM.10 = rowMeans(ARR.MM.10)





#####################################
#           Plotting:eps=0.1        #
#####################################

Ylimits = c(0.5,2)
# Ylimits.5 = c(0.5, max(arr.LS/arr.tau, arr.MM/arr.tau, arr.M/arr.tau)) #,arr.LTS/arr.tau, arr.S/arr.tau))
# Ylimits.2 = c(0.5, max(arr.LS.2/arr.tau.2, arr.MM.2/arr.tau.2, arr.M.2/arr.tau.2)) #,arr.LTS/arr.tau, arr.S/arr.tau))
# Ylimits.10 = c(0.5, max(arr.LS.10/arr.tau.10, arr.MM.10/arr.tau.10, arr.M.10/arr.tau.10)) #,arr.LTS/arr.tau, arr.S/arr.tau))


#jpeg("eps0.1.jpg", width = 800, height = 640, quality = 100)
par(mfrow = c(2,3))


plot(P, arr.tau.2/arr.tau.2, ylim = Ylimits, pch = 20, xlab = 'p', ylab = "", 
     main = TeX(r'($\epsilon = 0.1,$ $ N = 200$)'), cex.lab = 1.4, cex.axis = 1.2, cex.main = 1.4)
title(ylab="Ave(Residual Scale Estimate) (Tau Normalized)", line=2.6, cex.lab=1.2)
abline(h = 1, lwd = 2, col = "black")

lines(P, arr.LS.2/arr.tau.2, type = 'b', lty = 2, pch = 0)
lines(P, arr.M.2/arr.tau.2, type = "b", pch = 5,lty = 2)
lines(P, arr.MM.2/arr.tau.2, type = "b", pch = 6, lty = 2)
legend("topleft",c("Fast-Tau", "M", "MM", "LS"), lty=c(1, 2, 2, 2), pch = c(20, 5, 6, 0), cex=1.0)


plot(P, arr.tau/arr.tau, ylim = Ylimits, pch = 20, xlab = 'p', ylab = '',  
     main = TeX(r'($\epsilon = 0.1,$ $ N = 500$)'), cex.lab = 1.2, cex.axis = 1.2, cex.main = 1.4)
abline(h = 1, lwd = 2, col = "black")

lines(P, arr.LS/arr.tau, type = 'b', lty = 2, pch = 0)
lines(P, arr.M/arr.tau, type = "b", pch = 5,lty = 2)
lines(P, arr.MM/arr.tau, type = "b", pch = 6, lty = 2)
#legend("bottomright",c("Fast-Tau", "M", "MM", "LS"), lty=c(1, 2, 2, 2), pch = c(20, 5, 6, 0), cex=0.8)


plot(P, arr.tau.10/arr.tau.10, ylim = Ylimits, pch = 20, xlab = 'p',  ylab = '',
     main = TeX(r'($\epsilon = 0.1,$ $ N = 1500$)'), cex.lab = 1.2, cex.axis = 1.2, cex.main = 1.4)
abline(h = 1, lwd = 2, col = "black")

lines(P, arr.LS.10/arr.tau.10, type = 'b', lty = 2, pch = 0)
lines(P, arr.M.10/arr.tau.10, type = "b", pch = 5,lty = 2)
lines(P, arr.MM.10/arr.tau.10, type = "b", pch = 6, lty = 2)
#legend("bottomright",c("Fast-Tau", "M", "MM", "LS"), lty=c(1, 2, 2, 2), pch = c(20, 5, 6, 0), cex=0.8)
#dev.off()


#S and MM are coming to be the same 
#LTS and M are soming to be the same

#########################################
#--------eps = 0.2----------------------#
#########################################

#averaging over the 6 sample types:eps=0.2

arr.tau.c2= rowMeans(ARR.tau)
arr.LS.c2 = rowMeans(ARR.LS)
arr.M.c2 = rowMeans(ARR.M)
arr.MM.c2 = rowMeans(ARR.MM)

arr.tau.c2.2 = rowMeans(ARR.tau.2)
arr.LS.c2.2 = rowMeans(ARR.LS.2)
arr.M.c2.2 = rowMeans(ARR.M.2)
arr.MM.c2.2 = rowMeans(ARR.MM.2)

arr.tau.c2.10 = rowMeans(ARR.tau.10)
arr.LS.c2.10 = rowMeans(ARR.LS.10)
arr.M.c2.10 = rowMeans(ARR.M.10)
arr.MM.c2.10 = rowMeans(ARR.MM.10)





#####################################
#           Plotting:eps=0.2        #
#####################################

Ylimits = c(0.5,2)
# 
# Ylimits.5 = c(0.5, max(arr.LS/arr.tau, arr.MM/arr.tau, arr.M/arr.tau)) #,arr.LTS/arr.tau, arr.S/arr.tau))
# Ylimits.2 = c(0.5, max(arr.LS.2/arr.tau.2, arr.MM.2/arr.tau.2, arr.M.2/arr.tau.2)) #,arr.LTS/arr.tau, arr.S/arr.tau))
# Ylimits.10 = c(0.5, max(arr.LS.10/arr.tau.10, arr.MM.10/arr.tau.10, arr.M.10/arr.tau.10)) #,arr.LTS/arr.tau, arr.S/arr.tau))


#jpeg("eps0.1.jpg", width = 800, height = 640, quality = 100)
par(mfrow = c(1,3))


plot(P, arr.tau.c2.2/arr.tau.c2.2, ylim = Ylimits, pch = 20, xlab = 'p', ylab = "", 
     main = TeX(r'($\epsilon = 0.2,$ $ N = 200$)'), cex.lab = 1.4, cex.axis = 1.2, cex.main = 1.4)
title(ylab="Ave(Residual Scale Estimate) (Tau Normalized)", line=2.6, cex.lab=1.2)
abline(h = 1, lwd = 2, col = "black")

lines(P, arr.LS.c2.2/arr.tau.c2.2, type = 'b', lty = 2, pch = 0)
lines(P, arr.M.c2.2/arr.tau.c2.2, type = "b", pch = 5,lty = 2)
lines(P, arr.MM.c2.2/arr.tau.c2.2, type = "b", pch = 6, lty = 2)
#legend("bottomright",c("Fast-Tau", "M", "MM", "LS"), lty=c(1, 2, 2, 2), pch = c(20, 5, 6, 0), cex=1.4)


plot(P, arr.tau.c2/arr.tau.c2, ylim = Ylimits, pch = 20, xlab = 'p', ylab = '',  
     main = TeX(r'($\epsilon = 0.2,$ $ N = 500$)'), cex.lab = 1.2, cex.axis = 1.2, cex.main = 1.4)
abline(h = 1, lwd = 2, col = "black")

lines(P, arr.LS.c2/arr.tau.c2, type = 'b', lty = 2, pch = 0)
lines(P, arr.M.c2/arr.tau.c2, type = "b", pch = 5,lty = 2)
lines(P, arr.MM.c2/arr.tau.c2, type = "b", pch = 6, lty = 2)
#legend("bottomright",c("Fast-Tau", "M", "MM", "LS"), lty=c(1, 2, 2, 2), pch = c(20, 5, 6, 0), cex=0.8)


plot(P, arr.tau.c2.10/arr.tau.c2.10, ylim = Ylimits, pch = 20, xlab = 'p',  ylab = '',
     main = TeX(r'($\epsilon = 0.2,$ $ N = 1000$)'), cex.lab = 1.2, cex.axis = 1.2, cex.main = 1.4)
abline(h = 1, lwd = 2, col = "black")

lines(P, arr.LS.c2.10/arr.tau.c2.10, type = 'b', lty = 2, pch = 0)
lines(P, arr.M.c2.10/arr.tau.c2.10, type = "b", pch = 5,lty = 2)
lines(P, arr.MM.c2.10/arr.tau.c2.10, type = "b", pch = 6, lty = 2)
#legend("bottomright",c("Fast-Tau", "M", "MM", "LS"), lty=c(1, 2, 2, 2), pch = c(20, 5, 6, 0), cex=0.8)
#dev.off()

