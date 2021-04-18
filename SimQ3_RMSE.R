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
#ptm <- proc.time()
iter = 20
M = c(100,120,140, 160, 180,200)
P = c(2,5,10,20)


#-------------------------------------------#
# Case 1 - eps = 0.1, k = 2, t = 5, N = 1500 #
#-------------------------------------------#

N2 = 1500

#matrices for all (p,m) combinations
ARR.tau = matrix(0.0, nrow = length(P), ncol = length(M))         #Final
ARR.M = matrix(0.0, nrow = length(P), ncol = length(M))
ARR.LS = matrix(0.0, nrow = length(P), ncol = length(M))
ARR.MM = matrix(0.0, nrow = length(P), ncol = length(M))

for (m in 1:length(M)){
  for (p in 1:length(P)){
    
    sig.tau = matrix(0.0, nrow = P[p], ncol = iter)
    sig.M = matrix(0.0, nrow = P[p], ncol = iter)
    sig.LS = matrix(0.0, nrow = P[p], ncol = iter)
    sig.MM = matrix(0.0, nrow = P[p], ncol = iter)
    #sig.LTS = matrix(0.0, nrow = P[p], ncol = iter)
    #sig.S = matrix(0.0, nrow = P[p], ncol = iter)
    
    
    for (i in 1:iter){
      if ((i%%50==0)){print(i)}  
      model = makeData(P[p], 0.2, M[m])
      y = model$y
      X1 = model$X
      
      #bringing the column of 1's to the front
      X1 = X1[,c(P[p], 1:P[p]-1)]
      
      #X without the 1's column
      X1. = as.matrix(X1[,-1])
      
      TauModel = FastTau(x=X1, y=y, N=N2, kk=2, tt=5, approximate=0)
      sig.tau[,i] = TauModel$beta 
      
      LSModel <- lm(y ~ X1.)
      sig.LS[,i] = LSModel$coefficients
      
      #LTSModel <- ltsReg(y ~ X1.) #intercept automatically added
      #sig.LTS[i] = LTSModel$scale
      
      #m1 <- lmrob..M..fit(X1, y, beta.initial = coef(LTSModel), scale = LTSModel$scale,
      #control = lmrob.control(copute.rd=T)) #dataset should have an intercept column
      m1 = rlm(X1, y, method = "M", maxit = 10000)
      sig.M[,i] = m1$coefficients
      
      #Ops <- lmrob.control(copute.rd=T)
      #Ops$max.it <- N2
      #MMModel <- lmrob(y ~ X1., control = Ops)
      MMModel = rlm(X1, y, method = "MM", maxit = 10000)
      sig.MM[,i] = MMModel$coefficients
      
      
    }
    
    ARR.tau[p,m] = sqrt(mean(apply(sig.tau, 2, function(x) (sum(x^2)))))
    ARR.M[p,m] = sqrt(mean(apply(sig.M, 2, function(x) (sum(x^2)))))
    ARR.MM[p,m] = sqrt(mean(apply(sig.MM, 2, function(x) (sum(x^2)))))
    ARR.LS[p,m] = sqrt(mean(apply(sig.LS, 2, function(x) (sum(x^2)))))
    
  }
  print(m)
}


#averaging over the 6 sample types

arr.tau = rowMeans(ARR.tau)
arr.LS = rowMeans(ARR.LS)
arr.M = rowMeans(ARR.M)
arr.MM = rowMeans(ARR.MM)

arr.tau.2 = rowMeans(ARR.tau)
arr.LS.2 = rowMeans(ARR.LS)
arr.M.2 = rowMeans(ARR.M)
arr.MM.2 = rowMeans(ARR.MM)


#####################################
#           Plotting                #
#####################################


Ylimits.5 = c(0.0, 2)
Ylimits.2 = c(0.5, max(arr.LS.2/arr.tau.2, arr.MM.2/arr.tau.2, arr.M.2/arr.tau.2)) #,arr.LTS/arr.tau, arr.S/arr.tau))
Ylimits.10 = c(0.5, max(arr.LS.10/arr.tau.10, arr.MM.10/arr.tau.10, arr.M.10/arr.tau.10)) #,arr.LTS/arr.tau, arr.S/arr.tau))


#jpeg("eps0.1.jpg", width = 800, height = 640, quality = 100)
par(mfrow = c(1,2))


plot(P, arr.tau, ylim = c(0,3), pch = 20, xlab = 'p', ylab = 'ROOT MEAN SQUARED ERROR', type = 'o',lty = 1, lwd = 2, 
     main = TeX(r'($\epsilon = 0.1,$ $ N = 1500$)'), cex.lab = 1.2, cex.axis = 1.2, cex.main = 1.4)
#abline(h = 1, lwd = 2, col = "black")

lines(P, arr.LS, ylim = Ylimits.5, type = 'b', lty = 2, pch = 0)
lines(P, arr.M, ylim = Ylimits.5, type = "b",lty = 3, pch = 5)
lines(P, arr.MM, ylim = Ylimits.5, type = "b", lty = 4, pch = 6)
#legend("bottomright",c("Fast-Tau", "M", "MM", "LS"), lty=c(1, 3, 4, 2), pch = c(20, 5, 6, 0), cex=1.0)


plot(P, arr.tau.2, ylim = c(0, 3), pch = 20, xlab = 'p', ylab = '', type = 'o',lty = 1, lwd = 2, 
     main = TeX(r'($\epsilon = 0.2,$ $ N = 1500$)'), cex.lab = 1.2, cex.axis = 1.2, cex.main = 1.4)
#abline(h = 1, lwd = 2, col = "black")

lines(P, arr.LS.2, ylim = Ylimits.5, type = 'b', lty = 2, pch = 0)
lines(P, arr.M.2, ylim = Ylimits.5, type = "b",lty = 3, pch = 5)
lines(P, arr.MM.2, ylim = Ylimits.5, type = "b", lty = 4, pch = 6)
legend("bottomright",c("Fast-Tau", "M", "MM", "LS"), lty=c(1, 3, 4, 2), pch = c(20, 5, 6, 0), cex=1.0)








#----------------------------------------------------------------------------




# 
# iter = 2
# 
# beta = matrix(0.0, ncol=iter, nrow = 5)
# for (i in 1:iter){
#   if ((i%%50==0)){print(i)}  
#   model = makeData(5, 0.1, 100)
#   y = model$y
#   X1 = model$X
#   
#   #bringing the column of 1's to the front
#   X1 = X1[,c(5, 1:4)]
#   
#   #X without the 1's column
#   X1. = as.matrix(X1[,-1])
#   
#   TauModel = FastTau(x=X1, y=y, N=N2, kk=2, tt=5, approximate=0)
#   #sig.tau[i] = TauModel$scale 
#   beta[,i] = TauModel$beta
#   
#   
# }
# 
# s = 0
# for(i in 1:2){
#   a = sum(beta[,i]^2)
#   s = a + s 
# }
# sqrt(s/iter)
# 
# sqrt(mean(beta^2))
# beta0 = rep(0,5)
# sqrt((sum((beta[,1] - beta0)^2) + sum((beta[,2] - beta0)^2))/2)
# sqrt(mean(apply(beta, 2, function(x) (sum(x^2)))))
