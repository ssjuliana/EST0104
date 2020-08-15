# a)

require(mvnfast)

m1 = c(0,0,0)
m2 = c(3,3,3)
n = 40
m = 3

set.seed(123)

x <- rmixn(n = n, mu = rbind(m1,m2), sigma = list(diag(m), diag(m)), 
           w = c(0.5, 0.5), retInd = TRUE)

require(ggplot2)
require(tidyr)

x.df = gather(as.data.frame(x), v, x, V1:V3, factor_key = TRUE)

# os dados 
ggplot(data = x.df) + 
  geom_density(aes(x = x)) +
  facet_wrap( ~ v, ncol = 3) + 
  theme_bw()


# b) simulated annealing para m1 e m2

# a função vai poder colocar pesos nas normais 
# como default vou usar que cada um tem peso 0.5,
# igual a simulação dos dados 

require(mvtnorm)

f = function(m, x, p = c(0.5, 0.5)){ # log da verossimilhança
  m1 = m[1:3]
  m2 = m[4:6]
  
  p1 = p[1]
  p2 = p[2]
  
  sum(log(p1*dmvnorm(x = x , mean = m1) + p2*dmvnorm(x = x, mean = m2)))
}

simulated.annealing <- function(x, m0){
  
  n.iter = 500
  
  teta = matrix(0, nrow = n.iter, ncol = length(m0))
  teta[1,] <- m0
  
  f.eval = accept = rep(0, n.iter)
  f.eval[1] = f(teta[1,], x)
  
  
  for(i in 2:n.iter){
    
    temp = 1/log(1+i) 
    d = 0.5 * sqrt(temp)
    e = runif(length(m0),-d,d) 
    
    te1 = teta[(i-1),] + e
    u = runif(1)
    
    prob = min(exp((f(te1,x) - f(teta[(i-1),],x))/temp), 1)
    
    teta[i,] = (u <= prob)*te1 + (u > prob)*teta[(i-1),]
    f.eval[i] = f(teta[i,],x)
    
    accept[i] = (u <= prob)
  }
  
  pos = which.max(f.eval)
  teta.ot = teta[pos,]
  
  list("arg.opt" = teta.ot, "mean.accept" = mean(accept[-1]))
  
}

m0 = c(-1, -1, -1, 5, 5, 5)
simulated.annealing(x, m0 = m0)

# c) algoritmo EM 

# d) optim do R

optim(par = c(0, 0, 0, 3, 3, 3), fn = f, x = x, control = list(fnscale = -1))$par


