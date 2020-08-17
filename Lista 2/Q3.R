# a)

require(mvnfast)

m1 = c(0,0,0)
m2 = c(3,3,3)
n = 40
m = 3

set.seed(123)
m1 = c(0,0,0)
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
  teta[1,] = m0
  
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
  
  list("arg.opt" = teta.ot, "mean.accept" = mean(accept[-1]), 'it' = pos)
  
}

m0 = c(-1, -2, 1, 4, 7, 5)
simulated.annealing(x, m0 = m0)

# c) algoritmo EM 

EM.alg <- function(x, m0, eps = 1e-04){
  
  mu1Inicial = m0[1:3]
  mu2Inicial = m0[4:6]
  
  cc = 1
  
  conta = 0
  
  while(cc > eps){
    
    # Etapa E
    E = dmvnorm(x,mu1Inicial)/(dmvnorm(x,mu1Inicial) + 3*dmvnorm(x,mu2Inicial))
    
    # Etapa M
    mu1Par = colSums(x*E)/sum(E)
    mu2Par = colSums(x*(1 - E))/sum(1 - E)
    
    cc1 = (mu1Inicial- mu1Par)^2
    cc2 = (mu2Inicial - mu2Par)^2
    cc = mean(cc1 + cc2)
    mu1Inicial = mu1Par
    mu2Inicial = mu2Par
    
    conta = conta + 1
  }
  
  list('arg.opt' = c(mu1Par, mu2Par), 'E' = E, 'iter' = conta)
  
}

mix.true = attr(x, "index")
mix.EM = ifelse(round(EM.alg(x, m0)$E, 1) == 0, 2, 1)

table(mix.true, mix.EM)

# d) optim do R

optim(par = c(1, 0.5, -0.5, 4, 2.5, 2), fn = f, x = x, control = list(fnscale = -1))

# e) comparar

m0 = c(-1, -2, 1, 4, 7, 5)

set.seed(123)

SA = simulated.annealing(x, m0 = m0)
EM = EM.alg(x, m0)
OPT = optim(par = c(1, 0.5, -0.5, 4, 2.5, 2), fn = f, x = x, method = 'SANN', control = list(fnscale = -1))

round(SA$arg.opt, 2)
round(EM$arg.opt, 2)
round(OPT$par, 2)

SA$it
EM$iter
OPT$counts


# Trocar o valor inicial 

m0 = c(2, 2, 2, 2, 2, 2)

set.seed(123)

SA = simulated.annealing(x, m0 = m0)
EM = EM.alg(x, m0)
OPT = optim(par = m0, fn = f, x = x, method = 'SANN', control = list(fnscale = -1))

round(SA$arg.opt, 2)
round(EM$arg.opt, 2)
round(OPT$par, 2)

# nenhum dos três algoritmos possuem resultados satisfatórios

# f)

m2 = c(1,1,0)

set.seed(123)
x2 <- rmixn(n = n, mu = rbind(m1,m2), sigma = list(diag(m), diag(m)), 
            w = c(0.5, 0.5), retInd = TRUE)

# attr(x2, "index")

m0 = c(-1, -2, 1, 4, 7, 5)

set.seed(123)

SA = simulated.annealing(x2, m0 = m0)
EM = EM.alg(x2, m0)
OPT = optim(par = m0, fn = f, x = x2, method = 'SANN', control = list(fnscale = -1))


# valores verdadeiros 
c(m1, m2)

# valores estimados 
round(SA$arg.opt, 2)
round(EM$arg.opt, 2)
round(OPT$par, 2)

SA$it
EM$iter
OPT$counts

# podemos ver que o SA e o optim possuem comportamentos similares neste caso 
# eles não conseguem estimar corretamente os vetores de média 
# enquanto que o EM consegue. Foi o que melhor se saiu nesse caso. 

# trocar valor inicial
m0 = c(2, 2, 2, 2, 2, 2)

set.seed(123)

SA = simulated.annealing(x2, m0 = m0)
EM = EM.alg(x2, m0)
OPT = optim(par = m0, fn = f, x = x2, method = 'SANN', control = list(fnscale = -1))

round(SA$arg.opt, 2)
round(EM$arg.opt, 2)
round(OPT$par, 2)

# dessa vez todos os tres métodos falham na estimação dos vetores 
# de médias. 

# g) dados genéticos 

EM.alg <- function(x, m0, eps = 1e-04){
  
  aux1 = length(m0)/2
  
  mu1Inicial = m0[1:aux1]
  mu2Inicial = m0[(aux1 + 1):length(m0)]
  
  cc = 1
  conta = 0
  
  while(cc > eps){
    
    # Etapa E
    E = dmvnorm(x,mu1Inicial)/(dmvnorm(x,mu1Inicial) + 3*dmvnorm(x,mu2Inicial))
    
    # Etapa M
    mu1Par = colSums(x*E)/sum(E)
    mu2Par = colSums(x*(1 - E))/sum(1 - E)
    
    cc1 = (mu1Inicial- mu1Par)^2
    cc2 = (mu2Inicial - mu2Par)^2
    cc = mean(cc1 + cc2)
    mu1Inicial = mu1Par
    mu2Inicial = mu2Par
    
    conta = conta + 1
  }
  
  list('arg.opt' = c(mu1Par, mu2Par), 'E' = E, 'iter' = conta)
  
}
