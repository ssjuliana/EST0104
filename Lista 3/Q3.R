# a ) 

X = function(rep){
  
  X.rep = vector()
  
  for (i in 1:rep){
    soma = 0 
    cont = 0 
    
    while(!soma > 1){
      u = runif(1)
      soma = soma + u
      cont = cont + 1
    }
    X.rep[i] = cont
  }
  return(X.rep)
}

# a) E(X)

# Utilizando simulação de Monte Carlo com n = 10.000
n = 10000

x = X(n)

mean(x) # media 
sqrt(1/n^2*sum((x - mean(x))^2)) # sd do erro do monte carlo


# b) P(X ≥ 10)

x.prob = function(k, param, rep){
  
  b = matrix(NA, nrow = rep, ncol = k)
  for (i in 1:k) b[,i] = rbeta(rep, shape1 = param[1], shape2 = param[2])
  h = apply(b, 1, sum) <= 1 #ind quais das betas tem soma <= 1
  f = 1 # U(0,1)
  g = apply(dbeta(b, param[1], param[2]), 1, prod)
  return(h*f/g)
  
}

x = x.prob(k = 9, param = c(1.5, 12), rep = n)

mean(x) # media 
sqrt(1/n^2*sum((x - mean(x))^2)) # sd do erro do monte carlo

# c) P(X = 10)

x.prob.equal = function(k, param, rep){
  
  b = matrix(NA, nrow = rep, ncol = k)
  for (i in 1:k) b[,i] = rbeta(rep, shape1 = param[1], shape2 = param[2])
  h = apply(b[,-k], 1, sum) <= 1 & apply(b, 1, sum) > 1 #ind quais das betas tem soma <= 1
  f = 1 # U(0,1)
  g = apply(dbeta(b, param[1], param[2]), 1, prod)
  return(h*f/g)
  
}

x = x.prob.equal(k = 10, param = c(1.5, 12), rep = n)

mean(x) # media 
sqrt(1/n^2*sum((x - mean(x))^2)) # sd do erro do monte carlo
