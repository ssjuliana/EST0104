## Entrando a matriz X de dados 
# a entrada ij indica o numero de gols que i marcou em j

# DADOS BRASILEIRÃO 2018 

source('BR2018Q4.R')

n = ncol(X)    # numero de times

## Definir valores iniciais para O e D
d = rep(1,n)
o = rep(1,n)


block_relaxation <- function(o, d, X, epslon = 0.001){
  
  conta = 1 # iteraçoes 
  
  d_old = rep(0,ncol(X))
  o_old = rep(0,ncol(X))   
  
  O = matrix(o, ncol = 1)
  D = matrix(d, ncol = 1)
  
  while(sum((d-d_old)^2) + sum((o-o_old)^2) > epslon){ #criterio de convergencia
    o_old = o
    o = log(rowSums(X)/sum(exp(-d)))   #registra o valor antigo, e atualiza todos os o`s
    d_old = d
    d = -log(colSums(X)/sum(exp(o)))  #registra o valor antigo, e atualiza todos  d`s
    
    
    #Salva os valores na matriz  
    conta = conta + 1
    
    O = cbind(O, o)
    D = cbind(D, d)
  }
  
  it = conta - 1
  
  rownames(O) = rownames(D) = paste("x", 1:nrow(O), sep = "") #mudar para o nome do time
  colnames(O) = colnames(D) = paste("i", 0:it, sep = "")

  
  list("Defensivo" = list("arg.ot" = o, "seq.ot" = O), 
       "Ofensivo" = list("arg.ot" = d, "seq.ot" = D), 
       "iterações" = it)  
  
}

bl.r = block_relaxation(o, d, X)


# versão optim do R

# log-lik 
logL = function(theta, X){
  temp = 0
  
  for(i in 1:nrow(X)){
    for(j in 1:ncol(X)){
      temp = temp - exp(theta[i] - theta[n+j]) + X[i,j]*(theta[i] - theta[n+j]) 
      }
  }
  
  temp
}

theta0 = c(o, d)
opt = optim(theta0,logL, X = X, control = list(fnscale = -1))

# versão newton-raphson 


nr.optim = function(x0, f1, f2, epslon = 0.0001,...){
  
  cc = 1
  conta = 0
  x0 = matrix(x0, ncol = 1)
  z = x0
  
  while(cc > epslon && conta < 1000){
    x1 = x0 - solve(f2(x0, X), tol = 1e-200)%*%(f1(x0, X))
    
    cc = t(f1(x1, X))%*%f1(x1, X) # criteiro de convergencia
    x0 = x1
    conta = conta + 1
    z = cbind(z,x1)
    
  }
  
  rownames(z) = paste("x", 1:nrow(x0), sep = "")
  colnames(z) = paste("i", 0:conta, sep = "")
  list("arg.ot" = x0, "seq.ot" = z, "iterações" = conta)  
}

grad <- function(theta, X){
  
  A = matrix(NA, ncol = 1, nrow = 2*nrow(X))

  for(i in 1:nrow(X)){
    temp = 0
    for(j in 1:ncol(X)){
      temp = temp - exp(theta[i] - theta[n+j]) + X[i,j]
    }
    A[i,] <- temp
  }
  
  for(j in 1:nrow(X)){
    temp = 0
    for(i in 1:ncol(X)){
      temp = temp + exp(theta[i] - theta[n+j]) - X[i,j]
    }
    A[n+j,] <- temp
  }
  
  return(A)

}

hessian <- function(theta, X){
  
  H <- matrix(0, ncol = 2*ncol(X), nrow = 2*nrow(X))
  
  # diagonal 
  for(i in 1:nrow(X)){
    temp = 0
    for(j in 1:ncol(X)){
      temp = temp - exp(theta[i] - theta[n+j]) 
    }
    H[i,i] <- temp
  }
  
  for(j in 1:nrow(X)){
    temp = 0
    for(i in 1:ncol(X)){
      temp = temp - exp(theta[i] - theta[n+j])
    }
    H[n+j,n+j] <- temp
  }
  
  # o_i d_k 
  
  for(i in 1:nrow(X)){
    for(j in 1:ncol(X)){
      H[i,ncol(X)+j] = exp(theta[i] - theta[n+j]) 
    }
  }
  
  #d_j o_k
  
  for(j in 1:nrow(X)){
    for(i in 1:ncol(X)){
      H[ncol(X)+j, i] <- exp(theta[i] - theta[n+j])
    }
  }

  
  return(H)
}

# valor inicial os par do block
theta0 = c(block_relaxation(o, d, X)$Defensivo$arg.ot, block_relaxation(o, d, X)$Ofensivo$arg.ot)
nr = nr.optim(theta0, grad, hessian)

results = cbind(cbind(c(bl.r$Defensivo$arg.ot, bl.r$Ofensivo$arg.ot)), 
                 cbind(opt$par), nr$arg.ot)

colnames(results) = c('BR', 'NM', 'NR')
results

require(ggplot2)
require(tidyr)

r.df = gather(data.frame(time = colnames(X), 
                         par = c(rep('Ofensivo', 20), rep('Defensivo', 20)), 
                         results), 
              metodo, arg.ot, BR:NR, factor_key = TRUE)

ggplot(data = r.df, mapping = aes(x = time, y = arg.ot, color = metodo, group = metodo)) + 
  geom_line(linetype = 'dashed') +
  geom_point() + 
  facet_wrap( ~ par, ncol = 2) + 
  labs(x = 'Time', y = 'Estimativa', color = 'Método')
