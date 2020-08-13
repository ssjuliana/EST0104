## Entrando a matriz X de dados 
# a entrada ij indica o numero de gols que i marcou em j

n = 10    # numero de times

X = matrix(c(0,2,5,6,3,1,3,3,1,2,
           1,0,2,1,0,0,1,1,6,0,
           2,1,0,1,0,0,2,0,1,5,
           0,2,2,0,1,3,4,1,1,1,
           2,2,3,1,0,4,4,2,3,1,
           2,1,0,1,0,0,1,0,1,3,
           1,0,2,1,0,0,0,1,2,5,
           1,0,2,1,0,2,1,0,0,0,
           0,2,2,2,0,1,3,3,0,2,
           0,0,2,1,0,1,1,0,0,0), byrow=T,
         ncol=10, nrow=10)


## Definir valores iniciais para O e D
d = rep(1,10)
o = rep(1,10)


block_relaxation <- function(o, d, X, epslon = 0.001){
  
  conta = 1 # iteraçoes 
  
  d_old = rep(0,10)
  o_old = rep(0,10)   
  
  O = matrix(o, ncol = 1)
  D = matrix(d, ncol = 1)
  
  while(sum((d-d_old)^2) + sum((o-o_old)^2) > epslon){ #criterio de convergencia
    o_old = o
    o = log(rowSums(X)/sum(exp(-d)))   #registra o valor antigo, e atualiza todos os o`s
    d_old = dos
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

block_relaxation(o, d, X)


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
optim(theta0,logL, X = X)

# versão newton-raphson 

x0 = c(theta0)

nr.optim = function(x0, f1, f2, epslon = 0.0001,...){
  
  cc = 1
  conta = 0
  x0 = matrix(x0, ncol = 1)
  z = x0
  
  while(cc > epslon && conta < 1000){
    x1 = x0 - solve(f2(x0, X))%*%(f1(x0, X))
    
    cc = sum((x1 - x0)^2) # criteiro de convergencia
    x0 = x1
    conta = conta + 1
    z = cbind(z,x1)
    
  }
  
  rownames(z) = paste("x", 1:nrow(x0), sep = "")
  colnames(z) = paste("i", 0:conta, sep = "")
  list("arg.ot" = x0, "seq.ot" = z, "iterações" = conta)  
}

nr.optim(theta0, grad, hessian)

# nao está finalizado