#----------------------------------------------------------------------------
# a) simulação dos dados 

set.seed(1234)
n = 200 # tamanho da amostra
x = runif(200, 2, 40) # U(2, 40)
beta = c(60,-0.05) # valor verdadeiro dos coeficientes
e = rnorm(n, sd = 2)

y = beta[1]*exp(beta[2]*x) + e # y = b0e^b1X + e


require(ggplot2)
ggplot() + geom_point(aes(x,y)) + theme_bw() # plotagem dos dados

#----------------------------------------------------------------------------
# b) Newton-Raphson 

nr.optim = function(x0, f1, f2, epslon = 0.0001,...){
  
  cc = 1
  conta = 0
  x0 = matrix(x0, ncol = 1)
  z = x0
  
  while(cc > epslon && conta < 1000){
    x1 = x0 - solve(f2(x0, x, y))%*%(f1(x0, x, y))
    
    cc = t(f1(x1, x, y))%*%f1(x1, x, y) # criteiro de convergencia
    x0 = x1
    conta = conta + 1
    z = cbind(z,x1)
    
  }
  
  rownames(z) = paste("x", 1:nrow(x0), sep = "")
  colnames(z) = paste("i", 1:(conta + 1), sep = "")
  list("arg.ot" = x0, "seq.ot" = z, "iterações" = conta)  
}

# aplicar newton-raphson no problema

# grandiente
f1 = function(b,x,y){
  rbind(-2*sum((y - b[1]*exp(b[2]*x))*exp(b[2]*x)), 
        -2*sum((y - b[1]*exp(b[2]*x))*x*b[1]*exp(b[2]*x))
  )
}

# hessiana
f2 = function(b,x,y){
  A = matrix(0, ncol = 2,nrow = 2)
  A[1,1] = sum(-2*exp(b[2]*x)*(-exp(b[2]*x)))
  A[2,2] = sum(2*(-b[1]*x*exp(b[2]*x))*(-b[1]*x*exp(b[2]*x))) + sum(2*(y-b[1]*exp(b[2]*x))*(-b[1]*x^2*exp(b[2]*x)))
  A[2,1] = A[1,2] = sum(2*(-b[1]*x*exp(b[2]*x))*(-exp(b[2]*x))) + sum(2*(y-b[1]*exp(b[2]*x))*(-x*exp(b[2]*x)))
  
  
  return(A)
}

# o NR só funciona se o chute inicial for perto do valor verdadeiro, então 
# como valor inicial sera utilizado os valores estimados da regressão ajustada
# com a transformação log(Y) = log(b0) + b1x

reg = lm(log(y) ~ x)

# valor inicial: 
x0 = c(exp(reg$coefficients[1]), reg$coefficients[2])

nr.optim(x0, f1, f2, epslon = 0.001)

#----------------------------------------------------------------------------
# c) Gauss-Newton:

gn.optim = function(x0, J, r, epslon = 0.0001,...){
  cc = 1
  conta = 0
  x0 = matrix(x0, ncol = 1)
  z = x0
  
  while(cc > epslon && conta < 1000){
    x1 = x0 - solve(t(J(x0,x,y))%*%J(x0,x,y))%*%t(J(x0,x,y))%*%r(x0,x,y)
    
    cc = sum((x1 - x0)^2) # criteiro de convergencia
    x0 = x1
    conta = conta + 1
    z = cbind(z,x1)
  }
  
  rownames(z) = paste("x", 1:nrow(x0), sep = "")
  colnames(z) = paste("i", 1:(conta + 1), sep = "")
  list("arg.ot" = x0, "seq.ot" = z, "iterações" = conta)  
}


# residuos:
r <- function(b, x, y){
  y - b[1]*exp(b[2]*x)
}

# jacobiano:
J <- function(b, x, y){
  cbind(-exp(b[2]*x), -b[1]*x*exp(b[2]*x))
}

gn.optim(x0 = rbind(10, 0), J, r, epslon = 0.001)
