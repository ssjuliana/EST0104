X = c(6.2, 5.1, 7.6, 2.5, 3.5, 9.4, 4.1, 6.3, 3.0, 0.8)
Y = c(6.9, 5.1, 7.5, 11.1, 10.9, 4.2, 10.5, 6.8, 12.3, 14.3)

# a) Utilize o bootstrap para construir um intervalo de confian¸ca 95% para β1 e β2.

# grafico 
plot(X, Y, pch = 19)

reg = lm(Y ~ X)

b1 = reg$coefficients[1]
b2 = reg$coefficients[2]


bootstrap = function(B, X, Y){
  
  b1 = b2 = vector()
  data = data.frame(X,Y)
  
  for (i in 1:B){
    
    aux = sample(1:length(X), length(X), replace = TRUE)
    dados = data[aux,]
    
    reg = lm(Y ~ X, data = dados)
    
    b1[i] = reg$coefficients[1]
    b2[i] = reg$coefficients[2]
  }
  
  return(cbind(b1, b2))
  
}

dados = bootstrap(B = 10000, X, Y)
ic = apply(dados, 2, quantile, c(0.025, 0.975))


par(mfrow = c(1,2))
hist(dados[,1], main = expression(paste("Estimates for ", beta[0])), xlab = expression(hat(beta)[0]))
abline(v = b1, col = 'blue', lty = 2)
abline(v = ic[1], col = 'red', lty = 2)
abline(v = ic[2], col = 'red', lty = 2)
text(5, 4000, expression(paste(beta[0], ' [IC95%] = ')), cex = 0.6)
text(11, 4000, paste0(round(b1, 1), ' [', round(ic[1],2), ',', round(ic[2],2), ']'), cex = 0.6)

hist(dados[,2], main = expression(paste("Estimates for ", beta[1])), xlab = expression(hat(beta)[1]))
abline(v = b2, col = 'blue', lty = 2)
abline(v = ic[3], col = 'red', lty = 2)
abline(v = ic[4], col = 'red', lty = 2)
text(-3, 4000, expression(paste(beta[1], ' [IC95%] = ')), cex = 0.6)
text(-1.8, 4000, paste0(round(b2, 1), ' [', round(ic[3],2), ',', round(ic[4],2), ']'), cex = 0.6)


# b) utilize o bootstrap parametrico 


boot.parar = function(B, X, sigma, b_1, b_2){
  
  s = sigma
  b1 = b2 = vector()
  
  for (i in 1:B){
    
    y = round(b_1 + b_2*X + rnorm(length(X), mean = 0, sd = s), 1)
    dados = data.frame(X,y)
    
    reg = lm(y ~ X, data = dados)
    
    b1[i] = reg$coefficients[1]
    b2[i] = reg$coefficients[2]
    
    s = sigma(reg)
  }
  
  return(cbind(b1, b2))

}

dados = boot.parar(B = 10000, X, s = sigma(reg), b_1 = b1, b_2 = b2)
ic = apply(dados, 2, quantile, c(0.025, 0.975))


par(mfrow = c(1,2))
hist(dados[,1], main = expression(paste("Estimates for ", beta[0])), xlab = expression(hat(beta)[0]))
abline(v = b1, col = 'blue', lty = 2)
abline(v = ic[1], col = 'red', lty = 2)
abline(v = ic[2], col = 'red', lty = 2)
text(12, 4000, expression(paste(beta[0], ' [IC95%] = ')), cex = 0.6)
text(18, 4000, paste0(round(b1, 1), ' [', round(ic[1],2), ',', round(ic[2],2), ']'), cex = 0.6)

hist(dados[,2], main = expression(paste("Estimates for ", beta[1])), xlab = expression(hat(beta)[1]))
abline(v = b2, col = 'blue', lty = 2)
abline(v = ic[3], col = 'red', lty = 2)
abline(v = ic[4], col = 'red', lty = 2)
text(-3, 4000, expression(paste(beta[1], ' [IC95%] = ')), cex = 0.6)
text(-2, 4000, paste0(round(b2, 1), ' [', round(ic[3],2), ',', round(ic[4],2), ']'), cex = 0.6)


# c) teste de permutação

# temos que fazer as replicações sobre H0, ou seja 
# B1 = 0, portanto => Y = b0 + e 

permut.coef <- function(X,Y) {
  
  y <- sample(Y, replace = FALSE)
  fit = lm(y ~ X)
  coef = fit$coefficients[2]
  
  return(coef)
}

B = 5000
permut <- replicate(B, permut.coef(X,Y))

pval <- mean(abs(permut) >= abs(b2))


hist(permut, main = paste0('Histogram of the Permutation \nwith ', B, ' replicates'), xlab = expression(hat(beta)[1]))
abline(v = b2, col = 'red', lty = 2)
text(-0.87, 400, paste0('p-value: ', pval), cex = 0.6)