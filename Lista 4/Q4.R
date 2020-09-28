##################################################################################
###############  Modelo de Regressão
##################################################################################

## Definição dos parametros da priori
m0 = 0; m1 = 0; t0 = 1/200; t1 = 1/200; a = 1; b = 0.2

#Simulação dos dados
set.seed(13)
x = runif(50, 1, 9)
y = 5 - 0.5*x + rnorm(50)

# a) Qual a posteriori para esse modelo?

posteriori = function(v){
  ans = -Inf
  
  if(v[3]>0){
    ans = sum(dnorm(y, v[1] + v[2]*x, sqrt(1/v[3]), log = TRUE) + 
                dnorm(v[1], 0, sqrt(1/t0), log = TRUE) + dnorm(v[2], 0, sqrt(1/t1), log = TRUE) + 
                dgamma(v[3],a, b, log = TRUE))
  }
  
  ans
}


# b) Metropoliis Hastings com proposta qij ~ N(i, s^2)

# proposta
proposta = function(teta_old, var = 0.1){ 
  p = sample(c(1,2,3), 1)
  
  c(ifelse(1 == p, rnorm(1, mean = teta_old[1], sd = sqrt(var)), teta_old[1]), #beta0
    ifelse(2 == p, rnorm(1, mean = teta_old[2], sd = sqrt(var)), teta_old[2]), #beta1
    abs(ifelse(3 == p, rnorm(1, mean = teta_old[3], sd = sqrt(var)), teta_old[3]))) #precisao
}

prob_aceita = function(teta_old, teta_new){ # prob de aceitacao - como estamos trabalhando com o 
  min(exp(posteriori(teta_new) - posteriori(teta_old)), 1)} #log temos que trabalhar com a exp

# algoritmo
n_MCMC = 10000  #n de it do MCMC
MCMC_teta = matrix(NA, nrow = n_MCMC, ncol = 3)
MCMC_teta[1,] = c(0,0,1) #valor inicial 


for(it in 1:(n_MCMC-1)){    ### inicia o algoritmo
  
  teta_old = MCMC_teta[it,]
  teta_new = proposta(teta_old) ## cria uma proposta
  
  if(runif(1) < prob_aceita(teta_old,teta_new)){
    MCMC_teta[it+1,] = teta_new #com prob prob_aceita registra o valor da proposta no MCMC 
  } else {
    MCMC_teta[it+1,] = teta_old #caso contrario repete o valor antigo 
  }

}

require("coda")

X = mcmc(MCMC_teta)
plot(X)
summary(X)
(1 - rejectionRate(X)) # taxa de aceitacao
effectiveSize(X) # tamanho efetivo

# c) Avaliar o efeito de s^2 no modelo e encontrar um valor de s^2 ótimo

valor.var = seq(0.1, 10, 0.1) # 100 valores para testar
var.otimo = matrix(NA, nrow = length(valor.var), ncol = 6) # salvar os resultados


for(i in 1:length(valor.var)){
  
  n_MCMC = 5000  #n de it do MCMC
  MCMC_teta = matrix(NA, nrow = n_MCMC, ncol = 3)
  MCMC_teta[1,] = c(0,0,1) #valor inicial 
  
  
  for(it in 1:(n_MCMC-1)){    ### inicia o algoritmo
    
    teta_old = MCMC_teta[it,]
    teta_new = proposta(teta_old, var = valor.var[i]) ## cria uma proposta
    
    if(runif(1) < prob_aceita(teta_old,teta_new)){
      MCMC_teta[it+1,] = teta_new #com prob prob_aceita registra o valor da proposta no MCMC 
    } else {
      MCMC_teta[it+1,] = teta_old #caso contrario repete o valor antigo 
      }
    
  }
  
  Xa = mcmc(MCMC_teta)
  
  var.otimo[i, 1:3] = (1 - rejectionRate(Xa)) # taxa de aceitacao
  var.otimo[i, 4:6] = effectiveSize(Xa) # tamanho efetivo
  
}

colnames(var.otimo) = c('tx.aceitac_b0', 'tx.aceitac_b1', 'tx.aceitac_prec', 
                     'tam.efetivo_b0', 'tam.efetivo_b1', 'tam.efetivo_prec')


par(mfrow = c(1,2))
plot(valor.var, var.otimo[,1], col = 'red', type = 'l', ylim = c(0, 0.20), 
     ylab = 'taxa de aceitação', xlab = 'variância da proposta')
lines(valor.var, var.otimo[,2], col = 'blue')
lines(valor.var, var.otimo[,3], col = 'green')
lines(valor.var, apply(var.otimo[,1:3], 1, mean), col = 'gray', lty = 2)
legend('topright', 
       legend = c(expression(paste(beta[0])), 
                  expression(paste(beta[1])), 
                  expression(paste(phi)), 'valor médio'), 
       col = c("red", "blue", "green", "gray"), lty = c(1,1,1,2), cex = 0.8)
abline(v = valor.var[which.max(apply(var.otimo[,1:3], 1, mean))], lty = 2)


plot(valor.var, var.otimo[,4], col = 'red', type = 'l', ylim = c(0,300),
     ylab = 'tamanho efetivo', xlab = 'variância da proposta')
lines(valor.var, var.otimo[,5], col = 'blue')
lines(valor.var, var.otimo[,6], col = 'green')
lines(valor.var, apply(var.otimo[,1:3], 1, mean), col = 'gray', lty = 2)
legend('topright', 
       legend = c(expression(paste(beta[0])), 
                  expression(paste(beta[1])), 
                  expression(paste(phi)), 'valor médio'), 
       col = c("red", "blue", "green", "gray"), lty = c(1,1,1,2), cex = 0.8)
abline(v = valor.var[which.max(apply(var.otimo[,4:6], 1, mean))], lty = 2)

# d) estimativas media posteriori do modelo

# aqui uma observacao: precisa rodar novamente a letra b) antes de 
# rodar as estimativas

res = summary(X)
res$statistics[,1]

plot(X)

# e) Encontre as distribui¸c˜oes condicionais dos trˆes parˆametros relevantes do modelo, 
# e monte o algoritmo Gibbs Sampler

n_MCMC = 10000
MCMC_teta = matrix(ncol = n_MCMC, nrow = 3)     
MCMC_teta[,1] = c(0,0,1)

post = vector()               # guardar os valores de posteriori
post[1] = posteriori(MCMC_teta[,1])


x_2 = sum(x^2)  #para economizar cálculo
n = length(x)


for(it in 1:(n_MCMC-1)){    ### inicia o algoritmo
  
  teta_old = MCMC_teta[,it]
  teta_new = teta_old
  
  if(it%%3 == 1){
    s = 1/(t0 + teta_old[3]*n)
    m = teta_old[3]*sum(y-teta_old[2]*x)*s
    teta_new[1] = rnorm(1, m, sqrt(s))
    
  }else if(it%%3 == 2){
    s = 1/(t1 + teta_old[3]*x_2)
    m = teta_old[3]*sum((y-teta_old[1])*x)*s
    teta_new[2] = rnorm(1, m, sqrt(s))
    
  }else{
    beta = b + sum((y - teta_old[1] - teta_old[2]*x)^2)/2
    teta_new[3] = rgamma(1, shape = a + n/2, rate = beta)    
  }
  
  MCMC_teta[,it+1] = teta_new   
  post[it+1] = posteriori(teta_new)
}


X1 = mcmc(t(rbind(MCMC_teta,post)))


# f) Comparação dos dois metodos

res1 = summary(X1)


tab = rbind(c(res$statistics[,1], (1 - rejectionRate(X)) , effectiveSize(X)),
      c(res1$statistics[1:3,1], (1 - rejectionRate(X1)[1:3]) , effectiveSize(X1)[1:3]))

require(xtable)

xtable(tab)
