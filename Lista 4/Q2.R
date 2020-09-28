n_MCMC = 5000

funcao_p = function(x){  
  return(exp(-(x[1]^2/2)-(x[2]^2/2)))
}

proposta = function(x){
  x1 = x[1] + sample(-1:1, 1) 
  x2 = x[2] + sample(-1:1, 1) 
  
  return(c(x1,x2))     
}

prob_aceita = function(x_old,x_new){
  min(funcao_p(x_new)/funcao_p(x_old), 1)
}

MCMC_d = matrix(ncol = n_MCMC, nrow = 2)
MCMC_d[,1] = c(0,0)

it = 1
while(it < (n_MCMC-1)){ 
  
    x_old = MCMC_d[,it]
    x_new = x_old
    x_new = proposta(x_old) 
    
    if(runif(1) < prob_aceita(x_old,x_new)){
      MCMC_d[,it+1] = x_new 
    }else{
      MCMC_d[,it+1] = x_old
    }
    it = it + 1
}

rowMeans(MCMC_d, na.rm = TRUE)
apply(MCMC_d, 1, var, na.rm = TRUE)

# Observar os resultados
require(coda)
rownames(MCMC_d) = c('x1', 'x2')
resultado = mcmc(na.omit(t(MCMC_d)))

summary(resultado)
plot(resultado)

d = data.frame(na.omit(t(MCMC_d)))

require(plotly)
p <- plot_ly(d, x = ~x1, y = ~x2)

add_histogram2dcontour(p)

effectiveSize(resultado)
1-rejectionRate(resultado)
