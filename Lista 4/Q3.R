n_MCMC = 2500
n_MCMC1 = 2501

funcao_p = function(x){  ## função proporcional a densidade que queremos amostrar (aqui X é o vetor de dados)
  (x[1]+x[2]+x[3])^(5)
}

proposta = function(xi){    # aqui xi eh o elemento que queremos amostrar 
  candidatos = c(0:9)
  candidatos = candidatos[-which(candidatos == xi)]#lista posíveis mudanças
  sample(candidatos,1)         # escolhe um elemento da lista de candidatos com dist uniforme
}

#   p_proposta = 1/9   - todas as peopostas tem a mesma prob. 
prob_aceita = function(x_old,x_new){
  min(funcao_p(x_new)/funcao_p(x_old), 1)
}

MCMC_d = matrix(ncol = n_MCMC, nrow = 3)     #declara a matriz que conterá os valores dos dígitos ao longo do MCMC
#cada coluna eh uma ieração do MCMC

p = vector() #declara vetor que vai guardar resultado da função densidade para cada observação

MCMC_d[,1] = c(0,0,0) #escolhe valores iniciais

it = 1
while(it < (n_MCMC-1)){  ### inicia o algoritmo
  
  for (i in 1:3){         ### para cada um dos dígitos do número
    x_old = MCMC_d[,it]
    x_new = x_old
    x_new[i] = proposta(x_old[i])                ## cria uma proposta
    
    if(runif(1) < prob_aceita(x_old,x_new)){
      MCMC_d[,it+1] = x_new                         #com probabilidade prob_aceita registra o valor da proposta no MCMV 
    }else{
      MCMC_d[,it+1] = x_old                         #caso contrario repete o valor antigo 
    }
    it = it + 1
  }
}


#######################################   Observar os resultados

plot(MCMC_d[1,-c(1:n_MCMC/10)],type="l")
hist(MCMC_d[1,-c(1:n_MCMC/10)])

#################  Estimar E(X)
X = MCMC_d[1,]*100 + MCMC_d[2,]*10 + MCMC_d[3,]  #montando os numeros em função dos algarismos

plot(X)

hist(X, nclass = 100)                         ## Olhar a distribuição de X
hist(X, nclass = 500)

mean(X, na.rm = TRUE)
sd(X, na.rm = TRUE)


## Vizualização dinamica
plot(MCMC_d[1,1:30], MCMC_d[2,1:30], xlim = c(0,9), ylim = c(0,9), xlab = "unidade", ylab = "dezena", type = "b")
plot(MCMC_d[1,1:300], MCMC_d[2,1:300], xlim = c(0,9), ylim = c(0,9), xlab = "unidade", ylab = "dezena", type = "b")
plot(MCMC_d[1,1:3000], MCMC_d[2,1:3000], xlim = c(0,9), ylim = c(0,9), xlab = "unidade", ylab = "dezena" ,type="b")

require(coda)
resultado = mcmc(na.omit(t(MCMC_d)))

# modificada da aula

funcao_p = function(x){  ## função proporcional a densidade que queremos amostrar (aqui X é o vetor de dados)
  (x[1] + x[2] + x[3] + x[4] + x[5])^(5)
}

proposta = function(xi){    # aqui xi eh o elemento que queremos amostrar 
  candidatos = c(0:9)
  candidatos = candidatos[-which(candidatos == xi)]#lista posíveis mudanças
  sample(candidatos,1)         # escolhe um elemento da lista de candidatos com dist uniforme
}

#   p_proposta = 1/9   - todas as peopostas tem a mesma prob. 
prob_aceita = function(x_old,x_new){
  min(funcao_p(x_new)/funcao_p(x_old), 1)
}

MCMC_d = matrix(ncol = n_MCMC1, nrow = 5)     #declara a matriz que conterá os valores dos dígitos ao longo do MCMC
#cada coluna eh uma ieração do MCMC

p = vector() #declara vetor que vai guardar resultado da função densidade para cada observação

MCMC_d[,1] = c(0,0,0,0,0) #escolhe valores iniciais

it = 1

aux = n_MCMC1-1
while(it < aux){  ### inicia o algoritmo
  
  for (i in 1:5){         ### para cada um dos dígitos do número
    x_old = MCMC_d[,it]
    x_new = x_old
    x_new[i] = proposta(x_old[i])                ## cria uma proposta
    
    if(runif(1) < prob_aceita(x_old,x_new)){
      MCMC_d[,it+1] = x_new                         #com probabilidade prob_aceita registra o valor da proposta no MCMV 
    }else{
      MCMC_d[,it+1] = x_old                         #caso contrario repete o valor antigo 
    }
    it = it + 1
  }
}

#######################################   Observar os resultados

plot(MCMC_d[1,-c(1:n_MCMC1/10)],type="l")
hist(MCMC_d[1,-c(1:n_MCMC1/10)])

#################  Estimar E(X)
X = MCMC_d[1,]*10000 + MCMC_d[2,]*1000 + MCMC_d[3,]*100 + MCMC_d[4,]*10 + MCMC_d[5,]

plot(X)

hist(X, nclass = 100)                         ## Olhar a distribuição de X
hist(X, nclass = 500)

mean(X, na.rm = TRUE)
sd(X, na.rm = TRUE)

require(coda)
resultado2 = mcmc(na.omit(t(MCMC_d)))

effectiveSize(resultado)
effectiveSize(resultado2)