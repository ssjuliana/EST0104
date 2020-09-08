# questão 1 

# Indique o algoritmo para construir geradores de numeros aleatorios (partindo de n´umeros com 
# dist uniforme [0,1]) para

n = 1000
u = runif(n)

# a) PARETO

a = 3
b = 5

x = a/(1-u)^(1/b) 

require(EnvStats)

hist(x, freq = FALSE, main = paste0('Pareto Distribution with \n location = ', a, ' and shape = ', b))
lines(sort(x), dpareto(sort(x), location = a, shape = b), col = 'red', lty = 2)
legend(8, 0.5, legend = c("Simulation", "Pareto dist."), col = c("black", "red"), lty = 1:2, cex = 0.8)
text(11, 0.2, paste0('ks test p-value: ', round(ks.test(x,"ppareto", location = a, shape = b)$p.value, 4)))

# b) Gumbel padrão 

x = -log(-log(u))

require(evd)

hist(x, freq = FALSE, main = 'Standard Gumbel Distribution', xlim = c(-4, 10), ylim = c(0, 0.5))
lines(sort(x), dgumbel(sort(x)), col = 'red', lty = 2)
legend(3, 0.3, legend = c("Simulation", "Std. Gumbel dist."), col = c("black", "red"), lty = 1:2, cex = 0.8)
text(6.8, 0.1, paste0('ks test p-value: ', round(ks.test(x,"pgumbel")$p.value, 4)))

# c) Distribuição F 

dist.f = function(d1, d2){
  # utilizando box-muller para gerar as normal padrão 
  
  qui1 = matrix(NA, ncol = d1, nrow = n)
  for(i in 1:d1){
    u1 = runif(n/2); u2 = runif(n/2)
    x1 = sqrt(-2*log(u1))*cos(2*pi*u2); x2 = sqrt(-2*log(u1))*sin(2*pi*u2)
    
    qui1[,i] = c(x1,x2)^2 # normal padrão  
  }
  
  qui2 = matrix(NA, ncol = d2, nrow = n)
  for(j in 1:d2){
    u1 = runif(n/2); u2 = runif(n/2)
    x1 = sqrt(-2*log(u1))*cos(2*pi*u2); x2 = sqrt(-2*log(u1))*sin(2*pi*u2)
    
    qui2[,j] = c(x1,x2)^2 # normal padrão  
  }
  
  X = apply(qui1, 1, sum) #qui df = d1 
  Y = apply(qui2, 1, sum) #qui df = d2 
  
  distF = (X/d1)/(Y/d2)

  return(distF)
}

d1 = 5
d2 = 10

x = dist.f(d1, d2)

hist(x, freq = FALSE, main = paste0('F Distribution with \ndf1 = ', d1, ' and df2 = ', d2), ylim = c(0, 0.7))
lines(sort(x), df(sort(x), df1 = d1, df2 = d2), col = 'red', lty = 2)
legend(6, 0.4, legend = c("Simulation", "F dist."), col = c("black", "red"), lty = 1:2, cex = 0.8)
text(9, 0.1, paste0('ks test p-value: ', round(ks.test(x,"pf", df1 = d1, df2 = d2)$p.value, 4)))


# Binomial Negativa

conta = function(v.assum,x){
  z = v.assum
  nz = length(z)
  res = rep(0,nz)
  names(res) = paste(z)
  for(i in 1:nz){
    res[i] = sum(x==z[i])
  }
  return(res)
}

#teste de aderência xquadrado
aderencia = function(freq.obs, freq.esp){
  k = length(freq.obs)
  est.teste = sum(((freq.obs - freq.esp)^2)/freq.esp)
  valor.p = pchisq(q = est.teste, df= (k - 1),lower.tail = FALSE)
  list("Estatistica.teste" = est.teste, "valor.p" = valor.p)
}

#gerador de distribuições discretas
gna.disc3 = function(x, prob.fmp, n){
  y = rep(0,n)
  probs = cumsum(prob.fmp)
  i = 0
  while(i < n){
    u = runif(1)
    pos = sum(u > probs) + 1
    if(pos <= length(x)){
      i = i+1
      y[i] = x[pos]
    }
  }
  y
}


# massa da binomial negativa
fmp.bineg = function(x,r,p){
  choose(x + r - 1, r - 1)*(p^r)*((1 - p)^(x))
} 

lim = 0.999985
p = 0.3 #escolhendo p 
r = 10 # escolhendo r

k = 0
cc = 0

while(cc < lim){
  cc = sum(fmp.bineg(0:k, r, p))
  k = k + 1
}

fda.bineg = cumsum(fmp.bineg(0:k, r, p))
fda.bineg = fda.bineg[fda.bineg <= lim]

z = 0:length(fda.bineg)
ns = 1000
mger = gna.disc3(z,fmp.bineg(z,r,p),ns)
rger = rnbinom(ns,r,p)


#comparando o gerador criado com o gerador do R
mtab = conta(z,mger)
rtab = conta(z,rger)

#teste de aderência do qui quadrado
freq.e = fmp.bineg(z,r,p)*ns
freq.o = mtab

# resultado
plot(mtab, type = "h", main = paste0('Negative Binomial (r = ', r, ', p = ',p, ')'), ylab = '')
lines(rtab, type = "h", ylab = '', col = 'red', lty = 2)
legend(40, 30, legend = c("Simulation", "Neg. Binon dist."), col = c("black", "red"), lty = 1:2, cex = 0.8)
text(60, 10, paste0('X² test p-value: ', round(aderencia(freq.o,freq.e)$valor.p, 4)))


# fmp da distribuição bivariada 

fmp.bivar = function(x,y,n,p){
  choose(n,x)*choose(n,y)*((x^y*(n-x)^(n-y)*p^x*(1-p)^(n-x))/(n^n)) 
}


p = 0.3 #escolhendo p 
n = 10 #escolhendo n

# nesse caso, como x e y vao ate n 
# entao nao precisamos definir um limite sup pra alcançar
# e sim teremos as combinações de 0 até n tanto em x quanto em y

comb = expand.grid(0:n,0:n); colnames(comb) = c("x","y")

fmp.prob = fmp.bivar(comb[1:nrow(comb),'x'], comb[1:nrow(comb),'y'], n, p)
fda.bivar = cumsum(fmp.prob)

z = 0:length(fda.bivar)
ns = 1000
mger = gna.disc3(z,fmp.bivar(z,z,n,p),ns)
mtab = conta(z,mger)