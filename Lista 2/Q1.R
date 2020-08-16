f.obj <- function(z, mu, x){
  f = sum(z * (x - mu[1]))^2 + sum((1- z)*(x - mu[2])^2)
  
  return(f)
}

kmeans.func <- function(z,x){
  
  conta = 1
  z = z0
  
  z_old <- rep(0, times = nrow(x))
  
  Z = matrix(z)
  M = matrix(c(0,0))
  
  while(!identical(z, z_old)){ 
    
    mu = c(mean(x[which(z == 1)]), mean(x[which(z == 0)])) 
    z_old <- z  
  
    for(i in 1:nrow(x)) z[i] <- ifelse(f.obj(z = 1, mu, x[i]) < f.obj(z = 0, mu, x[i]), 1, 0)
    
    conta = conta + 1
    
    Z = cbind(Z,z)
    M = cbind(M,mu)
  }

  Z = ifelse(Z == 0, 2, Z)
  z = ifelse(z == 0, 2, z)
  
  it = conta - 1
  rownames(Z) = paste("x", 1:nrow(x), sep = "")
  colnames(Z) = paste("i", 0:it, sep = "")
  
  M <- M[,-1]
  
  rownames(M) <- paste("mu", 1:2, sep = "")
  colnames(M) <- paste("i", 1:it, sep = "")
  
  list("grupo" = list("z.arg.ot" = z, "z.seq.ot" = Z), 
       "media" = list("mu.arg.ot" = M, "mu.seq.ot" = mu),
       "it" = it)  
}


set.seed(123)
x = cbind(c(rnorm(45, mean = 5), rnorm(55, mean = 8)))

require(ggplot2)
ggplot() + geom_density(aes(x)) + theme_bw()

z0 = sample(c(0,1), length(x), replace = TRUE) # valor inicial
kmeans.func(z0, x)

