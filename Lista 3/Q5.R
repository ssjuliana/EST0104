# gerador dos dados 
norm.generator <- function(sigma, n = 1){
  
  u1 = runif(n); u2 = runif(n)
  z = sqrt(-2*log(u1))*cos(2*pi*u2)
  x = z*sigma
  
  return(x)
}

logret = numeric(100)
gr = numeric(100)

sigma = c(2,4,6,8) # valores de sigma escolhidos

for(t in 1:100){
  u = runif(1)
  aux = cut(u, breaks = c(0, 0.25, 0.5, 0.75, Inf), labels = c(1, 2, 3, 4))
  
  gr[t] = aux
  logret[t] = norm.generator(sigma = sigma[as.numeric(aux)])
}

df = data.frame(logret, gr, t = 1:100)
df$gr = as.factor(df$gr)

plot(df$t, df$logret, col = df$gr, type = 'b', xlab = 'time t', ylab = 'logret',
     main = 'Data generated from a normal distribution \n with 4 distinct sigmas')
legend('bottomright', legend = levels(df$gr), col = 1:4, cex = 0.8, pch = 1, bty = 'n', horiz = TRUE)

# b) 

mod <- depmix(logret ~ 1, family = gaussian(), nstates = 4, data = df)
set.seed(1)
fm2 <- fit(mod, verbose = FALSE)

summary(fm2)
print(fm2)

# Classification (inference task)
probs <- posterior(fm2)             # Compute probability of being in each state
head(probs)
rowSums(head(probs)[,2:5])          # Check that probabilities sum to 1

table(paste0('i', df$gr), paste0('s', probs$state))
