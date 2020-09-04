# gerador

X = function(a,b,u) a + b*u

n = 500
u = runif(n)

# primeiramente vamos considerar valores arbitrários para 

a.arb = 2 
b.arb = 4

# U(a - 1/2 ; a + 1/2)

a = a.arb - (1/2)
b = 1 

# aplicando 
x = X(a,b,u)

# comparando para ver se se adequa
hist(x, freq = FALSE, main = paste0('U(', a.arb - 1/2, ',', a.arb + 1/2,')'), xlim = c(a.arb - 1/2, a.arb + 1))
lines(sort(x), dunif(sort(x), a.arb - 0.5, a.arb + 0.5), col = 'red', lty = 2)
legend(a.arb + 0.55, 0.5, legend = c("Simulation", "Uniform"), col = c("black", "red"), lty = 1:2, cex = 0.8)
text(a.arb + 0.8, 1.0, paste0('ks test \np-value: ', round(ks.test(x,"punif", a.arb - 0.5, a.arb + 0.5)$p.value, 4)))

# U(0 ; b)

a = 0
b = b.arb 

# aplicando 
x = X(a,b,u)

# comparando para ver se se adequa
hist(x, freq = FALSE, main = paste0('U(', 0, ',', b.arb,')'), xlim = c(0, b.arb + 2))
lines(sort(x), dunif(sort(x), 0, b.arb), col = 'red', lty = 2)
legend(b.arb + 0.1, 0.15, legend = c("Simulation", "Uniform"), col = c("black", "red"), lty = 1:2, cex = 0.8)
text(b.arb + 1.1, 0.25, paste0('ks test \np-value: ', round(ks.test(x,"punif", 0, b.arb)$p.value, 4)))


# U(a ; b)

a = a.arb
b = b.arb - a.arb 

# aplicando 
x = X(a,b,u)

# comparando para ver se se adequa
hist(x, freq = FALSE, main = paste0('U(', a.arb, ',', b.arb,')'), xlim = c(a.arb, b.arb + 1))
lines(sort(x), dunif(sort(x), a.arb, b.arb), col = 'red', lty = 2)
legend(b.arb+0.1, 0.2, legend = c("Simulation", "Uniform"), col = c("black", "red"), lty = 1:2, cex = 0.8)
text(b.arb+0.6, 0.5, paste0('ks test \np-value: ', round(ks.test(x,"punif", a.arb, b.arb)$p.value, 4)))
