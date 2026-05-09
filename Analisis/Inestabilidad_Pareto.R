library(evir)

Nrep <- 3
part <- seq(0.01, 1, by = 0.01)  
n <- 100000
q<- 0.9 #percentil para establecer el umbral
p <- 0.999 #percentil para calcular después

xi_teo <- 0.75




#umbral a tomar: percentil q teórico.

u_teo <- (1/xi_teo)*((1-q)^(-xi_teo)-1)

cuant_teo <- (1/xi_teo)*((1-p)^(-xi_teo)-1)

colores <- c("black", "blue", "red") 
#si se aumenta Nrep aumentar los colores 

plot(NULL,xlim = range(part),ylim = c(0, 50),
     xlab = "Proporción de muestra",
     ylab = "Error relativo (%)",
     main = paste("Error rel (p =", p, ", xi =", xi_teo, ")"))

for (k in 1:Nrep) {
  data_pareto <- evir::rgpd(n, xi_teo)
  error_rel <- numeric(length(part))
  data_pareto <- evir::rgpd(n, xi_teo)
  error_rel <- numeric(length(part))
  for (i in seq_along(part)) {
    vectp <- head(data_pareto, floor(part[i] * n))
    Nu <- sum(vectp > u_teo)
    N  <- length(vectp)
    fit <- evir::gpd(vectp, threshold = u_teo)
    xi_hat    <- fit$par.ests[1]
    sigma_hat <- fit$par.ests[2]
    cuant_hat <- u_teo +
      (sigma_hat / xi_hat) *((N * (1 - p) / Nu)^(-xi_hat) - 1)
    error_rel[i] <-  100*abs(cuant_hat -cuant_teo)/ cuant_teo
  }
  lines(part, error_rel, col = colores[k], lwd = 1)
}

grid(nx = NA,ny = NULL,lty = 2, col = "gray", lwd = .5)

