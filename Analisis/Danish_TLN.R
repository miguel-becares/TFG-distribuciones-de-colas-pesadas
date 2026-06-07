library(evir)
library(evd)

source("../Funciones/TLN_funciones.r")

data(danish)

hist(log(danish) ,breaks =50,
     xlab="Millones de DKK",ylab="Frecuencia")

u_TLN <- exp(0.7)

danishtr <-danish[danish>u_TLN]
ldanishtr <-log(danishtr)
ldanish <- log(danish)
1-length(danishtr)/length(danish)

init <- c(mean(ldanishtr),sd(ldanishtr))

fit <- optim(par = init,fn = loglik_truncnorm,x = ldanishtr,
             u = log(u_TLN), 
             method = "L-BFGS-B",lower = c(-Inf, 1e-6),hessian = TRUE)

fit$par

mu_TLN    <- fit$par[1]
sigma_TLN <- fit$par[2]

qq_truncnorm(ldanish,log(u_TLN), mu=mu_TLN,sigma=sigma_TLN,r=1,
             xlim=c(0.7,6),ylim=c(0.7,6),
             main=paste("QQ-plot Normal truncada"))


compare_tn_cdf(ldanishtr,log(u_TLN),mu=mu_TLN,sigma=sigma_TLN,
               xlim=c(log(u_TLN),6),ylim=c(0,1))


compare_perc_tln(danish,u = u_TLN, mu = mu_TLN, sigma = sigma_TLN, np=80,
                 xlim=c(.9,.99),ylim=c(0,40))

q_tln_cola(c(.9,.95,.99),danish,u = u_TLN, mu = mu_TLN, sigma = sigma_TLN)
quantile(danish,c(.9,.95,.99))



