library(truncnorm)


source("../Funciones/DPG_funciones.r")
source("../Funciones/TLN_funciones.r")


ping_data <- read.csv("https://raw.githubusercontent.com/miguel-becares/TFG-distribuciones-de-colas-pesadas/main/Datos/ping.csv", stringsAsFactors = FALSE)

ping <- ping_data$tiempo_ms


#histogramas en escala logaritmica
hist(log(ping),breaks=50,ylab="Frecuencia",xlab="log(ping)",
     main="Histograma de log(ping)")
hist(log(ping),breaks=50,ylim=c(0,500),ylab="Frecuencia",xlab="log(ping)",
     main="Histograma de la cola de log(ping)"  )


#Prueba con Pareto generalizada.
hillfit <-evir::hill(ping,option="xi",ci=0,start=8,end=0.2*length(ping))
#evir::meplot(ping,omit=1,labels=TRUE,main="Función de exceso medio")
abline(v=4000)
k_HILL=4000

xi_HILL<-hillfit$y[which(hillfit$x==k_HILL)]
u_HILL <-sort(ping,decreasing=TRUE)[k_HILL]
pingtr <- ping[ping>=u_HILL]-u_HILL
beta_HILL<-mle_beta(pingtr,xi_HILL)

#estudio de los resultados
gpd_qqplot(pingtr,xi=xi_HILL,beta=beta_HILL,
           main="QQ-plot ping Pareto generalizada")

gpd_qqplot_escala_log(pingtr,xi=xi_HILL,beta=beta_HILL, 
                      xlim=c(1,3000),ylim = c(1,3000),
                      main="QQ-plot ping Pareto generalizada")

compare_gpd_cdf(pingtr,xi=xi_HILL,beta=beta_HILL,
                main=paste("CDF empírica vs teórica método de EMV u = ",u_HILL))

compare_gpd_cdf_escala_log(pingtr,xi=xi_HILL,beta=beta_HILL,
                main=paste("CDF empírica vs teórica método de EMV u = ",u_HILL))

compare_perc_gpd(pingtr,u_HILL,xi=xi_HILL,beta=beta_HILL,
                 np = 25,xlim=c(.99,.9999),ylim=c(100,4000))

#comparar cuantiles
q_gpd_cola(c(.9,.99,.999,.9999),pingtr,u = u_HILL,xi=xi_HILL,beta=beta_HILL)
quantile(ping,c(.9,.99,.999,.9999))


#Ajuste inadecuado, probamos con Log-normal Truncada

#Selección de umbral
hist(ping,breaks=2000,ylim=c(0,2000),xlim = c(12,200),ylab="Frecuencia",
     main="Histograma de ping ")
abline(v=44,col="red",lwd = 2)


u_TLN <- 44

lping <- log(ping)
pingtr <-ping[ping>u_TLN]
lpingtr <-log(pingtr)

n<-length(ping)
Nu<- length(pingtr)
1-Nu/n


init <- c(mean(lpingtr), sd(lpingtr))
fit <- optim(par = init,fn = loglik_truncnorm,x = lpingtr,
             u = log(u_TLN), #log(u) pq has hecho la transformación Y=log(X)
             method = "L-BFGS-B",lower = c(-Inf, 1e-6),hessian = TRUE)
fit$par
mu_TLN    <- fit$par[1]
sigma_TLN <- fit$par[2]

#estudio de los resultados
qq_truncnorm(lping,log(u_TLN), mu=mu_TLN,sigma=sigma_TLN,r=10,
             xlim=c(3.5,9),ylim=c(3,9),
             main=paste("QQ-plot Normal truncada"))

compare_tn_cdf(lpingtr,log(u_TLN),mu=mu_TLN,sigma=sigma_TLN,
               xlim=c(log(u_TLN),8),ylim=c(0,1))



compare_perc_tln(ping,u = u_TLN, mu = mu_TLN, sigma = sigma_TLN, np=40,
                 xlim=c(1-length(pingtr)/length(ping),.999),ylim=c(40,400))


#comparar cuantiles
q_tln_cola(c(.975,.99,.999),ping,u = u_TLN, mu = mu_TLN, sigma = sigma_TLN)
quantile(ping,c(.975,.99,.999))

