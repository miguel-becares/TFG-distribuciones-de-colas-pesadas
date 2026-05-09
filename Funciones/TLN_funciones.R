library(truncnorm)

#Funciones

loglik_truncnorm <- function(par, x, u){
  mu <- par[1]
  sigma <- par[2]
  
  if (!is.finite(mu) || !is.finite(sigma) || sigma <= 0)
    return(1e12)
  
  z <- (u - mu) / sigma
  
  # evitar colas numéricamente imposibles
  if (z > 8) return(1e12)
  
  k <- length(x)
  
  ll <- 
    -k * log(sigma) -
    sum((x - mu)^2) / (2 * sigma^2) -
    k * pnorm(z, log.p = TRUE, lower.tail = FALSE)
  
  if (!is.finite(ll)) return(1e12)
  
  -ll
}

compare_tn_cdf <- function(data, u, mu, sigma,
                           xlim=c(min(data),max(data)),ylim =x(0,1),
                           main = "Comparación CDF: Empírica vs Normal Truncada") {
  #Hay que meter log(ping) y log(u) como parámetros
  sdata <- sort(data)
  F_emp <- ecdf(sdata)
  xg <- seq(xlim[1],xlim[2], length.out =500)
  
  plot(sdata, F_emp(sdata),pch = 16, cex = 0.7,,col = "black",xlab = "x",
       ylab = "CDF",xlim=xlim,ylim = ylim,main = main)
  lines(xg,truncnorm::ptruncnorm(xg,a = u,b = Inf,
                                 mean = mu, sd = sigma),,col = "red",lwd = 2)
}

qq_truncnorm <-function(data,u,mu,sigma, r=ceiling(length(data)/2000),
                        main="QQ-plot Normal truncada en log(u)",
                        xlim=c(u,max(data)),ylim=c(3,8),
                        col_teo = "red",col_emp="black"){
  #Hay que meter log(ping) y log(u) como parámetros
  #r es para  representar solo uno de cada r puntos.
  
  datatr <- data[data>=u]
  sdatatr <- sort(datatr)
  n<-length(datatr)
  p <- (1:n - 0.5)/n
  
  idx <- seq(1,n, by = r)
  cteo <- truncnorm::qtruncnorm(p, a = u,b = Inf,mean = mu,sd = sigma)
  
  plot(cteo[idx],sdatatr[idx],
       xlim=xlim,ylim=ylim, pch = 19, col = col_emp,
       xlab="Cuantiles teóricos",ylab="Cuantiles empíricos",main=main)
  
  abline(0,1,col=col_teo,lwd=2)
}

q_tln_cola <- function(p,data,u,mu,sigma){
  n <- length(data)
  datatr <- data[data>u]
  Nu<-length(datatr)
  Fu<- plnorm(u,meanlog = mu,sdlog=sigma)
  s<- Fu+(1-Fu)*n/Nu*(p-1+Nu/n)
  qlnorm(s, meanlog = mu,sdlog=sigma)
}

compare_perc_tln<- function(data,u,mu,sigma,
                            xlim=c(.95,1),ylim=c(0,5000),np = 50,
                            main="Comparar percentil TLN"){
  n <- length(data)
  sdata <- sort(data)
  
  datatr <- data[data>u]
  sdatatr <- sort(datatr)
  Nu<-length(datatr)
  Fu<- plnorm(u,meanlog = mu,sdlog=sigma)
  xg <- seq(1-Nu/n, 0.99999,length.out = 2000)
  xg2<- seq(1-Nu/n, 0.99999,length.out = np)
  Fn <- ecdf(sdata)
  perc_tln<- function(t){
    s<- Fu+(1-Fu)*n/Nu*(t-1+Nu/n)
    qlnorm(s, meanlog = mu,sdlog=sigma)} 
  plot(xg2, quantile(data,xg2), xlim=xlim, ylim=ylim, main = main,
       xlab= "Percentil", ylab = "Cuantil", col = "black")
  lines(xg, perc_tln(xg),col = "red",lwd = 1.5)
}



