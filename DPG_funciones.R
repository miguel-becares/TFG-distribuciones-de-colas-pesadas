library(evir)

#Funciones

mle_beta <- function(x, xi) {
  
  Nu <- length(x)
  beta_init <- sd(x)
  loglik <- function(beta) {
    if (beta <= 0) return(-Inf)
    term <- 1 + xi * x / beta
    if (any(term <= 0)) return(-Inf)
    ll <- -Nu * log(beta) - (1 + 1/xi) * sum(log(term))
    return(ll)
  }
  
  opt <- optim(
    par = beta_init,
    fn = loglik,
    method = "L-BFGS-B",
    lower = 1e-12,
    control = list(fnscale = -1)
  )
  return(opt$par)
}

gpd_qqplot <-function(x, xi, beta,
                      main = "GPD QQ-plot",pch = 1, 
                      col_teo = "red",col_emp="black"){
  x_sort <- sort(x)
  Nu<- length(x)
  
  p<- (1:Nu)/(Nu+1)
  
  c_teor <- (beta/xi)*((1-p)^(-xi)-1)
  
  plot(x_sort, c_teor,col=col_emp,pch=pch,main=main,
       xlab="Cuantiles teóricos",ylab="Cuantiles empíricos")
  abline(0,1,col=col_teo,lwd=2)
}

gpd_qqplot_escala_log <-function(x, xi, beta,
                           main = "GPD QQ-plot",pch = 1, 
                           xlim = c(1,200),ylim = c(1,200),
                           col_teo = "red",col_emp="black"){
  x_sort <- sort(x)
  Nu<- length(x)
  
  p<- (1:Nu)/(Nu+1)
  c_teor <- (beta/xi)*((1-p)^(-xi)-1)
  
  
  
  plot(1+x_sort,1+ c_teor,col=col_emp,pch=pch,main=main,
       xlab="Cuantiles teóricos escala log",
       ylab="Cuantiles empíricos escala log"
       ,log="xy",xlim=xlim,ylim=ylim)
  abline(0,1,col=col_teo,lwd=2)
}

compare_gpd_cdf <- function(x,xi,beta,
                            xlim=c(0,1.1*max(x)),
                            col_emp = "black", col_gpd = "red",
                            main="CDF Empírica vs CDF Pareto generalizada"){
  Nu<- length(x)
  ecdf_x<- sort(x)
  ecdf_y<- (1:Nu)/Nu
  xg<-seq(xlim[1],xlim[2],length.out=500)
  
  pdpg <-function(z){
    1-(1+xi*z/beta)^(-1/xi)}
  
  plot(ecdf_x,ecdf_y, pch = 16, cex = 0.7,col=col_emp,
       xlim=xlim,ylim=c(0,1),main=main,xlab="Excedencias",ylab="Probabilidad acumulada")
  lines(xg,pdpg(xg),col = col_gpd,lwd=2)
}

compare_gpd_cdf_escala_log <- function(x,xi,beta,
                                 xlim=c(0,1.1*max(x)),
                                 col_emp = "black", col_gpd = "red",
                                 main="CDF Empírica vs CDF Pareto generalizada",
                                 ylim=c(0,1)){
  Nu<- length(x)
  ecdf_x<- sort(x)
  ecdf_y<- (1:Nu)/Nu
  xg<-seq(xlim[1],max(x)*1.2,length.out=2000)
  
  pdpg <-function(z){
    1-(1+xi*z/beta)^(-1/xi)}
  
  plot(1+ecdf_x,ecdf_y, pch = 19, cex = .7,col=col_emp
       ,ylim=c(0,1),main=main,
       xlab="Excedencias escala log",ylab="Probabilidad acumulada",
       log="x")
  lines(1+xg,pdpg(xg),col = col_gpd,lwd=2)
}

compare_gpd_pdf <- function(x, xi, beta,breaks = "FD",
                            xlim=c(min(x),1.1*max(x)),
                            main="Densidad empírica vs GPD") {
  Nu <- length(x)
  xg <- seq(xlim[1],xlim[2], length.out = 500)
  dgpd <- function(t) {
    (1 / beta) * (1 + xi * t / beta)^(-1 / xi - 1)}
  
  hist(x,breaks = breaks,freq=FALSE,col = "grey80",            
       border = "white",xlim = xlim,ylim = c(0,max(dgpd(xg))),
       xlab = "Excedencias",ylab = "Densidad",main = main)
  
  lines(xg, dgpd(xg),col = "red",lwd = 1.5)
}

compare_perc_gpd <- function(data,u,xi,beta,
                             xlim=c(0.99,1),ylim=c(0,300),np =50,
                             main="Cuantil teórico vs empírico"){
  datas <- sort(data)
  Fn <-ecdf(data)
  n<-length(data)
  Nu<-length(data[data>u])
  pg<-seq(xlim[1],.9999 ,length.out=500)
  pg2<-seq(xlim[1],.9999 ,length.out=np)
  percfun<-function(t){
    u+beta/xi*((n/Nu*(1-t))^(-xi)-1)}
  
  plot(pg2,quantile(data,pg2),xlim=xlim,ylim=ylim,
       xlab="Percentil",ylab="Cuantil",main=main)
  
  
  lines(pg,percfun(pg),col="red",lwd=1.5)
}

q_gpd_cola <- function(p,data,u,xi,beta){
  n<-length(data)
  Nu<-length(data[data>u])
  u+beta/xi*((n/Nu*(1-p))^(-xi)-1)
}

compare_gpd_sobre_cola <- function(x,u, xi, beta,breaks = "FD",
                                      xlim=c(min(x),1.1*max(x)),
                                      ylim= c(0,0.1),
                                      main="Densidad empírica vs GPD") {
  n <- length(x)
  Nu<- length(x[x>u])
  beta_u <- beta*((Nu/n)^(xi))
  nu_u <- u-beta/xi*(1-(Nu/n)^(xi))
  xg <- seq(nu_u,xlim[2], length.out = 500)
  print(beta_u)
  print(nu_u)
  dgpd <- function(t) {
    (1 / beta_u) * (1 + xi * (t-nu_u) / beta_u)^(-1 / xi - 1)}
  
  hist(x,breaks = breaks,freq=FALSE,col = "grey80",            
       border = "white",xlim = xlim,ylim = ylim,
       xlab = "Excedencias",ylab = "Densidad",main = main)
  
  lines(xg, dgpd(xg),col = "red",lwd = 1.5)
  abline(v= u,lwd = 2,col ="blue")
}


