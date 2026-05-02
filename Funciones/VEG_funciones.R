library(evd)


#Funciones
gev_qqplot <- function(x, mu, sigma, xi,
                       main = "QQ-plot distribución valores extremos generalizada",
                       xlab = "Cuantiles teóricos",
                       ylab = "Cuantiles empíricos",
                       pch = 1, col = "black") {
  
  x_sorted <- sort(x)
  n <- length(x_sorted)
  p <- (1:n) / (n + 1)
  q_theo <- mu + (sigma / xi) * ( (-log(p))^(-xi) - 1 )
  plot(q_theo, x_sorted,main = main,xlab = xlab,
       ylab = ylab,pch = pch,col = col)
  abline(0, 1, col = "red", lwd = 2)
}

compare_gev_cdf_frechet <- function(x, mu, sigma, xi,
                            main = "CDF empírica vs teórica",
                            xlim = c(min(x),max(x)),
                            col_emp = "black", col_gev = "red") {
  
  if (xi <= 0) {
    stop("Esta función es solo para xi > 0")
  }
  xmin <- mu - sigma / xi
  q99 <- mu + (sigma / xi) * ((-log(0.99))^(-xi) - 1)
  
  if (is.null(xlim)) { xlim <- c(xmin, q99) }
  n <- length(x)
  
  ecdf_x <- sort(x)
  ecdf_y <- (1:n) / n
  
  xg <- seq(xlim[1], xlim[2], length.out = 500)
  pgev <- function(z) {
    s <- 1 + xi * (z - mu) / sigma
    return(exp(-s^(-1 / xi)))
  }
  
  plot(ecdf_x, ecdf_y, pch = 16, cex = 0.7, col = col_emp,
       xlim = xlim, ylim = c(0,1),
       xlab = "x", ylab = "Probabilidad acumulada", main = main)
  lines(xg, pgev(xg), col = col_gev, lwd = 2)
}

compare_gev_cdf_weibull <- function(x, mu, sigma, xi,
                                    main = "CDF empírica vs teórica",
                                    xlim = NULL,
                                    col_emp = "black", col_gev = "red") {
  
  if (xi >= 0) {
    stop("Esta función es solo para xi < 0")
  }
  
  xmax <- mu - sigma / xi
  n <- length(x)
  x_sorted <- sort(x)
  ecdf_y <- (1:n) / n
  
  if (is.null(xlim)) {
    rango <- max(x_sorted) - min(x_sorted)
    margen <- 0.1 * rango  
    xlim <- c(min(x_sorted) - margen, max(x_sorted) + margen)
  }
  
  pgev_weibull <- function(z) {
    s <- 1 + xi * (z - mu) / sigma
    return(exp(-s^(-1 / xi)))
  }
  
  xg <- seq(xlim[1], xlim[2], length.out = 500)
  
  plot(x_sorted, ecdf_y,
       pch = 16, cex = 0.7, col = col_emp,
       xlim = xlim, ylim = c(0,1),
       xlab = "x", ylab = "Probabilidad acumulada",
       main = main)
  
  lines(xg, pgev_weibull(xg), col = col_gev, lwd = 2)
}


compare_gev_pdf <- function(x, mu, sigma, xi,breaks = "FD",
                                    main = "PDF empírica vs teórica (Weibull)",
                                    xlim = c(.8*min(x),1.2*max(x)),
                                    col_emp = "black",col_gev = "red") {
  n <- length(x)
  h <- hist(x, breaks = breaks, plot = TRUE, freq =FALSE,xlim = xlim)

  xg <- seq(xlim[1],xlim[2], length.out = 500)
  
  dgev_weibull <- function(z) {
    s <- 1 + xi * (z - mu) / sigma
    return((1 / sigma) *s^(-1/xi - 1) *exp(-s^(-1/xi)))
  }

  lines(xg, dgev_weibull(xg), col = col_gev, lwd = 2)
}


return_lev_grafica <- function(n,mu,sigma,xi, main="Gráfica nivel de retorno"){
  T<- seq(from=2,to=n,by=1)
  cm<-1.05*evd::qgev(1-1/n,loc = mu, scale = sigma, shape = xi)
  R <- 10
  idx <- seq(from=R-min(T)+1, to=n, by = R)
  plot(T[idx],
       evd::qgev(1-1/T[idx],loc = mu, scale = sigma, shape = xi)
       ,xlim=c(0,n),
       ylim=c(0.9*evd::qgev(1-1/2,loc = mu, scale = sigma, shape = xi),cm),
       col="red",xlab="Años",ylab="Nivel de retorno",
       main="Niveles de retorno")
  lines(T,evd::qgev(1-1/T,loc = mu, scale = sigma, shape = xi),
        col = "black",lwd = 1)
}

