library(evd)


#Funciones


return_lev_grafica_ini <- function(n,mu,sigma,xi, main="Gráfica nivel de retorno"){
  T<- seq(from=2,to=n,by=1)
  cuantil<- mu+(sigma/xi)*((-log(1-1/T))^(-xi)-1)
  cm<-1.05*(mu+(sigma/xi)*((-log(1-1/n))^(-xi)-1))
  R <- 10
  idx <- seq(from=R-min(T)+1, to=n, by = R)
  plot(T[idx],cuantil[idx],xlim=c(0,n),ylim=c(0.9*min(cuantil),cm),col="red",
       xlab="Años",ylab="Nivel de retorno",
       main="Niveles de retorno")
  lines(T,cuantil,col = "black",lwd = 1)
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

return_lev_grafica2 <- function(data, mu, sigma, xi, step = 5,
                                main = "GEV teórica vs cuantiles empíricos") {
  
  n <- length(data)
  
  T_theo <- seq(2, n, length.out = 500)
  
  q_theo <- mu + (sigma / xi) * ((-log(1 - 1 / T_theo))^(-xi) - 1)
  
  plot(T_theo, q_theo,type = "l",
       col = "black",lwd = 2,xlab = "Periodo de retorno (años)",
       ylab = "Nivel de retorno",main = main)
  
  T_emp <- seq(5, n, by = step)
  
  probs <- 1 - 1 / T_emp
  
  q_emp <- quantile(data, probs = probs, type = 7)
  
  points(T_emp, q_emp,col = "red",pch = 19)
  
}

gev_qqplot <- function(x, mu, sigma, xi,
                       main = "QQ-plot distribución valores extremos generalizada",
                       xlab = "Cuantiles teóricos",
                       ylab = "Cuantiles empíricos",
                       pch = 1, col = "black") {
  
  
  # Ordenar datos
  x_sorted <- sort(x)
  n <- length(x_sorted)
  
  # Probabilidades empíricas
  p <- (1:n) / (n + 1)
  #p <- (1:n - 0.5)/n
  # Cuantiles teóricos GEV 
  q_theo <- mu + (sigma / xi) * ( (-log(p))^(-xi) - 1 )
  
  # QQ-plot
  plot(q_theo, x_sorted,main = main,xlab = xlab,
       ylab = ylab,pch = pch,col = col)
  
  abline(0, 1, col = "red", lwd = 2)
}

compare_gev_pdf <- function(x, mu, sigma, xi,breaks = "FD",
                            main ="PDF empírica vs teórica",
                            xlim = c(min(x),max(x)),ylim = c(0,0.01)) {
  # soporte teórico
  xmin <- mu - sigma / xi
  
  # cuantil 0.99
  q99 <- mu + (sigma / xi) * ((-log(0.99))^(-xi) - 1)
  
  # xlim por defecto
  if (is.null(xlim)) {xlim <- c(xmin, q99)}
  
  # filtrar datos fuera de soporte
  t <- 1 + xi * (x - mu) / sigma
  x <- x[t > 0]
  n <- length(x)
  
  # histograma empírico (sin plot)
  h <- hist(x,breaks = breaks,plot = FALSE)
  
  mids <- h$mids
  widths <- diff(h$breaks)
  dens_emp <- h$counts / (n * widths)
  
  # grid GEV
  xg <- seq(0, 300, length.out = 500)
  
  dgev <- function(z) {
    s <- 1 + xi * (z - mu) / sigma
    (1 / sigma) *s^(-1 / xi - 1) *exp(-s^(-1 / xi))
  }
  
  plot(mids, dens_emp,
       pch = 16,cex = .7,col = "black",xlim = xlim,ylim=ylim,
       xlab = "x",ylab = "density",main = main)
  
  lines(xg, dgev(xg),col = "red",lwd = 2)
}

compare_gev_cdf_frechet <- function(x, mu, sigma, xi,
                            main = "CDF empírica vs teórica",
                            xlim = c(min(x),max(x)),
                            col_emp = "black", col_gev = "red") {
  
  # soporte teórico
  xmin <- mu - sigma / xi
  # cuantil 0.99
  q99 <- mu + (sigma / xi) * ((-log(0.99))^(-xi) - 1)
  
  # xlim por defecto
  if (is.null(xlim)) { xlim <- c(xmin, q99) }
  
  # filtrar datos fuera de soporte
  t <- 1 + xi * (x - mu) / sigma
  x <- x[t > 0]
  n <- length(x)
  
  # CDF empírica
  ecdf_x <- sort(x)
  ecdf_y <- (1:n) / n
  
  # grid GEV
  xg <- seq(xlim[1], xlim[2], length.out = 500)
  pgev <- function(z) {
    s <- 1 + xi * (z - mu) / sigma
    exp(-s^(-1 / xi))
  }
  
  # plot
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
  
  # --- soporte ---
  xmax <- mu - sigma / xi
  
  # --- filtrar datos fuera de soporte ---
  t <- 1 + xi * (x - mu) / sigma
  x <- x[t > 0]
  n <- length(x)
  
  # --- ordenar datos ---
  x_sorted <- sort(x)
  
  # --- CDF empírica ---
  ecdf_y <- (1:n) / n
  
  # --- xlim con margen ---
  if (is.null(xlim)) {
    xmin_data <- min(x_sorted)
    xmax_data <- max(x_sorted)
    rango <- xmax_data - xmin_data
    
    margen <- 0.1 * rango  # 5% a cada lado
    xlim <- c(xmin_data - margen, xmax_data + margen)
  }
  
  # --- función GEV ---
  pgev_weibull <- function(z) {
    s <- 1 + xi * (z - mu) / sigma
    
    out <- numeric(length(z))
    
    valid <- s > 0
    out[valid] <- exp(-s[valid]^(-1 / xi))
    
    out[s <= 0] <- 1
    
    return(out)
  }
  
  xg <- seq(xlim[1], xlim[2], length.out = 500)
  
  plot(x_sorted, ecdf_y,
       pch = 16, cex = 0.7, col = col_emp,
       xlim = xlim, ylim = c(0,1),
       xlab = "x", ylab = "Probabilidad acumulada",
       main = main)
  
  lines(xg, pgev_weibull(xg), col = col_gev, lwd = 2)
  
}
