library(evd)
library(evir)

source("../Funciones/VEG_funciones.r")

# Datos
missi_anual <- read.csv("https://raw.githubusercontent.com/miguel-becares/TFG-distribuciones-de-colas-pesadas/main/Datos/Mississippi_max_anual.csv", stringsAsFactors = FALSE)

mis_max_an <- missi_anual$flow_m3s

# Histograma
hist(mis_max_an,breaks =18,main = "Histograma de máximos anuales",
     xlab = "Caudal en metros cúbicos por segundo",
     ylab = "Fracuencia")

#Ajuste VEG
mis_fit <-evir::gev(mis_max_an)
mis_fit$par.ests
xi   <- mis_fit$par.ests[1]
sigma <- mis_fit$par.ests[2]
mu   <- mis_fit$par.ests[3]

#extremo superior del soporte
mu-sigma/xi

#estudio de los resultados
gev_qqplot(mis_max_an,mu=mu,sigma=sigma,xi=xi)

compare_gev_cdf_weibull(mis_max_an,mu=mu,sigma=sigma,xi=xi)

compare_gev_pdf(mis_max_an,mu=mu,sigma=sigma,xi=xi)

return_lev_grafica(n=length(mis_max_an),mu,sigma,xi)

#niveles de retorno a 5,10,25 años:
evd::qgev(1-1/5,loc = mu, scale = sigma, shape = xi)
evd::qgev(1-1/10,loc = mu, scale = sigma, shape = xi)
evd::qgev(1-1/25,loc = mu, scale = sigma, shape = xi)

#nivel de retorno asociado a la máxima observación
s<-evd::pgev(max(mis_max_an),loc = mu, scale = sigma, shape = xi)
1/(1-s)
