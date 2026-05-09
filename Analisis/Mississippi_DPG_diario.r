library(evd)
library(evir)

source("../Funciones/DPG_funciones.r")

missi_daily <- read.csv("https://raw.githubusercontent.com/miguel-becares/TFG-distribuciones-de-colas-pesadas/main/Datos/Mississippi_daily.csv", stringsAsFactors = FALSE)

missi_day_dis <- missi_daily$discharge

hist(missi_day_dis,breaks = 40, main="Histograma de datos diarios",
     xlab = "Caudal en metros cúbicos por segundo",ylab = "Frecuencia")

#ver que no haya Na
which(is.na(missi_day_dis))

#1º función media de excesos
evir::meplot(missi_day_dis,omit=100,labels=TRUE,pch = 18,
             main="Función de exceso medio")
abline(v=15000,col = "red",lwd=2)

u_EMV<-15000

#1.2 Ajuste por EMV
missitr_EMV<-missi_day_dis[missi_day_dis>u_EMV] -u_EMV
Nu<-length(missitr_EMV)
n <-length(missi_day_dis)

fit <-evir::gpd(missi_day_dis,threshold=u_EMV,method="ml")
fit$par.ests
xi_EMV<-fit$par.ests[1]
beta_EMV<-fit$par.ests[2]
u_EMV - beta_EMV/xi_EMV


#estudio de los resultados
gpd_qqplot(missitr_EMV,xi=xi_EMV,beta=beta_EMV,
           main="QQ-plot método de EMV u = 10.000")

compare_gpd_cdf(missitr_EMV,xi=xi_EMV,beta=beta_EMV,
                main=paste("CDF empírica vs teórica método de EMV u = ",u_EMV))

compare_gpd_pdf(missitr_EMV,xi=xi_EMV,beta=beta_EMV,
                main="Densidad empírica vs GDP método de EMV")

compare_perc_gpd(missi_day_dis,u_EMV,xi=xi_EMV,beta=beta_EMV,
                 np = 25,xlim=c(.99,1),ylim=c(15000,30000),
                 main =paste("Cuantil teórico vs empírico u = ",u_EMV))

#comparar percentiles
q_gpd_cola(c(.9,.99,.999),missi_day_dis,u = u_EMV,xi=xi_EMV,beta=beta_EMV)
quantile(missi_day_dis,c(.9,.99,.999))

#Pareto generalizada en la cola 
compare_gpd_sobre_cola(missi_day_dis,u=u_EMV,xi=xi_EMV,
                       beta=beta_EMV,breaks=200,xlim= c(0,u_EMV+10000),
                       ylim=c(0,0.0002),main = "")

compare_gpd_sobre_cola(missi_day_dis,u=u_EMV,xi=xi_EMV,
                       beta=beta_EMV,breaks=200,
                       xlim= c(5000,max(missi_day_dis)),
                       ylim=c(0,0.00008),main = "")

