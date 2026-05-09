library(evir)

source("../Funciones/DPG_funciones.r")


#datos
data(danish)

#análisis explorativo
hist(log(danish) ,breaks =50,
     xlab="Millones de DKK",ylab="Frecuencia")

#1º Método de estimación: Función media de excesos y EMV

u_EMV<-0.999

danishtr_EMV<-danish[danish>u_EMV] -u_EMV
n<-length(danish)
k_EMV<-length(danishtr_EMV)
1-k_EMV/length(danish)

fit <-evir::gpd(danish,threshold=u_EMV,method="ml")
fit$par.ests
xi_EMV<-fit$par.ests[1]
beta_EMV<-fit$par.ests[2]

#estudio de los resultados
gpd_qqplot(danishtr_EMV,xi=xi_EMV,beta=beta_EMV,
           main="QQ-plot método de EMV")

gpd_qqplot_escala_log(danishtr_EMV,xi=xi_EMV,beta=beta_EMV,
                      main="QQ-plot método de EMV escala logarítmica",
                      xlim=c(1,250),ylim=c(1,250))
compare_gpd_cdf(danishtr_EMV,xi=xi_EMV,beta=beta_EMV,
                main="CDF empírica vs teórica método de EMV")

compare_gpd_cdf_escala_log(danishtr_EMV,xi=xi_EMV,beta=beta_EMV,
                           main="CDF empírica vs teórica método de EMV",
                           xlim=c(0,6))

compare_gpd_pdf(danishtr_EMV,xi=xi_EMV,beta=beta_EMV,
                main="Densidad empírica vs GDP método de EMV")

compare_perc_gpd(danish,u_EMV,xi=xi_EMV,beta=beta_EMV,
                 np = 25,xlim=c(.99,1))

#comparar percentiles

q_gpd_cola(c(.9,.99,.999),danish,u = u_EMV,xi=xi_EMV,beta=beta_EMV)
quantile(danish,c(.9,.99,.999))

#2º método: Hill y EMV

hillfit <-evir::hill(danish,option="xi",ci=0,start=8,end=length(danish))
abline(v=1350)
abline(v=1450)
abline(v=1400,col = "red")
k=1400
n<- length(danish)
1-k/n

xi_HILL=hillfit$y[which(hillfit$x==k)]
u_HILL=sort(danish,decreasing=TRUE)[k]
danishtr_HILL <- danish[danish>u_HILL]-u_HILL
beta_HILL<-mle_beta(danishtr_HILL,xi_HILL)

#estudio de los resultados
gpd_qqplot(danishtr_HILL,xi=xi_HILL,beta=beta_HILL,main="QQ-plot método de HILL")


gpd_qqplot_escala_log(danishtr_HILL,xi=xi_HILL,beta=beta_HILL,
                      main="QQ-plot método de HILL escala logarítmica",
                      xlim=c(1,330),ylim=c(1,330))

compare_gpd_cdf_escala_log(danishtr_HILL,xi=xi_HILL,beta=beta_HILL,
                           main="CDF empírica vs teórica método de HILL",
                           xlim=c(0,6))

compare_perc_gpd(danish,u_HILL,xi=xi_HILL,beta=beta_HILL,
                 np = 25,xlim=c(.99,.9999),ylim=c(0,400))

#comparar percentiles
q_gpd_cola(c(.9,.99,.999),danish,u = u_HILL,xi=xi_HILL,beta=beta_HILL)
quantile(danish,c(.9,.99,.999))


#Obtener xi mediante la pendiente de la regresión lineal
#en la función de exceso medio
umbrales<-MEF$threshold[which((MEF$threshold>0.99) & (MEF$threshold<=100))]
excesos_medios<-MEF$meanExcess[which((MEF$threshold>0.99) & (MEF$threshold<=100))]

reg_lin<-lm(excesos_medios~umbrales)
reg_lin$coefficients
m<-reg_lin$coefficients[2]
xi_pendiente <- m/(1+m)

evir::meplot(danish,xlim=c(1,100),omit=1,labels=TRUE,
             main="Función de exceso medio con regresión lineal")
abline(reg_lin$coefficients[1],reg_lin$coefficients[2],
       col="red",lwd=1.5)




