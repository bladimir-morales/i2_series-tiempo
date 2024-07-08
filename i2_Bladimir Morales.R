rm(list=ls())

source('./Interrogatorios/i2/casa/functions.R')
library(tidyverse)
library(forecast)
library(tseries)
library(urca)
library(xtable)

bd <- read.table('./Interrogatorios/i2/casa/smo66b_raw.txt', header = T)
names(bd)

#write.table(bd, './Interrogatorios/i2/casa/INAM002C_raw.txt', row.names = F)

(start <- min(bd$age))
max(bd$age)
# Ajuste ----
X <- bd$smo66b_raw
X <- ts(X, start = start, frequency = 1)
n <- length(X)
t <- 1:n
plot(X, xlab = 'Tiempo', ylab = 'mm', col = 'black', lwd = 2,
     cex.lab = 5, cex.axis = 2)

pdf('./Interrogatorios/i2/casa/fig/serie.pdf', width = 8, height = 6)
autoplot(X, series="X")+
  ggtitle("")+
  geom_line(size=1.2, color = '#2B7DFF')+
  xlab("Año") + ylab("mm")+
  theme(legend.position="none")+
  theme(axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))
dev.off()

# Prueba de Dickey-Fuller Aumentada (ADF)
# H0: La serie no es estacionaria (tiene una raíz unitaria)
# H1: La serie es estacionaria
adf.test(X)
adf.test(X, k = 15)
# Prueba de KPSS
# H0: La serie es estacionaria.
# H1: La serie no es estacionaria.
kpss.test(X)

ndiffs(X, test = "kpss")
ndiffs(X, test = "adf")
ndiffs(X, test = "pp")

pdf('./Interrogatorios/i2/casa/fig/acfX.pdf', width = 8, height = 6)
acf(X, lag.max = 25, main = '', cex.axis=1.5, cex.lab= 1.7,
    xlab = 'Rezagos', ylim= c(-.3,1.1))
dev.off()
pdf('./Interrogatorios/i2/casa/fig/pacfX.pdf', width = 8, height = 6)
pacf(X, lag.max = 25, main = '', cex.axis=1.5, cex.lab= 1.7,
     xlab = 'Rezagos', ylab = 'PACF', ylim= c(-.3,1.1))
dev.off()

# M1 - AR(3)----
fit <- Arima(X - mean(X), order = c(3,0,0), include.mean = F)
summary(fit)
e <- fit$residuals
hatX <- fitted(fit)
## Significancia 
Z0 <- coef(fit)/sqrt(diag(fit$var.coef))
2*(1-pnorm(abs(Z0)))# ar2 no significativo

## Estacionariedad 
# Deben tener decaimiento exponencial para ser estacionarias
pdf('./Interrogatorios/i2/casa/fig/acf_e1.pdf', width = 8, height = 6)
acf(e,main = '', cex.axis=1.5, cex.lab= 1.7,
    xlab = 'Rezagos', ylim= c(-.3,1.1))
dev.off()
pdf('./Interrogatorios/i2/casa/fig/pacf_e1.pdf', width = 8, height = 6)
pacf(e,main = '', cex.axis=1.5, cex.lab= 1.7,
     xlab = 'Rezagos', ylim= c(-.3,1.1), ylab = 'PACF')
dev.off()
# Prueba de Dickey-Fuller Aumentada (ADF)
# H0: La serie no es estacionaria (tiene una raíz unitaria)
# H1: La serie es estacionaria
adf.test(e)
adf.test(e, k = 10)
ndiffs(e, test = 'adf')
# Prueba de KPSS
# H0: La serie es estacionaria.
# H1: La serie no es estacionaria.
kpss.test(e)
ndiffs(e, test = 'kpss')

## Invertivibilidad
plot(fit)
# Obtener las raíces del polinomio AR
roots <- polyroot(c(1, -phi))
# Verificar si los módulos de las raíces son mayores que uno
stable <- all(Mod(roots) > 1)
print(stable)#True es estable
## Blancura 
LSTS::Box.Ljung.Test(e, lag = 20, main = '')

## Homocedasticidad
plot(e)
plot(hatX, e, pch = 19, cex = .8,
     xlab = expression(hat(Y)))
abline(h = 0, col = 'red', lty = 2)
# H0: Homocedasticidad
# H1: Heterocedasticidad
lmtest::bptest(e ~ time(e))
lmtest::bptest(e ~ hatX)

## Normalidad
# H0: X es normal
# H1: X no es normal
nortest::lillie.test(e)
qqnorm(e, pch = 16, cex = .8)
qqline(e, col = 'red', lwd = 2)

# M2 - ARMA(3,3) ----
fit <- Arima(X - mean(X), order = c(3,0,3), include.mean = F)
summary(fit)
e <- fit$residuals
hatX <- fitted(fit)
## Significancia
phi <- coef(fit)[coef(fit)!=0]
Z0 <- phi/sqrt(diag(fit$var.coef))
2*(1-pnorm(abs(Z0)))

## Estacionariedad
# Deben tener decaimiento exponencial para ser estacionarias
acf(e)
pacf(e)
# Prueba de Dickey-Fuller Aumentada (ADF)
# H0: La serie no es estacionaria (tiene una raíz unitaria)
# H1: La serie es estacionaria
adf.test(e)
adf.test(e, k = 10)
ndiffs(e, test = 'adf')
# Prueba de KPSS
# H0: La serie es estacionaria.
# H1: La serie no es estacionaria.
kpss.test(e)
ndiffs(e, test = 'kpss')

## Invertivibilidad
plot(fit)
# Obtener las raíces del polinomio AR
roots <- polyroot(c(1, -phi))
# Verificar si los módulos de las raíces son mayores que uno
stable <- all(Mod(roots) > 1)
print(stable)#True es estable
## Blancura
LSTS::Box.Ljung.Test(e, lag = 20, main = '')

## Homocedasticidad
plot(e)
plot(hatX, e, pch = 19, cex = .8,
     xlab = expression(hat(Y)))
abline(h = 0, col = 'red', lty = 2)
# H0: Homocedasticidad
# H1: Heterocedasticidad
lmtest::bptest(e ~ time(e))
lmtest::bptest(e ~ hatX)

## Normalidad
# H0: X es normal
# H1: X no es normal
nortest::lillie.test(e)
qqnorm(e, pch = 16, cex = .8)
qqline(e, col = 'red', lwd = 2)


# Transformación por heterocedasticidad ----
(lambda <- BoxCox.lambda(X, method = 'loglik'))
(lambda <- BoxCox.lambda(X, method = 'guerrero'))

MASS::boxcox(X ~ time(X))
Y <- BoxCox(X, lambda = lambda)
#Y <- log(X)
plot(Y-mean(Y))

ndiffs(Y, test= 'kpss')
ndiffs(Y, test = 'adf')
ndiffs(Y, test = 'pp')

pdf('./Interrogatorios/i2/casa/fig/acfY.pdf', width = 8, height = 6)
acf(Y, lag.max = 25, main = '', cex.axis=1.5, cex.lab= 1.7,
    xlab = 'Rezagos', ylim= c(-.3,1.1))
dev.off()
pdf('./Interrogatorios/i2/casa/fig/pacfY.pdf', width = 8, height = 6)
pacf(Y, lag.max = 25, main = '', cex.axis=1.5, cex.lab= 1.7,
     xlab = 'Rezagos', ylab = 'PACF', ylim= c(-.3,1.1))
dev.off()
#LSTS::ts.diag(Y)

# M3- AR(3) ----
fit <- Arima(Y-mean(Y), order = c(3,0,0), include.mean = F)
summary(fit)
e <- fit$residuals
hatY <- fitted(fit)
## Significancia
Z0 <- coef(fit)/sqrt(diag(fit$var.coef))
pvalor <- round(2*(1-pnorm(abs(Z0))),4) # ar2 no significativo


(tabla <- data.frame('Parámetro'= c('phi_1', 'phi_2', 'phi_3'),
           'Estimador' = round(phi,3),
           'Z_0' = round(Z0,3),
           'p-valor'= pvalor))

## Estacionariedad
# Deben tener decaimiento exponencial para ser estacionarias
acf(e)
pacf(e)
# Prueba de Dickey-Fuller Aumentada (ADF)
# H0: La serie no es estacionaria (tiene una raíz unitaria)
# H1: La serie es estacionaria
adf.test(e)
adf.test(e, k = 10)
ndiffs(e, test = 'adf')
# Prueba de KPSS
# H0: La serie es estacionaria.
# H1: La serie no es estacionaria.
kpss.test(e)
ndiffs(e, test = 'kpss')

## Invertivibilidad
plot(fit)
# Obtener las raíces del polinomio AR
roots <- polyroot(c(1, -phi))
# Verificar si los módulos de las raíces son mayores que uno
stable <- all(Mod(roots) > 1)
print(stable)#True es estable
## Blancura
LSTS::Box.Ljung.Test(e, lag = 20, main = '')

LSTS::ts.diag(e)

## Homocedasticidad
plot(e)
plot(hatX, e, pch = 19, cex = .8,
     xlab = expression(hat(Y)))
abline(h = 0, col = 'red', lty = 2)
# H0: Homocedasticidad
# H1: Heterocedasticidad
lmtest::bptest(e ~ time(e))
lmtest::bptest(e ~ hatX)

## Normalidad
# H0: X es normal
# H1: X no es normal
nortest::lillie.test(e)
qqnorm(e, pch = 16, cex = .8)
qqline(e, col = 'red', lwd = 2)


# Modelo 3 fijo 
fixed = c(NA,0,NA)
fit <- Arima(Y-mean(Y), order = c(3,0,0), include.mean = F, fixed = fixed)
summary(fit)
e <- fit$residuals
hatY <- fitted(fit)
## Significancia
phi = coef(fit)[coef(fit)!=0]
Z0 <- phi/sqrt(diag(fit$var.coef))
pvalor <- round(2*(1-pnorm(abs(Z0))),4) # ar2 no significativo

(tabla1 <- data.frame('Estimador' = round(coef(fit),3),
                     'Z_0' = round(c(Z0[1],NA,Z0[2]),3),
                     'p-valor'= c(pvalor[1],NA,pvalor[2])))

xtable(cbind(tabla,tabla1))

# M4 - ARMA (3,3) ----
fixed <- c(NA, 0, NA)
fit <- Arima(Y-mean(Y), order = c(3,0,0), include.mean = F, fixed = fixed)
summary(fit)
e <- fit$residuals
hatY <- fitted(fit)
## Significancia
phi <- coef(fit)[coef(fit)!=0]
Z0 <- phi/sqrt(diag(fit$var.coef))
round(2*(1-pnorm(abs(Z0))),3)

## Estacionariedad
# Deben tener decaimiento exponencial para ser estacionarias
acf(e, lag.max = 15)
pacf(e, lag.max = 15)
# Prueba de Dickey-Fuller Aumentada (ADF)
# H0: La serie no es estacionaria (tiene una raíz unitaria)
# H1: La serie es estacionaria
adf.test(e)
adf.test(e, k = 10)
# Prueba de KPSS
# H0: La serie es estacionaria.
# H1: La serie no es estacionaria.
kpss.test(e)

## Invertivibilidad
plot(fit)
# Obtener las raíces del polinomio AR
roots <- polyroot(c(1, phi))
# Verificar si los módulos de las raíces son mayores que uno
stable <- all(Mod(roots) > 1)
print(stable)#True es estable
## Blancura
LSTS::Box.Ljung.Test(e, lag = 20, main = '')

## Homocedasticidad
plot(e)
plot(hatY, e, pch = 19, cex = .8,
     xlab = expression(hat(Y)))
abline(h = 0, col = 'red', lty = 2)
# H0: Homocedasticidad
# H1: Heterocedasticidad
lmtest::bptest(e ~ time(e))
lmtest::bptest(e ~ hatY)

## Normalidad
# H0: X es normal
# H1: X no es normal
nortest::lillie.test(e)
qqnorm(e, pch = 16, cex = .8)
qqline(e, col = 'red', lwd = 2)

# Estimación ----
## Durbin levinson ----
# AR(3) ----
fit_dl.loglik <- nlminb(start = c(.1,0,.1), objective = dl.loglik, 
                        serie = Y - mean(Y), p = 3, q = 0)
(phi_dl <- fit_dl.loglik$par)
fit_dl <- forecast::Arima(Y-mean(Y), order = c(3,0,0), fixed = phi_dl, include.mean = F)
summary(fit_dl)

# ACF
pdf('./Interrogatorios/i2/casa/fig/acf_dl_ar3.pdf', width = 8, height = 6)
acf(Y, lag.max = 20, main = '', cex.axis=1.5, cex.lab= 1.7,
    xlab = 'Rezagos', ylim= c(-.3,1.1))
ACF <- ARMAacf(ar = phi_dl, lag.max = 20)
lines(ACF ~ c(0:20), type = "p", col = "red", pch =20)
w <- Bartlett(ar = phi_dl, lag.max = 19)
LI <- ACF - qnorm(0.975) * sqrt(c(0,w)/n)
LS <- ACF + qnorm(0.975) * sqrt(c(0,w)/n)
lines(LI ~ c(0:20), lty = 2, col = "red")
lines(LS ~ c(0:20), lty = 2, col = "red")
dev.off()
## Periodograma 
aux <- LSTS::periodogram(Y, plot = F)
pdf('./Interrogatorios/i2/casa/fig/periodograma_dl_ar3.pdf', width = 8, height = 6)
plot(aux$periodogram ~ aux$lambda, type = "l", 
     xlab = expression(lambda), ylab = '',
     main = '', cex.axis=1.5, cex.lab= 1.5, col = 'gray', lwd = 3)
lines(LSTS::spectral.density(ar = phi_dl, sd = sqrt(.3517), 
                             lambda = aux$lambda) ~ aux$lambda, col = 'red',
      lwd = 2)
legend('topright', legend = c('Periodograma','D.Espec. estimada'), 
       col = c('gray','red'), lwd = 2)
dev.off()
## Ajuste en serie
pdf('./Interrogatorios/i2/casa/fig/ajuste_ar3.pdf', width = 8, height = 6)
plot(Y-mean(Y), col = 'gray', xlab = 'Tiempo', ylab = expression(Y-mu), lwd = 2,
     cex.lab = 5, cex.axis = 2)
lines(fit_dl$fitted, col = 'red', lwd = 2)
legend('bottomleft', legend = c('Y centrada', 'Ajuste'),
       col = c('gray','red'), lty = 1, lwd = 2)
dev.off()


## Análisis de residuos ----
fit <- Arima(Y-mean(Y), order = c(3,0,0), include.mean = F)
summary(fit)
## Significancia
Z0 <- phi_dl/sqrt(diag(fit$var.coef))
pvalor <- round(2*(1-pnorm(abs(Z0))),4) # ar2 no significativo

(tabla <- data.frame('Parámetro'= c('phi_1', 'phi_2', 'phi_3'),
                     'Estimador' = round(phi_dl,3),
                     'Z_0' = round(Z0,3),
                     'p-valor'= pvalor))

xtable(tabla)

e <- fit_dl$residuals
Yest <- fit_dl$fitted
## Estacionariedad
# Deben tener decaimiento exponencial para ser estacionarias
pdf('./Interrogatorios/i2/casa/fig/acf_res_ar3.pdf', width = 8, height = 6)
acf(e, lag.max = 25, main = '', cex.axis=1.5, cex.lab= 1.7,
    xlab = 'Rezagos', ylim= c(-.3,1.1))
dev.off()
pacf(e)
# Prueba de Dickey-Fuller Aumentada (ADF)
# H0: La serie no es estacionaria (tiene una raíz unitaria)
# H1: La serie es estacionaria
adf.test(e)
adf.test(e, k = 10)
ndiffs(e, test = 'adf')
# Prueba de KPSS
# H0: La serie es estacionaria.
# H1: La serie no es estacionaria.
kpss.test(e)
ndiffs(e, test = 'kpss')

## Invertivibilidad
pdf('./Interrogatorios/i2/casa/fig/inv_ar3.pdf', width = 8, height = 6)
plot(fit, main = '', xlab = 'Real', ylab = 'Imaginarios', cex.lab = 1.5)
dev.off()
# Obtener las raíces del polinomio AR
roots <- polyroot(c(1, -phi))
# Verificar si los módulos de las raíces son mayores que uno
stable <- all(Mod(roots) > 1)
print(stable)#True es estable
## Blancura
pdf('./Interrogatorios/i2/casa/fig/bl_ar3.pdf', width = 8, height = 6)
LSTS::Box.Ljung.Test(e, lag = 20, main = ' ')
dev.off()

## Homocedasticidad
pdf('./Interrogatorios/i2/casa/fig/res_ar3.pdf', width = 8, height = 6)
plot(e, xlab = 'Tiempo', ylab = 'Residuos', col = 'black', lwd = 2,
     cex.lab = 5, cex.axis = 2)
dev.off()
plot(hatX, e, pch = 19, cex = .8,
     xlab = expression(hat(Y)), cex.axis = 1.5, cex.lab = 1.5)
abline(h = 0, col = 'red', lty = 2)
# H0: Homocedasticidad
# H1: Heterocedasticidad
lmtest::bptest(e ~ time(e))
lmtest::bptest(e ~ hatX)

## Normalidad
# H0: X es normal
# H1: X no es normal
nortest::lillie.test(e)
pdf('./Interrogatorios/i2/casa/fig/norm_ar3.pdf', width = 8, height = 6)
qqnorm(e, pch = 16, cex = .8, xlab = 'Teóricos', ylab = 'Empíricos',
       cex.lab=1.5, cex.axis=1.5, main = '')
qqline(e, col = 'red', lwd = 2)
dev.off()

# ARMA(3,3) ----
fit_dl.loglik <- nlminb(start = rep(.1,6), objective = dl.loglik, 
                        serie = Y - mean(Y), p = 3, q = 3)
(phi_dl <- fit_dl.loglik$par)
fit_dl <- forecast::Arima(Y-mean(Y), order = c(3,0,3), fixed = phi_dl, include.mean = F)

# ACF
pdf('./Interrogatorios/i2/casa/fig/acf_dl_arma33.pdf', width = 8, height = 6)
acf(Y, lag.max = 20, main = '', cex.axis=1.5, cex.lab= 1.7,
    xlab = 'Rezagos', ylim= c(-.3,1.1))
ACF <- ARMAacf(ar = phi_dl[1:3], ma =phi_dl[4:6], lag.max = 20)
lines(ACF ~ c(0:20), type = "p", col = "red", pch =20)
w <- Bartlett(ar = phi_dl, lag.max = 20)
LI <- ACF - qnorm(0.975) * sqrt(c(0,abs(w))/n)
LS <- ACF + qnorm(0.975) * sqrt(c(0,w)/n)
lines(LI ~ c(0:20), lty = 2, col = "red")
lines(LS ~ c(0:20), lty = 2, col = "red")
dev.off()
## Periodograma 
aux <- LSTS::periodogram(Y, plot = F)
pdf('./Interrogatorios/i2/casa/fig/periodograma_dl_arma33.pdf', width = 8, height = 6)
plot(aux$periodogram ~ aux$lambda, type = "l", 
     xlab = expression(lambda), ylab = '',
     main = '', cex.axis=1.5, cex.lab= 1.5, col = 'gray', lwd = 3)
lines(LSTS::spectral.density(ar = phi_dl[1:3], ma =phi_dl[4:6], sd = sqrt(.3517), 
                             lambda = aux$lambda) ~ aux$lambda, col = 'red',
      lwd = 2)
legend('topright', legend = c('Periodograma','D.Espec. estimada'), 
       col = c('gray','red'), lwd = 2)
dev.off()
## Ajuste en serie
pdf('./Interrogatorios/i2/casa/fig/ajuste_arma33.pdf', width = 8, height = 6)
plot(Y-mean(Y), col = 'gray', xlab = 'Tiempo', ylab = expression(Y-mu), lwd = 2,
     cex.lab = 5, cex.axis = 2)
lines(fit_dl$fitted, col = 'red', lwd = 2)
legend('bottomleft', legend = c('Y centrada', 'Ajuste'),
       col = c('gray','red'), lty = 1, lwd = 2)
dev.off()

## Análisis de residuos ----
fit <- Arima(Y-mean(Y), order = c(3,0,3), include.mean = F)
summary(fit)
## Significancia
Z0 <- phi_dl/sqrt(diag(fit$var.coef))
pvalor <- round(2*(1-pnorm(abs(Z0))),3) 

(tabla <- data.frame('Parámetro'= c('phi_1', 'phi_2', 'phi_3','theta_1','theta_2','theta_3'),
                     'Estimador' = round(phi_dl,3),
                     'Z_0' = round(Z0,3),
                     'p-valor'= pvalor))
xtable(tabla)

e <- fit_dl$residuals
Yest <- fit_dl$fitted
## Estacionariedad
# Deben tener decaimiento exponencial para ser estacionarias
pdf('./Interrogatorios/i2/casa/fig/acf_res_arma33.pdf', width = 8, height = 6)
acf(e, lag.max = 25, main = '', cex.axis=1.5, cex.lab= 1.7,
    xlab = 'Rezagos', ylim= c(-.3,1.1))
dev.off()
pacf(e)
# Prueba de Dickey-Fuller Aumentada (ADF)
# H0: La serie no es estacionaria (tiene una raíz unitaria)
# H1: La serie es estacionaria
adf.test(e)
adf.test(e, k = 10)
ndiffs(e, test = 'adf')
# Prueba de KPSS
# H0: La serie es estacionaria.
# H1: La serie no es estacionaria.
kpss.test(e)
ndiffs(e, test = 'kpss')

## Invertivibilidad
pdf('./Interrogatorios/i2/casa/fig/inv_arma33.pdf', width = 8, height = 6)
plot(fit, main = '', xlab = 'Real', ylab = 'Imaginarios', cex.lab = 1.5)
dev.off()
# Obtener las raíces del polinomio AR
roots <- polyroot(c(1, -phi))
# Verificar si los módulos de las raíces son mayores que uno
stable <- all(Mod(roots) > 1)
print(stable)#True es estable
## Blancura
pdf('./Interrogatorios/i2/casa/fig/bl_arma33.pdf', width = 8, height = 6)
LSTS::Box.Ljung.Test(e, lag = 20, main = ' ')
dev.off()

## Homocedasticidad
pdf('./Interrogatorios/i2/casa/fig/res_arma33.pdf', width = 8, height = 6)
plot(e, xlab = 'Tiempo', ylab = 'Residuos', col = 'black', lwd = 2,
     cex.lab = 5, cex.axis = 2)
dev.off()
plot(hatX, e, pch = 19, cex = .8,
     xlab = expression(hat(Y)), cex.axis = 1.5, cex.lab = 1.5)
abline(h = 0, col = 'red', lty = 2)
# H0: Homocedasticidad
# H1: Heterocedasticidad
lmtest::bptest(e ~ time(e))
lmtest::bptest(e ~ hatX)

## Normalidad
# H0: X es normal
# H1: X no es normal
nortest::lillie.test(e)
pdf('./Interrogatorios/i2/casa/fig/norm_arma3.pdf', width = 8, height = 6)
qqnorm(e, pch = 16, cex = .8, xlab = 'Teóricos', ylab = 'Empíricos',
       cex.lab=1.5, cex.axis=1.5, main = '')
qqline(e, col = 'red', lwd = 2)
dev.off()



# Whittle ----
## AR (3) ----
phi_w <- nlminb(objective = w.loglik, serie = Y-mean(Y), 
                start = c(.1,0,.1,1), p = 3, q = 0)$par
phi_w
fit_w <- forecast::Arima(Y-mean(Y), order = c(3,0,0), fixed = phi_w[1:3], 
                         include.mean = F)

fit <- forecast::Arima(Y-mean(Y), order = c(3,0,0),include.mean = F)
summary(fit_w)
names(summary(fit_w))

e <- fit_w$residuals
hatY <- fitted(fit_w)
## Significancia
Z0 <- coef(fit_w)/sqrt(diag(fit$var.coef))
pvalor <- round(2*(1-pnorm(abs(Z0))),3) # ar2 no significativo


(tabla <- data.frame('Parámetro'= c('phi_1', 'phi_2', 'phi_3','sigma'),
                     'Estimador' = round(phi_w,3),
                     'Z_0' = c(round(Z0,3),NA),
                     'p-valor'= c(pvalor,NA))
  )
xtable(tabla)
# ACF
acf(Y, lag.max = 20)
ACF <- ARMAacf(ar = phi_w[1:3], lag.max = 20)
lines(ACF ~ c(0:20), type = "p", col = "red", pch =20)
w <- Bartlett(ar = phi, lag.max = 10)
LI <- ACF - qnorm(0.975) * sqrt(c(0,w)/n)
LS <- ACF + qnorm(0.975) * sqrt(c(0,w)/n)
lines(LI ~ c(0:20), lty = 2, col = "red")
lines(LS ~ c(0:20), lty = 2, col = "red")

## Periodograma 
aux <- LSTS::periodogram(Y, plot = F)
plot(aux$periodogram ~ aux$lambda, type = "l")
lines(LSTS::spectral.density(ar = phi_w[1:3], sd = phi_w[4], 
                             lambda = aux$lambda) ~ aux$lambda, col = 'red')

## Ajuste en serie
plot(Y-mean(Y), col = 'gray')
lines(fit_w$fitted, col = 'red')
  
# Comparación de modelos ----
Y_naive <- lag(as.numeric(Y), 1)
Y_naive[1] <- 0

me_dl <- ME(Y, fit_dl$fitted)
rmse_dl <- RMSE(Y, fit_dl$fitted)
mae_dl <- MAE(Y, fit_dl$fitted)
mpe_dl <- MPE(Y, fit_dl$fitted)
mape_dl <- MAPE(Y, fit_dl$fitted)
mase_dl <- MASE(as.numeric(Y), fit_dl$fitted, Y_naive)

me_w <- ME(Y, fit_w$fitted)
rmse_w <- RMSE(Y, fit_w$fitted)
mae_w <- MAE(Y, fit_w$fitted)
mpe_w <- MPE(Y, fit_w$fitted)
mape_w <- MAPE(Y, fit_w$fitted)
mase_w <- MASE(as.numeric(Y), fit_w$fitted, Y_naive)

tabla <- data.frame('Medida' = c('ME','RMSE','MAE','MPE','MAPE','MASE'),
           'Durbin-Levinson' = round(c(me_dl, rmse_dl, mae_dl, mpe_dl, mape_dl, mase_dl),3),
           'Whittle' = round(c(me_w, rmse_w, mae_w, mpe_w, mape_w, mase_w),3)
           )

xtable(tabla)

######################################################
modelo <- fit
Y <- ts(c(X, 0.90), start = start(X))
fit2 <- forecast::Arima(Y-mean(X), order = c(3,0,0), fixed = modelo$coef, include.mean = F)
pre <- forecast::forecast(fit2, h = 10)
pre$x <- pre$x+mean(X)
pre$mean <- pre$mean+mean(X)
pre$upper <- pre$upper+mean(X)
pre$lower <- pre$lower+mean(X)
plot(pre)

c(rep(.2,5),rep(.5,10),rep(.6,5))
