source('Data_gene.R')
source('Update_theta.R')
source('SE.R')
###################################################################
######################### Data Generation #########################
###################################################################
############ First block size with no error: 10^5
############ Second block size with error: 10^5
############ Last two thetas have error-in-variables
n1 <- n2 <- 10^5 
theta <- c(1,0.5,0.5,0.2,0.2)
umu = 0; corru = 0.1
SigmaU <- matrix(corru, 2, 2) + diag(1-corru, 2)
d = 4; corr = 0.1; q = d - 2
mu = rep(0,d); Sigma <- matrix(corr, d, d) + diag(1-corr, d); model.error = 1

Dmat = Datageneration(n1, n2, mu, umu, Sigma, SigmaU, model.error, q = q, theta)


######################### Update estimates ########################
#Coef.ZX
#C.coef.ZX
#Cum.C.coef.WS
Result = OnlineMeasearment(Dmat$y1, Dmat$y2, Dmat$Z_1, Dmat$W_1, Dmat$Z_2, Dmat$X_2, Dmat$X_1)

########################### Caculate SE ###########################
SE= SEcase(n1, n2, Result$NewV, Result$C.V, Result$ZW, Dmat$W_1, Dmat$y1, Dmat$X_2, 
                   Result$Coef.ZX, Result$C.coef.ZX, Result$Cum.C.coef.WS, Result$C.Resi, Result$Resi_zx) 
