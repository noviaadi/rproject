library(copula)
library(fCopulae)
library(QRM)
library(MASS)

#(a) Memplot Penutupan harga saham
MIT <- read.csv("MIT.csv",header=TRUE) 
MIT.Close <- MIT$Close
plot(MIT.Close,type="l") 

HMC <- read.csv("HMC.csv",header=TRUE) 
HMC.Close <- HMC$Close
plot(HMC.Close,type="l") 

#(b)Menghitung logreturn
RM <- diff(log(MIT.Close),lag=1) # menghitung log return
head(RM) # melihat enam data pertama
tail(RM) # melihat enam data terakhir
RH <- diff(log(HMC.Close),lag=1) # menghitung log return
head(RH) # melihat enam data pertama
tail(RH) # melihat enam data terakhir
RMH <- cbind(RM,RH) # mengabung nilai log return NVS & PFE
RMH <- as.matrix(RMH) #
n <- nrow(RMH)
head(RMH) # mengecek data 
tail(RMH)


splom2(RMH[,1:2], cex=0.4, col.mat="blue")
U <- pobs(RMH[,1:2]) #amatan semu
splom2(U, cex=0.4,)
library(scatterplot3d)
scatterplot3d(RM,RH, color = "blue",pch = ".")

hist(RMH,breaks="Scott",probability=TRUE,col="#f2b722") # breaks Scott

#(c)
library(fBasics)
basicStats(RMH)

#(d)
rop <- cor(RMH, method="pearson") 
rop
rop <- rop[1,2]
rop
tau <- Kendall(RMH) # bisa juga cor(returns,method="kendall")
tau
tau <- tau[1,2] 
tau
ros <- cor(RMH, method="spearman") 
ros
ros <- ros[1,2]
ros

#no 1.3

# Estimasi Marginal

est.normal1 <- fitdistr(RM, "normal") # distribusi normal
est.normal1

AIC(est.normal1) #Nilai AIC
BIC(est.normal1) #Nilai BIC

est.normal2 <- fitdistr(RH, "normal")
est.normal2

AIC(est.normal2) #Nilai AIC 
BIC(est.normal2) #Nilai BIC 

## Margin Distribusi Student *t*

#est.t1 <- fitdistr(RM, "t") # distribusi t
#est.t1

#AIC(est.t1) #Nilai AIC
#BIC(est.t1) #Nilai BIC 

#est.t2 <- fitdistr(RT, "t")
#est.t2

#AIC(est.t2) #Nilai AIC
#BIC(est.t2) #Nilai BIC
par.Gumbel <- 1/(1-tau) #parameter kopula Gumbel
par.Gumbel
par.Clayton<-(2*tau)/(1-tau) #parameter kopula Clayton
par.Clayton

# Menghitung Amatan Semu (*Pseudo Observations*)
U <- pobs(RMH)
head(U)
tail(U)
# Kopula normal (Gauss)
NormCop <- normalCopula(rop, dim = 2, dispstr = "un")
NormCop # 
NormCopEst <- fitCopula(NormCop,U, method="mpl")
NormCopEst
# Kopula t
tCop <- tCopula(rop, dim = 2, dispstr="un", df.fixed=FALSE)
tCop
tCopEst <- fitCopula(tCop, U, method="mpl", estimate.variance=TRUE)
tCopEst
# Kopula GUmbel
GumbCop <- gumbelCopula(par.Gumbel, dim =2)
GumbCop
GumbCopEst <- fitCopula(GumbCop,U, method="mpl")
GumbCopEst
# Kopula Clayton
ClaytonCop <- claytonCopula(param =par.Clayton, dim = 2)
ClaytonCop
ClaytonCopEst <- fitCopula(ClaytonCop,U, method="mpl")
ClaytonCopEst
# Kopula Frank
FrankCop <- frankCopula(param = NA_real_, dim = 2)
FrankCop
FrankCopEst <- fitCopula(FrankCop,U,method="mpl")
FrankCopEst

set.seed(2106) # set.seed agar bisa direproduksi
r<-100000

## Simulasi Kopula Gauss

normcop.obj <- normalCopula(0.1737, dim = 2, dispstr = "un")
sim.return.normcop<-rCopula(r,normcop.obj)# Nilai simulasi

## Simulasi Kopula Student *t*

tcop.obj <- tCopula(0.1694, dim = 2, dispstr = "un", df=23.5537)
sim.return.tcop<-rCopula(r,tcop.obj) # Simulated values


## Simulasi Kopula Gumbel

##  simulasi kopula Gumbel}
Gumbelcop.obj<- gumbelCopula(1.093, dim =2)    
sim.return.Gumbelcop <- rCopula(r, Gumbelcop.obj)


## Simulasi Kopula Clayton

#  simulasi kopula Clayton}
Claytoncop.obj<- claytonCopula(0.2358, dim = 2)
sim.return.ClaytonCop<-rCopula(r,Claytoncop.obj)# Nilai Simulasi


## Simulasi kopula Frank
#  simulasi kopula Frank}
Frankcop.obj<- frankCopula(param = 1.087 , dim = 2)
sim.retun.Frankcop<-rCopula(r,Frankcop.obj)# Nilai Simulasi

# Simulasi Bobot
#A : 0.7 ;0,3
wa <- cbind(0.7,0.3) # bobot masing-masing risiko
#B
wb<- cbind(0.5,0.5)

#Kopula Gauss Distribusi Normal
#  kopula Gauss margin normal skenario A
sim.RM <- qnorm(sim.return.normcop[,1], mean=est.normal1$estimate[[1]], sd=est.normal1$estimate[[2]])
sim.RH <- qnorm(sim.return.normcop[,2], mean=est.normal2$estimate[[1]], sd=est.normal2$estimate[[2]])
MC.data <- cbind(sim.RM, sim.RH)
MC.Lsim.A <- -(MC.data%*%t(wa))
MC.Lsim.B <- -(MC.data%*%t(wb))

#Menghitung VaR dan TVaR
alfa <- c(0.95, 0.99,0.995) # nilai-nilai var dan tvar yang akan dihitung
#VaR
quantile(MC.Lsim.A, alfa) #skenario A
quantile(MC.Lsim.B, alfa) #skenario B
#TVaR skenario A
mean(MC.Lsim.A[MC.Lsim.A > quantile(MC.Lsim.A, alfa[1])])
mean(MC.Lsim.A[MC.Lsim.A > quantile(MC.Lsim.A, alfa[2])])
mean(MC.Lsim.A[MC.Lsim.A > quantile(MC.Lsim.A, alfa[3])])
cbind(mean(MC.Lsim.A[MC.Lsim.A > quantile(MC.Lsim.A, alfa[1])]),mean(MC.Lsim.A[MC.Lsim.A > quantile(MC.Lsim.A, alfa[2])]),mean(MC.Lsim.A[MC.Lsim.A > quantile(MC.Lsim.A, alfa[3])]))
#TVaR skenario B
mean(MC.Lsim.B[MC.Lsim.B > quantile(MC.Lsim.B, alfa[1])])
mean(MC.Lsim.B[MC.Lsim.B > quantile(MC.Lsim.B, alfa[2])])
mean(MC.Lsim.B[MC.Lsim.B > quantile(MC.Lsim.B, alfa[3])])
cbind(mean(MC.Lsim.B[MC.Lsim.B > quantile(MC.Lsim.B, alfa[1])]),mean(MC.Lsim.B[MC.Lsim.B > quantile(MC.Lsim.B, alfa[2])]),mean(MC.Lsim.B[MC.Lsim.B > quantile(MC.Lsim.B, alfa[3])]))

## Kopula Student *t* Margin Distribusi Normal

#  kopula t margin normal Skenario A
sim.RMst <- qnorm(sim.return.tcop[,1], mean=est.normal1$estimate[[1]], sd=est.normal1$estimate[[2]])
sim.RHst <- qnorm(sim.return.tcop[,2], mean=est.normal2$estimate[[1]], sd=est.normal2$estimate[[2]])
MC.data.st<-cbind(sim.RMst, sim.RHst)
MC.Lsim.stA <- -(MC.data.st%*%t(wa)) 
MC.Lsim.stB <- -(MC.data.st%*%t(wb)) 
#VaR
quantile(MC.Lsim.stA, alfa)
quantile(MC.Lsim.stB, alfa)
#TVaR Skenario A
mean(MC.Lsim.stA[MC.Lsim.stA > quantile(MC.Lsim.stA, alfa[1])])
mean(MC.Lsim.stA[MC.Lsim.stA > quantile(MC.Lsim.stA, alfa[2])])
mean(MC.Lsim.stA[MC.Lsim.stA > quantile(MC.Lsim.stA, alfa[3])])
cbind(mean(MC.Lsim.stA[MC.Lsim.stA > quantile(MC.Lsim.stA, alfa[1])]),mean(MC.Lsim.stA[MC.Lsim.stA > quantile(MC.Lsim.stA, alfa[2])]),mean(MC.Lsim.stA[MC.Lsim.stA > quantile(MC.Lsim.stA, alfa[3])]))
#TVaR Skenario B
mean(MC.Lsim.stB[MC.Lsim.stB > quantile(MC.Lsim.stB, alfa[1])])
mean(MC.Lsim.stB[MC.Lsim.stB > quantile(MC.Lsim.stB, alfa[2])])
mean(MC.Lsim.stA[MC.Lsim.stB > quantile(MC.Lsim.stB, alfa[3])])
cbind(mean(MC.Lsim.stB[MC.Lsim.stB > quantile(MC.Lsim.stB, alfa[1])]),mean(MC.Lsim.stB[MC.Lsim.stB > quantile(MC.Lsim.stB, alfa[2])]),mean(MC.Lsim.stB[MC.Lsim.stB > quantile(MC.Lsim.stB, alfa[3])]))

## Kopula Gumbel Margin Distribusi Normal

#  kopula Gumbel margin normal}
sim.RMg <- qnorm(sim.return.Gumbelcop[,1], mean=est.normal1$estimate[[1]], sd=est.normal1$estimate[[2]])
sim.RHg <- qnorm(sim.return.Gumbelcop[,2], mean=est.normal2$estimate[[1]], sd=est.normal2$estimate[[2]])
MC.data.gu <-cbind(sim.RMg, sim.RHg)
MC.Lsim.gA <- -(MC.data.gu%*%t(wa))
MC.Lsim.gB <- -(MC.data.gu%*%t(wb))
#VaR
quantile(MC.Lsim.gA, alfa)
quantile(MC.Lsim.gB, alfa)
#TVaR skenario A
mean(MC.Lsim.gA[MC.Lsim.gA > quantile(MC.Lsim.gA, alfa[1])])
mean(MC.Lsim.gA[MC.Lsim.gA > quantile(MC.Lsim.gA, alfa[2])])
mean(MC.Lsim.gA[MC.Lsim.gA > quantile(MC.Lsim.gA, alfa[3])])
cbind(mean(MC.Lsim.gA[MC.Lsim.gA> quantile(MC.Lsim.gA, alfa[1])]),mean(MC.Lsim.gA[MC.Lsim.gA > quantile(MC.Lsim.gA, alfa[2])]),mean(MC.Lsim.gA[MC.Lsim.gA > quantile(MC.Lsim.gA, alfa[3])]))
#TVaR skenario B
mean(MC.Lsim.gB[MC.Lsim.gB > quantile(MC.Lsim.gB, alfa[1])])
mean(MC.Lsim.gB[MC.Lsim.gB > quantile(MC.Lsim.gB, alfa[2])])
mean(MC.Lsim.gB[MC.Lsim.gB > quantile(MC.Lsim.gB, alfa[3])])
cbind(mean(MC.Lsim.gB[MC.Lsim.gB> quantile(MC.Lsim.gB, alfa[1])]),mean(MC.Lsim.gB[MC.Lsim.gB > quantile(MC.Lsim.gB, alfa[2])]),mean(MC.Lsim.gB[MC.Lsim.gB > quantile(MC.Lsim.gB, alfa[3])]))

## Kopula Clayton Margin Distribusi Normal

#  kopula Clayton margin normal}
sim.RMc <- qnorm(sim.return.ClaytonCop[,1], mean=est.normal1$estimate[[1]], sd=est.normal1$estimate[[2]])
sim.RHc <- qnorm(sim.return.ClaytonCop[,2], mean=est.normal2$estimate[[1]], sd=est.normal2$estimate[[2]])
MC.data.c<- cbind(sim.returns1, sim.returns2)
MC.Lsim.cA <- -(MC.data.c%*%t(wa))
MC.Lsim.cB <- -(MC.data.c%*%t(wb))
#VaR
quantile(MC.Lsim.cA, alfa)
quantile(MC.Lsim.cB, alfa)
#TVaR skenario A
mean(MC.Lsim.cA[MC.Lsim.cA > quantile(MC.Lsim.cA, alfa[1])])
mean(MC.Lsim.cA[MC.Lsim.cA > quantile(MC.Lsim.cA, alfa[2])])
mean(MC.Lsim.cA[MC.Lsim.cA > quantile(MC.Lsim.cA, alfa[3])])
cbind(mean(MC.Lsim.cA[MC.Lsim.cA > quantile(MC.Lsim.cA, alfa[1])]),mean(MC.Lsim.cA[MC.Lsim.cB > quantile(MC.Lsim.cA, alfa[2])]),mean(MC.Lsim.cA[MC.Lsim.cA > quantile(MC.Lsim.cA, alfa[3])]))
#TVaR skenario B
mean(MC.Lsim.cB[MC.Lsim.cB > quantile(MC.Lsim.cB, alfa[1])])
mean(MC.Lsim.cB[MC.Lsim.cB > quantile(MC.Lsim.cB, alfa[2])])
mean(MC.Lsim.cA[MC.Lsim.cB > quantile(MC.Lsim.cB, alfa[3])])
cbind(mean(MC.Lsim.cB[MC.Lsim.cB > quantile(MC.Lsim.cB, alfa[1])]),mean(MC.Lsim.cB[MC.Lsim.cB > quantile(MC.Lsim.cB, alfa[2])]),mean(MC.Lsim.cB[MC.Lsim.cB > quantile(MC.Lsim.cB, alfa[3])]))

#  kopula Frank margin normal}
sim.RMf <- qnorm(sim.retun.Frankcop[,1], mean=est.normal1$estimate[[1]], sd=est.normal1$estimate[[2]])
sim.RHf <- qnorm(sim.retun.Frankcop[,2], mean=est.normal2$estimate[[1]], sd=est.normal2$estimate[[2]])
MC.data.f<-cbind(sim.RMf, sim.RHf)
MC.Lsim.fA <- -(MC.data.f%*%t(wa))
MC.Lsim.fB <- -(MC.data.f%*%t(wb))
#VaR 
quantile(MC.Lsim.fA, alfa) 
quantile(MC.Lsim.fB, alfa) 
#TVaR skenario A
mean(MC.Lsim.fA[MC.Lsim.fA > quantile(MC.Lsim.fA, alfa[1])])
mean(MC.Lsim.fA[MC.Lsim.fA > quantile(MC.Lsim.fA, alfa[2])])
mean(MC.Lsim.fA[MC.Lsim.fA > quantile(MC.Lsim.fA, alfa[3])])
cbind(mean(MC.Lsim.fA[MC.Lsim.fA > quantile(MC.Lsim.fA, alfa[1])]),mean(MC.Lsim.fA[MC.Lsim.fA > quantile(MC.Lsim.fA, alfa[2])]),mean(MC.Lsim.fA[MC.Lsim.fA > quantile(MC.Lsim.fA, alfa[3])]))
#TVaR skenario B
mean(MC.Lsim.fB[MC.Lsim.fB > quantile(MC.Lsim.fB, alfa[1])])
mean(MC.Lsim.fB[MC.Lsim.fB > quantile(MC.Lsim.fB, alfa[2])])
mean(MC.Lsim.fB[MC.Lsim.fB > quantile(MC.Lsim.fB, alfa[3])])
cbind(mean(MC.Lsim.fB[MC.Lsim.fB > quantile(MC.Lsim.fB, alfa[1])]),mean(MC.Lsim.fB[MC.Lsim.fB > quantile(MC.Lsim.fB, alfa[2])]),mean(MC.Lsim.fB[MC.Lsim.fB > quantile(MC.Lsim.fB, alfa[3])]))
