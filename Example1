#install.packages("ggplot2")
library(ggplot2)
#install.packages("invgamma")
library(invgamma)
#install.packages("ContourFunctions")
library(ContourFunctions)
#install.packages("rgl")
library(rgl)
#install.packages("latex2exp")
library(latex2exp)
#install.packages("pracma")
library(pracma)
#install.packages("EnvStats")
library(EnvStats)
library(grid)
#install.packages("gridExtra")
library(gridExtra)


#### Example 1 ####
#rm(list = ls())

# values to plot the graph of the posterior
thetas <- seq(0.01, 4, .01)

# values under H0
theta00 <- c(0.3,0.6,0.9, 1.2,1.5,1.8,2.1,2.4) # H0
theta0=2.4

# model y: Exp(1/theta), i.e. E[y]=theta
# prior Jeffreys'
# posterior is an inverse gamma

# scen A : n =6
set.seed(1)
yA=rexp(6,rate =1/1.2)+.089327
sum(yA)
mean(yA)
alpha_A <-nA<- 6      # n
delta_A <- 7.2    # t_n

# scen B ; n=12
yB=c(yA, yA)
alpha_B <-nB<- 6*2      # n
delta_B <- 7.2*2    # t_n= sum(y)

# scen C ; n=20
set.seed(1)
yC=rexp(20,1/1.2)-0.1072
mean(yC)
alpha_C <-nC<- 20      # n
delta_C <- sum(yC)  # t_n

# scen D ; n=40
yD=c(yC, yC)
alpha_D <-nD<- 40     # n
delta_D <- sum(yD)   # t_n


# Posterior distributions
g_post_A <- dinvgamma(x = thetas, shape = alpha_A, rate = delta_A)
g_post_B <- dinvgamma(x = thetas, shape = alpha_B, rate = delta_B)
g_post_C <- dinvgamma(x = thetas, shape = alpha_C, rate = delta_C)
g_post_D <- dinvgamma(x = thetas, shape = alpha_D, rate = delta_D)

# Medians
m1_A <- qinvgamma(0.5, shape = alpha_A, rate = delta_A)
m1_B <- qinvgamma(0.5, shape = alpha_B, rate = delta_B)
m1_C <- qinvgamma(0.5, shape = alpha_C, rate = delta_C)
m1_D <- qinvgamma(0.5, shape = alpha_D, rate = delta_D)

# Single plots
data_A <- data.frame(thetas=thetas, g_posterior=g_post_A)
data_B <- data.frame(thetas=thetas, g_posterior=g_post_B)
data_C <- data.frame(thetas=thetas, g_posterior=g_post_C)
data_D <- data.frame(thetas=thetas, g_posterior=g_post_D)


# Discrepancy measure
delta_H = function(m1, theta_H, theta, alpha, delta){
if(m1 < theta_H){
    out = 1-2*(1-pinvgamma(theta_H, shape = alpha, rate = delta))
  }
  else {
    out = 1-2*(pinvgamma(theta_H, shape = alpha, rate = delta))
  }
}

delta_H_A=sapply(theta00, function (theta0)  delta_H(m1_A, theta_H = theta0,
                     theta = theta, alpha = alpha_A, delta = delta_A))
delta_H_B=sapply(theta00, function (theta0)   delta_H(m1_B, theta_H = theta0,
                     theta = theta, alpha = alpha_B, delta = delta_B))
delta_H_C=sapply(theta00, function (theta0)  delta_H(m1_C, theta_H = theta0,
                     theta = theta, alpha = alpha_C, delta = delta_C))
delta_H_D=sapply(theta00, function (theta0)  delta_H(m1_D, theta_H = theta0,
                        theta = theta, alpha = alpha_D, delta = delta_D))


 
round(delta_H_A,3)
round(delta_H_B,3)
round(delta_H_C,3)
round(delta_H_D,3)



### approximations of the BDM 
# log likelihood
llik=function(x,y) sum(dexp(y, rate = 1/x, log = T))
# log posterior
lpost=function(x,y) sum(dexp(y, rate = 1/x, log = T))+log(1/x)
# score
score= function (x,y) -length(y)/x+ sum(y)/x^2
# third derivative
third_deriv <- function(x,alpha,delta) {
  result=  -(alpha)/x^3+6*delta/x^4
  return(result)
}
 
# functions to implement the method by Zhou et al 2024
library(Matrix)

# Define zeta functions
zeta1 <- function(kappa) { dnorm(kappa) / pnorm(kappa) }
zeta2 <- function(kappa) {
  z1 <- zeta1(kappa)
  -z1 * (kappa + z1)
}
zeta3 <- function(kappa) {
  z1 <- zeta1(kappa)
  z2 <- zeta2(kappa)
  -z2 * (kappa + 2 * z1)
}

# Line search equation for kappa
solve_kappa <- function(R) {
  f <- function(kappa) {
    z1 <- zeta1(kappa)
    z2 <- zeta2(kappa)
    z3 <- zeta3(kappa)
    (kappa * z3^(2/3) - R * z1) * (z3^(2/3) + R * z2) + R^2 * z1 * z2
  }
  uniroot(f, lower = -10, upper = 10)$root
}

# Derivative Matching (DM) main function
derivative_matching <- function(m, J, t) {
  p <- length(m)
  
  # Step 1: Compute R
  t_1_3 <- sign(t) * abs(t)^(1/3)
  R <- sum(t_1_3 * solve(J, t_1_3))  # equivalent to (t^(1/3))' J^{-1} (t^(1/3))
  
  # Step 2: Solve for kappa
  kappa_star <- solve_kappa(R)
  
  # Step 3: Compute d*, Sigma*, mu*
  d_star <- sign(t) * abs(t / zeta3(kappa_star))^(1/3)
  Sigma_star <- solve(J + zeta2(kappa_star) * (d_star %*% t(d_star)))
  mu_star <- m - zeta1(kappa_star) * Sigma_star %*% d_star
  
  list(mu = mu_star, Sigma = Sigma_star, d = d_star)
}

##  case A
optA= optim(fn = function (x) -lpost(x, yA),par = 1, hessian = T)    
theta_tildeA=optA$par
j_theta_tildeA=optA$hessian

optAL= optim(fn = function (x) -llik(x, yA),par = 1, hessian = T)    
theta_hatA=optAL$par
j_theta_hatA=optAL$hessian

# I order
(BDM_IO_A=sapply(theta00, function (theta0) 
  2*(pnorm(abs(theta0-theta_hatA), sd=1/sqrt(j_theta_hatA)))-1))
# high order
rp<-function(theta, thetahat) sign(thetahat- theta)*sqrt(2*(llik(thetahat, yA)-llik(theta, yA)))
qb<- function(theta, thetahat) {
  score(theta, y = yA)* 1/sqrt(j_theta_hatA)*(thetahat/theta)^-1
}
rb<-function(theta, thetahat=theta_hatA) rp(theta, thetahat)+1/rp(theta, thetahat)*
  log(qb(theta, thetahat)/rp(theta, thetahat ))
rbV=Vectorize(rb, "theta")
 
BDM_HO_A=sapply(theta00, function (theta0) 
  2*pnorm(abs(rbV(theta0, theta_hatA)))-1)
 
# SKS
h0= sqrt(nA)*(theta0-theta_tildeA )
j_theta_tildeA=optA$hessian
omega_tildeA=nA/j_theta_tildeA
l3A=third_deriv(theta_tildeA,alpha_A,delta_A)
(BDM_SN_A=sapply(theta00, function (theta0) 2*(pnorm(sqrt(nA)* abs(theta0-theta_tildeA )/sqrt(omega_tildeA)))-
   2*sign(theta0-theta_tildeA )*l3A/(6*nA^(3/2))*omega_tildeA^(3/2)*
    dnorm(sqrt(nA)* (theta0-theta_tildeA ),0, sqrt(omega_tildeA))*
    ((sqrt(nA)*(theta0-theta_tildeA ))^2/omega_tildeA+2))-1)

(BDM_SN_Al=sapply(theta00, function (theta0) (pnorm(sqrt(nA)*(theta0-theta_tildeA )/sqrt(omega_tildeA)))-
                   l3A/(6*nA^(3/2))*omega_tildeA^(3/2)*
                   dnorm(sqrt(nA)* (theta0-theta_tildeA ),0, sqrt(omega_tildeA))*
                   ((sqrt(nA)*(theta0-theta_tildeA ))^2/omega_tildeA+2)))

(BDM_SN_Au=sapply(theta00, function (theta0) 1-pnorm(sqrt(nA)* (theta0-theta_tildeA )/sqrt(omega_tildeA))+
                   l3A/(6*nA^(3/2))*omega_tildeA^(3/2)*
                   dnorm(sqrt(nA)* (theta0-theta_tildeA ),0, sqrt(omega_tildeA))*
                   ((sqrt(nA)*(theta0-theta_tildeA ))^2/omega_tildeA+2)))

BDM_SN_A=sapply(1:length(BDM_SN_Au),function(k) 2*min(abs(BDM_SN_Au[k]-0.5),abs(BDM_SN_Al[k]-0.5)))
 
snpdf=function (x) 2*dnorm((x-theta_tildeA),0, 1/sqrt(j_theta_tildeA))*
  pnorm((x-theta_tildeA)^3*sqrt(2*pi)*l3A/(12))
snpdfV=Vectorize(snpdf)
BDM_SN_num_A=sapply(theta00, function (theta0) 
  2*(max((integrate(snpdfV, lower=theta0, upper=Inf)$value),
                                                      (integrate(snpdfV, upper=theta0, lower=-Inf)$value))))-1
 

DM=derivative_matching(optA$par,optA$hessian, l3A)
BDM_SN_A_DM=sapply(theta00, function (theta0) 
  2*min(sn::psn(theta0,xi=DM$mu,alpha = DM$d,omega = sqrt(DM$Sigma)),
      1-sn::psn(theta0,xi=DM$mu,alpha = DM$d,omega = sqrt(DM$Sigma))))


# PLOT POSTERIOR
sn_post_A =  2*dnorm((thetas-theta_tildeA),0,1/sqrt(j_theta_tildeA))*
  pnorm((thetas-theta_tildeA)^3*l3A*sqrt(2*pi)/(12))
data_A_SN <- data.frame(thetas=thetas, sn_posterior=sn_post_A)
p1ASN <-  ggplot(data_A_SN, aes(x=thetas,y=sn_posterior))+
  geom_line(size = 0.4)+
  geom_vline(xintercept=m1_A,size = 0.4,col="blue3", linetype="twodash")+
  # geom_vline(xintercept=theta0,size = 0.4,col="red3")+
  labs(title="", x=TeX("$\\theta"),
       y=TeX("SNormal"))+
  scale_x_continuous(expand = c(0, 0), breaks = seq(0,4,0.5)) + 
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,2,0.5), limits = c(0,2)) +
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5))

pdfn<- function (x) dnorm( (x-theta_hatA),sd = 1/sqrt(j_theta_hatA))
n_post_A=pdfn(thetas)
data_A_N <- data.frame(thetas=thetas, n_posterior=n_post_A)
p1AN <-  ggplot(data_A_N, aes(x=thetas,y=n_posterior))+
  geom_line(size = 0.4)+
  geom_vline(xintercept=m1_B,size = 0.4,col="blue3", linetype="twodash")+
  #  geom_vline(xintercept=theta0,size = 0.4,col="red3")+
  labs(title="", x=TeX("$\\theta"),
       y=TeX("$Normal-IO"))+
  scale_x_continuous(expand = c(0, 0), breaks = seq(0,4,0.5)) + 
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,2,0.5), limits = c(0,2)) +
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5))

# avoid numerical problems
library(likelihoodAsy)
lik  <- function(theta, data)
{
   lik=sum(dexp(data$y, rate = 1/theta[1], log = T))+
     dnorm(theta[2],0,50, log=T)
  return(lik)
}
# We need a function for simulating a data set
genData  <- function(theta, data) 
{
  out <- data
   
  Y=rexp(length(data$y),   rate = 1/theta[1])
  #out=list()
  out$y <-Y
  return(out)
  
}
data=list(y=yA)

# SN DERIVATIVE MATCHING
SNDM_POST=sapply(thetas, function(x) sn::dsn(x,xi=DM$mu,alpha = DM$d,omega = sqrt(DM$Sigma)))
data_A_SN_DM=data_A
data_A_SN_DM$g_posterior=SNDM_POST

# high order
thetasL=seq(from=min(thetas), to=max(thetas), length=40000)
RSTARS=sapply(thetasL, function (theta0) 
   ( (rbV(theta0, theta_hatA))))
RSTARS[is.nan(RSTARS)]=0


data_A_HO=data_A
resampled=sample(thetasL,400000,dnorm(qnorm((pnorm(RSTARS)))),replace=T)
d=density(resampled,from = min(thetas), to=max(thetas),n = 400)
d$x==thetas
data_A_HO$g_posterior=d$y
data_A_N$group <- "IOrder"
data_A_SN$group <- "SKS"
data_A_SN_DM$group <- "SN"
data_A_HO$group <- "HOrder"
data_A$group <- "IGamma"
colnames(data_A)
colnames(data_A_SN)[2]="g_posterior"
colnames(data_A_N)[2]="g_posterior"
combined_data <- rbind(data_A_N,data_A_SN, data_A_HO, data_A, data_A_SN_DM)

# Create the overimposed plot
p_all <- ggplot(combined_data, aes(x = thetas, y = g_posterior, color = group)) +
  geom_line(size = 0.6) +
  geom_vline(xintercept = m1_A, size = 0.4, col = "blue3", linetype = "twodash") +
  labs(title = "n=6",
       x = TeX("$\\theta$"),
       y = TeX("$\\pi(\\theta | y)$"),
       color = " ") + 
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 4, 0.5)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 2, 0.5), limits = c(0, 2.2)) +
  scale_color_manual(
    values = c(
      "HOrder" = "#E55",   # Orange, distinct for colorblind users
      "IGamma" = "#69B126",   # Sky blue, distinguishable for colorblind users
      "IOrder" = "#E71",   # Teal, works well for colorblind users
      "SKS" = "#F9E442",     # Yellow, works for colorblind users
      "SN" = "#0092B2"  )) +     # Dark blue, also colorblind-friendly    ))+
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

p_all



##  case B
optB= optim(fn = function (x) -lpost(x, yB),par = 1, hessian = T)    
#equal to  mean(y)
theta_tildeB=optB$par
j_theta_tildeB=optB$hessian
optBL= optim(fn = function (x) -llik(x, yB),par = 1, hessian = T)    
#equal to  mean(y)
theta_hatB=optBL$par
j_theta_hatB=optBL$hessian
# I order
(BDM_IO_B=sapply(theta00, 
                 function(theta0) 2*(pnorm(abs((theta0-theta_hatB)),
                      sd = 1/sqrt(j_theta_hatB)))-1))

# high order
rp<-function(theta, thetahat) sign(thetahat- theta)*sqrt(2*(llik(thetahat, yB)-llik(theta, yB)))

qb<- function(theta, thetahat) {
  score(theta, y = yB)* 1/sqrt(j_theta_hatB)*(thetahat/theta)^-1
}

score(theta_tildeB, y = yB)
rb<-function(theta, thetahat=theta_hatB) rp(theta, thetahat)+1/rp(theta, thetahat)*
  log(qb(theta, thetahat)/rp(theta, thetahat ))
rbV=Vectorize(rb, "theta") 
integrate(f = function (x) dnorm(rbV(x, theta_hatB)),0, Inf)
BDM_HO_B=sapply(theta00, function (theta0) 
  2*pnorm(abs(rbV(theta0, theta_hatB)))-1)
BDM_HO_B
 
# SKS
h0= sqrt(nB)*(theta0-theta_tildeB )
omega_tildeB=nB/j_theta_tildeB
l3B=third_deriv(theta_tildeB,alpha_B,delta_B)
(BDM_SN_B=sapply(theta00, function (theta0) 2*(pnorm(h0/sqrt(omega_tildeB)))-
    2*l3B/(6*nB^(3/2))*omega_tildeB^(3/2)*
    dnorm(h0,0, sqrt(omega_tildeB))*
    (h0^2/omega_tildeB+2)-1))

(BDM_SN_B=sapply(theta00, function (theta0) 2*(pnorm(sqrt(nB)* abs(theta0-theta_tildeB )/
                                                       sqrt(omega_tildeB)))-sign((theta0-theta_tildeB ))*
                   2*l3A/(6*nB^(3/2))*omega_tildeB^(3/2)*
                   dnorm(sqrt(nB)* (theta0-theta_tildeB ),0, sqrt(omega_tildeB))*
                   ((sqrt(nB)*(theta0-theta_tildeB ))^2/omega_tildeB+2))-1)

#(BDM_SN_B=2*(pnorm(sqrt(nB)*(theta0-theta_tildeB),0,sqrt(nB)/sqrt(j_theta_tildeB)))-
 #   2*l3B/(6*nB^(3/2))*(nB/j_theta_tildeB)^(3/2)*
  #  dnorm(sqrt(nB)*(theta0-theta_tildeB),0, sqrt(nB)/sqrt(j_theta_tildeB))*
   # (nB*(theta0-theta_tildeB)^2/(nB/j_theta_tildeB)+2)-1)


snpdf=function (x) 2*dnorm((x-theta_tildeB),0, 1/sqrt(j_theta_tildeB))*
  pnorm((x-theta_tildeB)^3*sqrt(2*pi)*l3B/(12))
snpdfV=Vectorize(snpdf)
#check that it integrates to 1
integrate(snpdf, lower=0, upper=8)$value
integrate(snpdf, lower=0, upper=Inf)$value
BDM_SN_num_B=sapply(theta00, function (theta0) 2*(max((integrate(snpdfV, lower=theta0, upper=Inf)$value),
  (integrate(snpdfV, upper=theta0, lower=-Inf)$value))))-1
round(BDM_SN_num_B,3)

# SN
DM=derivative_matching(optB$par,optB$hessian, l3B)

BDM_SN_B_DM=sapply(theta00, function (theta0) 
  2*min(sn::psn(theta0,xi=DM$mu,alpha = DM$d,omega = sqrt(DM$Sigma)),
        1-sn::psn(theta0,xi=DM$mu,alpha = DM$d,omega = sqrt(DM$Sigma))))


 
# SKS
sn_post_B = sn_post_B = 2*dnorm((thetas-theta_tildeB),0,1/sqrt(j_theta_tildeB))*
  pnorm((thetas-theta_tildeB)^3*l3B*sqrt(2*pi)/(12))

data_B_SN <- data.frame(thetas=thetas, sn_posterior=sn_post_B)
p1BSN <-  ggplot(data_B_SN, aes(x=thetas,y=sn_posterior))+
  geom_line(size = 0.4)+
  geom_vline(xintercept=m1_B,size = 0.4,col="blue3", linetype="twodash")+
  # geom_vline(xintercept=theta0,size = 0.4,col="red3")+
  labs(title="", x=TeX("$\\theta"),
       y=TeX("SNormal"))+
  scale_x_continuous(expand = c(0, 0), breaks = seq(0,4,0.5)) + 
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,2,0.5), limits = c(0,2)) +
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5))

pdfn<- function (x) dnorm( (x-theta_hatB),sd = 1/sqrt(j_theta_hatB))
n_post_B=pdfn(thetas)
data_B_N <- data.frame(thetas=thetas, n_posterior=n_post_B)
p1BN <-  ggplot(data_B_N, aes(x=thetas,y=n_posterior))+
  geom_line(size = 0.4)+
  geom_vline(xintercept=m1_B,size = 0.4,col="blue3", linetype="twodash")+
  # geom_vline(xintercept=theta0,size = 0.4,col="red3")+
  labs(title="", x=TeX("$\\theta"),
       y=TeX("Normal-IO"))+
  scale_x_continuous(expand = c(0, 0), breaks = seq(0,4,0.5)) + 
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,2,0.5), limits = c(0,2)) +
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5))


thetasL=seq(from=min(thetas), to=max(thetas), length=40000)
RSTARS=sapply(thetasL, function (theta0) 
  ( (rbV(theta0, theta_hatB))))
RSTARS[is.nan(RSTARS)]=0

# SN DERIVATIVE MATCHING
SNDM_POST=sapply(thetas, function(x) sn::dsn(x,xi=DM$mu,alpha = DM$d,omega = sqrt(DM$Sigma)))
data_B_SN_DM=data_B
data_B_SN_DM$g_posterior=SNDM_POST

data_B_HO=data_B
resampled=sample(thetasL,400000,dnorm(qnorm((pnorm(RSTARS)))),replace=T)
d=density(resampled,from = min(thetas), to=max(thetas),n = 400)
d$x==thetas
data_B_HO$g_posterior=d$y
data_B_N$group <- "IOrder"
data_B_SN$group <- "SKS"
data_B_SN_DM$group <- "SN"
data_B_HO$group <- "HOrder"
data_B$group <- "IGamma"
colnames(data_B)
colnames(data_B_SN)[2]="g_posterior"
colnames(data_B_N)[2]="g_posterior"
combined_data <- rbind(data_B_N,data_B_SN, data_B_HO, data_B, data_B_SN_DM)

# Create the overimposed plot
p_all <- ggplot(combined_data, aes(x = thetas, y = g_posterior, color = group)) +
  geom_line(size = 0.6) +
  geom_vline(xintercept = m1_A, size = 0.4, col = "blue3", linetype = "twodash") +
  labs(title = "n=12",
       x = TeX("$\\theta$"),
       y = TeX("$\\pi(\\theta | y)$"),
       color = " ") +
  scale_color_manual(
    values = c(
      "HOrder" = "#E55",   # Orange, distinct for colorblind users
      "IGamma" = "#69B126",   # Sky blue, distinguishable for colorblind users
      "IOrder" = "#E71",   # Teal, works well for colorblind users
      "SKS" = "#F9E442",     # Yellow, works for colorblind users
      "SN" = "#0092B2"  )) +     # Dark blue, also colorblind-friendly    ))+
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 4, 0.5)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 2, 0.5), limits = c(0, 2.2)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

p_all
# results BDM 
library(xtable)

results6=rbind(
  BDM_HO_A,
  BDM_IO_A, delta_H_A,
 
  BDM_SN_A,
 BDM_SN_num_A,
 1-BDM_SN_A_DM
 )
results6
results12=rbind(
  BDM_HO_B,
  
  delta_H_B, BDM_IO_B,
  BDM_SN_B,
  BDM_SN_num_B,1-BDM_SN_B_DM
)
results6[is.nan(results6)]=0
results6[1<(results6)]=1
 
matplot((t(results6)), type="l",axes=F, col=c("#E55","#69B126","#E71", "#F9E442",8, "#0092B2" ),lty=1, ylab="BDM",
        main="n=6")
axis(1, 1:8,theta00)
axis(2,seq(0,1, by=0.1))
legend("bottomright", col=c("#E55","#69B126","#E71", "#F9E442",8, "#0092B2" ),
       lty=1, 
       bty ="n",
       legend = c("HOrder","IGamma","IOrder",
                                          "SKS", "SKS-num", "SN"), cex=0.7)


results12[is.nan(results12)]=0
results12[1<(results12)]=1
matplot((t(results12)), type="l",axes=F,  col=c("#E55","#69B126","#E71", "#F9E442",8, "#0092B2" ),lty=1, ylab="BDM",
        main="n=12")
axis(1, 1:8,theta00)
axis(2,seq(0,1, by=0.1))
legend("bottomright",  col=c("#E55","#69B126","#E71", "#F9E442",8, "#0092B2" ),lty=1,bty="n", legend = c("HOrder","IGamma","IOrder",
                                                                "SKS", "SKS-num", "SN"), cex=0.7)

matplot((t(results12)), type="l",axes=F, col=5:1,lty=5:1, ylab="BDM",
        main="n=12")
axis(1, 1:8,theta00)
axis(2,seq(0,1, by=0.1))
legend("bottomright", col=c(2,5,3,8,"violet",8),lty=1,bty="n", legend = c("HOrder","IGamma","IOrder",
                                                                          "SKS", "SKS-num"), cex=0.7)

results=cbind(rbind(
    BDM_IO_A,
  BDM_HO_A,
  BDM_SN_A,
  BDM_SN_num_A,
  delta_H_A),
  rbind(
    BDM_IO_B,
    BDM_HO_B,
    BDM_SN_B,
    BDM_SN_num_B,
    delta_H_B))
 
xtable::xtable(results)
round(cbind(rbind(
  BDM_IO_A,
  BDM_HO_A,
  BDM_SN_A,
  BDM_SN_num_A,
  delta_H_A),
  rbind(
    BDM_IO_B,
    BDM_HO_B,
    BDM_SN_B,
    BDM_SN_num_B,
    delta_H_B)),3)



########## 

# C and D

##  case C
optC= optim(fn = function (x) -lpost(x, yC),par = 1, hessian = T)    
theta_tildeC=optC$par
j_theta_tildeC=optC$hessian

optCL= optim(fn = function (x) -llik(x, yC),par = 1, hessian = T)    
theta_hatC=optCL$par
j_theta_hatC=optCL$hessian

# I order
(BDM_IO_C=sapply(theta00, function (theta0) 
  2*(pnorm(abs(theta0-theta_hatC), sd=1/sqrt(j_theta_hatC)))-1))
# high order
rp<-function(theta, thetahat) sign(thetahat- theta)*sqrt(2*(llik(thetahat, yC)-llik(theta, yC)))
qb<- function(theta, thetahat) {
  score(theta, y = yC)* 1/sqrt(j_theta_hatC)*(thetahat/theta)^-1
}
rb<-function(theta, thetahat=theta_hatC) rp(theta, thetahat)+1/rp(theta, thetahat)*
  log(qb(theta, thetahat)/rp(theta, thetahat ))
rbV=Vectorize(rb, "theta")

BDM_HO_C=sapply(theta00, function (theta0) 
  2*pnorm(abs(rbV(theta0, theta_hatC)))-1)

# SN 
h0= sqrt(nC)*(theta0-theta_tildeC )
j_theta_tildeC=optC$hessian
omega_tildeC=nC/j_theta_tildeC
l3C=third_deriv(theta_tildeC,alpha_C,delta_C)
(BDM_SN_C=sapply(theta00, function (theta0) 2*(pnorm(sqrt(nC)* abs(theta0-theta_tildeC )/sqrt(omega_tildeC)))-
                   2*sign(theta0-theta_tildeC )*l3C/(6*nC^(3/2))*omega_tildeC^(3/2)*
                   dnorm(sqrt(nC)* (theta0-theta_tildeC ),0, sqrt(omega_tildeC))*
                   ((sqrt(nC)*(theta0-theta_tildeC ))^2/omega_tildeC+2))-1)

(BDM_SN_Cl=sapply(theta00, function (theta0) (pnorm(sqrt(nC)*(theta0-theta_tildeC )/sqrt(omega_tildeC)))-
                    l3C/(6*nC^(3/2))*omega_tildeC^(3/2)*
                    dnorm(sqrt(nC)* (theta0-theta_tildeC ),0, sqrt(omega_tildeC))*
                    ((sqrt(nC)*(theta0-theta_tildeC ))^2/omega_tildeC+2)))

(BDM_SN_Cu=sapply(theta00, function (theta0) 1-pnorm(sqrt(nC)* (theta0-theta_tildeC )/sqrt(omega_tildeC))+
                    l3C/(6*nC^(3/2))*omega_tildeC^(3/2)*
                    dnorm(sqrt(nC)* (theta0-theta_tildeC ),0, sqrt(omega_tildeC))*
                    ((sqrt(nC)*(theta0-theta_tildeC ))^2/omega_tildeC+2)))

BDM_SN_C=sapply(1:length(BDM_SN_Cu),function(k) 2*min(abs(BDM_SN_Cu[k]-0.5),abs(BDM_SN_Cl[k]-0.5)))

snpdf=function (x) 2*dnorm((x-theta_tildeC),0, 1/sqrt(j_theta_tildeC))*
  pnorm((x-theta_tildeC)^3*sqrt(2*pi)*l3C/(12))
snpdfV=Vectorize(snpdf)
BDM_SN_num_C=sapply(theta00, function (theta0) 
  2*(max((integrate(snpdfV, lower=theta0, upper=Inf)$value),
         (integrate(snpdfV, upper=theta0, lower=-Inf)$value))))-1



# 
DM=derivative_matching(optC$par,optC$hessian, l3C)
BDM_SN_C_DM=sapply(theta00, function (theta0) 
  2*min(sn::psn(theta0,xi=DM$mu,alpha = DM$d,omega = sqrt(DM$Sigma)),
        1-sn::psn(theta0,xi=DM$mu,alpha = DM$d,omega = sqrt(DM$Sigma))))

# skew approx
sn_post_C = sn_post_B = 2*dnorm((thetas-theta_tildeC),0,1/sqrt(j_theta_tildeC))*
  pnorm((thetas-theta_tildeC)^3*l3C*sqrt(2*pi)/(12))

data_C_SN <- data.frame(thetas=thetas, sn_posterior=sn_post_C)
p1CSN <-  ggplot(data_C_SN, aes(x=thetas,y=sn_posterior))+
  geom_line(size = 0.4)+
  geom_vline(xintercept=m1_B,size = 0.4,col="blue3", linetype="twodash")+
  # geom_vline(xintercept=theta0,size = 0.4,col="red3")+
  labs(title="", x=TeX("$\\theta"),
       y=TeX("SNormal"))+
  scale_x_continuous(expand = c(0, 0), breaks = seq(0,4,0.5)) + 
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,2,0.5), limits = c(0,2)) +
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5))

pdfn<- function (x) dnorm( (x-theta_hatC),sd = 1/sqrt(j_theta_hatC))
n_post_C=pdfn(thetas)
data_C_N <- data.frame(thetas=thetas, n_posterior=n_post_C)
p1CN <-  ggplot(data_C_N, aes(x=thetas,y=n_posterior))+
  geom_line(size = 0.4)+
  geom_vline(xintercept=m1_B,size = 0.4,col="blue3", linetype="twodash")+
  # geom_vline(xintercept=theta0,size = 0.4,col="red3")+
  labs(title="", x=TeX("$\\theta"),
       y=TeX("Normal-IO"))+
  scale_color_manual(
    values = c(
      "HOrder" = "#E55",   # Orange, distinct for colorblind users
      "IGamma" = "#69B126",   # Sky blue, distinguishable for colorblind users
      "IOrder" = "#E71",   # Teal, works well for colorblind users
      "SKS" = "#F9E442",     # Yellow, works for colorblind users
      "SN" = "#0092B2"  )) +     # Dark blue, also colorblind-friendly    ))+
  
  scale_x_continuous(expand = c(0, 0), breaks = seq(0,4,0.5)) + 
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,2,0.5), limits = c(0,2)) +
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5))


thetasL=seq(from=min(thetas), to=max(thetas), length=40000)
RSTARS=sapply(thetasL, function (theta0) 
  ( (rbV(theta0, theta_hatC))))
RSTARS[is.nan(RSTARS)]=0


SNDM_POST=sapply(thetas, function(x) sn::dsn(x,xi=DM$mu,alpha = DM$d,omega = sqrt(DM$Sigma)))
data_C_SN_DM=data_C
data_C_SN_DM$g_posterior=SNDM_POST

data_C_HO=data_C
resampled=sample(thetasL,400000,dnorm(qnorm((pnorm(RSTARS)))),replace=T)
d=density(resampled,from = min(thetas), to=max(thetas),n = 400)
d$x==thetas
data_C_HO$g_posterior=d$y
data_C_N$group <- "IOrder"
data_C_SN$group <- "SKS"
data_C_HO$group <- "HOrder"
data_C$group <- "IGamma"
colnames(data_C)
colnames(data_C_SN)[2]="g_posterior"
colnames(data_C_N)[2]="g_posterior"
combined_data <- rbind(data_C_N,data_C_SN, data_C_HO, data_C, data_C_SN_DM)

# Create the overimposed plot
p_all <- ggplot(combined_data, aes(x = thetas, y = g_posterior, color = group)) +
  geom_line(size = 0.6) +
  geom_vline(xintercept = m1_A, size = 0.4, col = "blue3", linetype = "twodash") +
  labs(title = "n=20",
       x = TeX("$\\theta$"),
       y = TeX("$\\pi(\\theta | y)$"),
       color = " ") +
  scale_color_manual(
    values = c(
      "HOrder" = "#E55",   # Orange, distinct for colorblind users
      "IGamma" = "#69B126",   # Sky blue, distinguishable for colorblind users
      "IOrder" = "#E71",   # Teal, works well for colorblind users
      "SKS" = "#F9E442",     # Yellow, works for colorblind users
      "SN" = "#0092B2"  )) +     # Dark blue, also colorblind-friendly    ))+
  
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 4, 0.5)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 2, 0.5), limits = c(0, 2.2)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

p_all

##### D


##  case D
optD= optim(fn = function (x) -lpost(x, yD),par = 1, hessian = T)    
theta_tildeD=optD$par
j_theta_tildeD=optD$hessian

optDL= optim(fn = function (x) -llik(x, yD),par = 1, hessian = T)    
theta_hatD=optDL$par
j_theta_hatD=optDL$hessian

# I order
(BDM_IO_D=sapply(theta00, function (theta0) 
  2*(pnorm(abs(theta0-theta_hatD), sd=1/sqrt(j_theta_hatD)))-1))
# high order
rp<-function(theta, thetahat) sign(thetahat- theta)*sqrt(2*(llik(thetahat, yD)-llik(theta, yD)))
qb<- function(theta, thetahat) {
  score(theta, y = yD)* 1/sqrt(j_theta_hatD)*(thetahat/theta)^-1
}
rb<-function(theta, thetahat=theta_hatD) rp(theta, thetahat)+1/rp(theta, thetahat)*
  log(qb(theta, thetahat)/rp(theta, thetahat ))
rbV=Vectorize(rb, "theta")

BDM_HO_D=sapply(theta00, function (theta0) 
  2*pnorm(abs(rbV(theta0, theta_hatD)))-1)

# SN 
h0= sqrt(nD)*(theta0-theta_tildeD )
j_theta_tildeD=optD$hessian
omega_tildeD=nD/j_theta_tildeD
l3D=third_deriv(theta_tildeD,alpha_D,delta_D)
(BDM_SN_D=sapply(theta00, function (theta0) 2*(pnorm(sqrt(nD)* abs(theta0-theta_tildeD )/sqrt(omega_tildeD)))-
                   2*sign(theta0-theta_tildeD )*l3D/(6*nD^(3/2))*omega_tildeD^(3/2)*
                   dnorm(sqrt(nD)* (theta0-theta_tildeD ),0, sqrt(omega_tildeD))*
                   ((sqrt(nD)*(theta0-theta_tildeD ))^2/omega_tildeD+2))-1)

(BDM_SN_Dl=sapply(theta00, function (theta0) (pnorm(sqrt(nD)*(theta0-theta_tildeD )/sqrt(omega_tildeD)))-
                    l3D/(6*nD^(3/2))*omega_tildeD^(3/2)*
                    dnorm(sqrt(nD)* (theta0-theta_tildeD ),0, sqrt(omega_tildeD))*
                    ((sqrt(nD)*(theta0-theta_tildeD ))^2/omega_tildeD+2)))

(BDM_SN_Du=sapply(theta00, function (theta0) 1-pnorm(sqrt(nD)* (theta0-theta_tildeD )/sqrt(omega_tildeD))+
                    l3D/(6*nD^(3/2))*omega_tildeD^(3/2)*
                    dnorm(sqrt(nD)* (theta0-theta_tildeD ),0, sqrt(omega_tildeD))*
                    ((sqrt(nD)*(theta0-theta_tildeD ))^2/omega_tildeD+2)))

BDM_SN_D=sapply(1:length(BDM_SN_Du),function(k) 2*min(abs(BDM_SN_Du[k]-0.5),abs(BDM_SN_Dl[k]-0.5)))

snpdf=function (x) 2*dnorm((x-theta_tildeD),0, 1/sqrt(j_theta_tildeD))*
  pnorm((x-theta_tildeD)^3*sqrt(2*pi)*l3D/(12))
snpdfV=Vectorize(snpdf)
BDM_SN_num_D=sapply(theta00, function (theta0) 
  2*(max((integrate(snpdfV, lower=theta0, upper=Inf)$value),
         (integrate(snpdfV, upper=theta0, lower=-Inf)$value))))-1


 
# skew approx
DM=derivative_matching(optD$par,optD$hessian, l3D)
BDM_SN_D_DM=sapply(theta00, function (theta0) 
  2*min(sn::psn(theta0,xi=DM$mu,alpha = DM$d,omega = sqrt(DM$Sigma)),
        1-sn::psn(theta0,xi=DM$mu,alpha = DM$d,omega = sqrt(DM$Sigma))))

sn_post_D = sn_post_D = 2*dnorm((thetas-theta_tildeD),0,1/sqrt(j_theta_tildeD))*
  pnorm((thetas-theta_tildeD)^3*l3B*sqrt(2*pi)/(12))

data_D_SN <- data.frame(thetas=thetas, sn_posterior=sn_post_D)
p1DSN <-  ggplot(data_D_SN, aes(x=thetas,y=sn_posterior))+
  geom_line(size = 0.4)+
  geom_vline(xintercept=m1_D,size = 0.4,col="blue3", linetype="twodash")+
  # geom_vline(xintercept=theta0,size = 0.4,col="red3")+
  labs(title="", x=TeX("$\\theta"),
       y=TeX("SNormal"))+
  scale_color_manual(
    values = c(
      "HOrder" = "#E55",   # Orange, distinct for colorblind users
      "IGamma" = "#69B126",   # Sky blue, distinguishable for colorblind users
      "IOrder" = "#E71",   # Teal, works well for colorblind users
      "SKS" = "#F9E442",     # Yellow, works for colorblind users
      "SN" = "#0092B2"  )) +     # Dark blue, also colorblind-friendly    ))+
  
  scale_x_continuous(expand = c(0, 0), breaks = seq(0,4,0.5)) + 
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,2,0.5), limits = c(0,2)) +
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5))

pdfn<- function (x) dnorm( (x-theta_hatD),sd = 1/sqrt(j_theta_hatD))
n_post_D=pdfn(thetas)
data_D_N <- data.frame(thetas=thetas, n_posterior=n_post_D)
p1DN <-  ggplot(data_D_N, aes(x=thetas,y=n_posterior))+
  geom_line(size = 0.4)+
  geom_vline(xintercept=m1_D,size = 0.4,col="blue3", linetype="twodash")+
  # geom_vline(xintercept=theta0,size = 0.4,col="red3")+
  labs(title="", x=TeX("$\\theta"),
       y=TeX("Normal-IO"))+
  scale_x_continuous(expand = c(0, 0), breaks = seq(0,4,0.5)) + 
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,2,0.5), limits = c(0,2)) +
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5))


thetasL=seq(from=min(thetas), to=max(thetas), length=40000)
RSTARS=sapply(thetasL, function (theta0) 
  ( (rbV(theta0, theta_hatD))))
RSTARS[is.nan(RSTARS)]=0



SNDM_POST=sapply(thetas, function(x) sn::dsn(x,xi=DM$mu,alpha = DM$d,omega = sqrt(DM$Sigma)))
data_D_SN_DM=data_D
data_D_SN_DM$g_posterior=SNDM_POST

data_D_HO=data_D
resampled=sample(thetasL,400000,dnorm(qnorm((pnorm(RSTARS)))),replace=T)
d=density(resampled,from = min(thetas), to=max(thetas),n = 400)
d$x==thetas
data_D_HO$g_posterior=d$y
data_D_N$group <- "IOrder"
data_D_SN$group <- "SKS"
data_D_SN_DM$group <- "SN"
data_D_HO$group <- "HOrder"
data_D$group <- "IGamma"
colnames(data_D)
colnames(data_D_SN)[2]="g_posterior"
colnames(data_D_N)[2]="g_posterior"
combined_data <- rbind(data_D_N,data_D_SN, data_D_HO, data_D, data_D_SN_DM)

# Create the overimposed plot
p_all <- ggplot(combined_data, aes(x = thetas, y = g_posterior, color = group)) +
  geom_line(size = 0.6) +
  geom_vline(xintercept = m1_A, size = 0.4, col = "blue3", linetype = "twodash") +
  labs(title = "n=40",
       x = TeX("$\\theta$"),
       y = TeX("$\\pi(\\theta | y)$"),
       color = " ") +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 4, 0.5)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 2, 0.5), limits = c(0, 2.9)) +
  theme_classic() +
  scale_color_manual(
    values = c(
      "HOrder" = "#E55",   # Orange, distinct for colorblind users
      "IGamma" = "#69B126",   # Sky blue, distinguishable for colorblind users
      "IOrder" = "#E71",   # Teal, works well for colorblind users
      "SKS" = "#F9E442",     # Yellow, works for colorblind users
      "SN" = "#0092B2"  )) +     # Dark blue, also colorblind-friendly    ))+
  theme(plot.title = element_text(hjust = 0.5))

p_all
results20=rbind(
  BDM_HO_C,
  delta_H_C, BDM_IO_C,
  BDM_SN_C,
  BDM_SN_num_C, 
  1-BDM_SN_C_DM
)
results20[is.nan(results20)]=0
results20[1<(results20)]=1
matplot((t(results20)), type="l",axes=F,  col=c("#E55","#69B126","#E71", "#F9E442",8, "#0092B2" ),lty=1, ylab="BDM",
        main="n=20")
axis(1, 1:8,theta00)
axis(2,seq(0,1, by=0.1))
legend("bottomright",  col=c("#E55","#69B126","#E71", "#F9E442",8, "#0092B2" ),lty=1,bty="n", legend = c("HOrder","IGamma","IOrder",
                                                                                                         "SKS", "SKS-num", "SN"), cex=0.7)

 
round(1-BDM_SN_D_DM,2)

results40=rbind(
  BDM_HO_D,
  delta_H_D, BDM_IO_D,
  BDM_SN_D,
  BDM_SN_num_D, 1-BDM_SN_D_DM)
 
results40[is.nan(results40)]=0
results40[1<(results40)]=1
matplot((t(results40)), type="l",axes=F,  col=c("#E55","#69B126","#E71", "#F9E442",8, "#0092B2" ),lty=1, ylab="BDM",
        main="n=40")
axis(1, 1:8,theta00)
axis(2,seq(0,1, by=0.1))
legend("bottomright",  col=c("#E55","#69B126","#E71", "#F9E442",8, "#0092B2" ),lty=1,bty="n", legend = c("HOrder","IGamma","IOrder",
                                                                                                         "SKS", "SKS-num", "SN"), cex=0.7)

matplot((t(results40)), type="l",axes=F, col=c(2,3,5,"violet",8),lty=1, ylab="BDM",
        main="n=40")
axis(1, 1:8,theta00)
axis(2,seq(0,1, by=0.1))
legend("bottomright", col=c(2,3,5,"violet",8),lty=1, bty ="n",
       legend = c("HOrder","IGamma","IOrder",
                  "SKS", "SKS-num"), cex=0.7)


results=cbind(rbind(
  
  BDM_IO_A,
  BDM_HO_A,
  BDM_SN_A,
  BDM_SN_num_A,
  delta_H_A),
  rbind(
    BDM_IO_B,
    BDM_HO_B,
    BDM_SN_B,
    BDM_SN_num_B,
    delta_H_B))

xtable::xtable(rbind(results6, results12, results20, results40))
round(cbind(rbind(
  BDM_IO_A,
  BDM_HO_A,
  BDM_SN_A,
  BDM_SN_num_A,
  delta_H_A),
  rbind(
    BDM_IO_B,
    BDM_HO_B,
    BDM_SN_B,
    BDM_SN_num_B,
    delta_H_B)),3)



########## 

# C and D



####### PLOT of posterior
p1A <- ggplot(data_A, aes(x=thetas,y=g_posterior))+
  geom_line(size = 0.4)+
  geom_vline(xintercept=m1_A,size = 0.4,col="blue3", linetype="twodash")+
  #geom_vline(xintercept=theta0,size = 0.4,col="red3")+
  labs(title="n = 6", x=TeX("$\\theta"),
       y=TeX("$IG(\\theta | n,t_n)"))+
  scale_x_continuous(expand = c(0, 0), breaks = seq(0,4,0.5)) + 
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,2,0.5), limits = c(0,2)) +
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5))

p1B <- ggplot(data_B, aes(x=thetas,y=g_posterior))+
  geom_line(size = 0.4)+
  geom_vline(xintercept=m1_A,size = 0.4,col="blue3", linetype="twodash")+
  # geom_vline(xintercept=theta0,size = 0.4,col="red3")+
  labs(title="n = 12", x=TeX("$\\theta"),
       y=TeX("$IG(\\theta | n,t_n)"))+
  scale_x_continuous(expand = c(0, 0), breaks = seq(0,4,0.5)) + 
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,2,0.5), limits = c(0,2)) +
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5))

p1C <- ggplot(data_C, aes(x=thetas,y=g_posterior))+
  geom_line(size = 0.4)+
  geom_vline(xintercept=m1_A,size = 0.4,col="blue3", linetype="twodash")+
  #  geom_vline(xintercept=theta0,size = 0.4,col="red3")+
  labs(title="n = 20", x=TeX("$\\theta"),
       y=TeX("$IG(\\theta | n,t_n)"))+
  scale_x_continuous(expand = c(0, 0), breaks = seq(0,4,0.5)) + 
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,2,0.5), limits = c(0,2)) +
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5))

p1D <- ggplot(data_D, aes(x=thetas,y=g_posterior))+
  geom_line(size = 0.4)+
  geom_vline(xintercept=m1_A,size = 0.4,col="blue3", linetype="twodash")+
  #  geom_vline(xintercept=theta0,size = 0.4,col="red3")+
  labs(title="n = 40", x=TeX("$\\theta"),
       y=TeX("$IG(\\theta | n,t_n)"))+
  scale_x_continuous(expand = c(0, 0), breaks = seq(0,4,0.5)) + 
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,2,0.5), limits = c(0,2.2)) +
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5))

library(ggplot2)
library(latex2exp)


# Labels: median and H
m1 <- TeX("$m_1$")
H <- TeX("$\\theta_H$")



sn_post_A =  2*dnorm((thetas-theta_tildeA),0,1/sqrt(j_theta_tildeA))*
  pnorm((thetas-theta_tildeA)^3*l3A*sqrt(2*pi)/(12))
data_A_SN <- data.frame(thetas=thetas, sn_posterior=sn_post_A)
p1ASN <-  ggplot(data_A_SN, aes(x=thetas,y=sn_posterior))+
  geom_line(size = 0.4)+
  geom_vline(xintercept=m1_A,size = 0.4,col="blue3", linetype="twodash")+
 # geom_vline(xintercept=theta0,size = 0.4,col="red3")+
  labs(title="", x=TeX("$\\theta"),
       y=TeX("SNormal"))+
  scale_x_continuous(expand = c(0, 0), breaks = seq(0,4,0.5)) + 
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,2,0.5), limits = c(0,2)) +
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5))

pdfn<- function (x) dnorm( (x-theta_hatA),sd = 1/sqrt(j_theta_hatA))
n_post_A=pdfn(thetas)
data_A_N <- data.frame(thetas=thetas, n_posterior=n_post_A)
p1AN <-  ggplot(data_A_N, aes(x=thetas,y=n_posterior))+
  geom_line(size = 0.4)+
  geom_vline(xintercept=m1_B,size = 0.4,col="blue3", linetype="twodash")+
#  geom_vline(xintercept=theta0,size = 0.4,col="red3")+
  labs(title="", x=TeX("$\\theta"),
       y=TeX("$Normal-IO"))+
  scale_x_continuous(expand = c(0, 0), breaks = seq(0,4,0.5)) + 
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,2,0.5), limits = c(0,2)) +
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5))

# skew approx
sn_post_B = sn_post_B = 2*dnorm((thetas-theta_tildeB),0,1/sqrt(j_theta_tildeB))*
  pnorm((thetas-theta_tildeB)^3*l3B*sqrt(2*pi)/(12))

data_B_SN <- data.frame(thetas=thetas, sn_posterior=sn_post_B)
p1BSN <-  ggplot(data_B_SN, aes(x=thetas,y=sn_posterior))+
  geom_line(size = 0.4)+
  geom_vline(xintercept=m1_B,size = 0.4,col="blue3", linetype="twodash")+
 # geom_vline(xintercept=theta0,size = 0.4,col="red3")+
  labs(title="", x=TeX("$\\theta"),
       y=TeX("SNormal"))+
  scale_x_continuous(expand = c(0, 0), breaks = seq(0,4,0.5)) + 
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,2,0.5), limits = c(0,2)) +
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5))

pdfn<- function (x) dnorm( (x-theta_hatB),sd = 1/sqrt(j_theta_hatB))
n_post_B=pdfn(thetas)
data_B_N <- data.frame(thetas=thetas, n_posterior=n_post_B)
p1BN <-  ggplot(data_B_N, aes(x=thetas,y=n_posterior))+
  geom_line(size = 0.4)+
  geom_vline(xintercept=m1_B,size = 0.4,col="blue3", linetype="twodash")+
 # geom_vline(xintercept=theta0,size = 0.4,col="red3")+
  labs(title="", x=TeX("$\\theta"),
       y=TeX("Normal-IO"))+
  scale_x_continuous(expand = c(0, 0), breaks = seq(0,4,0.5)) + 
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,2,0.5), limits = c(0,2)) +
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5))



grid.arrange( p1A,  p1B,p1ASN, p1BSN,  p1AN, p1BN,ncol = 2)


sn_post_C =  2*dnorm((thetas-theta_tildeC),0,1/sqrt(j_theta_tildeC))*
  pnorm((thetas-theta_tildeC)^3*l3C*sqrt(2*pi)/(12))
data_C_SN <- data.frame(thetas=thetas, sn_posterior=sn_post_C)
p1CSN <-  ggplot(data_C_SN, aes(x=thetas,y=sn_posterior))+
  geom_line(size = 0.4)+
  geom_vline(xintercept=m1_C,size = 0.4,col="blue3", linetype="twodash")+
  # geom_vline(xintercept=theta0,size = 0.4,col="red3")+
  labs(title="", x=TeX("$\\theta"),
       y=TeX("SNormal"))+
  scale_x_continuous(expand = c(0, 0), breaks = seq(0,4,0.5)) + 
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,2,0.5), limits = c(0,2.2)) +
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5))

pdfn<- function (x) dnorm( (x-theta_hatC),sd = 1/sqrt(j_theta_hatC))
n_post_C=pdfn(thetas)
data_C_N <- data.frame(thetas=thetas, n_posterior=n_post_C)
p1CN <-  ggplot(data_C_N, aes(x=thetas,y=n_posterior))+
  geom_line(size = 0.4)+
  geom_vline(xintercept=m1_B,size = 0.4,col="blue3", linetype="twodash")+
  #  geom_vline(xintercept=theta0,size = 0.4,col="red3")+
  labs(title="", x=TeX("$\\theta"),
       y=TeX("$Normal-IO"))+
  scale_x_continuous(expand = c(0, 0), breaks = seq(0,4,0.5)) + 
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,2,0.5), limits = c(0,2.2)) +
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5))



sn_post_D =  2*dnorm((thetas-theta_tildeD),0,1/sqrt(j_theta_tildeD))*
  pnorm((thetas-theta_tildeD)^3*l3D*sqrt(2*pi)/(12))
data_D_SN <- data.frame(thetas=thetas, sn_posterior=sn_post_D)
p1DSN <-  ggplot(data_D_SN, aes(x=thetas,y=sn_posterior))+
  geom_line(size = 0.4)+
  geom_vline(xintercept=m1_D,size = 0.4,col="blue3", linetype="twodash")+
  # geom_vline(xintercept=theta0,size = 0.4,col="red3")+
  labs(title="", x=TeX("$\\theta"),
       y=TeX("SNormal"))+
  scale_x_continuous(expand = c(0, 0), breaks = seq(0,4,0.5)) + 
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,2,0.5), limits = c(0,2.2)) +
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5))

pdfn<- function (x) dnorm( (x-theta_hatD),sd = 1/sqrt(j_theta_hatD))
n_post_D=pdfn(thetas)
data_D_N <- data.frame(thetas=thetas, n_posterior=n_post_D)
p1DN <-  ggplot(data_D_N, aes(x=thetas,y=n_posterior))+
  geom_line(size = 0.4)+
  geom_vline(xintercept=m1_B,size = 0.4,col="blue3", linetype="twodash")+
  #  geom_vline(xintercept=theta0,size = 0.4,col="red3")+
  labs(title="", x=TeX("$\\theta"),
       y=TeX("$Normal-IO"))+
  scale_x_continuous(expand = c(0, 0), breaks = seq(0,4,0.5)) + 
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,2,0.5), limits = c(0,2.2)) +
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5))


grid.arrange( p1C,  p1D,p1CSN, p1DSN,  p1CN, p1DN,ncol = 2)


######
