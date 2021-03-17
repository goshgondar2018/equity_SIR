library(deSolve)
library(EnvStats)
library(reshape2)
library(ggpubr)
library(tidyverse)
library(invgamma)
library(patchwork)

# define parameters
CBR_tot <- (31.8/1000)
TFR_1<-6.4; TFR_2<-5.6; TFR_3<-4.9; TFR_4<-4.3; TFR_5<-2.6

CBR_1<-((5*CBR_tot)/(1+(TFR_2/TFR_1)+(TFR_3/TFR_1)+(TFR_4/TFR_1)+(TFR_5/TFR_1)))
CBR_2<-CBR_1*(TFR_2/TFR_1); CBR_3<-CBR_1*(TFR_3/TFR_1); CBR_4<-CBR_1*(TFR_4/TFR_1); CBR_5<-CBR_1*(TFR_5/TFR_1)

nu1<-CBR_1/365; nu2<-CBR_2/365; nu3<-CBR_3/365; nu4<-CBR_4/365; nu5<-CBR_5/365

CFR_1<- 0.0218; CFR_2<-0.0189; CFR_3<-0.0160; CFR_4<-0.0131; CFR_5<-0.0102 

p1<-0.432; p2<-0.499; p3<-0.544; p4<-0.586; p5<-0.743
vacc_efficacy<-0.85
mean_p=mean(c(p1,p2,p3,p4,p5))

mu1<-nu1; mu2<-nu2; mu3<-nu3; mu4<-nu4; mu5<-nu5

N1 <- 15000; N2 <- 15000; N3 <- 15000; N4 <- 15000; N5 <- 15000 
N_total <- N1+N2+N3+N4+N5

# ESIR model
dx.dt.ESIR <- function(t, y, parms) { #MAKE SURE YOU UPDATE THE NUMBER OF PEOPLE PER INCOME QUINTILE
  # AJ: first line only true when nu and mu are equal in quintile 
  dS1 <- (1-parms["p1"])*parms["nu1"]*(y["S1"]+y["I1"]+y["R1"])+
    - parms["beta11"]*(y["S1"]*y["I1"])+                     
    - parms["beta21"]*(y["S1"]*y["I2"])+
    - parms["beta31"]*(y["S1"]*y["I3"])+
    - parms["beta41"]*(y["S1"]*y["I4"])+
    - parms["beta51"]*(y["S1"]*y["I5"])+
    - parms["mu1"]*y["S1"]
  
  dI1 <-   parms["beta11"]*(y["S1"]*y["I1"])+
    + parms["beta21"]*(y["S1"]*y["I2"])+
    + parms["beta31"]*(y["S1"]*y["I3"])+
    + parms["beta41"]*(y["S1"]*y["I4"])+
    + parms["beta51"]*(y["S1"]*y["I5"])+
    - parms["gamma"]*y["I1"]+ 
    - (parms["mu1"]*y["I1"])
  
  dR1 <- (parms["p1"])*parms["nu1"]*(y["S1"]+y["I1"]+y["R1"])+
    + parms["gamma"]*(1-parms["d1"])*y["I1"]+
    - parms["mu1"]*y["R1"]
  
  dD1 <- parms["gamma"]*(parms["d1"])*y["I1"]
  
  dS2 <-  (1-parms["p2"])*parms["nu2"]*(y["S2"]+y["I2"]+y["R2"])+
    - parms["beta12"]*(y["S2"]*y["I1"])+
    - parms["beta22"]*(y["S2"]*y["I2"])+
    - parms["beta32"]*(y["S2"]*y["I3"])+
    - parms["beta42"]*(y["S2"]*y["I4"])+
    - parms["beta52"]*(y["S2"]*y["I5"])+
    - parms["mu2"]*y["S2"]
  
  dI2 <- parms["beta12"]*(y["S2"]*y["I1"])+ 
    + parms["beta22"]*(y["S2"]*y["I2"])+
    + parms["beta32"]*(y["S2"]*y["I3"])+ 
    + parms["beta42"]*(y["S2"]*y["I4"])+
    + parms["beta52"]*(y["S2"]*y["I5"])+
    - parms["gamma"]*y["I2"]+
    - (parms["mu2"]*y["I2"])
  
  dR2 <- (parms["p2"])*parms["nu2"]*(y["S2"]+y["I2"]+y["R2"])+
    + parms["gamma"]*(1-parms["d2"])*y["I2"]+
    - parms["mu2"]*y["R2"]
  
  dD2 <- parms["gamma"]*(parms["d2"])*y["I2"]
  
  dS3 <-   (1-parms["p3"])*parms["nu3"]*(y["S3"]+y["I3"]+y["R3"])+
    - parms["beta13"]*(y["S3"]*y["I1"])+
    - parms["beta23"]*(y["S3"]*y["I2"])+
    - parms["beta33"]*(y["S3"]*y["I3"])+
    - parms["beta43"]*(y["S3"]*y["I4"])+
    - parms["beta53"]*(y["S3"]*y["I5"])+
    - parms["mu3"]*y["S3"]
  
  dI3 <- parms["beta13"]*(y["S3"]*y["I1"])+ 
    + parms["beta23"]*(y["S3"]*y["I2"])+
    + parms["beta33"]*(y["S3"]*y["I3"])+
    + parms["beta43"]*(y["S3"]*y["I4"])+
    + parms["beta53"]*(y["S3"]*y["I5"])+
    - parms["gamma"]*y["I3"]+ 
    - (parms["mu3"]*y["I3"])
  
  dR3 <- (parms["p3"])*parms["nu3"]*(y["S3"]+y["I3"]+y["R3"])+
    + parms["gamma"]*(1-parms["d3"])*y["I3"]+
    - parms["mu3"]*y["R3"]
  
  dD3 <- parms["gamma"]*(parms["d3"])*y["I3"]
  
  dS4 <-  (1-parms["p4"])*parms["nu4"]*(y["S4"]+y["I4"]+y["R4"])+
    - parms["beta14"]*(y["S4"]*y["I1"])+
    - parms["beta24"]*(y["S4"]*y["I2"])+
    - parms["beta34"]*(y["S4"]*y["I3"])+
    - parms["beta44"]*(y["S4"]*y["I4"])+
    - parms["beta54"]*(y["S4"]*y["I5"])+
    - parms["mu4"]*y["S4"]
  
  dI4 <- parms["beta14"]*(y["S4"]*y["I1"])+ 
    + parms["beta24"]*(y["S4"]*y["I2"])+
    + parms["beta34"]*(y["S4"]*y["I3"])+ 
    + parms["beta44"]*(y["S4"]*y["I4"])+
    + parms["beta54"]*(y["S4"]*y["I5"])+
    - parms["gamma"]*y["I4"]+ 
    - (parms["mu4"]*y["I4"])
  
  dR4 <- (parms["p4"])*parms["nu4"]*(y["S4"]+y["I4"]+y["R4"])+
    + parms["gamma"]*(1-parms["d4"])*y["I4"]+
    - parms["mu4"]*y["R4"]
  
  dD4 <- parms["gamma"]*(parms["d4"])*y["I4"]
  
  dS5 <-  (1-parms["p5"])*parms["nu5"]*(y["S5"]+y["I5"]+y["R5"])+
    - parms["beta15"]*(y["S5"]*y["I1"])+
    - parms["beta25"]*(y["S5"]*y["I2"])+
    - parms["beta35"]*(y["S5"]*y["I3"])+
    - parms["beta45"]*(y["S5"]*y["I4"])+
    - parms["beta55"]*(y["S5"]*y["I5"])+
    - parms["mu5"]*y["S5"]
  
  dI5 <-parms["beta15"]*(y["S5"]*y["I1"])+
    + parms["beta25"]*(y["S5"]*y["I2"])+
    + parms["beta35"]*(y["S5"]*y["I3"])+ 
    + parms["beta45"]*(y["S5"]*y["I4"])+
    + parms["beta55"]*(y["S5"]*y["I5"])+
    - parms["gamma"]*y["I5"]+ 
    - (parms["mu5"]*y["I5"])
  
  dR5 <- (parms["p5"])*parms["nu5"]*(y["S5"]+y["I5"]+y["R5"])+
    + parms["gamma"]*(1-parms["d5"])*y["I5"]+
    - parms["mu5"]*y["R5"]
  
  dD5 <- parms["gamma"]*(parms["d5"])*y["I5"]
  
  return(list(c(dS1, dI1, dR1, dD1, dS2, dI2, dR2, dD2, dS3, dI3, dR3, dD3, dS4, dI4, dR4, dD4,
                dS5, dI5, dR5, dD5)))
  
}

# set time sequence, 400 days
times_vector <- seq(from=1, to=400, by=1)

# define function to run ESIR model for each transmission matrix scenario 
ESIR_transmission_matrix_scenarios<-function(transmission_matrix_scenario){
  n=1000 #1000
  df_all=data.frame(matrix(nrow=n*400,ncol=42))
  df_all[,2]<-rep(seq(1,400),n)
  colnames(df_all)<-c("iteration","days","I1","I2","I3","I4","I5",
                      "S1_stand","S2_stand","S3_stand","S4_stand","S5_stand",
                      "I1_stand","I2_stand","I3_stand","I4_stand","I5_stand",
                      "R1_stand","R2_stand","R3_stand","R4_stand","R5_stand",
                      "D1_stand","D2_stand","D3_stand","D4_stand","D5_stand",
                      "current_pop_1","current_pop_2","current_pop_3","current_pop_4","current_pop_5",
                      "D1","D2","D3","D4","D5","R1","R2","R3","R4","R5")
  n2=100
  
  for (i in 1:n){
    
    # sample from a uniform distribution for each quintile-specific birth rate parameter
    nu1_sample=runif(1,min=0.5*nu1,max=1.5*nu1)
    nu2_sample=runif(1,min=0.5*nu2,max=1.5*nu2) 
    nu3_sample=runif(1,min=0.5*nu3,max=1.5*nu3)
    nu4_sample=runif(1,min=0.5*nu4,max=1.5*nu4)
    nu5_sample=runif(1,min=0.5*nu5,max=1.5*nu5)
    nu_mean_sample=mean(nu1_sample,nu2_sample,nu3_sample,nu4_sample,nu5_sample)
    
    mu1_sample=nu1_sample
    mu2_sample=nu2_sample
    mu3_sample=nu3_sample
    mu4_sample=nu4_sample
    mu5_sample=nu5_sample
    mu_mean_sample=mean(mu1_sample,mu2_sample,mu3_sample,mu4_sample,mu5_sample)
    
    # sample from a uniform distribution for each CFR parameter
    CFR_1_sample=runif(1,min=0.5*CFR_1,max=1.5*CFR_1)
    CFR_2_sample=runif(1,min=0.5*CFR_2,max=1.5*CFR_2)
    CFR_3_sample=runif(1,min=0.5*CFR_3,max=1.5*CFR_3)
    CFR_4_sample=runif(1,min=0.5*CFR_4,max=1.5*CFR_4)
    CFR_5_sample=runif(1,min=0.5*CFR_5,max=1.5*CFR_5)
    CFR_mean_sample=mean(CFR_1_sample,CFR_2_sample,CFR_3_sample,CFR_4_sample,CFR_5_sample)
    
    # sample from a gamma distribution for the recovery rate
    gamma_sample=rinvgamma(1,shape=15)
    
    # sample from a uniform distribution for the vaccination coverage rate parameters
    p1_sample=rtri(1,min=0.5*p1,mode=p1,max=1.5*p1)
    p2_sample=rtri(1,min=0.5*p2,mode=p2,max=1.5*p2)
    p3_sample=rtri(1,min=0.5*p3,mode=p3,max=1.5*p3)
    p4_sample=rtri(1,min=0.5*p4,mode=p4,max=1.5*p4)
    p5_sample=rtri(1,min=0.5*p5,mode=p5,max=1.5*p5)
    
    df_all[(400*(i-1)+1):(400*(i-1)+400),1]<-rep(i,400)
    
    # sample transmission rate terms
    if (transmission_matrix_scenario == 'half within quintile'| transmission_matrix_scenario=='high within quintile'|transmission_matrix_scenario=='low within quintile'){
      if (transmission_matrix_scenario=='half within quintile'){
        all_R0=runif(5,min=10,max=22)
        Rii=all_R0/2
        Rii=sort(Rii,decreasing=TRUE)
        sum_Rij=all_R0/2
      }
      else if (transmission_matrix_scenario=='high within quintile') {
        all_R0=runif(5,min=10,max=22)
        Rii=0.9*all_R0
        Rii=sort(Rii,decreasing=TRUE)
        sum_Rij=0.1*all_R0
      }
      else {
        all_R0=runif(5,min=10,max=22)
        Rii=0.1*all_R0
        Rii=sort(Rii,decreasing=TRUE)
        sum_Rij=0.9*all_R0
      }
      
      sum_Rij=sort(sum_Rij,decreasing=TRUE)
      sum_Rij_rounded=round(sum_Rij,1) #round(sum_Rij,0)
      sum_Rij_rounded=sort(sum_Rij_rounded,decreasing=TRUE)
      
      # draw four random values for Rij using a gamma dist. 
      # with mean equal to the previously sampled sum of Rij terms for quintile 1 / n
      each_value_R1j_df=data.frame(matrix(nrow=n2,ncol=4))
      for (j in 1:n2){
        each_value=runif(4,1,sum_Rij[1]/2-1)
        each_value_R1j_df[j,]=each_value
      }
      
      # calculate sums for each 4 Rij sampled combination
      # select the combinations that yield the (rounded) sum of the previously sampled sum of Rij terms for quintile 1
      # that is, is most closely consistent with the previously designated sum of R1j terms
      each_value_R1j_df$sum=rowSums(each_value_R1j_df)
      each_value_R1j_df$sum_rounded=round(each_value_R1j_df$sum,1) #round(each_value_R1j_df$sum,0)
      R1j=each_value_R1j_df[which(each_value_R1j_df$sum_rounded==sum_Rij_rounded[1]),1:4]
      R1j_sorted=t(apply(R1j,1,sort))
      
      # identify the optimal combination of R1j terms as that which has the smallest distance of the sum of R1j terms 
      # from the previously designated sum of R1j terms (if there is more than 1 combination in R1j_sorted)
      # define each of the optimal R1j terms from this selected combination
      if(length(R1j_sorted)>0){
        R1j_optimal_index=rowSums(R1j_sorted)[which(rowSums(R1j_sorted)-sum_Rij[1]==min(rowSums(R1j_sorted)-sum_Rij[1]))][1]
        if (nrow(R1j_sorted)>1){
          R12_optimal=R1j_sorted[,4][names(R1j_optimal_index)][[1]]
          R13_optimal=R1j_sorted[,3][names(R1j_optimal_index)][[1]]
          R14_optimal=R1j_sorted[,2][names(R1j_optimal_index)][[1]]
          R15_optimal=R1j_sorted[,1][names(R1j_optimal_index)][[1]]
        }
        else {
          R12_optimal=R1j_sorted[,4]
          R13_optimal=R1j_sorted[,3]
          R14_optimal=R1j_sorted[,2]
          R15_optimal=R1j_sorted[,1]
        }
      }
      
      # draw three random values for R1j using a gamma dist. 
      # with mean equal to the previously sampled sum of Rij terms for quintile 1 / n - previously defined optimal R12
      # (noting the symmetry R12=R21)
      each_value_R2j_df=data.frame(matrix(nrow=n2,ncol=3))
      for (j in 1:n2){
        if(sum_Rij[2]-R12_optimal>0) {
          each_value=runif(3,min=1,max=((2/3)*(sum_Rij[2]-R12_optimal))-1)
          each_value_R2j_df[j,1:3]=each_value
        }
      }
      
      each_value_R2j_df$R21=R12_optimal
      each_value_R2j_df$sum<-rowSums(each_value_R2j_df)
      each_value_R2j_df$sum_rounded=round(each_value_R2j_df$sum,1) #round(each_value_R2j_df$sum,0)
      R2j=each_value_R2j_df[which(each_value_R2j_df$sum_rounded==sum_Rij_rounded[2]),1:4]
      R2j_sorted=t(apply(R2j,1,sort))
      
      if(length(R2j_sorted)>0){
        R2j_optimal_index=rowSums(R2j_sorted)[which(rowSums(R2j_sorted)-sum_Rij[2]==min(rowSums(R2j_sorted)-sum_Rij[2]))][1]
        if (nrow(R2j_sorted)>1){
          R23_optimal=R2j_sorted[,3][names(R2j_optimal_index)][[1]]
          R24_optimal=R2j_sorted[,2][names(R2j_optimal_index)][[1]]
          R25_optimal=R2j_sorted[,1][names(R2j_optimal_index)][[1]]
        }
        else {
          R23_optimal=R2j_sorted[,3]
          R24_optimal=R2j_sorted[,2]
          R25_optimal=R2j_sorted[,1]
        }
      }
      
      # draw two random values for R3j using a gamma dist. 
      # with mean equal to the previously sampled sum of Rij terms for quintile 1 / n - previously defined optimal R23 and R13
      # (noting the symmetry R23=R32 and R13=R31)
      each_value_R3j_df=data.frame(matrix(nrow=n2,ncol=2))
      for (j in 1:n2){
        if(sum_Rij[3]-R23_optimal-R13_optimal>0) {
          each_value=runif(2,min=1,max=(sum_Rij[3]-R23_optimal-R13_optimal)-1)
          each_value_R3j_df[j,1:2]=each_value
        }
      }
      
      each_value_R3j_df$R31=R13_optimal
      each_value_R3j_df$R32=R23_optimal
      
      each_value_R3j_df$sum<-rowSums(each_value_R3j_df)
      each_value_R3j_df$sum_rounded=round(each_value_R3j_df$sum,1) #round(each_value_R3j_df$sum,1)
      R3j=each_value_R3j_df[which(each_value_R3j_df$sum_rounded==sum_Rij_rounded[3]),1:4]
      R3j_sorted=t(apply(R3j,1,sort))
      
      if(length(R3j_sorted)>0){
        R3j_optimal_index=rowSums(R3j_sorted)[which(rowSums(R3j_sorted)-sum_Rij[3]==min(rowSums(R3j_sorted)-sum_Rij[3]))][1]
        if (nrow(R3j_sorted)>1){
          R34_optimal=R3j_sorted[,2][names(R3j_optimal_index)][[1]]
          R35_optimal=R3j_sorted[,1][names(R3j_optimal_index)][[1]]
        }
        else {
          R34_optimal=R3j_sorted[,2]
          R35_optimal=R3j_sorted[,1]
        }
      }
      
      # draw two random values for R4j using a gamma dist. 
      # with mean equal to the previously sampled sum of Rij terms for quintile 1 / n - previously defined optimal R34 and R24 and R14
      # (noting the symmetry R34=R43 and R24=R42 and R14=R41)
      each_value_R4j_df=data.frame(matrix(nrow=n2,ncol=1))
      for (j in 1:n2){
        if(sum_Rij[4]-R34_optimal-R24_optimal-R14_optimal>0) {
          each_value=runif(1,min=1,max=(2*(sum_Rij[4]-R34_optimal-R24_optimal-R14_optimal))-1)
          each_value_R4j_df[j,1]=each_value
        }
      }
      
      each_value_R4j_df$R41=R14_optimal
      each_value_R4j_df$R42=R24_optimal
      each_value_R4j_df$R43=R34_optimal
      
      each_value_R4j_df$sum<-rowSums(each_value_R4j_df)
      each_value_R4j_df$sum_rounded=round(each_value_R4j_df$sum,1) #round(each_value_R4j_df$sum,1)
      
      R4j=each_value_R4j_df[which(each_value_R4j_df$sum_rounded==sum_Rij_rounded[4]),1:4]
      R4j_sorted=t(apply(R4j,1,sort))
      
      if(length(R4j_sorted)>0){
        R4j_optimal_index=rowSums(R4j_sorted)[which(rowSums(R4j_sorted)-sum_Rij[4]==min(rowSums(R4j_sorted)-sum_Rij[4]))][1]
        if (nrow(R4j_sorted)>1){
          R45_optimal=R4j_sorted[,1][names(R4j_optimal_index)][[1]]
        }
        else{
          R45_optimal=R4j_sorted[,1]
        }
      }
      
      betas<-matrix(c(beta11=(Rii[1]*(gamma_sample))/N1, 
                      beta21=(R12_optimal*(gamma_sample))/N1, 
                      beta31=(R13_optimal*(gamma_sample))/N1, 
                      beta41=(R14_optimal*(gamma_sample))/N1, 
                      beta51=(R15_optimal*(gamma_sample))/N1,
                      beta12=(R12_optimal*(gamma_sample))/N2, 
                      beta22=(Rii[2]*(gamma_sample))/N2,
                      beta32=(R23_optimal*(gamma_sample))/N2, 
                      beta42=(R24_optimal*(gamma_sample))/N2, 
                      beta52=(R25_optimal*(gamma_sample))/N2,
                      beta13=(R13_optimal*(gamma_sample))/N3,
                      beta23=(R23_optimal*(gamma_sample))/N3, 
                      beta33=(Rii[3]*(gamma_sample))/N3, 
                      beta43=(R34_optimal*(gamma_sample))/N3, 
                      beta53=(R35_optimal*(gamma_sample))/N3,
                      beta14=(R14_optimal*(gamma_sample))/N4, 
                      beta24=(R24_optimal*(gamma_sample))/N4, 
                      beta34=(R34_optimal*(gamma_sample))/N4, 
                      beta44=(Rii[4]*(gamma_sample))/N4, 
                      beta54=(R45_optimal*(gamma_sample))/N4,
                      beta15=(R15_optimal*(gamma_sample))/N5, 
                      beta25=(R25_optimal*(gamma_sample))/N5, 
                      beta35=(R35_optimal*(gamma_sample))/N5, 
                      beta45=(R45_optimal*(gamma_sample))/N5, 
                      beta55=(Rii[5]*(gamma_sample))/N5),
                    nrow=5,ncol=5,byrow=TRUE)
      
    }
    else if (transmission_matrix_scenario=='Mexico contact survey'){
      
      f<-0.03
      
      L11 <- 1.25; L12 <- 0.9; L13 <- 0.8; L14 <- 0.75; L15 <- 0.76
      L22 <- 1.16; L23 <- 1.05; L24 <- 0.90; L25 <- 0.88 
      L33 <- 1.25; L34 <- 1.25; L35 <- 1.1; L44 <- 1.49 
      L45 <- 1.4; L55 <- 2.3 
      
      beta <- 16*gamma_sample
      
      num_c11 <- L11*beta*(N1^2+N2^2+N3^2+N4^2+N5^2+
                             2*N1*N2+2*N1*N3+2*N1*N4+2*N1*N5+
                             2*N2*N3+2*N2*N4+2*N2*N5+
                             2*N3*N4+2*N3*N5+
                             2*N4*N5)
      denom_c11 <- L11*N1^2*f + N2^2*f*L22+N3^2*f*L33+N4^2*f*L44+N5^2*f*L55+
        N1*N2*(f+f)*L12+N1*N3*(f+f)*L13+N1*N4*(f+f)*L14+N1*N5*(f+f)*L15+
        N2*N3*(f+f)*L23+N2*N4*(f+f)*L24+N2*N5*(f+f)*L25+
        N3*N4*(f+f)*L34+N3*N5*(f+f)*L35+
        N4*N5*(f+f)*L45
      
      c11 <- num_c11/denom_c11; c12 <- c11*(L12/L11); c13 <- c11*(L13/L11); c14 <- c11*(L14/L11); c15 <- c11*(L15/L11)
      c22 <- c11*(L22/L11); c23 <- c11*(L23/L11); c24 <- c11*(L24/L11); c25 <- c11*(L25/L11)
      c33 <- c11*(L33/L11); c44 <- c11*(L44/L11); c45 <- c11*(L45/L11)
      c34 <- c11*(L34/L11); c35 <- c11*(L35/L11) ; c55 <- c11*(L55/L11)
      
      beta11 <- (f*c11)/N1; beta12 <- (f*c12)/N2; beta13 <- (f*c13)/N3
      beta14 <- (f*c14)/N4; beta15 <- (f*c15)/N5
      
      beta21 <- (f*c12)/N1; beta22 <- (f*c22)/N2; beta23 <- (f*c23)/N3
      beta24 <- (f*c24)/N4; beta25 <- (f*c25)/N5
      
      beta31 <- (f*c13)/N1; beta32 <- (f*c23)/N2; beta33 <- (f*c33)/N3
      beta34 <- (f*c34)/N4; beta35 <- (f*c35)/N5
      
      beta41 <- (f*c14)/N1; beta42 <- (f*c24)/N2; beta43 <- (f*c34)/N3
      beta44 <- (f*c44)/N4; beta45 <- (f*c45)/N5
      
      beta51 <- (f*c15)/N1; beta52 <- (f*c25)/N2; beta53 <- (f*c35)/N3
      beta54 <- (f*c45)/N4; beta55 <- (f*c55)/N5
      
      betas<-matrix(c(beta11=beta11, beta21=beta21, beta31=beta31, 
                      beta41=beta41, beta51=beta51,
                      beta12=beta12, beta22=beta22, beta32=beta32, 
                      beta42=beta42, beta52=beta52,
                      beta13=beta13, beta23=beta23, beta33=beta33, 
                      beta43=beta43, beta53=beta53,
                      beta14=beta14, beta24=beta24, beta34=beta34, 
                      beta44=beta44, beta54=beta54,
                      beta15=beta15, beta25=beta25, beta35=beta35, 
                      beta45=beta45, beta55=beta55),
                    nrow=5,ncol=5)
      
    }
    
    else {
      R0_sample=runif(5,min=10,max=22)
      R0_3=sort(R0_sample,decreasing=TRUE)[3]
      betas<-matrix(c(beta11=(R0_3*(gamma_sample))/N_total,beta21=beta11,beta31=beta11,beta41=beta11,beta51=beta11,
                        beta12=beta11,beta22=beta11,beta32=beta11,beta42=beta11,beta52=beta11,beta13=beta11,
                        beta23=beta11,beta33=beta11,beta43=beta11,beta53=beta11,beta14=beta11,beta24=beta11,
                        beta34=beta11,beta44=beta11,beta54=beta11,beta15=beta11,beta25=beta11,beta35=beta11,
                        beta45=beta11,beta55=beta11),
                    nrow=5,ncol=5,byrow=TRUE)
    }
    
    
    # set the baseline parameters
    parms_scenario<-c(beta11=betas[1,1], beta21=betas[1,2], beta31=betas[1,3], beta41=betas[1,4], beta51=betas[1,5],
                      beta12=betas[2,1], beta22=betas[2,2], beta32=betas[2,3], beta42=betas[2,4], beta52=betas[2,5],
                      beta13=betas[3,1], beta23=betas[3,2], beta33=betas[3,3], beta43=betas[3,4], beta53=betas[3,5],
                      beta14=betas[4,1], beta24=betas[4,2], beta34=betas[4,3], beta44=betas[4,4], beta54=betas[4,5],
                      beta15=betas[5,1], beta25=betas[5,2], beta35=betas[5,3], beta45=betas[5,4], beta55=betas[5,5],
                      nu=c(nu1_sample,nu2_sample,nu3_sample,nu4_sample,nu5_sample), 
                      mu=c(mu1_sample,mu2_sample,mu3_sample,mu4_sample,mu5_sample),
                      d=c(CFR_1_sample,CFR_2_sample,CFR_3_sample,CFR_4_sample,CFR_5_sample), 
                      gamma= gamma_sample,
                      p=c(vacc_efficacy*p1_sample,vacc_efficacy*p2_sample,vacc_efficacy*p3_sample,vacc_efficacy*p4_sample,vacc_efficacy*p5_sample))
    
    ## Set initial conditions for S, I, R and N for each income quintile
    START.S1 <- (N1*(1-(p1*vacc_efficacy)))-1; START.I1 <- 1; START.R1 <- (N1*p1*vacc_efficacy); START.D1 <- 0
    
    START.S2 <- (N2*(1-(p2*vacc_efficacy)))-1; START.I2 <- 1; START.R2 <- (N2*p2*vacc_efficacy); START.D2 <- 0
    
    START.S3 <- (N3*(1-(p3*vacc_efficacy)))-1; START.I3 <- 1; START.R3 <- (N3*p3*vacc_efficacy); START.D3 <- 0
    
    START.S4 <- (N4*(1-(p4*vacc_efficacy)))-1; START.I4 <- 1; START.R4 <- (N4*p4*vacc_efficacy); START.D4 <- 0
    
    START.S5 <- (N5*(1-(p5*vacc_efficacy)))-1; START.I5 <- 1; START.R5 <- (N5*p5*vacc_efficacy); START.D5 <- 0
    
    START.S_ALL=START.S1+START.S2+START.S3+START.S4+START.S5
    START.I_ALL=START.I1+START.I2+START.I3+START.I4+START.I5
    START.R_ALL=START.R1+START.R2+START.R3+START.R4+START.R5
    START.D_ALL=START.D1+START.D2+START.D3+START.D4+START.D5
    
    if (transmission_matrix_scenario=='homogeneous-base case'){
      ESIR.output <- lsoda(y=c(S1=START.S_ALL, I1=START.I_ALL, R1=START.R_ALL, D1=START.D_ALL,
                               S2=START.S_ALL, I2=START.I_ALL, R2=START.R_ALL, D2=START.D_ALL,
                               S3=START.S_ALL, I3=START.I_ALL, R3=START.R_ALL, D3=START.D_ALL,
                               S4=START.S_ALL, I4=START.I_ALL, R4=START.R_ALL, D4=START.D_ALL,
                               S5=START.S_ALL, I5=START.I_ALL, R5=START.R_ALL, D5=START.D_ALL),
                           times=times_vector, 
                           func=dx.dt.ESIR, 
                           parms=parms_scenario)
    }
    else {
      ESIR.output <- lsoda(y=c(S1=START.S1, I1=START.I1, R1=START.R1, D1=START.D1,
                               S2=START.S2, I2=START.I2, R2=START.R2, D2=START.D2,
                               S3=START.S3, I3=START.I3, R3=START.R3, D3=START.D3,
                               S4=START.S4, I4=START.I4, R4=START.R4, D4=START.D4,
                               S5=START.S5, I5=START.I5, R5=START.R5, D5=START.D5), 
                           times=times_vector, 
                           func=dx.dt.ESIR, 
                           parms=parms_scenario)
    }
    # run ESIR models
    
    
    df_all[(400*(i-1)+1):(400*(i-1)+400),3:7]=ESIR.output[,c("I1","I2","I3","I4","I5")]
    
    # standardize size of each quintile by compartment to proportion of current population for each quantile
    ESIR.output_quintile1=data.frame(ESIR.output[,c("S1","I1","R1","D1")])
    ESIR.output_quintile1$pop_current=rowSums(ESIR.output_quintile1)
    ESIR.output_quintile2=data.frame(ESIR.output[,c("S2","I2","R2","D2")])
    ESIR.output_quintile2$pop_current=rowSums(ESIR.output_quintile2)
    ESIR.output_quintile3=data.frame(ESIR.output[,c("S3","I3","R3","D3")])
    ESIR.output_quintile3$pop_current=rowSums(ESIR.output_quintile3)
    ESIR.output_quintile4=data.frame(ESIR.output[,c("S4","I4","R4","D4")])
    ESIR.output_quintile4$pop_current=rowSums(ESIR.output_quintile4)
    ESIR.output_quintile5=data.frame(ESIR.output[,c("S5","I5","R5","D5")])
    ESIR.output_quintile5$pop_current=rowSums(ESIR.output_quintile5)
    
    
    df_all[(400*(i-1)+1):(400*(i-1)+400),8:42]=c(ESIR.output[,"S1"]/ESIR.output_quintile1$pop_current,
                                                 ESIR.output[,"S2"]/ESIR.output_quintile2$pop_current,
                                                 ESIR.output[,"S3"]/ESIR.output_quintile3$pop_current,
                                                 ESIR.output[,"S4"]/ESIR.output_quintile4$pop_current,
                                                 ESIR.output[,"S5"]/ESIR.output_quintile5$pop_current,
                                                 ESIR.output[,"I1"]/ESIR.output_quintile1$pop_current,
                                                 ESIR.output[,"I2"]/ESIR.output_quintile2$pop_current,
                                                 ESIR.output[,"I3"]/ESIR.output_quintile3$pop_current,
                                                 ESIR.output[,"I4"]/ESIR.output_quintile4$pop_current,
                                                 ESIR.output[,"I5"]/ESIR.output_quintile5$pop_current,
                                                 ESIR.output[,"R1"]/ESIR.output_quintile1$pop_current,
                                                 ESIR.output[,"R2"]/ESIR.output_quintile2$pop_current,
                                                 ESIR.output[,"R3"]/ESIR.output_quintile3$pop_current,
                                                 ESIR.output[,"R4"]/ESIR.output_quintile4$pop_current,
                                                 ESIR.output[,"R5"]/ESIR.output_quintile5$pop_current,
                                                 ESIR.output[,"D1"]/ESIR.output_quintile1$pop_current,
                                                 ESIR.output[,"D2"]/ESIR.output_quintile2$pop_current,
                                                 ESIR.output[,"D3"]/ESIR.output_quintile3$pop_current,
                                                 ESIR.output[,"D4"]/ESIR.output_quintile4$pop_current,
                                                 ESIR.output[,"D5"]/ESIR.output_quintile5$pop_current,
                                                 ESIR.output_quintile1$pop_current,
                                                 ESIR.output_quintile2$pop_current,
                                                 ESIR.output_quintile3$pop_current,
                                                 ESIR.output_quintile4$pop_current,
                                                 ESIR.output_quintile5$pop_current,
                                                 ESIR.output[,"D1"], ESIR.output[,"D2"], ESIR.output[,"D3"],
                                                 ESIR.output[,"D4"], ESIR.output[,"D5"], ESIR.output[,"R1"],
                                                 ESIR.output[,"R2"], ESIR.output[,"R3"], ESIR.output[,"R4"],
                                                 ESIR.output[,"R5"])
  }
  return(list(df_all, betas))
}

ESIR_vaccination_scenarios<-function(vaccination_scenario_list){
  n=1000 #1000
  
  n2=100
  
  all_parameter_output=data.frame(matrix(nrow=48,ncol=6000))
  for (i in 1:n){
    all_R0=runif(5,min=10,max=22)
    Rii=all_R0/2
    Rii=sort(Rii,decreasing=TRUE)
    sum_Rij=all_R0/2
    
    sum_Rij=sort(sum_Rij,decreasing=TRUE)
    sum_Rij_rounded=round(sum_Rij,1) #round(sum_Rij,0)
    sum_Rij_rounded=sort(sum_Rij_rounded,decreasing=TRUE)
    
    # draw four random values for Rij using a gamma dist. 
    # with mean equal to the previously sampled sum of Rij terms for quintile 1 / n
    each_value_R1j_df=data.frame(matrix(nrow=n2,ncol=4))
    for (j in 1:n2){
      each_value=runif(4,1,sum_Rij[1]/2-1)
      each_value_R1j_df[j,]=each_value
    }
    
    # calculate sums for each 4 Rij sampled combination
    # select the combinations that yield the (rounded) sum of the previously sampled sum of Rij terms for quintile 1
    # that is, is most closely consistent with the previously designated sum of R1j terms
    each_value_R1j_df$sum=rowSums(each_value_R1j_df)
    each_value_R1j_df$sum_rounded=round(each_value_R1j_df$sum,1) #round(each_value_R1j_df$sum,0)
    R1j=each_value_R1j_df[which(each_value_R1j_df$sum_rounded==sum_Rij_rounded[1]),1:4]
    R1j_sorted=t(apply(R1j,1,sort))
    
    # identify the optimal combination of R1j terms as that which has the smallest distance of the sum of R1j terms 
    # from the previously designated sum of R1j terms (if there is more than 1 combination in R1j_sorted)
    # define each of the optimal R1j terms from this selected combination
    if(length(R1j_sorted)>0){
      R1j_optimal_index=rowSums(R1j_sorted)[which(rowSums(R1j_sorted)-sum_Rij[1]==min(rowSums(R1j_sorted)-sum_Rij[1]))][1]
      if (nrow(R1j_sorted)>1){
        R12_optimal=R1j_sorted[,4][names(R1j_optimal_index)][[1]]
        R13_optimal=R1j_sorted[,3][names(R1j_optimal_index)][[1]]
        R14_optimal=R1j_sorted[,2][names(R1j_optimal_index)][[1]]
        R15_optimal=R1j_sorted[,1][names(R1j_optimal_index)][[1]]
      }
      else {
        R12_optimal=R1j_sorted[,4]
        R13_optimal=R1j_sorted[,3]
        R14_optimal=R1j_sorted[,2]
        R15_optimal=R1j_sorted[,1]
      }
    }
    
    # draw three random values for R1j using a gamma dist. 
    # with mean equal to the previously sampled sum of Rij terms for quintile 1 / n - previously defined optimal R12
    # (noting the symmetry R12=R21)
    each_value_R2j_df=data.frame(matrix(nrow=n2,ncol=3))
    for (j in 1:n2){
      if(sum_Rij[2]-R12_optimal>0) {
        each_value=runif(3,min=1,max=((2/3)*(sum_Rij[2]-R12_optimal))-1)
        each_value_R2j_df[j,1:3]=each_value
      }
    }
    
    each_value_R2j_df$R21=R12_optimal
    each_value_R2j_df$sum<-rowSums(each_value_R2j_df)
    each_value_R2j_df$sum_rounded=round(each_value_R2j_df$sum,1) #round(each_value_R2j_df$sum,0)
    R2j=each_value_R2j_df[which(each_value_R2j_df$sum_rounded==sum_Rij_rounded[2]),1:4]
    R2j_sorted=t(apply(R2j,1,sort))
    
    if(length(R2j_sorted)>0){
      R2j_optimal_index=rowSums(R2j_sorted)[which(rowSums(R2j_sorted)-sum_Rij[2]==min(rowSums(R2j_sorted)-sum_Rij[2]))][1]
      if (nrow(R2j_sorted)>1){
        R23_optimal=R2j_sorted[,3][names(R2j_optimal_index)][[1]]
        R24_optimal=R2j_sorted[,2][names(R2j_optimal_index)][[1]]
        R25_optimal=R2j_sorted[,1][names(R2j_optimal_index)][[1]]
      }
      else {
        R23_optimal=R2j_sorted[,3]
        R24_optimal=R2j_sorted[,2]
        R25_optimal=R2j_sorted[,1]
      }
    }
    
    # draw two random values for R3j using a gamma dist. 
    # with mean equal to the previously sampled sum of Rij terms for quintile 1 / n - previously defined optimal R23 and R13
    # (noting the symmetry R23=R32 and R13=R31)
    each_value_R3j_df=data.frame(matrix(nrow=n2,ncol=2))
    for (j in 1:n2){
      if(sum_Rij[3]-R23_optimal-R13_optimal>0) {
        each_value=runif(2,min=1,max=(sum_Rij[3]-R23_optimal-R13_optimal)-1)
        each_value_R3j_df[j,1:2]=each_value
      }
    }
    
    each_value_R3j_df$R31=R13_optimal
    each_value_R3j_df$R32=R23_optimal
    
    each_value_R3j_df$sum<-rowSums(each_value_R3j_df)
    each_value_R3j_df$sum_rounded=round(each_value_R3j_df$sum,1) #round(each_value_R3j_df$sum,1)
    R3j=each_value_R3j_df[which(each_value_R3j_df$sum_rounded==sum_Rij_rounded[3]),1:4]
    R3j_sorted=t(apply(R3j,1,sort))
    
    if(length(R3j_sorted)>0){
      R3j_optimal_index=rowSums(R3j_sorted)[which(rowSums(R3j_sorted)-sum_Rij[3]==min(rowSums(R3j_sorted)-sum_Rij[3]))][1]
      if (nrow(R3j_sorted)>1){
        R34_optimal=R3j_sorted[,2][names(R3j_optimal_index)][[1]]
        R35_optimal=R3j_sorted[,1][names(R3j_optimal_index)][[1]]
      }
      else {
        R34_optimal=R3j_sorted[,2]
        R35_optimal=R3j_sorted[,1]
      }
    }
    
    # draw two random values for R4j using a gamma dist. 
    # with mean equal to the previously sampled sum of Rij terms for quintile 1 / n - previously defined optimal R34 and R24 and R14
    # (noting the symmetry R34=R43 and R24=R42 and R14=R41)
    each_value_R4j_df=data.frame(matrix(nrow=n2,ncol=1))
    for (j in 1:n2){
      if(sum_Rij[4]-R34_optimal-R24_optimal-R14_optimal>0) {
        each_value=runif(1,min=1,max=(2*(sum_Rij[4]-R34_optimal-R24_optimal-R14_optimal))-1)
        each_value_R4j_df[j,1]=each_value
      }
    }
    
    each_value_R4j_df$R41=R14_optimal
    each_value_R4j_df$R42=R24_optimal
    each_value_R4j_df$R43=R34_optimal
    
    each_value_R4j_df$sum<-rowSums(each_value_R4j_df)
    each_value_R4j_df$sum_rounded=round(each_value_R4j_df$sum,1) #round(each_value_R4j_df$sum,1)
    
    R4j=each_value_R4j_df[which(each_value_R4j_df$sum_rounded==sum_Rij_rounded[4]),1:4]
    R4j_sorted=t(apply(R4j,1,sort))
    
    if(length(R4j_sorted)>0){
      R4j_optimal_index=rowSums(R4j_sorted)[which(rowSums(R4j_sorted)-sum_Rij[4]==min(rowSums(R4j_sorted)-sum_Rij[4]))][1]
      if (nrow(R4j_sorted)>1){
        R45_optimal=R4j_sorted[,1][names(R4j_optimal_index)][[1]]
      }
      else{
        R45_optimal=R4j_sorted[,1]
      }
    }
    
    # sample from a uniform distribution for each quintile-specific birth rate parameter
    nu1_sample=runif(1,min=0.5*nu1,max=1.5*nu1)#rtri(1,min=0,max=0.1/365,mode=nu1)
    nu2_sample=runif(1,min=0.5*nu2,max=1.5*nu2) #rtri(1,min=0,max=0.1/365,mode=nu2)
    nu3_sample=runif(1,min=0.5*nu3,max=1.5*nu3)#rtri(1,min=0,max=0.1/365,mode=nu3)
    nu4_sample=runif(1,min=0.5*nu4,max=1.5*nu4)#rtri(1,min=0,max=0.1/365,mode=nu4)
    nu5_sample=runif(1,min=0.5*nu5,max=1.5*nu5)#rtri(1,min=0,max=0.1/365,mode=nu5)
    nu_mean_sample=mean(nu1_sample,nu2_sample,nu3_sample,nu4_sample,nu5_sample)
    
    mu1_sample=nu1_sample
    mu2_sample=nu2_sample
    mu3_sample=nu3_sample
    mu4_sample=nu4_sample
    mu5_sample=nu5_sample
    mu_mean_sample=mean(mu1_sample,mu2_sample,mu3_sample,mu4_sample,mu5_sample)
    
    # sample from a triangular distribution for each CFR parameter
    CFR_1_sample=runif(1,min=0.5*CFR_1,max=1.5*CFR_1)#rtri(1,min=0,max=0.1,mode=CFR_1)
    CFR_2_sample=runif(1,min=0.5*CFR_2,max=1.5*CFR_2)#rtri(1,min=0,max=0.1,mode=CFR_2)
    CFR_3_sample=runif(1,min=0.5*CFR_3,max=1.5*CFR_3)#rtri(1,min=0,max=0.1,mode=CFR_3)
    CFR_4_sample=runif(1,min=0.5*CFR_4,max=1.5*CFR_4)#rtri(1,min=0,max=0.1,mode=CFR_4)
    CFR_5_sample=runif(1,min=0.5*CFR_5,max=1.5*CFR_5)#rtri(1,min=0,max=0.1,mode=CFR_5)
    CFR_mean_sample=mean(CFR_1_sample,CFR_2_sample,CFR_3_sample,CFR_4_sample,CFR_5_sample)
    
    # sample from a gamma distribution for the recovery rate
    gamma_sample=rinvgamma(1,shape=15)
    
    
    # define the beta matrix as (Rij*(gamma+CFR_i+(1-CFR_i)*mu_i))/N_i
    
    R0_3=sort(all_R0,decreasing=TRUE)[3]
    
    N_total=N1+N2+N3+N4+N5
    betas_homogeneous<-matrix(c(beta11=(R0_3*(gamma_sample))/N_total, 
                                beta21=(R0_3*(gamma_sample))/N_total,
                                beta31=(R0_3*(gamma_sample))/N_total, 
                                beta41=(R0_3*(gamma_sample))/N_total,
                                beta51=(R0_3*(gamma_sample))/N_total,
                                beta12=(R0_3*(gamma_sample))/N_total,
                                beta22=(R0_3*(gamma_sample))/N_total,
                                beta32=(R0_3*(gamma_sample))/N_total,
                                beta42=(R0_3*(gamma_sample))/N_total,
                                beta52=(R0_3*(gamma_sample))/N_total,
                                beta13=(R0_3*(gamma_sample))/N_total,
                                beta23=(R0_3*(gamma_sample))/N_total,
                                beta33=(R0_3*(gamma_sample))/N_total,
                                beta43=(R0_3*(gamma_sample))/N_total,
                                beta53=(R0_3*(gamma_sample))/N_total,
                                beta14=(R0_3*(gamma_sample))/N_total, 
                                beta24=(R0_3*(gamma_sample))/N_total,
                                beta34=(R0_3*(gamma_sample))/N_total,
                                beta44=(R0_3*(gamma_sample))/N_total,
                                beta54=(R0_3*(gamma_sample))/N_total,
                                beta15=(R0_3*(gamma_sample))/N_total,
                                beta25=(R0_3*(gamma_sample))/N_total,
                                beta35=(R0_3*(gamma_sample))/N_total,
                                beta45=(R0_3*(gamma_sample))/N_total,
                                beta55=(R0_3*(gamma_sample))/N_total),
                              nrow=5,ncol=5,byrow=TRUE)
    
    # set the baseline parameters
    parms_baseline_all=data.frame(matrix(nrow=47,ncol=7))
    for (vaccination_scenario in vaccination_scenario_list){
      if (vaccination_scenario == 'base case'){
        
        parms_baseline<-c(beta11=betas[1,1], beta21=betas[1,2], beta31=betas[1,3], beta41=betas[1,4], beta51=betas[1,5],
                          beta12=betas[2,1], beta22=betas[2,2], beta32=betas[2,3], beta42=betas[2,4], beta52=betas[2,5],
                          beta13=betas[3,1], beta23=betas[3,2], beta33=betas[3,3], beta43=betas[3,4], beta53=betas[3,5],
                          beta14=betas[4,1], beta24=betas[4,2], beta34=betas[4,3], beta44=betas[4,4], beta54=betas[4,5],
                          beta15=betas[5,1], beta25=betas[5,2], beta35=betas[5,3], beta45=betas[5,4], beta55=betas[5,5],
                          nu=c(nu1_sample,nu2_sample,nu3_sample,nu4_sample,nu5_sample), 
                          mu=c(mu1_sample,mu2_sample,mu3_sample,mu4_sample,mu5_sample),
                          d=c(CFR_1_sample,CFR_2_sample,CFR_3_sample,CFR_4_sample,CFR_5_sample), 
                          gamma= gamma_sample,
                          p=c(vacc_efficacy*p1,vacc_efficacy*p2,vacc_efficacy*p3,vacc_efficacy*p4,vacc_efficacy*p5),
                          scenario=vaccination_scenario)
      }
      else if (vaccination_scenario == 'mean vaccination coverage'){
        betas<-matrix(c(beta11=(Rii[1]*(gamma_sample))/N1, 
                        beta21=(R12_optimal*(gamma_sample))/N1, 
                        beta31=(R13_optimal*(gamma_sample))/N1, 
                        beta41=(R14_optimal*(gamma_sample))/N1, 
                        beta51=(R15_optimal*(gamma_sample))/N1,
                        beta12=(R12_optimal*(gamma_sample))/N2, 
                        beta22=(Rii[2]*(gamma_sample))/N2,
                        beta32=(R23_optimal*(gamma_sample))/N2, 
                        beta42=(R24_optimal*(gamma_sample))/N2, 
                        beta52=(R25_optimal*(gamma_sample))/N2,
                        beta13=(R13_optimal*(gamma_sample))/N3,
                        beta23=(R23_optimal*(gamma_sample))/N3, 
                        beta33=(Rii[3]*(gamma_sample))/N3, 
                        beta43=(R34_optimal*(gamma_sample))/N3, 
                        beta53=(R35_optimal*(gamma_sample))/N3,
                        beta14=(R14_optimal*(gamma_sample))/N4, 
                        beta24=(R24_optimal*(gamma_sample))/N4, 
                        beta34=(R34_optimal*(gamma_sample))/N4, 
                        beta44=(Rii[4]*(gamma_sample))/N4, 
                        beta54=(R45_optimal*(gamma_sample))/N4,
                        beta15=(R15_optimal*(gamma_sample))/N5, 
                        beta25=(R25_optimal*(gamma_sample))/N5, 
                        beta35=(R35_optimal*(gamma_sample))/N5, 
                        beta45=(R45_optimal*(gamma_sample))/N5, 
                        beta55=(Rii[5]*(gamma_sample))/N5),
                      nrow=5,ncol=5,byrow=TRUE)
        
        parms_baseline<-c(beta11=betas[1,1], beta21=betas[1,2], beta31=betas[1,3], beta41=betas[1,4], beta51=betas[1,5],
                          beta12=betas[2,1], beta22=betas[2,2], beta32=betas[2,3], beta42=betas[2,4], beta52=betas[2,5],
                          beta13=betas[3,1], beta23=betas[3,2], beta33=betas[3,3], beta43=betas[3,4], beta53=betas[3,5],
                          beta14=betas[4,1], beta24=betas[4,2], beta34=betas[4,3], beta44=betas[4,4], beta54=betas[4,5],
                          beta15=betas[5,1], beta25=betas[5,2], beta35=betas[5,3], beta45=betas[5,4], beta55=betas[5,5],
                          nu=c(nu1_sample,nu2_sample,nu3_sample,nu4_sample,nu5_sample), 
                          mu=c(mu1_sample,mu2_sample,mu3_sample,mu4_sample,mu5_sample),
                          d=c(CFR_1_sample,CFR_2_sample,CFR_3_sample,CFR_4_sample,CFR_5_sample), 
                          gamma= gamma_sample,
                          p=c(vacc_efficacy*mean_p,vacc_efficacy*mean_p,vacc_efficacy*mean_p,vacc_efficacy*mean_p,vacc_efficacy*mean_p),
                          scenario=vaccination_scenario)
      }
      else if (vaccination_scenario == '50% relative increase'){
        parms_baseline<-c(beta11=betas[1,1], beta21=betas[1,2], beta31=betas[1,3], beta41=betas[1,4], beta51=betas[1,5],
                          beta12=betas[2,1], beta22=betas[2,2], beta32=betas[2,3], beta42=betas[2,4], beta52=betas[2,5],
                          beta13=betas[3,1], beta23=betas[3,2], beta33=betas[3,3], beta43=betas[3,4], beta53=betas[3,5],
                          beta14=betas[4,1], beta24=betas[4,2], beta34=betas[4,3], beta44=betas[4,4], beta54=betas[4,5],
                          beta15=betas[5,1], beta25=betas[5,2], beta35=betas[5,3], beta45=betas[5,4], beta55=betas[5,5],
                          nu=c(nu1_sample,nu2_sample,nu3_sample,nu4_sample,nu5_sample), 
                          mu=c(mu1_sample,mu2_sample,mu3_sample,mu4_sample,mu5_sample),
                          d=c(CFR_1_sample,CFR_2_sample,CFR_3_sample,CFR_4_sample,CFR_5_sample), 
                          gamma= gamma_sample,
                          p=c(vacc_efficacy*1.5*p1,vacc_efficacy*1.5*p2,vacc_efficacy*1.5*p3,vacc_efficacy*1.5*p4,vacc_efficacy*1.5*p5),
                          scenario=vaccination_scenario)
      }
      else if (vaccination_scenario == 'vaccination highest quintile'){
        parms_baseline<-c(beta11=betas[1,1], beta21=betas[1,2], beta31=betas[1,3], beta41=betas[1,4], beta51=betas[1,5],
                          beta12=betas[2,1], beta22=betas[2,2], beta32=betas[2,3], beta42=betas[2,4], beta52=betas[2,5],
                          beta13=betas[3,1], beta23=betas[3,2], beta33=betas[3,3], beta43=betas[3,4], beta53=betas[3,5],
                          beta14=betas[4,1], beta24=betas[4,2], beta34=betas[4,3], beta44=betas[4,4], beta54=betas[4,5],
                          beta15=betas[5,1], beta25=betas[5,2], beta35=betas[5,3], beta45=betas[5,4], beta55=betas[5,5],
                          nu=c(nu1_sample,nu2_sample,nu3_sample,nu4_sample,nu5_sample), 
                          mu=c(mu1_sample,mu2_sample,mu3_sample,mu4_sample,mu5_sample),
                          d=c(CFR_1_sample,CFR_2_sample,CFR_3_sample,CFR_4_sample,CFR_5_sample), 
                          gamma= gamma_sample,
                          p=c(vacc_efficacy*p5,vacc_efficacy*p5,vacc_efficacy*p5,vacc_efficacy*p5,vacc_efficacy*p5),
                          scenario=vaccination_scenario)
      }
      else if (vaccination_scenario=='full coverage'){
        parms_baseline<-c(beta11=betas[1,1], beta21=betas[1,2], beta31=betas[1,3], beta41=betas[1,4], beta51=betas[1,5],
                          beta12=betas[2,1], beta22=betas[2,2], beta32=betas[2,3], beta42=betas[2,4], beta52=betas[2,5],
                          beta13=betas[3,1], beta23=betas[3,2], beta33=betas[3,3], beta43=betas[3,4], beta53=betas[3,5],
                          beta14=betas[4,1], beta24=betas[4,2], beta34=betas[4,3], beta44=betas[4,4], beta54=betas[4,5],
                          beta15=betas[5,1], beta25=betas[5,2], beta35=betas[5,3], beta45=betas[5,4], beta55=betas[5,5],
                          nu=c(nu1_sample,nu2_sample,nu3_sample,nu4_sample,nu5_sample), 
                          mu=c(mu1_sample,mu2_sample,mu3_sample,mu4_sample,mu5_sample),
                          d=c(CFR_1_sample,CFR_2_sample,CFR_3_sample,CFR_4_sample,CFR_5_sample), 
                          gamma= gamma_sample,
                          p=c(vacc_efficacy*rep(1,5)),
                          scenario=vaccination_scenario)
      }
      else if (vaccination_scenario=='homogeneous-base case'){
        parms_baseline<-c(beta11=betas_homogeneous[1,1], beta21=betas_homogeneous[1,2], beta31=betas_homogeneous[1,3], beta41=betas_homogeneous[1,4], beta51=betas_homogeneous[1,5],
                          beta12=betas_homogeneous[2,1], beta22=betas_homogeneous[2,2], beta32=betas_homogeneous[2,3], beta42=betas_homogeneous[2,4], beta52=betas_homogeneous[2,5],
                          beta13=betas_homogeneous[3,1], beta23=betas_homogeneous[3,2], beta33=betas_homogeneous[3,3], beta43=betas_homogeneous[3,4], beta53=betas_homogeneous[3,5],
                          beta14=betas_homogeneous[4,1], beta24=betas_homogeneous[4,2], beta34=betas_homogeneous[4,3], beta44=betas_homogeneous[4,4], beta54=betas_homogeneous[4,5],
                          beta15=betas_homogeneous[5,1], beta25=betas_homogeneous[5,2], beta35=betas_homogeneous[5,3], beta45=betas_homogeneous[5,4], beta55=betas_homogeneous[5,5],
                          nu=c(nu1_sample,nu2_sample,nu3_sample,nu4_sample,nu5_sample), 
                          mu=c(mu1_sample,mu2_sample,mu3_sample,mu4_sample,mu5_sample),
                          d=c(CFR_1_sample,CFR_2_sample,CFR_3_sample,CFR_4_sample,CFR_5_sample), 
                          gamma= gamma_sample,
                          p=c(p1,p2,p3,p4,p5),
                          scenario='homogeneous-base case')
      }
      else {
        parms_baseline<-c(beta11=betas[1,1], beta21=betas[1,2], beta31=betas[1,3], beta41=betas[1,4], beta51=betas[1,5],
                          beta12=betas[2,1], beta22=betas[2,2], beta32=betas[2,3], beta42=betas[2,4], beta52=betas[2,5],
                          beta13=betas[3,1], beta23=betas[3,2], beta33=betas[3,3], beta43=betas[3,4], beta53=betas[3,5],
                          beta14=betas[4,1], beta24=betas[4,2], beta34=betas[4,3], beta44=betas[4,4], beta54=betas[4,5],
                          beta15=betas[5,1], beta25=betas[5,2], beta35=betas[5,3], beta45=betas[5,4], beta55=betas[5,5],
                          nu=c(nu1_sample,nu2_sample,nu3_sample,nu4_sample,nu5_sample), 
                          mu=c(mu1_sample,mu2_sample,mu3_sample,mu4_sample,mu5_sample),
                          d=c(CFR_1_sample,CFR_2_sample,CFR_3_sample,CFR_4_sample,CFR_5_sample), 
                          gamma= gamma_sample,
                          p=c(vacc_efficacy*p1,vacc_efficacy*p2,vacc_efficacy*p3,vacc_efficacy*p4,vacc_efficacy*1.5*p5),
                          scenario='only top quintile')
      }
      # run ESIR models
      parms_baseline_all[,which(vaccination_scenario_list %in% vaccination_scenario)]=parms_baseline
    }
    all_parameter_output[1:47,(7*(i-1)+1):(7*i)]=parms_baseline_all
    all_parameter_output[48,(7*(i-1)+1):(7*i)]=i
  }
  return(all_parameter_output)
}

vaccination_scenario_list=c('base case', 'mean vaccination coverage', 
                            '50% relative increase', 'vaccination highest quintile', 'full coverage','homogeneous-base case', 'only top quintile')
set.seed(96)
params_vacc_scenarios=ESIR_vaccination_scenarios(vaccination_scenario_list)

ALL<-data.frame(matrix(nrow=17,ncol=ncol(params_vacc_scenarios)))
colnames(ALL)<-seq(1,ncol(params_vacc_scenarios))
rownames(ALL)<-c('prop_dead_quintile1','prop_dead_quintile2',
                 'prop_dead_quintile3','prop_dead_quintile4',
                 'prop_dead_quintile5',"pop_quintile1","pop_quintile2",
                 "pop_quintile3","pop_quintile4","pop_quintile5",
                 "num_dead_quintile1","num_dead_quintile2","num_dead_quintile3",
                 "num_dead_quintile4","num_dead_quintile5",
                 'scenario', 'iteration')

for(i in 1:ncol(params_vacc_scenarios)){
  scenario_curr=params_vacc_scenarios[47,i]
  iteration_curr=as.numeric(params_vacc_scenarios[48,i])
  parameters_curr=as.numeric(params_vacc_scenarios[1:46,i])
  names(parameters_curr)=c("beta11","beta21","beta31","beta41","beta51",
                           "beta12","beta22",'beta32', "beta42", "beta52",
                           "beta13", "beta23", "beta33", "beta43", "beta53",
                           "beta14", "beta24", "beta34", "beta44", "beta54",
                           "beta15", "beta25", "beta35", "beta45", "beta55",
                           "nu1", "nu2", "nu3", "nu4", "nu5",
                           "mu1", "mu2", "mu3", "mu4", "mu5",
                           "d1", "d2", "d3", "d4", "d5",
                           "gamma", "p1", "p2", "p3", "p4", "p5")
  
  if (scenario_curr=='homogeneous-base case'){
    ## Set initial conditions for S, I, R and N for each income quintile
    START.S1 <- (N1*(1-(p1*vacc_efficacy)))-1; START.I1 <- 1; START.R1 <- (N1*p1*vacc_efficacy); START.D1 <- 0
    
    START.S2 <- (N2*(1-(p2*vacc_efficacy)))-1; START.I2 <- 1; START.R2 <- (N2*p2*vacc_efficacy); START.D2 <- 0
    
    START.S3 <- (N3*(1-(p3*vacc_efficacy)))-1; START.I3 <- 1; START.R3 <- (N3*p3*vacc_efficacy); START.D3 <- 0
    
    START.S4 <- (N4*(1-(p4*vacc_efficacy)))-1; START.I4 <- 1; START.R4 <- (N4*p4*vacc_efficacy); START.D4 <- 0
    
    START.S5 <- (N5*(1-(p5*vacc_efficacy)))-1; START.I5 <- 1; START.R5 <- (N5*p5*vacc_efficacy); START.D5 <- 0
    
    START.S_ALL=START.S1+START.S2+START.S3+START.S4+START.S5
    START.I_ALL=START.I1+START.I2+START.I3+START.I4+START.I5
    START.R_ALL=START.R1+START.R2+START.R3+START.R4+START.R5
    START.D_ALL=START.D1+START.D2+START.D3+START.D4+START.D5
    
    ESIR.output <- lsoda(y=c(S1=START.S_ALL, I1=START.I_ALL, R1=START.R_ALL, D1=START.D_ALL,
                             S2=START.S_ALL, I2=START.I_ALL, R2=START.R_ALL, D2=START.D_ALL,
                             S3=START.S_ALL, I3=START.I_ALL, R3=START.R_ALL, D3=START.D_ALL,
                             S4=START.S_ALL, I4=START.I_ALL, R4=START.R_ALL, D4=START.D_ALL,
                             S5=START.S_ALL, I5=START.I_ALL, R5=START.R_ALL, D5=START.D_ALL),
                         times=times_vector, 
                         func=dx.dt.ESIR, 
                         parms=parameters_curr)
  }
  else if(scenario_curr=='base case') {
    START.S1 <- (N1*(1-(p1*vacc_efficacy)))-1; START.I1 <- 1; START.R1 <- (N1*p1*vacc_efficacy); START.D1 <- 0
    
    START.S2 <- (N2*(1-(p2*vacc_efficacy)))-1; START.I2 <- 1; START.R2 <- (N2*p2*vacc_efficacy); START.D2 <- 0
    
    START.S3 <- (N3*(1-(p3*vacc_efficacy)))-1; START.I3 <- 1; START.R3 <- (N3*p3*vacc_efficacy); START.D3 <- 0
    
    START.S4 <- (N4*(1-(p4*vacc_efficacy)))-1; START.I4 <- 1; START.R4 <- (N4*p4*vacc_efficacy); START.D4 <- 0
    
    START.S5 <- (N5*(1-(p5*vacc_efficacy)))-1; START.I5 <- 1; START.R5 <- (N5*p5*vacc_efficacy); START.D5 <- 0
    
    ESIR.output <- lsoda(y=c(S1=START.S1, I1=START.I1, R1=START.R1, D1=START.D1,
                             S2=START.S2, I2=START.I2, R2=START.R2, D2=START.D2,
                             S3=START.S3, I3=START.I3, R3=START.R3, D3=START.D3,
                             S4=START.S4, I4=START.I4, R4=START.R4, D4=START.D4,
                             S5=START.S5, I5=START.I5, R5=START.R5, D5=START.D5), 
                         times=times_vector, 
                         func=dx.dt.ESIR, 
                         parms=parameters_curr)
  }
  else if (scenario_curr=='mean vaccination coverage') {
    START.S1 <- (N1*(1-(mean_p*vacc_efficacy)))-1; START.I1 <- 1; 
    START.R1 <- (N1*mean_p*vacc_efficacy); START.D1 <- 0
    
    START.S2 <- (N2*(1-(mean_p*vacc_efficacy)))-1; START.I2 <- 1
    START.R2 <- (N2*mean_p*vacc_efficacy); START.D2 <- 0
    
    START.S3 <- (N3*(1-(mean_p*vacc_efficacy)))-1; START.I3 <- 1
    START.R3 <- (N3*mean_p*vacc_efficacy); START.D3 <- 0
    
    START.S4 <- (N4*(1-(mean_p*vacc_efficacy)))-1; START.I4 <- 1
    START.R4 <- (N4*mean_p*vacc_efficacy); START.D4 <- 0
    
    START.S5 <- (N5*(1-(mean_p*vacc_efficacy)))-1; START.I5 <- 1
    START.R5 <- (N5*mean_p*vacc_efficacy); START.D5 <- 0
    
    ESIR.output <- lsoda(y=c(S1=START.S1, I1=START.I1, R1=START.R1, D1=START.D1,
                             S2=START.S2, I2=START.I2, R2=START.R2, D2=START.D2,
                             S3=START.S3, I3=START.I3, R3=START.R3, D3=START.D3,
                             S4=START.S4, I4=START.I4, R4=START.R4, D4=START.D4,
                             S5=START.S5, I5=START.I5, R5=START.R5, D5=START.D5), 
                         times=times_vector, 
                         func=dx.dt.ESIR, 
                         parms=parameters_curr)
  }
  
  else if (scenario_curr=='50% relative increase') {
    START.S1 <- (N1*(1-(p1*1.5*vacc_efficacy)))-1; START.I1 <- 1; 
    START.R1 <- (N1*p1*1.5*vacc_efficacy); START.D1 <- 0
    
    START.S2 <- (N2*(1-(p2*1.5*vacc_efficacy)))-1; START.I2 <- 1
    START.R2 <- (N2*p2*1.5*vacc_efficacy); START.D2 <- 0
    
    START.S3 <- (N3*(1-(p3*1.5*vacc_efficacy)))-1; START.I3 <- 1
    START.R3 <- (N3*p3*1.5*vacc_efficacy); START.D3 <- 0
    
    START.S4 <- (N4*(1-(p4*1.5*vacc_efficacy)))-1; START.I4 <- 1
    START.R4 <- (N4*p4*1.5*vacc_efficacy); START.D4 <- 0
    
    START.S5 <- (N5*(1-(p5*1.5*vacc_efficacy)))-1; START.I5 <- 1
    START.R5 <- (N5*p5*1.5*vacc_efficacy); START.D5 <- 0
    
    ESIR.output <- lsoda(y=c(S1=START.S1, I1=START.I1, R1=START.R1, D1=START.D1,
                             S2=START.S2, I2=START.I2, R2=START.R2, D2=START.D2,
                             S3=START.S3, I3=START.I3, R3=START.R3, D3=START.D3,
                             S4=START.S4, I4=START.I4, R4=START.R4, D4=START.D4,
                             S5=START.S5, I5=START.I5, R5=START.R5, D5=START.D5), 
                         times=times_vector, 
                         func=dx.dt.ESIR, 
                         parms=parameters_curr)
  }
  
  else if (scenario_curr=='vaccination highest quintile') {
    START.S1 <- (N1*(1-(p5*vacc_efficacy)))-1; START.I1 <- 1; 
    START.R1 <- (N1*p5*vacc_efficacy); START.D1 <- 0
    
    START.S2 <- (N2*(1-(p5*vacc_efficacy)))-1; START.I2 <- 1
    START.R2 <- (N2*p5*vacc_efficacy); START.D2 <- 0
    
    START.S3 <- (N3*(1-(p5*vacc_efficacy)))-1; START.I3 <- 1
    START.R3 <- (N3*p5*vacc_efficacy); START.D3 <- 0
    
    START.S4 <- (N4*(1-(p5*vacc_efficacy)))-1; START.I4 <- 1
    START.R4 <- (N4*p5*vacc_efficacy); START.D4 <- 0
    
    START.S5 <- (N5*(1-(p5*vacc_efficacy)))-1; START.I5 <- 1
    START.R5 <- (N5*p5*vacc_efficacy); START.D5 <- 0
    
    ESIR.output <- lsoda(y=c(S1=START.S1, I1=START.I1, R1=START.R1, D1=START.D1,
                             S2=START.S2, I2=START.I2, R2=START.R2, D2=START.D2,
                             S3=START.S3, I3=START.I3, R3=START.R3, D3=START.D3,
                             S4=START.S4, I4=START.I4, R4=START.R4, D4=START.D4,
                             S5=START.S5, I5=START.I5, R5=START.R5, D5=START.D5), 
                         times=times_vector, 
                         func=dx.dt.ESIR, 
                         parms=parameters_curr)
  }
  
  else if (scenario_curr== 'full coverage'){
    
    START.S1 <- (N1*(1-(1*vacc_efficacy)))-1; START.I1 <- 1; 
    START.R1 <- (N1*1*vacc_efficacy); START.D1 <- 0
    
    START.S2 <- (N2*(1-(1*vacc_efficacy)))-1; START.I2 <- 1
    START.R2 <- (N2*1*vacc_efficacy); START.D2 <- 0
    
    START.S3 <- (N3*(1-(1*vacc_efficacy)))-1; START.I3 <- 1
    START.R3 <- (N3*1*vacc_efficacy); START.D3 <- 0
    
    START.S4 <- (N4*(1-(1*vacc_efficacy)))-1; START.I4 <- 1
    START.R4 <- (N4*1*vacc_efficacy); START.D4 <- 0
    
    START.S5 <- (N5*(1-(1*vacc_efficacy)))-1; START.I5 <- 1
    START.R5 <- (N5*1*vacc_efficacy); START.D5 <- 0
    
    ESIR.output <- lsoda(y=c(S1=START.S1, I1=START.I1, R1=START.R1, D1=START.D1,
                             S2=START.S2, I2=START.I2, R2=START.R2, D2=START.D2,
                             S3=START.S3, I3=START.I3, R3=START.R3, D3=START.D3,
                             S4=START.S4, I4=START.I4, R4=START.R4, D4=START.D4,
                             S5=START.S5, I5=START.I5, R5=START.R5, D5=START.D5), 
                         times=times_vector, 
                         func=dx.dt.ESIR, 
                         parms=parameters_curr)
  }
  
  else{
    START.S1 <- (N1*(1-(p1*vacc_efficacy)))-1; START.I1 <- 1; 
    START.R1 <- (N1*p1*vacc_efficacy); START.D1 <- 0
    
    START.S2 <- (N2*(1-(p2*vacc_efficacy)))-1; START.I2 <- 1
    START.R2 <- (N2*p2*vacc_efficacy); START.D2 <- 0
    
    START.S3 <- (N3*(1-(p3*vacc_efficacy)))-1; START.I3 <- 1
    START.R3 <- (N3*p3*vacc_efficacy); START.D3 <- 0
    
    START.S4 <- (N4*(1-(p4*vacc_efficacy)))-1; START.I4 <- 1
    START.R4 <- (N4*p4*vacc_efficacy); START.D4 <- 0
    
    START.S5 <- (N5*(1-(p5*1.5*vacc_efficacy)))-1; START.I5 <- 1
    START.R5 <- (N5*p5*1.5*vacc_efficacy); START.D5 <- 0
    
    ESIR.output <- lsoda(y=c(S1=START.S1, I1=START.I1, R1=START.R1, D1=START.D1,
                             S2=START.S2, I2=START.I2, R2=START.R2, D2=START.D2,
                             S3=START.S3, I3=START.I3, R3=START.R3, D3=START.D3,
                             S4=START.S4, I4=START.I4, R4=START.R4, D4=START.D4,
                             S5=START.S5, I5=START.I5, R5=START.R5, D5=START.D5), 
                         times=times_vector, 
                         func=dx.dt.ESIR, 
                         parms=parameters_curr)
  }
  
  ESIR.output_quintile1=data.frame(ESIR.output[,c("S1","I1","R1")])
  ESIR.output_quintile1$pop_current=rowSums(ESIR.output_quintile1)
  ESIR.output_quintile2=data.frame(ESIR.output[,c("S2","I2","R2")])
  ESIR.output_quintile2$pop_current=rowSums(ESIR.output_quintile2)
  ESIR.output_quintile3=data.frame(ESIR.output[,c("S3","I3","R3")])
  ESIR.output_quintile3$pop_current=rowSums(ESIR.output_quintile3)
  ESIR.output_quintile4=data.frame(ESIR.output[,c("S4","I4","R4")])
  ESIR.output_quintile4$pop_current=rowSums(ESIR.output_quintile4)
  ESIR.output_quintile5=data.frame(ESIR.output[,c("S5","I5","R5")])
  ESIR.output_quintile5$pop_current=rowSums(ESIR.output_quintile5)
  
  ESIR.output_prop_dead_quintile1=data.frame(ESIR.output)%>%
    summarise(prop_D1=D1[400]/ESIR.output_quintile1$pop_current[400])
  
  ESIR.output_prop_dead_quintile2=data.frame(ESIR.output)%>%
    summarise(prop_D2=D2[400]/ESIR.output_quintile2$pop_current[400])
  
  ESIR.output_prop_dead_quintile3=data.frame(ESIR.output)%>%
    summarise(prop_D3=D3[400]/ESIR.output_quintile3$pop_current[400])
  
  ESIR.output_prop_dead_quintile4=data.frame(ESIR.output)%>%
    summarise(prop_D4=D4[400]/ESIR.output_quintile4$pop_current[400])
  
  ESIR.output_prop_dead_quintile5=data.frame(ESIR.output)%>%
    summarise(prop_D5=D5[400]/ESIR.output_quintile5$pop_current[400])
  
  ALL[1:15,i]<-as.numeric(c(ESIR.output_prop_dead_quintile1$prop_D1,ESIR.output_prop_dead_quintile2$prop_D2,
                            ESIR.output_prop_dead_quintile3$prop_D3,ESIR.output_prop_dead_quintile4$prop_D4,
                            ESIR.output_prop_dead_quintile5$prop_D5,ESIR.output_quintile1$pop_current[400],
                            ESIR.output_quintile2$pop_current[400],ESIR.output_quintile3$pop_current[400],
                            ESIR.output_quintile4$pop_current[400],ESIR.output_quintile5$pop_current[400],
                            data.frame(ESIR.output)$D1[400],data.frame(ESIR.output)$D2[400],
                            data.frame(ESIR.output)$D3[400],data.frame(ESIR.output)$D4[400],
                            data.frame(ESIR.output)$D5[400]))
  ALL[16,i]<-scenario_curr
  ALL[17,i]<-iteration_curr 
}

# extract and summarize simulation output by group (means and 95% CI)
ALL_final=data.frame(t(ALL))
ALL_final<-ALL_final%>%
  mutate(prop_dead_quintile1=as.numeric(as.character(prop_dead_quintile1)))%>%
  mutate(prop_dead_quintile2=as.numeric(as.character(prop_dead_quintile2)))%>%
  mutate(prop_dead_quintile3=as.numeric(as.character(prop_dead_quintile3)))%>%
  mutate(prop_dead_quintile4=as.numeric(as.character(prop_dead_quintile4)))%>%
  mutate(prop_dead_quintile5=as.numeric(as.character(prop_dead_quintile5)))%>%
  mutate(pop_quintile1=as.numeric(as.character(pop_quintile1)))%>%
  mutate(pop_quintile2=as.numeric(as.character(pop_quintile2)))%>%
  mutate(pop_quintile3=as.numeric(as.character(pop_quintile3)))%>%
  mutate(pop_quintile4=as.numeric(as.character(pop_quintile4)))%>%
  mutate(pop_quintile5=as.numeric(as.character(pop_quintile5)))%>%
  mutate(num_dead_quintile1=as.numeric(as.character(num_dead_quintile1)))%>%
  mutate(num_dead_quintile2=as.numeric(as.character(num_dead_quintile2)))%>%
  mutate(num_dead_quintile3=as.numeric(as.character(num_dead_quintile3)))%>%
  mutate(num_dead_quintile4=as.numeric(as.character(num_dead_quintile4)))%>%
  mutate(num_dead_quintile5=as.numeric(as.character(num_dead_quintile5)))

## designate base case scenario as that with DHS-reported quintile-specific coverage rates
## evaluate differences in the #D under each vaccination scenario compared to base case
ALL_summz=ALL_final%>%
  group_by(iteration)%>%
  mutate(order=seq(1,n()))%>%
  mutate(base_case_quintile1=ifelse(order==1,prop_dead_quintile1,0))%>%
  mutate(base_case_quintile2=ifelse(order==1,prop_dead_quintile2,0))%>%
  mutate(base_case_quintile3=ifelse(order==1,prop_dead_quintile3,0))%>%
  mutate(base_case_quintile4=ifelse(order==1,prop_dead_quintile4,0))%>%
  mutate(base_case_quintile5=ifelse(order==1,prop_dead_quintile5,0))%>%
  mutate(diff_quintile1=ifelse(order==1,base_case_quintile1,(prop_dead_quintile1-base_case_quintile1[1])*pop_quintile1))%>%
  mutate(diff_quintile2=ifelse(order==1,base_case_quintile2,(prop_dead_quintile2-base_case_quintile2[1])*pop_quintile2))%>%
  mutate(diff_quintile3=ifelse(order==1,base_case_quintile3,(prop_dead_quintile3-base_case_quintile3[1])*pop_quintile3))%>%
  mutate(diff_quintile4=ifelse(order==1,base_case_quintile4,(prop_dead_quintile4-base_case_quintile4[1])*pop_quintile4))%>%
  mutate(diff_quintile5=ifelse(order==1,base_case_quintile5,(prop_dead_quintile5-base_case_quintile5[1])*pop_quintile5))

ALL_summz_mean=ALL_summz%>%
  group_by(scenario)%>%
  summarise(mean(diff_quintile1),mean(diff_quintile2),
            mean(diff_quintile3),mean(diff_quintile4),
            mean(diff_quintile5))

ALL_summz_LQ=ALL_summz%>%
  group_by(scenario)%>%
  summarise(quantile(diff_quintile1,c(0.025,0.975))[[1]],quantile(diff_quintile2,c(0.025,0.975))[[1]],
            quantile(diff_quintile3,c(0.025,0.975))[[1]],quantile(diff_quintile4,c(0.025,0.975))[[1]],
            quantile(diff_quintile5,c(0.025,0.975))[[1]])

ALL_summz_UQ=ALL_summz%>%
  group_by(scenario)%>%
  summarise(quantile(diff_quintile1,c(0.025,0.975))[[2]],quantile(diff_quintile2,c(0.025,0.975))[[2]],
            quantile(diff_quintile3,c(0.025,0.975))[[2]],quantile(diff_quintile4,c(0.025,0.975))[[2]],
            quantile(diff_quintile5,c(0.025,0.975))[[2]])

colnames(ALL_summz_mean)<-c("scenario","1","2","3","4","5","6")

colnames(ALL_summz_LQ)<-c("scenario","1","2","3","4","5","6")

colnames(ALL_summz_UQ)<-c("scenario","1","2","3","4","5","6")

ALL_summz_mean_final=gather(ALL_summz_mean,quintile,value,"1":"5")
ALL_summz_LQ_final=gather(ALL_summz_LQ,quintile,value,"1":"5")
ALL_summz_UQ_final=gather(ALL_summz_UQ,quintile,value,"1":"5")

ALL_scenarios_mean_LQ=merge(ALL_summz_mean_final,ALL_summz_LQ_final,by=c("quintile","scenario"))
ALL_scenarios_mean_LQ_UQ=merge(ALL_scenarios_mean_LQ,ALL_summz_UQ_final,by=c("quintile","scenario"))
colnames(ALL_scenarios_mean_LQ_UQ)[3:5]<-c("mean","LQ","UQ")
ALL_scenarios_mean_LQ_UQ=ALL_scenarios_mean_LQ_UQ[which(ALL_scenarios_mean_LQ_UQ$scenario!='base case'),]

write.csv(ALL_scenarios_mean_LQ_UQ,"./out/ALL_scenarios_mean_LQ_UQ_Feb12.csv",row.names=F)

