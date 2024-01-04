#Code to for sensitivity analysis for all models from Beck-Johnson et al.
#© 2023 Colorado State University 

library(adaptivetau) 

exposed5=1
exposed10=c(1,3)
exposed15=c(1,2,4)
exposed20=c(1,2,5)
exposed50=c(1,3,5,13)
exposed100=c(1,5,10,25)
exposed200=c(2,10,20,50)
exposed500=c(5,25,50,100)
exposed1000=c(10,50,100)
exposed2500=c(25,100)
exposed5000=c(50,100)
exposed10000=100

exposed<-list(exposed5, exposed10, exposed15, exposed20, exposed50, exposed100, exposed200, exposed500, exposed1000, exposed2500, exposed5000, exposed10000)

herdsize<-c(5,10,15,20,50,100,200,500,1000,2500,5000,10000)




#Sensitivity analysis

require(lhs)

np<-1000

lhs.sample<-maximinLHS(np,4)

#parameter ranges, max and min

#All model parameter range for beta

beta.min<-0.34
beta.max<-141.62

#Clinical First model parameter ranges

phi.min<-(1/2.94)
phi.max<-(1/5.9)

omega.min<-(1/0.16)
omega.max<-(1/1.3)

theta.min<-(1/0.3)
theta.max<-(1/4.8)


parms.set1<-cbind(beta<-lhs.sample[,1]*(beta.max-beta.min)+beta.min, lhs.sample[,2]*(phi.max-phi.min)+phi.min, lhs.sample[,3]*(omega.max-omega.min)+omega.min, lhs.sample[,4]*(theta.max-theta.min)+theta.min)


#Infectious First model parameter ranges, viremia in OPF parameters

nu.min<-(1/0.18)
nu.max<-(1/0.84)

rho.min<-(1/2.8)
rho.max<-(1/5.28)

mu.min<-(1/3.72)
mu.max<-(1/5.43)


parms.set2<-cbind(beta<-lhs.sample[,1]*(beta.max-beta.min)+beta.min, lhs.sample[,2]*(nu.max-nu.min)+nu.min, lhs.sample[,3]*(rho.max-rho.min)+rho.min, lhs.sample[,4]*(mu.max-mu.min)+mu.min)


#Infectious First model parameter ranges, viremia in blood parameters

nu.b.min<-(1/1.49)
nu.b.max<-(1/3.86)

rho.b.min<-(1/1.44)
rho.b.max<-(1/2.30)

mu.b.min<-(1/2.09)
mu.b.max<-(1/2.67)


parms.set2b<-cbind(beta<-lhs.sample[,1]*(beta.max-beta.min)+beta.min, lhs.sample[,2]*(nu.b.max-nu.b.min)+nu.b.min, lhs.sample[,3]*(rho.b.max-rho.b.min)+rho.b.min, lhs.sample[,4]*(mu.b.max-mu.b.min)+mu.b.min)


#Infectious First model parameter ranges, viremia in nasal fluid parameters

nu.n.min<-(1/1.69)
nu.n.max<-(1/4.45)

rho.n.min<-(1/1.34)
rho.n.max<-(1/1.75)

mu.n.min<-(1/3.30)
mu.n.max<-(1/8.417)


parms.set2n<-cbind(beta<-lhs.sample[,1]*(beta.max-beta.min)+beta.min, lhs.sample[,2]*(nu.n.max-nu.n.min)+nu.n.min, lhs.sample[,3]*(rho.n.max-rho.n.min)+rho.n.min, lhs.sample[,4]*(mu.n.max-mu.n.min)+mu.n.min)


#Base model parameter ranges
lhs.samplebase<-maximinLHS(np,3)
sigma.min<-(1/3.1)
sigma.max<-(1/7.2)

gamma.min<-(1/0.3)
gamma.max<-(1/4.8)


parms.set<-cbind(beta<-lhs.samplebase[,1]*(beta.max-beta.min)+beta.min, lhs.samplebase[,2]*(sigma.max-sigma.min)+sigma.min, lhs.samplebase[,3]*(gamma.max-gamma.min)+gamma.min)



#Density Dependent, clinical first model

transitions.v1 = list(
  c(S = -1, E = +1), # transmission
  c(E = -1, C = +1), #becoming clinical
  c(C = -1, I = +1), # becoming infectious
  c(I = -1, R = +1) # recovery
)

Var1DD <- function(x, p, t) {
  return(c(p[1]*x['S']*x['I'], #transmission rate
           p[2]*x['E'],#becoming clinical rate
           p[3]*x['C'], #becoming infectious rate
           p[4]*x['I']  )) #recovery rate
}


for(h in 1:length(herdsize)){
for(j in 1:length(exposed[[h]])){
  
  init.values.v1 = c(
    S = herdsize[h]-exposed[[h]][j], # susceptible animals--base herd size
    E = exposed[[h]][j],   # exposed animals (1% of herd)
    C = 0, #clinical not infectious
    I = 0,   # infectious
    R = 0)    # recovered
    

for(k in 1:dim(parms.set1)[1]){

numsims=1000

parms1 <- parms.set1[k,]

  out.v1.DD.c=matrix(nrow=101,ncol=1000)
  out.v1.DD.c<-as.data.frame(out.v1.DD.c)
  
  out.v1.DD.i=matrix(nrow=101,ncol=1000)
  out.v1.DD.i<-as.data.frame(out.v1.DD.i)
  
  set.seed(478394)
  
  for(i in 1:1000){
    sim=ssa.adaptivetau(init.values.v1, transitions.v1, Var1DD, parms1, tf=100)
    sim=as.data.frame(sim)
    idx <- diff(ceiling(sim$time)) == 1
    tempc<-sim[idx,4]
    tempi<-sim[idx,5]
    repc=c(tempc, rep(NA,length=(101-length(tempc))))
    repi=c(tempi, rep(NA,length=(101-length(tempi))))
    out.v1.DD.c[,i]=repc
    out.v1.DD.i[,i]=repi
  }
  
  out.v1.DD.c[is.na(out.v1.DD.c)] <- 0
  out.v1.DD.i[is.na(out.v1.DD.i)] <- 0
 
  seed<-rep(exposed[[h]][j],length=101)
  herd<-rep(herdsize[h],length=101)

	out.v1.DD.c=cbind(herd,seed,out.v1.DD.c)
	out.v1.DD.i=cbind(herd,seed,out.v1.DD.i)

  
  filenameC=paste("Var1DDC",herdsize[h],"S",exposed[[h]][j],"sens",k,"csv",sep=".")
  filenameI=paste("Var1DDI",herdsize[h],"S",exposed[[h]][j],"sens",k,"csv",sep=".")
  
  write.csv(out.v1.DD.c,file=filenameC)
  write.csv(out.v1.DD.i,file=filenameI)
  
}
}

}


#Frequency dependent, clinicial first model

Var1FD <- function(x, p, t) {
  return(c(p[1]*x['S']*x['I']/(x['S']+x['E']+x['C']+x['I']+x['R']), #transmission rate
           p[2]*x['E'],#becoming clinical rate
           p[3]*x['C'], #becoming infectious rate
           p[4]*x['I']  )) #recovery rate
}


for(h in 1:length(herdsize)){
for(j in 1:length(exposed[[h]])){
  
  init.values.v1 = c(
    S = herdsize[h]-exposed[[h]][j], # susceptible animals--base herd size
    E = exposed[[h]][j],   # exposed animals (1% of herd)
    C = 0, #clinical not infectious
    I = 0,   # infectious
    R = 0)    # recovered
  
for(k in 1:dim(parms.set1)[1]){

numsims=1000

parms1 <- parms.set1[k,]

   numsims=1000
  
  out.v1.FD.c=matrix(nrow=101,ncol=1000)
  out.v1.FD.c<-as.data.frame(out.v1.FD.c)
  
  out.v1.FD.i=matrix(nrow=101,ncol=1000)
  out.v1.FD.i<-as.data.frame(out.v1.FD.i)
  
  set.seed(478394)
  
  for(i in 1:1000){
    sim=ssa.adaptivetau(init.values.v1, transitions.v1, Var1FD, parms1, tf=100)
    sim=as.data.frame(sim)
    idx <- diff(ceiling(sim$time)) == 1
    tempc<-sim[idx,4]
    tempi<-sim[idx,5]
    repc=c(tempc, rep(NA,length=(101-length(tempc))))
    repi=c(tempi, rep(NA,length=(101-length(tempi))))
    out.v1.FD.c[,i]=repc
    out.v1.FD.i[,i]=repi
  }
  
  out.v1.FD.c[is.na(out.v1.FD.c)] <- 0
  out.v1.FD.i[is.na(out.v1.FD.i)] <- 0
 
  seed<-rep(exposed[[h]][j],length=101)
  herd<-rep(herdsize[h],length=101)

	out.v1.FD.c=cbind(herd,seed,out.v1.FD.c)
	out.v1.FD.i=cbind(herd,seed,out.v1.FD.i)

  
  filenameC=paste("Var1FDC",herdsize[h],"S",exposed[[h]][j],"sens",k,"csv",sep=".")
  filenameI=paste("Var1FDI",herdsize[h],"S",exposed[[h]][j],"sens",k,"csv",sep=".")
  
  write.csv(out.v1.FD.c,file=filenameC)
  write.csv(out.v1.FD.i,file=filenameI)
  
}
}
}








#Density Dependent, infectious first, OPF parameters


  transitions.v2 = list(
    c(S = -1, E = +1), # transmission
    c(E = -1, I = +1), #becoming infectious
    c(I = -1, C = +1), # becoming clincial
    c(C = -1, R = +1) # recovery
  )
  
  Var2DD <- function(x, p, t) {
    return(c(p[1]*x['S']*(x['I']+x['C']), #transmission rate
             p[2]*x['E'],#becoming infectious rate
             p[3]*x['I'], #becoming clinical rate
             p[4]*x['C']  )) #recovery rate
  }
  

for(h in 1:length(herdsize)){
for(j in 1:length(exposed[[h]])){
  
   init.values.v2 = c(
    S = herdsize[h]-exposed[[h]][j], # susceptible animals--base herd size
    E = exposed[[h]][j],   # exposed animals (1% of herd)
    I = 0, #infectious not clinical
    C = 0,   #clinical and  infectious
    R = 0)    # recovered
    
for(k in 1:dim(parms.set2)[1]){

numsims=1000

parms2 <- parms.set2[k,]


  out.v2.DD.c=matrix(nrow=101,ncol=1000)
  out.v2.DD.c<-as.data.frame(out.v2.DD.c)
  
  out.v2.DD.i=matrix(nrow=101,ncol=1000)
  out.v2.DD.i<-as.data.frame(out.v2.DD.i)
  
  set.seed(478394)
  
  for(i in 1:1000){
    sim=ssa.adaptivetau(init.values.v2, transitions.v2, Var2DD, parms2, tf=100)
    sim=as.data.frame(sim)
    idx <- diff(ceiling(sim$time)) == 1
    tempc<-sim[idx,4]
    tempi<-sim[idx,5]
    repc=c(tempc, rep(NA,length=(101-length(tempc))))
    repi=c(tempi, rep(NA,length=(101-length(tempi))))
    out.v2.DD.c[,i]=repc
    out.v2.DD.i[,i]=repi
  }
  
  out.v2.DD.c[is.na(out.v2.DD.c)] <- 0
  out.v2.DD.i[is.na(out.v2.DD.i)] <- 0
 
  seed<-rep(exposed[[h]][j],length=101)
  herd<-rep(herdsize[h],length=101)

	out.v2.DD.c=cbind(herd,seed,out.v2.DD.c)
	out.v2.DD.i=cbind(herd,seed,out.v2.DD.i)

  
  filenameC=paste("Var2oDDC",herdsize[h],"S",exposed[[h]][j],"sens",k,"csv",sep=".")
  filenameI=paste("Var2oDDI",herdsize[h],"S",exposed[[h]][j],"sens",k,"csv",sep=".")
  
  write.csv(out.v2.DD.c,file=filenameC)
  write.csv(out.v2.DD.i,file=filenameI)
  
}
}
}



  #Frequency Dependent, infectious first, OPF parameters
  
Var2FD <- function(x, p, t) {
  return(c(p[1]*x['S']*(x['I']+x['C'])/(x['S']+x['E']+x['I']+x['C']+x['R']), #transmission rate
           p[2]*x['E'],#becoming infectious rate
           p[3]*x['I'], #becoming clinical rate
           p[4]*x['C']  )) #recovery rate
}



for(h in 1:length(herdsize)){
for(j in 1:length(exposed[[h]])){
  
    init.values.v2 = c(
    S = herdsize[h]-exposed[[h]][j], # susceptible animals--base herd size
    E = exposed[[h]][j],   # exposed animals (1% of herd)
    I = 0, #infectious not clinical
    C = 0,   #clinical and  infectious
    R = 0)    # recovered
  
 for(k in 1:dim(parms.set2)[1]){

numsims=1000

parms2 <- parms.set2[k,]

 
  
  out.v2.FD.c=matrix(nrow=101,ncol=1000)
  out.v2.FD.c<-as.data.frame(out.v2.FD.c)
  
  out.v2.FD.i=matrix(nrow=101,ncol=1000)
  out.v2.FD.i<-as.data.frame(out.v2.FD.i)
  
  set.seed(478394)
  
  for(i in 1:1000){
    sim=ssa.adaptivetau(init.values.v2, transitions.v2, Var2FD, parms2, tf=100)
    sim=as.data.frame(sim)
    idx <- diff(ceiling(sim$time)) == 1
    tempc<-sim[idx,4]
    tempi<-sim[idx,5]
    repc=c(tempc, rep(NA,length=(101-length(tempc))))
    repi=c(tempi, rep(NA,length=(101-length(tempi))))
    out.v2.FD.c[,i]=repc
    out.v2.FD.i[,i]=repi
  }
  
  out.v2.FD.c[is.na(out.v2.FD.c)] <- 0
  out.v2.FD.i[is.na(out.v2.FD.i)] <- 0
 
  seed<-rep(exposed[[h]][j],length=101)
  herd<-rep(herdsize[h],length=101)

	out.v2.FD.c=cbind(herd,seed,out.v2.FD.c)
	out.v2.FD.i=cbind(herd,seed,out.v2.FD.i)

  
  filenameC=paste("Var2oFDC",herdsize[h],"S",exposed[[h]][j],"csv",sep=".")
  filenameI=paste("Var2oFDI",herdsize[h],"S",exposed[[h]][j],"csv",sep=".")
  
  write.csv(out.v2.FD.c,file=filenameC)
  write.csv(out.v2.FD.i,file=filenameI)
  
}
}
}


#Density Dependent, infectious first, blood parameters


init.values.v2 = c(
    S = 4, # susceptible animals--base herd size
    E = 1,   # exposed animals (1% of herd)
    I = 0, #infectious not clinical
    C = 0,   #clinical and  infectious
    R = 0)    # recovered
  
  
  transitions.v2 = list(
    c(S = -1, E = +1), # transmission
    c(E = -1, I = +1), #becoming infectious
    c(I = -1, C = +1), # becoming clinical
    c(C = -1, R = +1) # recovery
  )
  
  Var2DD <- function(x, p, t) {
    return(c(p[1]*x['S']*(x['I']+x['C']), #transmission rate
             p[2]*x['E'],#becoming infectious rate
             p[3]*x['I'], #becoming clinical rate
             p[4]*x['C']  )) #recovery rate
  }
  

for(h in 1:length(herdsize)){
for(j in 1:length(exposed[[h]])){
  
   init.values.v2 = c(
    S = herdsize[h]-exposed[[h]][j], # susceptible animals--base herd size
    E = exposed[[h]][j],   # exposed animals (1% of herd)
    I = 0, #infectious not clinical
    C = 0,   #clinical and  infectious
    R = 0)    # recovered
 

for(k in 1:dim(parms.set2b)[1]){

numsims=1000

parms2 <- parms.set2b[k,]

  
  out.v2.DD.c=matrix(nrow=101,ncol=1000)
  out.v2.DD.c<-as.data.frame(out.v2.DD.c)
  
  out.v2.DD.i=matrix(nrow=101,ncol=1000)
  out.v2.DD.i<-as.data.frame(out.v2.DD.i)
  
  set.seed(478394)
  
  for(i in 1:1000){
    sim=ssa.adaptivetau(init.values.v2, transitions.v2, Var2DD, parms2, tf=100)
    sim=as.data.frame(sim)
    idx <- diff(ceiling(sim$time)) == 1
    tempc<-sim[idx,4]
    tempi<-sim[idx,5]
    repc=c(tempc, rep(NA,length=(101-length(tempc))))
    repi=c(tempi, rep(NA,length=(101-length(tempi))))
    out.v2.DD.c[,i]=repc
    out.v2.DD.i[,i]=repi
  }
  
  out.v2.DD.c[is.na(out.v2.DD.c)] <- 0
  out.v2.DD.i[is.na(out.v2.DD.i)] <- 0
 
  seed<-rep(exposed[[h]][j],length=101)
  herd<-rep(herdsize[h],length=101)

	out.v2.DD.c=cbind(herd,seed,out.v2.DD.c)
	out.v2.DD.i=cbind(herd,seed,out.v2.DD.i)

  
  filenameC=paste("Var2bDDC",herdsize[h],"S",exposed[[h]][j],"sens",k,"csv",sep=".")
  filenameI=paste("Var2bDDI",herdsize[h],"S",exposed[[h]][j],"sens",k,"csv",sep=".")
  
  write.csv(out.v2.DD.c,file=filenameC)
  write.csv(out.v2.DD.i,file=filenameI)
  
}
}
}

  #Frequency Dependent, infectious first, blood parameters
  
Var2FD <- function(x, p, t) {
  return(c(p[1]*x['S']*(x['I']+x['C'])/(x['S']+x['E']+x['I']+x['C']+x['R']), #transmission rate
           p[2]*x['E'],#becoming infectious rate
           p[3]*x['I'], #becoming clinical rate
           p[4]*x['C']  )) #recovery rate
}



for(h in 1:length(herdsize)){
for(j in 1:length(exposed[[h]])){
  
    init.values.v2 = c(
    S = herdsize[h]-exposed[[h]][j], # susceptible animals--base herd size
    E = exposed[[h]][j],   # exposed animals (1% of herd)
    I = 0, #infectious not clinical
    C = 0,   #clinical and  infectious
    R = 0)    # recovered
  
for(k in 1:dim(parms.set2b)[1]){

numsims=1000

parms2 <- parms.set2b[k,]
  
  out.v2.FD.c=matrix(nrow=101,ncol=1000)
  out.v2.FD.c<-as.data.frame(out.v2.FD.c)
  
  out.v2.FD.i=matrix(nrow=101,ncol=1000)
  out.v2.FD.i<-as.data.frame(out.v2.FD.i)
  
  set.seed(478394)
  
  for(i in 1:1000){
    sim=ssa.adaptivetau(init.values.v2, transitions.v2, Var2FD, parms2, tf=100)
    sim=as.data.frame(sim)
    idx <- diff(ceiling(sim$time)) == 1
    tempc<-sim[idx,4]
    tempi<-sim[idx,5]
    repc=c(tempc, rep(NA,length=(101-length(tempc))))
    repi=c(tempi, rep(NA,length=(101-length(tempi))))
    out.v2.FD.c[,i]=repc
    out.v2.FD.i[,i]=repi
  }
  
  out.v2.FD.c[is.na(out.v2.FD.c)] <- 0
  out.v2.FD.i[is.na(out.v2.FD.i)] <- 0
 
  seed<-rep(exposed[[h]][j],length=101)
  herd<-rep(herdsize[h],length=101)

	out.v2.FD.c=cbind(herd,seed,out.v2.FD.c)
	out.v2.FD.i=cbind(herd,seed,out.v2.FD.i)

  
  filenameC=paste("Var2bFDC",herdsize[h],"S",exposed[[h]][j],"sens",k,"csv",sep=".")
  filenameI=paste("Var2bFDI",herdsize[h],"S",exposed[[h]][j],"sens",k,"csv",sep=".")
  
  write.csv(out.v2.FD.c,file=filenameC)
  write.csv(out.v2.FD.i,file=filenameI)
  
}
}
  }



#Density Dependent, infectious first, nasal fluid parameters


  transitions.v2 = list(
    c(S = -1, E = +1), # transmission
    c(E = -1, I = +1), #becoming infectious
    c(I = -1, C = +1), # becoming clinical
    c(C = -1, R = +1) # recovery
  )
  
  Var2DD <- function(x, p, t) {
    return(c(p[1]*x['S']*(x['I']+x['C']), #transmission rate
             p[2]*x['E'],#becoming infectious rate
             p[3]*x['I'], #becoming clinical rate
             p[4]*x['C']  )) #recovery rate
  }
  



for(h in 1:length(herdsize)){
for(j in 1:length(exposed[[h]])){
  
   init.values.v2 = c(
    S = herdsize[h]-exposed[[h]][j], # susceptible animals--base herd size
    E = exposed[[h]][j],   # exposed animals (1% of herd)
    I = 0, #infectious not clinical
    C = 0,   #clinical and  infectious
    R = 0)    # recovered

for(k in 1:dim(parms.set2n)[1]){

numsims=1000

parms2 <- parms.set2n[k,]

  out.v2.DD.c=matrix(nrow=101,ncol=1000)
  out.v2.DD.c<-as.data.frame(out.v2.DD.c)
  
  out.v2.DD.i=matrix(nrow=101,ncol=1000)
  out.v2.DD.i<-as.data.frame(out.v2.DD.i)
  
  set.seed(478394)
  
  for(i in 1:1000){
    sim=ssa.adaptivetau(init.values.v2, transitions.v2, Var2DD, parms2, tf=100)
    sim=as.data.frame(sim)
    idx <- diff(ceiling(sim$time)) == 1
    tempc<-sim[idx,4]
    tempi<-sim[idx,5]
    repc=c(tempc, rep(NA,length=(101-length(tempc))))
    repi=c(tempi, rep(NA,length=(101-length(tempi))))
    out.v2.DD.c[,i]=repc
    out.v2.DD.i[,i]=repi
  }
  
  out.v2.DD.c[is.na(out.v2.DD.c)] <- 0
  out.v2.DD.i[is.na(out.v2.DD.i)] <- 0
 
  seed<-rep(exposed[[h]][j],length=101)
  herd<-rep(herdsize[h],length=101)

	out.v2.DD.c=cbind(herd,seed,out.v2.DD.c)
	out.v2.DD.i=cbind(herd,seed,out.v2.DD.i)

  
  filenameC=paste("Var2nDDC",herdsize[h],"S",exposed[[h]][j],"sens",k,"csv",sep=".")
  filenameI=paste("Var2nDDI",herdsize[h],"S",exposed[[h]][j],"sens",k,"csv",sep=".")
  
  write.csv(out.v2.DD.c,file=filenameC)
  write.csv(out.v2.DD.i,file=filenameI)
  
}
}
}

  #Frequency Dependent, infectious first, nasal fluid parameters
  
Var2FD <- function(x, p, t) {
  return(c(p[1]*x['S']*(x['I']+x['C'])/(x['S']+x['E']+x['I']+x['C']+x['R']), #transmission rate
           p[2]*x['E'],#becoming infectious rate
           p[3]*x['I'], #becoming clinical rate
           p[4]*x['C']  )) #recovery rate
}



for(h in 1:length(herdsize)){
for(j in 1:length(exposed[[h]])){
  
    init.values.v2 = c(
    S = herdsize[h]-exposed[[h]][j], # susceptible animals--base herd size
    E = exposed[[h]][j],   # exposed animals (1% of herd)
    I = 0, #infectious not clinical
    C = 0,   #cincial and  infectious
    R = 0)    # recovered
  
for(k in 1:dim(parms.set2n)[1]){

numsims=1000

parms2 <- parms.set2n[k,]
 
  
  out.v2.FD.c=matrix(nrow=101,ncol=1000)
  out.v2.FD.c<-as.data.frame(out.v2.FD.c)
  
  out.v2.FD.i=matrix(nrow=101,ncol=1000)
  out.v2.FD.i<-as.data.frame(out.v2.FD.i)
  
  set.seed(478394)
  
  for(i in 1:1000){
    sim=ssa.adaptivetau(init.values.v2, transitions.v2, Var2FD, parms2, tf=100)
    sim=as.data.frame(sim)
    idx <- diff(ceiling(sim$time)) == 1
    tempc<-sim[idx,4]
    tempi<-sim[idx,5]
    repc=c(tempc, rep(NA,length=(101-length(tempc))))
    repi=c(tempi, rep(NA,length=(101-length(tempi))))
    out.v2.FD.c[,i]=repc
    out.v2.FD.i[,i]=repi
  }
  
  out.v2.FD.c[is.na(out.v2.FD.c)] <- 0
  out.v2.FD.i[is.na(out.v2.FD.i)] <- 0
 
  seed<-rep(exposed[[h]][j],length=101)
  herd<-rep(herdsize[h],length=101)

	out.v2.FD.c=cbind(herd,seed,out.v2.FD.c)
	out.v2.FD.i=cbind(herd,seed,out.v2.FD.i)

  
  filenameC=paste("Var2nFDC",herdsize[h],"S",exposed[[h]][j],"sens",k,"csv",sep=".")
  filenameI=paste("Var2nFDI",herdsize[h],"S",exposed[[h]][j],"sens",k,"csv",sep=".")
  
  write.csv(out.v2.FD.c,file=filenameC)
  write.csv(out.v2.FD.i,file=filenameI)
  
}
}
}


#Density Dependent, base

   transitions = list(
    c(S = -1, E = +1), # transmission
    c(E = -1, I = +1), # becoming infectious
    c(I = -1, R = +1) # recovery
  )
  BaseDD <- function(x, p, t) {
    return(c(p[1]*x['S']*x['I'], #transmission rate
             p[2]*x['E'],#becoming infectious rate
             p[3]*x['I']  )) #recovery rate
  }



for(h in 1:length(herdsize)){
for(j in 1:length(exposed[[h]])){
  
  init.values = c(
    S = herdsize[h]-exposed[[h]][j], # susceptible animals--base herd size
    E = exposed[[h]][j],   # exposed animals (1% of herd)
    I = 0,   # infectious
    R = 0)    # recovered 
  
for(k in 1:dim(parms.set)[1]){

numsims=1000

parms <- parms.set[k,]

  
  out.base.DD=matrix(nrow=101,ncol=1000)
  out.base.DD<-as.data.frame(out.base.DD)
  

  set.seed(478394)
  
  for(i in 1:1000){
     sim=ssa.adaptivetau(init.values, transitions, BaseDD, parms, tf=100)
    sim=as.data.frame(sim)
    idx <- diff(ceiling(sim$time)) == 1
    temp<-sim[idx,4]
    rep=c(temp, rep(NA,length=(101-length(temp))))
    out.base.DD[,i]=rep
  }
  
  out.base.DD[is.na(out.base.DD)] <- 0
 
  seed<-rep(exposed[[h]][j],length=101)
  herd<-rep(herdsize[h],length=101)

	out.base.DD=cbind(herd,seed,out.base.DD)

  
  filename=paste("BaseDDI",herdsize[h],"S",exposed[[h]][j],"sens",k,"csv",sep=".")
  
  write.csv(out.base.DD,file=filename)
  
}
}
}

  #Frequency Dependent, base
  
BaseFD <- function(x, p, t) {
  return(c(p[1]*x['S']*x['I']/(x['S']+x['E']+x['I']+x['R']), #transmission rate
           p[2]*x['E'],#becoming infectious rate
           p[3]*x['I']  )) #recovery rate
}

for(h in 1:length(herdsize)){
for(j in 1:length(exposed[[h]])){
   
  init.values = c(
    S = herdsize[h]-exposed[[h]][j], # susceptible animals--base herd size
    E = exposed[[h]][j],   # exposed animals (1% of herd)
    I = 0,   # infectious
    R = 0)    # recovered 
  
 for(k in 1:dim(parms.set)[1]){

numsims=1000

parms <- parms.set[k,]

  out.base.FD=matrix(nrow=101,ncol=1000)
  out.base.FD<-as.data.frame(out.base.FD)
  

  set.seed(478394)
  
  for(i in 1:1000){
     sim=ssa.adaptivetau(init.values, transitions, BaseFD, parms, tf=100)
    sim=as.data.frame(sim)
    idx <- diff(ceiling(sim$time)) == 1
    temp<-sim[idx,4]
    rep=c(temp, rep(NA,length=(101-length(temp))))
    out.base.FD[,i]=rep
  }
  
  out.base.FD[is.na(out.base.FD)] <- 0
 
  seed<-rep(exposed[[h]][j],length=101)
  herd<-rep(herdsize[h],length=101)

	out.base.DD=cbind(herd,seed,out.base.DD)

  
  filename=paste("BaseFDI",herdsize[h],"S",exposed[[h]][j],"sens",k,"csv",sep=".")
  
  write.csv(out.base.FD,file=filename)
  }
}
}
