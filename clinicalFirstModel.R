
#Code to run Clinical First model from Beck-Johnson et al.for all herd sizes and seed sizes

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



transitions.v1 = list(
  c(S = -1, E = +1), # transmission
  c(E = -1, C = +1), #becoming clinical
  c(C = -1, I = +1), # becoming infectious
  c(I = -1, R = +1) # recovery
)

#Density dependent Clinical First  model

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
    C = 0, #clincial not infectious
    I = 0,   # infectious
    R = 0)    # recovered
  
   numsims=1000
  
  parms1 <- c(beta=21.84, phi=(1/4.03), omega=(1/0.52), theta=(1/1.3))
  
  
  out.v1.DD.e=matrix(nrow=101,ncol=1000)
  out.v1.DD.e<-as.data.frame(out.v1.DD.e)
  
  out.v1.DD.c=matrix(nrow=101,ncol=1000)
  out.v1.DD.c<-as.data.frame(out.v1.DD.c)
  
  out.v1.DD.i=matrix(nrow=101,ncol=1000)
  out.v1.DD.i<-as.data.frame(out.v1.DD.i)
  
  set.seed(478394)
  
  for(i in 1:1000){
    sim=ssa.adaptivetau(init.values.v1, transitions.v1, Var1DD, parms1, tf=100)
    sim=as.data.frame(sim)
    idx <- diff(ceiling(sim$time)) == 1
    tempe<-sim[idx,3]
    tempc<-sim[idx,4]
    tempi<-sim[idx,5]
    repe=c(tempe, rep(NA,length=(101-length(tempc))))
    repc=c(tempc, rep(NA,length=(101-length(tempc))))
    repi=c(tempi, rep(NA,length=(101-length(tempi))))
    out.v1.DD.e[,i]=repe
    out.v1.DD.c[,i]=repc
    out.v1.DD.i[,i]=repi
  }
  
  out.v1.DD.e[is.na(out.v1.DD.e)] <- 0
  out.v1.DD.c[is.na(out.v1.DD.c)] <- 0
  out.v1.DD.i[is.na(out.v1.DD.i)] <- 0
 
  seed<-rep(exposed[[h]][j],length=101)
  herd<-rep(herdsize[h],length=101)

  out.v1.DD.e=cbind(herd,seed,out.v1.DD.e)
  out.v1.DD.c=cbind(herd,seed,out.v1.DD.c)
	out.v1.DD.i=cbind(herd,seed,out.v1.DD.i)

	filenameE=paste("Var1DDE",herdsize[h],"S",exposed[[h]][j],"csv",sep=".")
  filenameC=paste("Var1DDC",herdsize[h],"S",exposed[[h]][j],"csv",sep=".")
  filenameI=paste("Var1DDI",herdsize[h],"S",exposed[[h]][j],"csv",sep=".")
  
  write.csv(out.v1.DD.e,file=filenameE)
  write.csv(out.v1.DD.c,file=filenameC)
  write.csv(out.v1.DD.i,file=filenameI)
  
}
}


#Frequency dependent Clinical First  model


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
    C = 0, #clincial not infectious
    I = 0,   # infectious
    R = 0)    # recovered
  
   numsims=1000
  
  parms1 <- c(beta=21.84, phi=(1/4.03), omega=(1/0.52), theta=(1/1.3))
  
  
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

  
  filenameC=paste("Var1FDC",herdsize[h],"S",exposed[[h]][j],"csv",sep=".")
  filenameI=paste("Var1FDI",herdsize[h],"S",exposed[[h]][j],"csv",sep=".")
  
  write.csv(out.v1.FD.c,file=filenameC)
  write.csv(out.v1.FD.i,file=filenameI)
  
}
}


