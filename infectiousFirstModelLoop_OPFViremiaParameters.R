#Code to run the Infectious First model with parameters from viremia in oesophageal-pharyngeal fluid from Beck-Johnson et al.for all herd sizes and seed sizes
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



  transitions.v2 = list(
    c(S = -1, E = +1), # transmission
    c(E = -1, I = +1), #becoming infectious
    c(I = -1, C = +1), # becoming clincial
    c(C = -1, R = +1) # recovery
  )
  
  
  #Density dependent Infectious First  model
  
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
    C = 0,   #cincial and  infectious
    R = 0)    # recovered

   numsims=1000
  
parms2 <- c(beta=21.84, nu=(1/0.46), rho=(1/3.67), mu=(1/4.45))
  
  
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

  
  filenameC=paste("Var2oDDC",herdsize[h],"S",exposed[[h]][j],"csv",sep=".")
  filenameI=paste("Var2oDDI",herdsize[h],"S",exposed[[h]][j],"csv",sep=".")
  
  write.csv(out.v2.DD.c,file=filenameC)
  write.csv(out.v2.DD.i,file=filenameI)
  
}
}


  #Frequency dependent Infectious First  model
  
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
  
   numsims=1000
  
parms2 <- c(beta=21.84, nu=(1/0.46), rho=(1/3.67), mu=(1/4.45))
  
  
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
