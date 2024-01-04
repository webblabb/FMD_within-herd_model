
#Code to run the base model from Beck-Johnson et al. for all herd sizes and seed sizes
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



 transitions = list(
    c(S = -1, E = +1), # transmission
    c(E = -1, I = +1), # becoming infectious
    c(I = -1, R = +1) # recovery
  )
 
 #Density dependent base model
 
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
  
    numsims=1000
  
  parms <- c(beta=21.84, sigma=(1/4.55), gamma=(1/1.3))
  
  
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

  
  filename=paste("BaseDDI",herdsize[h],"S",exposed[[h]][j],"csv",sep=".")
  
  write.csv(out.base.DD,file=filename)
  
}
}


#Frequency dependent base model

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
  
    numsims=1000
  
  parms <- c(beta=21.84, sigma=(1/4.55), gamma=(1/1.3))
  
  
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

  
  filename=paste("BaseFDI",herdsize[h],"S",exposed[[h]][j],"csv",sep=".")
  
  write.csv(out.base.FD,file=filename)
  }
}
