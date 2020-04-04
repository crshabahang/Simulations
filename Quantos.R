#Cameron R. Shabahang
#Pricing of FX Quanto Options
#Function to establish robust linear congruential generator
RobustLCG = function(numran)
{
  #Seeds
  Im1<-2147483563;
  Im2<-2147483399;
  #Inverse of first seed
  Am<-1.0/Im1;
  #Index back
  Imm1<-Im1-1;
  Ia1<-40014;
  Ia2<-40692;
  Iq1<-53668;
  Iq2<-52774;
  Ir1<-12211;
  Ir2<-3791;
  #Generate random numbers between 1 and 32
  NTAB<-32;
  NDIV<-as.integer(1+Imm1/NTAB);
  EPS<-1.2e-7;
  RNMX<-1.-EPS;
  idum2<-123456789;
  idum<-1000;
  
  ran2<-0;
  iy<-0;
  iv<-rep(0,NTAB);
  random<-rep(0,numran);
  icount<-1;
  for(icount in 1:numran){
    if(idum<=0){
      idum<-max(-idum,1);
      idum2<-idum;
      j<-NTAB+8;
      while(j>0){
        k<-as.integer(idum/Iq1);
        idum<-Ia1*(idum-k*Iq1)-k*Ir1;
        if(idum<0){idum<-idum+Im1}
        if(j<=NTAB){iv[j]<-idum}
        j<j-1
      }
      iy<-iv[1]
    }
    k<-as.integer(idum/Iq1)
    idum<-Ia1*(idum-k*Iq1)-k*Ir1
    if(idum<0){idum=idum+Im1}
    k<-as.integer(idum2/Iq2)
    idum2<-Ia2*(idum2-k*Iq2)-k*Ir2
    if(idum2<0){idum2<-idum2+Im2}
    j<-as.integer(iy/NDIV)+1
    iy<-iv[j]-idum2
    iv[j]<-idum
    if(iy<1){iy<-iy+Imm1}
    ran2<-min(Am*iy,RNMX)
    random[icount]<-ran2
    icount<-icount+1
  }
  #Returns the random number 
  return(random)
}

#Generate Inverse Normal Distribution
#Beasley-Springer-Moro algorithm
InverseNormal <- function(u){
  
  a0=2.50662823884
  a1=-18.61500062529
  a2=41.39119773534
  a3=-25.44106049637
  b0=-8.47351093090
  b1=23.08336743743
  b2=-21.06224101826
  b3=3.13082909833
  c0=0.3374754822726147
  c1=0.9761690190917186
  c2=0.1607979714918209
  c3=0.0276438810333863
  c4=0.0038405729373609
  c5=0.0003951896511919
  c6=0.0000321767881768
  c7=0.0000002888167364
  c8=0.0000003960315187
  y =  u-0.5
  
  if (abs(y)<0.42)
  {
    r = y*y
    x = y*(((a3*r+a2)*r+a1)*r+a0)/((((b3*r+b2)*r+b1)*r+b0)*r+1)
    
  }
  else
  {
    r = u
    if (y>0){r = 1-u}
    r = log(-log(r))
    x = c0+r*(c1+r*(c2+r*(c3+r*(c4+r*(c5+r*(c6+r*(c7+r*c8)))))))
    if(y<0){x=-x}
    
  }
  return(x)
}
# Implement uniform random number generator (INPUT: seed and numran)
# FX rates quoted as foreign currency/USD
seed=1000
numran=500000
inputvar=c(seed,numran)
rand_uniform_c=RobustLCG(inputvar)
rand_norm1=InverseNormal(x)
seed=2000
numran=500000
inputvar=c(seed,numran)
rand_uniform_=RobustLCG(inputvar)
rand_norm2=InverseNormal(x)
rand_eps1=rand_norm1
rho_c=0.25
rand_eps2<-rho_c*rand_norm1+(sqrt(1-rho_c^2))*rand_norm2
cor(rand_eps1,rand_eps2)

# Price paths and average
# FX rates quotes
S0<-100
rUSD<-0.03
VolStock<-0.25
T=3/12
FX0<-0.8
rFX<-0.02
VolFX<-0.2
KFX<-S0*FX0
numpath<-180000
ST_Vals<-c(rep(0),numpath)
ST_Rets<-c(rep(0),numpath)
FX_Vals<-c(rep(0),numpath)
FX_Rets<-c(rep(0),numpath)
Quanto_Vals<-c(rep(0),numpath)
jcount<-1
while(jcount <= numpath) {
  ST<-S0*(exp((rUSD-0.5*(VolStock^2))*T+VolStock*sqrt(T)*rand_eps1[jcount]))
  ST_Vals[jcount]<-ST
  ST_Rets[jcount]<-log(ST_Vals[jcount]/S0)
  FX<-FX0*(exp((rFX-rUSD-0.5*(VolFX^2))*T+VolFX*sqrt(T)*rand_eps2[jcount]))
  FX_Vals[jcount]<-FX
  FX_Rets[jcount]<-log(FX_Vals[jcount]/FX0)
  Quanto<-max(ST*FX-KFX,0)
  Quanto_Vals[jcount]<-Quanto
  jcount=jcount+1
  #  cat("A",ST",ST_Vals[12,jcount],ST_Rets[jcount],"\n")
}
#
mean(Quanto_Vals*exp(-rFX*T))
sd(Quanto_Vals*exp(-rFX*T))
Stderr_Quanto<-sd(Quanto_Vals*exp(-rFX*T))/sqrt(numpath)
mean(Quanto_Vals*exp(-rFX*T))-Stderr_Quanto*qnorm(0.975)
mean(Quanto_Vals*exp(-rFX*T))+Stderr_Quanto*qnorm(0.975)
#
hist(FX_Rets, breaks=20,freq=F)
FXRet_Vol<-VolFX*sqrt(T)
rFXm=exp(rFX*T)-1
curve(dnorm(x, mean=rFXm, sd=FXRet_Vol),from=-4, to=4,add=TRUE,lwd=2)
#
hist(ST_Rets, breaks=20,freq=F)
Vol_Ret<-VolStock*sqrt(T)
rUSDm=exp(rUSD*T)-1
curve(dnorm(x, mean=rUSDm, sd=Vol_Ret),from=-4, to=4,add=TRUE,lwd=2)


