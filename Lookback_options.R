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
#Inverse transform sampling 
NormalDistGenerator = function(NumElements)
{
  #Implements the robest LCG function on a number of elements
  uniform = RobustLCG(NumElements)
  NormalDist = matrix(0,NumElements,1)
  for (i in(1:NumElements))
  {
    #Inverse normal of the uniform distribution yields the normal distribution
    NormalDist[i] = InverseNormal(uniform[i])
  }
  return(NormalDist)
}

#Price for one path. As per problem statement, the price at each interval along the path is noted.
#The maximum of these, is considered in the pay off as shown in line 122
PricePath = function(S0,K,sigma,r,tau,NInt,dZ,OptionType)
{
  
  Prices = matrix(0,NInt,1)
  S = S0
  dt = tau/NInt
  for (i in (1:NInt))
  {
    Prices[i] = S+S*((r-sigma)*dt+sigma*dZ[i]*sqrt(dt))
    S = Prices[i]
  }
  if (OptionType == 'a')
  {
    LBKPrice = max(0,(max(Prices)-K))  
  }
  if (OptionType == 'b')
  {
    LBKPrice = max(0,(K-max(Prices)))
  }
  
  return(LBKPrice)
}


#Look Back Options
LookBackOption = function(S0,K,sigma,r,tau,NInt,NPaths,OptionType)
{
  print('here')
  print(OptionType)
  #Generate #intervals x #paths random numbers 
  Z = NormalDistGenerator(NInt*NPaths)
  dZ = matrix(0,NPaths,NInt)
  #Stores the price trajectory for each path
  PriceHistory = matrix(0,NPaths,NInt)
  
  #Arrange into #path rows, with each row containing #Int elements
  for (i in (1:NPaths))
  {
    for (j in (1:NInt))
    {
      dZ[i,j] = Z[(i-1)*NInt + j]
    }
  }
  
  #Calculate Price for each path
  PriceByPath = matrix(0,NPaths,1)
  
  
  for (i in (1:NPaths))
  {
    
    PriceByPath[i] = exp(-r*tau)*(PricePath(S0,K,sigma,r,tau,NInt,dZ[i,],OptionType))
  }
  
  #This list named 'Output' is returned. 1st element is Mean, 2nd element is Variance
  Output = matrix(0,2,1)
  Output[1] = mean(PriceByPath)
  Output[2] = var(PriceByPath)
  return(Output)
}


#Option Type A with 1000 Paths
P1000A = LookBackOption(100,105,0.5,0.05,1,12,1000,'a')
P10000A = LookBackOption(100,105,0.5,0.05,1,12,10000,'a')
P100000A = LookBackOption(100,105,0.5,0.05,1,12,100000,'a')
print('For Option type A, we have')
cat('1) 1000 paths, Price :',P1000A[1],'Variance :',P1000A[2])
cat('2) 10000 paths, Price :',P10000A[1],'Variance :',P10000A[2])
cat('3) 100000 paths, Price :',P100000A[1],'Variance :',P100000A[2])

P1000B = LookBackOption(100,105,0.5,0.05,1,12,1000,'b')
P10000B = LookBackOption(100,105,0.5,0.05,1,12,10000,'b')
P100000B = LookBackOption(100,105,0.5,0.05,1,12,100000,'b')
print('For Option type B, we have')
cat('1) 1000 paths, Price :',P1000B[1],'Variance :',P1000B[2])
cat('2) 10000 paths, Price :',P10000B[1],'Variance :',P10000B[2])
cat('3) 10000 paths, Price :',P100000B[1],'Variance :',P100000B[2])
