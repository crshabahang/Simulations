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

#Price for one path. 
#The price at each interval along the path is noted.
#The maximum of the prices at each interval is the payoff
PricePath = function(S0,K,sigma,r,tau,inter,dW,type)
{
  #Two-dimensional data structure
  Prices = matrix(0,inter,1)
  S = S0
  #Interval is the time to expiry divided by the number of intervals
  dt = tau/inter
  for (i in (1:inter))
  {
    #Stochastic stock price process
    Prices[i] = S+S*((r)*dt+sigma*dW[i]*sqrt(dt))
    S = Prices[i]
  }
  if (type == 'a')
  {
    #Fixed lookback call option 
    LookBackPrice = max(0,(max(Prices)-K))  
  }
  if (type == 'b')
  {
    #Floating lookback call option
    LookBackPrice = max(0,(Prices[inter]-min(Prices)))
  }
  
  return(LookBackPrice)
}


#Lookback Options
LookBackOption = function(S0,K,sigma,r,tau,inter,paths,type)
{
  print('count')
  print(type)
  #Generate random numbers 
  W = NormalDistGenerator(inter*paths)
  dW = matrix(0,paths,inter)
  #Records the price history for each path
  PriceHistory = matrix(0,paths,inter)
  
  #Arrange into rows corresponding with the number of paths
  #Each row contains elements corresponding with the number of elements 
  for (i in (1:paths))
  {
    for (j in (1:inter))
    {
      dW[i,j] = W[(i-1)*inter + j]
    }
  }
  
  #Calculate price for each path via discounting
  intervalPrice = matrix(0,paths,1)
  
  for (i in (1:paths))
  {
    intervalPrice[i] = exp(-r*tau)*(PricePath(S0,K,sigma,r,tau,inter,dW[i,],type))
  }
  
  #Returns a list named 'output' with the first element being the mean and the second variance
  output = matrix(0,2,1)
  output[1] = mean(intervalPrice)
  output[2] = var(intervalPrice)
  return(output)
}


#Option Type A
P1000A = LookBackOption(100,105,0.5,0.05,1,12,1000,'a')
P10000A = LookBackOption(100,105,0.5,0.05,1,12,10000,'a')
P100000A = LookBackOption(100,105,0.5,0.05,1,12,100000,'a')
print('For the fixed lookback call option, the prices are: ')
cat('1) 1000 paths, Price :',P1000A[1],'Variance :',P1000A[2])
cat('2) 10000 paths, Price :',P10000A[1],'Variance :',P10000A[2])
cat('3) 100000 paths, Price :',P100000A[1],'Variance :',P100000A[2])

#95% confidence intervals

error = qnorm(0.975)*sqrt(P1000A[2])/sqrt(1000)
left = P1000A[1]-error
right = P1000A[1]+error
print('The 95% confidence interval for 1000 paths is:')
print(left)
print(right)

error = qnorm(0.975)*sqrt(P10000A[2])/sqrt(10000)
left = P10000A[1]-error
right = P10000A[1]+error
print('The 95% confidence interval for 10000 paths is:')
print(left)
print(right)

error = qnorm(0.975)*sqrt(P100000A[2])/sqrt(100000)
left = P100000A[1]-error
right = P100000A[1]+error
print('The 95% confidence interval for 100000 paths is:')
print(left)
print(right)

#Option Type B
P1000B = LookBackOption(100,105,0.5,0.05,1,12,1000,'b')
P10000B = LookBackOption(100,105,0.5,0.05,1,12,10000,'b')
P100000B = LookBackOption(100,105,0.5,0.05,1,12,100000,'b')
print('For the floating lookback call option, the prices are:')
cat('1) 1000 paths, Price :',P1000B[1],'Variance :',P1000B[2])
cat('2) 10000 paths, Price :',P10000B[1],'Variance :',P10000B[2])
cat('3) 10000 paths, Price :',P100000B[1],'Variance :',P100000B[2])

#95% confidence intervals

error = qnorm(0.975)*sqrt(P1000B[2])/sqrt(1000)
left = P1000B[1]-error
right = P1000B[1]+error
print('The 95% confidence interval for 1000 paths is:')
print(left)
print(right)

error = qnorm(0.975)*sqrt(P10000B[2])/sqrt(10000)
left = P10000B[1]-error
right = P10000B[1]+error
print('The 95% confidence interval for 10000 paths is:')
print(left)
print(right)

error = qnorm(0.975)*sqrt(P100000B[2])/sqrt(100000)
left = P100000B[1]-error
right = P100000B[1]+error
print('The 95% confidence interval for 100000 paths is:')
print(left)
print(right)

