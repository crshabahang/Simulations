#LCG function
lcg <- function(a,c,m,run.length,seed) {
  x <- rep(0,run.length)
  
  x[1] <- seed
  for (i in 1:(run.length-1)) {
    x[i+1] <- (a * x[i] + c) %% m
  }
  U <- x/m # scale all of the x's to
  # produce uniformly distributed
  # random numbers between [0,1)
  return(list(U=U))
}

#Beasley-Springer-Moro algorithm for approximating the inverse normal
#Input: u, a scalar or matrix with elements between 0 and 1
#Output: x, an approximation for the inverse normal at u
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
  y = u - 0.5
  
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
    if(y<0){x = -x}
    
  }
  return(x)
  
}


#Generate Unifrom Distribution
UDF<-lcg(6,7,23,101,-1000)

UDF <- as.vector(unlist(UDF))

uniform = matrix(0,100,1)

for (i in (1:100))
{
  uniform[i] = UDF[i+1]
}

uniform <- as.vector(unlist(uniform))

part1a = 1 + uniform*1000





#UDF for roulette
UDF<-lcg(6,7,23,101,-1000)

UDF <- as.vector(unlist(UDF))

uniform = matrix(0,100,1)

for (i in (1:100))
{
  uniform[i] = UDF[i+1]
}

#
roulette = matrix(0,10,1)

for (i in (1:10))
{
  roulette[i] = uniform[i*5]
  
}

#10 plays. My number is 32
roulette = floor(1 + 36*roulette)
payoff = -10

for (i in (1:10))
  
{
  if (roulette[i]==32)
  {
    payoff = payoff + 36
  }
}

roulette = matrix(0,100,1)

for (i in (1:100))
{
  roulette[i] = uniform[i]
  
}

#100 plays. My number is 32
roulette = floor(1 + 36*roulette)
payoff2 = -100
for (i in (1:100))
  
{
  if (roulette[i]==32)
  {
    payoff2 = payoff2 + 36
  }
}

roulette = matrix(0,100,1)

for (i in (1:100))
{
  roulette[i] = uniform[i]
  
}


#1000 plays
UDF<-lcg(6,7,23,1001,-1000)

UDF <- as.vector(unlist(UDF))

uniform = matrix(0,1000,1)

for (i in (1:1000))
{
  uniform[i] = UDF[i+1]
}

#
roulette = matrix(0,1000,1)

for (i in (1:1000))
{
  roulette[i] = uniform[i]
  
}

roulette = floor(1 + 36*roulette)

payoff3 = -1000
for (i in (1:1000))
  
{
  if (roulette[i]==32)
  {
    payoff3 = payoff3 + 36
  }
}

#SUnfair roulette

roulette = matrix(0,1000,1)

for (i in (1:1000))
{
  roulette[i] = uniform[i]
  
}

print(roulette)

roulette = floor(1 + 37*roulette)

print(roulette)


payoff4 = -1000
for (i in (1:1000))
  
{
  if (roulette[i]==32)
  {
    payoff4 = payoff4 + 36
  }
}

#Simulating 10,000 Draws from a normal distribution 
UDF<-lcg(6,7,23,10001,-1000)

UDF <- as.vector(unlist(UDF))

uniform = matrix(0,10000,1)

for (i in (1:10000))
{
  uniform[i] = UDF[i+1]
}


NormalDist = matrix(0,10000,1)

for (i in(1:10000))
{
  NormalDist[i] = InverseNormal(uniform[i])
}


#Simulating 100 Draws from a normal distribution 
UDF<-lcg(6,7,23,101,-1000)

UDF <- as.vector(unlist(UDF))

uniform = matrix(0,100,1)

for (i in (1:100))
{
  uniform[i] = UDF[i+1]
}


NormalDist = matrix(0,100,1)

for (i in(1:100))
{
  NormalDist[i] = InverseNormal(uniform[i])
}

hist(NormalDist,freq=FALSE,breaks=10)
