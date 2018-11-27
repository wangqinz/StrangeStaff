install.packages("runjags")
library(runjags)

#### Poisson

pois="model{

for (i in 1:k){
    for (j in 1:N){
        Y[i,j]~dpois(exp(bn[i]+u[j])) }
}

for (i in 1:k){
    b[i]~dnorm(0,1)
}

m=mean(b)

for (i in 1:k){
    bn[i]=b[i]-m
}

for(i in 1:N){
    u[i]~dunif(-1000,1000)
}

}"


#### multi

multi="model{

for (j in 1:N){
    Y[1:k,j]~dmulti(p,n)
}

for (i in 1:k){
    b[i]~dnorm(0,1)
}

m=mean(b)

for (i in 1:k){
    bn[i]=b[i]-m
}

for (i in 1:k){
    p[i]=exp(bn[i])
}


}"

N=10
n=100 p=c(0.1,0.5,0.4) k=length(p) Y=rmultinom(N,n,p)

r=run.jags(model=pois,monitor=c("bn"),data=list(Y=Y,k=k,N=N))
rr=run.jags(model=multi,monitor=c("bn"),data=list(Y=Y,k=k,N=N,n=n))

plot(r)
