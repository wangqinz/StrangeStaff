# Newton method (TWO PWAYS)

# PART 1: from STAT-S 611 class
## including an example

########
number<-1000
mean<-10
sd<-2

x <- rnorm(number, mean, sd)
n <- length(x)
mx <- mean(x)
SS <- (n - 1)*var(x)

#### closed form
mx
var(x) * (n - 1) / n

#### Good strating points
g <- log(5)
mu <- 9

#### diverging starting points
g <- log(1000)
mu <- 30

#### method loop
for (i in 1:20){
  loglike <- n/2 * g + 0.5*exp(-g) * (SS + n * (mx - mu)^2)
  cat('mu[',i,'] = ', mu, '\tsigma^2[',i,'] = ', exp(g), '\t loglike =' , -loglike, '\n', sep ='')
  emg <- exp(-g)
  grad <- c(n * (mu - mx) * emg, 
            n/2 - emg/2*(SS + n*(mx - mu)^2)
  )
  hess <- matrix(c(n * emg, -n * (mu - mx) * emg, -n * (mu - mx) * emg, emg/2*(SS + n*(mx - mu)^2)), ncol = 2)
  if (prod(eigen(hess)$values == 0)) {
    warning('Hessian matrix not positive definite in iteration ',i)
  }
  th <- c(mu, g) - solve(hess) %*% grad
  mu <- th[1]
  g <- th[2]
}

mu
exp(g)


########
#### Scoring
#### choose some different starting points

g <- log(4)
mu <- 10

#### previously diverging starting points
g <- log(1000)
mu <- 30

tol <- 1e-15
maxiter <- 100
iter <- 1
loglike <- -n/2 * g - 0.5*exp(-g) * (SS + n * (mx - mu)^2)

repeat{
  cat('mu[',iter,'] = ', mu, '\tsigma^2[',iter,'] = ', exp(g), '\t loglike =' , loglike, '\n', sep ='')
  emg <- exp(-g)
  grad <- c(n * (mu - mx) * emg, 
            n/2 - emg/2*(SS + n*(mx - mu)^2)
  )
  info <- 1/n * diag(c(1/emg, 2))
  th <- c(mu, g) - info %*% grad
  mu <- th[1]
  g <- th[2]
  loglike_new <- -n/2 * g - 0.5*exp(-g) * (SS + n * (mx - mu)^2)
  if (iter == maxiter) stop("maxiter reached!")
  if ( abs(loglike_new - loglike)/abs(loglike) < tol) {
    loglike <- loglike_new
    break
  }
  loglike <- loglike_new
  iter <- iter + 1
}

mu
exp(g)

#### closed form
mx
var(x) * (n - 1) / n

# PART 2: self-made, strange, ugly but useful code
#### provide an example








n <- 20

xbar <- 28.0

x2sum <- 16500

fun.1 <- function(alpha){
    n <- 20
    xbar <- 28.0
    x2sum <- 16500
    return(xbar*xbar*gamma(1+2/alpha)-1/n*x2sum*gamma(1+1/alpha)*gamma(1+1/alpha))
}

fun.2 <- function(alpha){
    n <- 20
    xbar <- 28.0
    x2sum <- 16500
    return(-xbar*xbar*2/alpha/alpha*digamma(1+2/alpha)*gamma(1+2/alpha)+x2sum/n*2/alpha/alpha*gamma(1+1/alpha)*gamma(1+1/alpha)*digamma(1+1/alpha))
}

alpha.0 <- 1.0

alpha.1 <- alpha.0 - fun.1(alpha.0)/fun.2(alpha.0)

eps <- 10^(-6)

while (abs(alpha.1-alpha.0) > eps){
    print(alpha.0)
    print(alpha.1)
    alpha.0 <- alpha.1
    alpha.1 <- alpha.0 - fun.1(alpha.0)/fun.2(alpha.0)
}

beta.est <- xbar/gamma(1+1/alpha.1)

#adding iteration time

alpha.0 <- 1.0

alpha.1 <- alpha.0 - fun.1(alpha.0)/fun.2(alpha.0)

eps <- 10^(-6)

i.iter <- 1

while (abs(alpha.1-alpha.0) > eps & (i.iter < 100)){
    print(alpha.0)
    print(alpha.1)
    print(i.iter)
    alpha.0 <- alpha.1
    alpha.1 <- alpha.0 - fun.1(alpha.0)/fun.2(alpha.0)
    i.iter <- i.iter+1
}

beta.est <- xbar/gamma(1+1/alpha.1)

# bisection

i.iter <- 1

alpha.0 <- 10^(-1)

alpha.1 <- 100

eps <- 10^(-6)

while (abs(alpha.1-alpha.0) > eps & (i.iter < 100)){
    print(alpha.0)
    print(alpha.1)
    print(i.iter)
    if (sign(fun.1(alpha.0)) == sign(fun.1(alpha.1))) { print("find new initial interval");break }
    mid.alpha <- (alpha.0+alpha.1)/2
    if (sign(fun.1(mid.alpha)) == sign(fun.1(alpha.0))){
        alpha.0 <- mid.alpha
    }
    else{
        alpha.1 <- mid.alpha
    }
    i.iter <- i.iter+1
}
