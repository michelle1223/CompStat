#Exer 4.2
x = data.frame(enc=0:16,
               freq=c(379,299,222,145,109,95,73,59,45,30,24,12,4,2,0,1,1))
N = sum(x$freq)
y = rep(x$enc,x$freq)


#For Poisson distribution, mean and variance should be approximately equal
show(c(mean(y),var(y)))
hist(y,breaks=-0.5+c(0:17),freq=F)
pois.hat = mean(y)
z = 0:16
for(i in 1:length(z))
  lines(c(z[i]+0.1,z[i]+0.1),c(0,dpois(z[i],pois.hat)),lwd=5,col=3)

legend("topright",c("Poisson"),lty=1,col=3)


#Log-likelihood function
loglik = function(alpha,beta,mu,lambda,x)
{
  l = 0
  for(i in 1:length(x$enc))
  {
    e = x$enc[i]
    n = x$freq[i]
    if(e==0)
      l = l + n*log(alpha+beta*exp(-mu)+(1-alpha-beta)*exp(-lambda))
    else
      l = l + n*log(beta*exp(-mu)*mu^e+
                    (1-alpha-beta)*exp(-lambda)*lambda^e)-log(gamma(e+1))
  }
  l
}

#-loglik for optimization, parameters in a vector and transformed
minusloglik = function(param,x)
{
  alpha = exp(param[1])/(1+exp(param[1]))    #Transforming to [0,1] interval
  beta = exp(param[2])/(1+exp(param[2]))     #Transforming to [0,1] interval
  mu = exp(param[3])                        #Transforming to positive value
  lambda = exp(param[4])                    #Transforming to positive value

  l = 0
  for(i in 1:length(x$enc))
  {
    e = x$enc[i]
    n = x$freq[i]
    if(e==0)
      l = l + n*log(alpha+beta*exp(-mu)+(1-alpha-beta)*exp(-lambda))
    else
      l = l + n*log(beta*exp(-mu)*mu^e+
                    (1-alpha-beta)*exp(-lambda)*lambda^e)-log(gamma(e+1))
  }
  show(c(alpha,beta,mu,lambda,l))
  -l
}
param.init = c(0,-1,0,1)
res = optim(param.init,minusloglik,x=x)

param = res$par
#Transforming back to original scale
alpha = exp(param[1])/(1+exp(param[1]))    #Transforming to [0,1] interval
beta = exp(param[2])/(1+exp(param[2]))     #Transforming to [0,1] interval
mu = exp(param[3])                        #Transforming to positive value
lambda = exp(param[4])                    #Transforming to positive value
print("Estimates from Newton's method")
show(c(alpha,beta,mu,lambda))

#Show results
hist(y,freq=F)
z = 0:16
prob = (beta*exp(-mu)*mu^z + 
       (1-alpha-beta)*exp(-lambda)*lambda^z)/gamma(z+1)
prob[1] = prob[1]+alpha
for(i in 1:length(z))
  lines(c(z[i]+0.5,z[i]+0.5),c(0,prob[i]),lwd=5,col=2)


#EM algorithm

alpha = 0.6
beta = 0.3
mu = 1
lambda = 10
i = 0:16
eps = 0.001
l = loglik(alpha,beta,mu,lambda,x)
more = TRUE
while(more)
{
 l.old = l
 pi = (beta*exp(-mu)*mu^i + (1-alpha-beta)*exp(-lambda)*lambda^i)
 pi[1] = pi[1]+alpha
 zstat0 = alpha/pi[1]
 tstat = beta*exp(-mu)*mu^i/pi
 pstat = (1-alpha-beta)*exp(-lambda)*lambda^i/pi
 alpha = x$freq[1]*zstat0/N
 beta = sum(x$freq*tstat)/N
 mu = sum(i*x$freq*tstat)/sum(x$freq*tstat)
 lambda = sum(i*x$freq*pstat)/sum(x$freq*pstat)
 param = c(log(alpha/(1-alpha)),log(beta/(1-beta)),log(mu),log(1-alpha-beta))
 l = loglik(alpha,beta,mu,lambda,x)
 more = abs(l-l.old)>eps
 show(c(alpha,beta,mu,lambda,l))
}
print("Estimates from the EM algorithm")
show(c(alpha,beta,mu,lambda))
#Not working properly:
#hist(y,freq=F)     #This gives a histogram where breaks do not match the categories


hist(y,breaks=-0.5+c(0:17),freq=F)
#hist(y,freq=F)
z = 0:16
prob = (beta*exp(-mu)*mu^z + (1-alpha-beta)*exp(-lambda)*lambda^z)/gamma(z+1)
prob[1] = prob[1]+alpha
for(i in 1:length(z))
  lines(c(z[i]-0.1,z[i]-0.1),c(0,prob[i]),lwd=5,col=2)

pois.hat = mean(y)
for(i in 1:length(z))
  lines(c(z[i]+0.1,z[i]+0.1),c(0,dpois(z[i],pois.hat)),lwd=5,col=3)

legend("topright",c("Poisson mixtures","Poisson"),lty=1,col=2:3)
