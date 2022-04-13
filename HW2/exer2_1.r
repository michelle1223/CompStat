#Exer2.1
x = c(1.77,-0.23,2.76,3.80,3.47,56.75,-1.34,4.24,-2.44,3.29,3.71,-2.40,4.53,-0.07,-1.05,-13.87,-2.53,-1.75,0.27,43.21)
loglik = function(theta,x=x)
{
  sum(-log(pi)-log(1+(x-theta)^2))
 }
 theta = seq(-100,100,length=1000)
 lv = theta
 for(i in 1:length(theta))
  lv[i] = loglik(theta[i],x)
 plot(theta,lv,type="l")
 
 lscore = function(theta,x=x)
 {
 	2*sum((x-theta)/(1+(x-theta)^2))
 }
 
 l2score = function(theta,x=x)
 {
 	-sum((1-(x-theta)^2)/(1+(x-theta)^2))
 }
 
 #Newton-Raphson
 theta.start = c(-11,-1,0,1.5,4,4.7,7,8,38,55)
 maxit = 10000
 theta.sol = rep(NA,length(theta.start))
 numit.sol = rep(NA,length(theta.start))
 lval.sol = rep(NA,length(theta.start))
 for(k in 1:length(theta.start))
 {
 	theta = theta.start[k]
 	numit = 0
 	more = TRUE
 	while(more)
 	{
	 	l = loglik(theta,x=x)
 		l1 = lscore(theta,x=x)
 		l2 =	 l2score(theta,x=x)
 		theta = theta - l1/l2
 		numit = numit+1
 		if(abs(l1)<1e-8 | numit ==maxit)
 			more = FALSE
 	}
 	theta.sol[k] = theta
 	lval.sol[k] = loglik(theta,x=x)
 	numit.sol[k] = numit
 }
 points(theta.sol,lval.sol)
 
 
 #Bisection
 a=-1
 b=1
 more = TRUE
 while(more)
 {
 	theta = (a+b)/2
 	l1a = lscore(a,x=x)
 	l1b = lscore(b,x=x)
 	l1theta = lscore(theta,x=x)
 	if(l1a*l1theta<0)
 		b = theta
 	else
 		a=theta
 	theta = (a+b)/2
 	more = (b-a) > 1e-8
 }
 points(theta,loglik(theta,x=x),pch=2,col=3)
 
 #Fixed-point
 G = function(theta,x=x,alpha=alpha)
 {
 	alpha*lscore(theta,x=x)+1
 }
 alpha= 1 
 alpha = 0.64
 alpha = 0.25
 theta = -1
 more = TRUE
 numit = 0
 maxit=100
 while(more)
 {
 	l1 = lscore(theta,x=x)
 	theta = theta+alpha*l1
 	numit = numit+1
 	more = (abs(l1) > 1e-8) & (numit < maxit)
 }
points(theta,loglik(theta,x=x),pch=3,col=4)

#Secant method
theta.prev = -3
theta = -3
#theta.prev = -3
#theta = -3
 	numit = 0
 	more = TRUE
 	while(more)
 	{
	 	l = loglik(theta,x=x)
 		l1 = lscore(theta,x=x)
 		l1.prev = lscore(theta.prev,x=x)
 		l2 = (l1-l1.prev)/(theta-theta.prev)
 		theta = theta - l1/l2
 		numit = numit+1
 		if(abs(l1)<1e-8 | numit ==maxit)
 			more = FALSE
 	}
points(theta,loglik(theta,x=x),pch=4,col=5)

