fGG=function(X,lambda,beta,alpha)
{
  ((lambda*beta)/(gamma(alpha)))*(lambda*X)^((alpha*beta)-1)*
    (exp(-(lambda*X)^(beta)))
}
fLBGG=function(X,lambda,beta,alpha)
{
  (lambda*beta)/(gamma(alpha+(1/beta)))*((lambda*X)^(alpha*beta))*
    (exp(-(lambda*X)^(beta)))
}
#### The Weighted Mixture Generalized Gamma Distribution (WMGGD)
########### pdf
dWMGG=function(X,lambda,beta,alpha){
  (lambda/(lambda+1))*fGG(X,lambda,beta,alpha)+(1/(lambda+1))*fLBGG(X,lambda,beta,alpha)
}
#dWMGG(X,lambda,beta,alpha)
############ cdf 
install.packages('zipfR')
library(zipfR)
pWMGG=function(X,lambda,beta,alpha){
  1-((lambda/(lambda+1))*Igamma(alpha, (lambda*X)^beta)/gamma(alpha))-((1/(lambda+1))*Igamma(alpha+(1/beta),(lambda*X)^beta)/gamma(alpha+(1/beta)))
}
#pWMGG(X,lambda,beta,alpha)
######### WMGG-random numbers generation procedure ### 1 time
fM=function(X,lambda,beta,alpha){
  dWMGG(X,lambda,beta,alpha)/dgamma(X, shape=beta, scale = lambda)
}

rWMGG=function(n,lambda,beta,alpha){
  M <-optimise(f=function(X){fM(X,lambda,beta,alpha)},interval=c(0,5),maximum=T)$objective
  X=NULL
  while(length(X)<n){
    y=rgamma(n*M, shape=beta, scale = lambda)
    u=runif(n*M,0,M)
    X=c(X,y[u<fM(y,lambda,beta,alpha)])
  }
  X=X[1:n]
  MGG<-matrix(X, ncol=1)
  return(MGG)
}
#X<-rWMGG(n,lambda,beta,alpha)
############# EM Algorithm
EM=function(X){
  ###Real
  n=length(X)
  flcomplete1<-function(beta)
  {
    
    (beta[3]*beta[2]+1)*n*log(beta[1])-n*log(beta[1]+1)+(ak-1)*n*log(gamma(beta[3]+(1/beta[2])))-n*ak*log(gamma(beta[3]))+n*log(beta[2])+(beta[3]*beta[2]-1)*sum(X)-(beta[1]^beta[2])*sum(X)
    
  } 
  ####################### find initial value#################################
  Xbar <-apply(X,2,mean)
  alpha.old <- 0.5/(log(Xbar)-mean(log(X)))
  beta.old <- Xbar/alpha.old
  lambda.old <- Xbar/(apply(X,2,var))
  
  cond2=1
  cond3=1
  cond4=1
  while((cond2>0.01)&(cond3>0.01)&(cond4>0.01)){
    Ez=(lambda.old/(lambda.old+1))*fGG(X,lambda.old,beta.old,alpha.old)/((1/(lambda.old+1))*fLBGG(X,lambda.old,beta.old,alpha.old)+(lambda.old/(lambda.old+1))*fGG(X,lambda.old,beta.old,alpha.old))
    ak=mean(Ez)
    jib <-nlminb(c(lambda.old,beta.old,alpha.old), flcomplete1,scale = 100)
    lambda.new<-jib$par[1]
    beta.new<-jib$par[2]
    alpha.new<-jib$par[3]
  

    cond2=abs(lambda.new-lambda.old)
    cond3=abs(beta.new-beta.old)
    cond4=abs(alpha.new-alpha.old)
    
    
    lambda.old=lambda.new
    beta.old=beta.new
    alpha.old=alpha.new
  }
  return(list(lambda=lambda.new,beta=beta.new,alpha=alpha.new))
}
#EM(X)
###### quantile function
qWMGG=function(lambda,beta,alpha){
  u=runif(1,0,1)
  WMGG=(1/(lambda+1))*(Igamma.inv(alpha,u*gamma(alpha),lower=TRUE))^(1/beta)+((1/((lambda+1)*lambda))*(Igamma.inv(alpha+(1/beta),u*gamma(alpha+(1/beta)),lower=TRUE))^(1/beta))
  return(WMGG)
}
##qWMGG(lambda,beta,alpha)
