\name{WMGGD}
\alias{WMGGD}
\title{The Weighted Mixture Generalized Gamma Distribution}
\usage{
dWMGG(X,lambda,beta,alpha) 
pWMGG(X,lambda,beta,alpha)
qWMGG(lambda,beta,alpha) 
rWMGG(n,lambda,beta,alpha)
EM(X)
}
\arguments{
\item{X}{vector of data.}
\item{n}{a number of observations.}
\item{lambda}{a value of the parameter lambda, where lambda>0.}
\item{beta}{a value of the parameter beta, where beta>0.}
\item{alpha}{a value of the parameter alpha, where alpha>0.}
}
\value{
\item{dWMGG}{gives the density of a weighted mixture generalized gamma distribution.}
\item{pWMGG}{gives the distribution function of a weighted mixture generalized gamma distribution.}
\item{qWMGG}{gives the quantile function of a weighted mixture generalized gamma distribution.}
\item{rWMGG}{generates random numbers of a weighted mixture generalized gamma distribution.}
\item{EM}{EM gives parameter estimation function a weighted mixture generalized gamma distribution by using the expectation-maximization algorithm (EM).}
}
\description{
Density, distribution function, quantile function, random generation and parameter estimation function of the weighted mixture generalized gamma distribution as Phaphan and Abdullahi (2022) introduced.
}
\details{
Phaphan and Abdullahi (2022) introduced a weighted mixture generalized gamma distribution, which was a combination of generalized gamma distribution and length biased generalized gamma distribution. There were three parameters in this distribution: lambda, beta and alpha. 
}
\examples{
#install.packages("devtools")
library(devtools)
#install_github("wikanda-phaphan/WMGGD")
library(WMGGD)

? WMGGD
lambda <-0.7
beta <-1
alpha <-3
n <- 50
X<-rWMGG(n,lambda,beta,alpha)

dWMGG(X,lambda,beta,alpha)
pWMGG(X,lambda,beta,alpha)
qWMGG(lambda,beta,alpha)
EM(X)

}
\references{
I. Abdullahi, W. Phaphan. Some properties of the new mixture of generalized gamma distribution, Lobachevskii Journal of Mathematics, Vol.43, No.9, 2111–2121, 2022.
}
