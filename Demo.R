#======================================================================================
# program: demo.R
# note   : code submission for paper "code submission for paper "Integrative Sparse 
#                                                                Partial Least Squares".
# purpose: illustrate how to implement the proposed method 
#=========================================================================================

# step (1) source "iSPLS_Functions" file to load pre-defined functions
rm(list=ls()) 
source("iSPLS_Functions")

# step (2) generate datasets
data.sim <- function(n,p,q,l){
  cov <- matrix(1,nr=p,nc=p)
  for(i in 1:p){
    for(j in  1:i){
      cov[i,j]<- 0.7^(abs(i-j))
      cov[j,i] <- 0.7^(abs(i-j))
    }
  }
  
  x <- mvrnorm(n,rep(0,p),cov)
  
  beta1 <- c(rep(2,5),rep(3,5),rep(0,90))
  y <- matrix(0,nr=n,nc=q)
  for(j in 1:q){
    y[,j] <- x%*%beta1 + rnorm(n,0.1)
  }
  
  colnames(x) <- paste("Set",l,"x",c(1:p),sep=" ")
  colnames(y) <- paste("Set",l,"y",c(1:q),sep=" ")
  data.com <- data.frame(y,x)
  return(as.matrix(data.com))
}

## generate datasets
L <- 3  # the number of datasets    
n <- 20  # the number of observations in each dataset
p <- 100  # the number of explanatory variables in each dataset 
q <- 5  # the number of dependent variables in each dataset 

D1 <- NULL
for(l in 1:L){
  data <- data.sim(n=20,l,p=100,q=5)
  D1 <- cbind(D1, data)
}

y1 <- D1[,c(1:5,106:110,211:215)]
x1 <- D1[,-c(1:5,106:110,211:215)]

rslt <- ispls(x1, y1, L=3, mu1=0.05, mu3=0.5, kappa=0.05, type1="homo", type2="m",
                  scale.x = TRUE, scale.y = TRUE, maxstep = 50, trace = FALSE) 

