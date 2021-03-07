library(mvtnorm)
library(MASS)

## kernel function
k_f=function(u){
  return(dnorm(u,0,1))}

k_f02=function(u){
  return((3/4)*(1-u^2)*ifelse(abs(u)<=1,1,0))}


### weight function
w_f=function(X,X_new,p,h){
  W=0
  for(i in 1:n){
    W[i]=prod(k_f(matrix(X[i,]-X_new,p,1)/h))}
  return(W)
}

w_f_b=function(X,X_new,p,h,B){
  w_B=0
  for(i in 1:n){
    w_list=0
    for(t in 1:d){
      w_list[t]=t(B[,t])%*%matrix((X[i,]-X_new)/h,p,1)
    }
    w_B[i]=prod(k_f(w_list))}
  return(w_B)
}


#### my data
n=50
p=10
X=matrix(0,n,p)
for(a in 1:p){
  X[,a] <- rnorm(n,0,1)
}
e<-rnorm(n,0,0.1)
Y=(X[,1]*(X[,1]+X[,2]+1))+0.5*e 
h=6*(n^(-1/7))

####### B true

B_true=cbind(c(1,0,0,0,0,0,0,0,0,0),c(0,1,0,0,0,0,0,0,0,0))
d=dim(B_true)[2]


# B initail
B_ini=matrix(0,p,d)
### some param in while
tol=0.000001

set.seed(123)
sigma <-diag(1,p)

X_new=rmvnorm(1, mean=c(0,0,0,0,0,0,0,0,0,0), sigma=sigma)

### MAVE FUNTION

MAVE=function(X,Y,X_new,h,B_ini,d,tol,W_f,w_f_b){
  
  
  
  n=dim(X)[1]
  p=dim(X)[2]
  distance=1
  B=B_ini
  Iter=0
  
  while(distance>tol && Iter<1000){
    
    # choose k=1,2,...,d
    Iter=Iter+1
    if(Iter%%d==0){ k=d }else { k=Iter%%d }
    
    # C0、F0、D0、E0
    C0=F0=matrix(0,p,1)
    D0=matrix(0,p,p)
    E0=0
    
    # choose  w  with B ?
    if(Iter==1){W=w_f(X,X_new,p,h)} else{W=w_f_b(X,X_new,p,h,B)}
    
    for(i in 1:n){
      C0=C0+W[i]*matrix(X[i,]-X_new,p,1)
      D0=D0+W[i]*matrix(X[i,]-X_new,p,1)%*%t(matrix(X[i,]-X_new,p,1))
      E0=E0+W[i]*Y[i]
      F0=F0+(W[i]*matrix(X[i,]-X_new,p,1)*Y[i])
    }
    
    
    # M
    M=matrix(NaN,(d+1),(d+1))
    M[1,1]=1
    M[1,2:(d+1)]=t(C0)%*%B
    M[2:(d+1),1]=t(B)%*%C0
    M[2:(d+1),2:(d+1)]=t(B)%*%D0%*%B
    
    
    # H
    H=matrix(NaN,(d+1),1)
    H[1,1]=E0
    H[2:(d+1),1]=t(B)%*%F0
    
    # acde ,a0,d0,c0,e0
    acde=ginv(M)%*%H
    a0=acde[1]
    d0=acde[2:(d+1)][k]
    c0_e0=matrix(acde[2:(d+1)][-k],d-1,1)
    
    # A
    A=matrix(NaN,(p+d-1),(p+d-1))
    B_k=matrix(B[,-k],p,(d-1))
    A[1:p,1:p]=(d^2)*D0
    A[1:p,(p+1):(p+d-1)]=B_k
    A[(p+1):(p+d-1),1:p]=t(B_k)
    A[(p+d-1),(p+d-1)]=0
    
    
    A_plus= ginv(A)
    ## L
    L=matrix(NaN,(p+d-1),1)
    L[1:(p+d-1-1),1]=F0-(a0*C0)-(D0%*%B_k%*%c0_e0)
    L[p+d-1,1]=0
    
    ## S, b_hat
    S=A_plus%*%L
    b_hat=matrix(S[1:(p+d-1-1),1],(p+d-1-1),1)
    B_new=B
    B_new[,k]=b_hat
    
    # DISTANCE
    # prod(B_new)*prod(B)!=0
    
    if ( Iter>2 ){
      B_new_proj=B_new%*%solve(t(B_new)%*%B_new)%*%t(B_new)
      B_proj=B%*%solve(t(B)%*%B)%*%t(B)
      distance=sqrt(sum((B_new_proj-B_proj)^2))}else{ distance=1 }
    
    #print(B_new)
    #print(B)
    #print(distance)
    #print(Iter)
    
    B=B_new
    #print(B)
  }
  
  result=list(B_hat=B,Iter=Iter)
  return(result)
  
}
mave=MAVE(X,Y,X_new,h,B_ini,d,tol,W_f,w_f_b)

######################################################################

#### Simulstion 500 times ##################################

n=50
B_true=cbind(c(1,0,0,0,0,0,0,0,0,0),c(0,1,0,0,0,0,0,0,0,0))
d=dim(B_true)[2]
p=dim(B_true)[1]
h=6*(n^(-1/7))


set.seed(123)
sigma <-diag(1,p)
X_new=rmvnorm(1, mean=c(0,0,0,0,0,0,0,0,0,0), sigma=sigma)

# B initail
B_ini=B_true
  #matrix(0,p,d)


### some param. in while
tol=10^(-6)
dis1=0
dis2=0
dis=0
record=list()

for(i in 1:500){
  
  # data
  print(i)
  X=matrix(0,n,p)
  for(a in 1:p){
    X[,a] <- rnorm(n,0,1)
  }
  
  e<-rnorm(n,0,0.1)
  Y=(X[,1]*(X[,1]+X[,2]+1))+0.5*e 
  p=dim(X)[2]
  
  tryCatch({
    mave=MAVE(X,Y,X_new,h,B_ini,d,tol,W_f,w_f_b)
    B_hat=mave$B_hat
    record[[i]]=B_hat
    
    # DISTANCE
    B_hat_proj=B_hat%*%solve(t( B_hat)%*% B_hat)%*%t(B_hat)
    B_true_proj=B_true%*%solve(t(B_true)%*%B_true)%*%t(B_true)
    dis[i]=sqrt(sum((B_hat_proj-B_true_proj)^2))
    
    I=diag(1,p,p)
    B1=B_hat[,1]
    B2=B_hat[,2]
    dis1[i]=sqrt(sum(((I-B_true%*%solve(t(B_true)%*%B_true)%*%t(B_true))%*%B1)^2))
    dis2[i]=sqrt(sum(((I-B_true%*%solve(t(B_true)%*%B_true)%*%t(B_true))%*%B2)^2))
    
    
  },error = function(err) {
    print(paste("ERROR:  ",err)) 
    print(i)
  })
}

mean(dis,na.rm = T)
mean(dis1,na.rm = T)
mean(dis2,na.rm = T)





