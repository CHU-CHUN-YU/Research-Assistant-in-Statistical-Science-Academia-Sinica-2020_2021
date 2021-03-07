install.packages ('mvtnorm')
library(mvtnorm)

#################### generate x dataset #######################

x_data01<-function(n){
sigma <-diag(1,5)
x <- rmvnorm(n, mean=c(0,0,0,0,0), sigma=sigma)
e<-rnorm(n,0,1)
result=list(x=x,e=e)
return(result)
} #N(0,1)
x_data02<-function(n){
  x=matrix(0,n,5)
  for (i in 1:5){
    x[,i]=runif(n, min = 0, max = 1)
  }
  e<-runif(n, min = 0, max = 1)
  result=list(x=x,e=e)
  return(result)
} #U(0,1)


################ y=g(x) ###################################

g_1<-function(x,e){
  y=((5+x[,1]+x[,2]+x[,3])^2)+e  
  return(y)
}

g_2<-function(x,e){
  
  y=(x[,1])/(0.5+(x[,2]+1.5)^2)+e
  return(y)
}


############ Sliced inverse regerssion #############################
##################### sir_function(x,y,H,k) ####################


sir_function <- function(x,y,H,k) {
  # 樣本平均
  sample_mean=colMeans(x)
  #樣本 變異數矩陣
  sample_var=var(x)
  # 標準化
  x=scale(x,center = TRUE,scale = TRUE)
  data=cbind(x,y)
   # y 由小至大排序，切H=10分
  n=nrow(data)
  data_sort=data[order(data[,ncol(x)+1]),]
  data_split=as.data.frame(cbind(data_sort,rep(1:H,each=nrow(data_sort)/H)))
  group=split(data_split,data_split[,ncol(data_sort)+1])
  
  # m_hat p_hat
  mean_matrix=matrix(0,H,ncol(x))
  p_hat=matrix(0,H,1)
  for (i in 1:H){
    g=as.data.frame(group[i])
    mean_matrix[i,]=colMeans(g[,1:ncol(x)])
    p_hat[i,]=nrow(g)/n
  }
  
#  weighted covariance matrix--- v hat
  
  V_hat=0
  for (i in 1:H){
    
    s=p_hat[i]*mean_matrix[i,]%*%t(mean_matrix[i,])
    V_hat=V_hat+s
  }

  #  find the eigenvalues and the eigenvectors for V 
  E=eigen(V_hat, only.values = FALSE, EISPACK = FALSE)

  # K=1   the K largest eigenvectors 
  eta_hat=E$vectors[,1:k]

  # sapmle_covariance^(-1/2)
  e=eigen(sample_var)  
  v=e$values
  q=e$vectors
  sample_var_sqrt=q%*%diag(1/sqrt(v))%*%t(q)
  
  # bhat
  b_hat=t(eta_hat)%*%sample_var_sqrt
  
  # x_new
  x_new=t(b_hat%*%t(x))
  
  result=list(b_hat=b_hat,eta_hat=eta_hat,E=E)
  return(result)
}


####### Evaluation -Use projection matrix ##################


distance<-function(b_hat,b_true,k){
  b_true=matrix(b_true,5,k)
  b_hat=matrix(b_hat,5,k)
  B_HAT=b_hat%*%solve(t(b_hat)%*%b_hat)%*%t(b_hat)
  B_TRUE=b_true%*%solve(t(b_true)%*%b_true)%*%t(b_true)
  diff=B_HAT-B_TRUE
  d=norm(as.matrix(diff),'2')
  return(d)
}





