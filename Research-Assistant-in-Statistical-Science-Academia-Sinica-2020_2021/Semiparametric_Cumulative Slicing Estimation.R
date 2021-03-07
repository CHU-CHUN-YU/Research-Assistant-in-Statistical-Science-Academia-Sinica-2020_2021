install.packages ('mvtnorm')

library(mvtnorm)

#################### generate x  #######################

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
g_3<-function(x,e){
  
  y=x[,1]*(x[,1]+x[,2])+e
}

############## Cumulative Slicing Estimation #######################
###############################  CUME_function(x,y,H,k) #############
# H=n
CUME_function <- function(x,y,H,k) {
  
  
  # 樣本平均
  sample_mean=colMeans(x)
  #樣本 變異數矩陣
  sample_var=var(x)
  # 標準化
  x=scale(x,center = TRUE,scale = TRUE)
  
  data=cbind(x,y)
  
  # y 由小至大排序，選擇切點yi
  n=nrow(data)
  data_sort=data[order(data[,ncol(x)+1]),]
  data_sort=as.data.frame(data_sort)
  
  # m_hat / mean matrix
  mean_matrix=matrix(0,nrow(x),ncol(x))
  
  for (i in 1:H){
    g=data_sort[1:i,1:ncol(x)]
    mean_matrix[i,]=colMeans(g)
  }
  
  # weigh function
  #  weight_function<-function(yi){
  #   value=1
  #   return(value) 
  #  }
  
  #  weighted covariance matrix--- M hat
  
  M_hat=0
  for (i in 1:H){
    s=mean_matrix[i,]%*%t(mean_matrix[i,])*1
    M_hat=M_hat+s
  }
  
  #  find the eigenvalues and the eigenvectors for V 
  E=eigen(M_hat, only.values = FALSE, EISPACK = FALSE)
  
  # K=1   the K largest eigenvectors 
  
  eta_hat=E$vectors[,1:k]
  
  #sapmle_covariance^(-1/2)
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



