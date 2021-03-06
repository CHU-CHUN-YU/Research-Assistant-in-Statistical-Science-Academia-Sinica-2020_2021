rm(list = ls())
########################## data02 ##############################
# r=(1,1.5)/ai=0/bi=0
# x=noraml(0,1)
# g_function (xi-xj)

g_function <- function(x1,x2) {
  diff=(x1-x2)
  return(diff)
}  ###  2 method
simulation_data<-function(L,n,p){
  
  num=L
  r_star=t(matrix(c(1,1.5),1,p))
  L=c(0,log(log(n)),sqrt(log(n)),log(n))
  
  
  ############ a_star (n-1)
  a_star=0
  for (i_2 in 1:n){
    a_star[i_2]=0
    
    #(n-i_2)*L[num]/(n-1)
  }
  ############# b_star (n-1)
  b_star=0
  for (i_3 in 1:n){
    b_star[i_3]=a_star[i_3]
  }
  
  #############  X
  X=matrix(0,n,p)
  for (i_1 in 1:p){
    X[,i_1]=rnorm(n,0,1)
  }
  
  ############# z
  zij=matrix(0,n*(n-1),p)
  for(k_1 in 1:p){
    z=matrix(NaN,n,n)
    for (i_4 in 1:n){
      for(j_1 in 1:n){
        if (i_4==j_1)
        {z[i_4,j_1]=NA}
        else
        {z[i_4,j_1]=g_function(X[i_4,k_1],X[j_1,k_1])}
      }
    }
    a1=as.vector(t(z))
    b1=a1[-which(is.na(a1))]
    zij[1:(n*(n-1)),k_1]=t(b1)
  }
  #zij
  
  ######### v1
  v1=matrix(0,n*(n-1),n)
  for(k_2 in 1:n){
    v=matrix(NaN,n,n)
    for (i_5 in 1:n){
      for(j_2 in 1:n){
        if (i_5==j_2)
        {v[i_5,j_2]=NA}
        else if(i_5==k_2)
        {v[i_5,j_2]=1}
        else
        {v[i_5,j_2]=0}
      }
    }
    a2=as.vector(t(v))
    b2=a2[-which(is.na(a2))]
    v1[1:(n*(n-1)),k_2]=t(b2)
  } 
  
  ######### v2
  v2=matrix(0,n*(n-1),n)
  for(k_3 in 1:n){
    v=matrix(NaN,n,n)
    for (i_6 in 1:n){
      for(j_3 in 1:n){
        if (i_6==j_3)
        {v[i_6,j_3]=NA}
        else if(j_3==k_3)
        {v[i_6,j_3]=1}
        else
        {v[i_6,j_3]=0}
      }
    }
    a3=as.vector(t(v))
    b3=a3[-which(is.na(a3))]
    v2[1:(n*(n-1)),k_3]=t(b3)
  }  
  
  ##### W
  w=cbind(zij,v1,v2)
  
  ##### theta
  theta=matrix(c(r_star,t(a_star),t(b_star)),p+2*(n),1)
  
  ### prob
  link_prob=exp(w%*%theta)/(1+exp(w%*%theta))
  
  
  #### convert to adj matrix
  value=link_prob
  m=diag(NaN,n,n)
  
  for (i in 1:n){
    for(j in 1:n){
      if (i==j)
      {m[i,j]=0}
      else
      {m[i,j]=value[1]}&&{value=value[-1]}
    }
  }
  
  adj=matrix(rbinom(n*n,1,as.vector(m)),n,n) #matrix
  
  
  ### vector
  
  adj_vector=diag(NaN,n,n)
  for(i in 1:n){
    for(j in 1:n){
      if (i!=j)
      {adj_vector[i,j]=adj[i,j]}
    }
  }
  
  adj_vector=as.vector(t(adj_vector))
  adj_vector=adj_vector[-which(is.na(adj_vector))]
  
  
  result=list(X=X,w=w,link_prob=link_prob,adj=adj,adj_vector=adj_vector,z=w[,1:p],theta=theta,v1=v1,v2=v2)
  return(result)
}


## quadruple_nodes
quadruple_nodes<-function(n,data){
  
  p=dim(data$X)[2]
  link=data$adj
  x <- seq(n)
  fir_2=t(combn(x,2))  # c n取2
  time=0
  det_list=list()
  z_det=0
  
  # 判斷所有組合的 z---z_det\det_list\times
  for (k in (1:dim(fir_2)[1])){
    a=fir_2[k,]
    x1=x[-a]
    sec_2=t(combn(x1,2)) 
    
    for (t in (1:dim(sec_2)[1])){
      b=sec_2[t,]
      det_matrix=matrix(0,2,2)
      
      for (i in 1:2){
        for(j in 1:2){
          det_matrix[i,j]=link[a[i],b[j]]
        }
      }
      z=((det_matrix[1,1]-det_matrix[1,2])-(det_matrix[2,1]-det_matrix[2,2]))/2
      if(z==1|z==-1){
        time=time+1
        det_list=append(det_list,list(cbind(matrix(a,2,1),matrix(b,2,1))))
        z_det[time]=z
      }
    }
  }
  
  ## z_list 
  z_list=list()
  for (i in 1:p){
    # z_list
    z=data$z[,i]
    z1=diag(0,n,n)
    for(i in 1:n){
      for ( j in 1:n){
        if(i==j){
          z1[i,j]=NaN
        } 
        else{
          z1[i,j]=z[1]
          z=z[-1]
        }
      }
    }
    z_list=append(z_list,list(z1))
  }
  
  ## 分別將各組合 z_i1j1, z_i2j2, z_i1j2, z_i2j1
  ## zij
  zij=matrix(0,length(det_list),4*p)
  for (k in (1:length(det_list))){
    a=matrix(det_list[[k]],2,2)[,1]
    b=matrix(det_list[[k]],2,2)[,2]
    
    for(i in 1:p){
      
      zp=as.matrix(z_list[[i]]) 
      zij[k,i]=zp[a[1],b[1]]
      zij[k,i+p]=zp[a[2],b[2]]
      zij[k,i+2*p]=zp[a[1],b[2]]
      zij[k,i+3*p]=zp[a[2],b[1]]
      
    }
  }
  
  result=list(time=time,det_list=det_list,z_det=z_det,z_list=z_list,zij=zij)
  return(result)
} 



## Sliced Average Variance Estimation ################
CSAVE<- function(x,y,k,p) {
  # 樣本平均
  sample_mean=colMeans(x)
  #樣本 變異數矩陣
  sample_var=var(x)
  
  # 標準化
  x=scale(x,center = TRUE,scale = TRUE)
  data=cbind(x,y)
  # y 由小至大排序，
  n=nrow(data)
  data_sort=data[order(data[,ncol(x)+1]),]
  
  # 1 and 0
  group_1=data[data[,ncol(data)]==1,]
  group_0=data[data[,ncol(data)]==0,]
  
  # sapmle_covariance^(-1/2)
  e=eigen(sample_var)  
  v=e$values
  q=e$vectors
  sample_var_sqrt=q%*%diag(1/sqrt(v))%*%solve(q)
  
  
  ##  f hat--p(y=1)
  f_hat=nrow(group_1)/n
  
  ## delta
  cov_1=var(group_1[,-(k*p+1)])
  cov_0=var(group_0[,-(k*p+1)])
  
  # I
  I= diag(k*p)
  #M
  M=(1-f_hat)*(I-cov_0)^2+(I-cov_1)^2
  
  #  find the eigenvalues and the eigenvectors for M 
  E=eigen(M)
  
  # K=1   the K largest eigenvectors 
  eta_hat=E$vectors[,1:k]
  
  # bhat
  b_hat=t(eta_hat)%*%sample_var_sqrt
  
  result=list(b_hat=b_hat,sample_var_sqrt=sample_var_sqrt,sample_var=sample_var,M=M,f_hat=f_hat)
  return(result)
}


####### Evaluation -Use projection matrix ##################
### distance
distance<-function(b_hat,b_true){
  
  B_HAT=b_hat%*%solve(t(b_hat)%*%b_hat)%*%t(b_hat)
  B_TRUE=b_true%*%solve(t(b_true)%*%b_true)%*%t(b_true)
  
  diff=B_HAT-B_TRUE
  d=norm(as.matrix(diff),'2')
  result=list(d=d,B_HAT=B_HAT,B_TRUE=B_TRUE)
  return(result)
}



