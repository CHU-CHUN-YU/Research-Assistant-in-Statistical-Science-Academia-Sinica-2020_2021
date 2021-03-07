
###############################################################
# L=c(0,log(log(n)),sqrt(log(n)),log(n))
###############################################################

### g_function 
g_function <- function(x1,x2) {
  diff=abs(x1-x2)
  return(diff)
}

###  network data 
p=2
n=8
simulation_data<-function(L,n,p){
  
  num=L
  r_star=t(matrix(c(1,1.5),1,p))
  L=c(0,log(log(n)),sqrt(log(n)),log(n))
  
  
  ############ a_star (n-1)
  a_star=0
  for (i_2 in 1:(n-1)){
    a_star[i_2]<-(n-i_2)*L[num]/(n-1)
  }
  ############# b_star (n-1)
  b_star=0
  for (i_3 in 1:n-1){
    b_star[i_3]=a_star[i_3]
  }
  
  #############  X
  X=matrix(0,n,p)
  for (i_1 in 1:n){
    X[i_1,]=rbeta(p,2,2,ncp = 0)
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
  
  v1=v1[,-n] #drop the last column
  
  
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
  v2=v2[,-n] #drop the last column
  
  ##### W
  w=cbind(zij,v1,v2)
  
  ##### theta
  theta=matrix(c(r_star,t(a_star),t(b_star)),p+2*(n-1),1)
  
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
data=simulation_data(2,n,p)

### quadruple of distinct nodes 

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
quadruple_nodes=quadruple_nodes(n,data)





