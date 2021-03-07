## kernel function
k_f=function(u){
  return(dnorm(u,0,1))}

## my data
n=500
X1=rnorm(n,mean=0,sd=2)
X2=rnorm(n,mean=0,sd=2)
X=cbind(X1,X2)
h=c(1,1)
e=rnorm(n,mean=0,sd=1)
Y=sin(X1)+cos(X2)+e  

## want to predict 
X_new=rbind(c(1,2.3))

#### Local ploynomial linear regression with  d>1
local_ploy_linear_reg=function(X_new,k_f,h,X,Y){
  
  n=dim(X)[1]
  p=dim(X)[2]
  n_new=dim(X_new)[1]
  r_x=0
  
  l_X_T=matrix(NaN,n_new,n)
  
  for(i in 1:n_new){
    
    #############################
    #   X_x   
    #   W_x  
    X_x=matrix(NaN,n,(p+1)) # X_x 
    X_x[,1]=1
    W_x=diag(NaN,n,n) # W_x 
    
    for(r in 1:n){
      v=0 
      for(c in 1:p){
        X_x[r,c+1]=X[r,c]-X_new[i,c]
        v[c]=k_f((X[r,c]-X_new[i,c])/h[c])
      }
    W_x[r,r]=prod(v)
     
    }  
    ##################################
    # 例外處理
    tryCatch({
      
      ### l(x).trans
      e_T=append(1,rep(0,p))
      l_X_T[i,]=e_T%*%solve(t(X_x)%*% W_x%*%X_x)%*%t(X_x)%*%W_x
      B_coeff=solve(t(X_x)%*% W_x%*%X_x)%*%t(X_x)%*%W_x%*%Y
      ### r(x)
      r_x[i]= l_X_T[i,]%*%Y
      
      ### L --smoothing matr
      
    },error = function(err) {
      print(paste("ERROR:  ",err)) 
      print(i)
    })
  }
  result=list(X=X,Y=Y,r_x=r_x,X_new=X_new,l_X_T=l_X_T,B_coeff=B_coeff)
  return(result)
}

## r(x)_hat
req=local_ploy_linear_reg(X_new,k_f,h,X,Y)

req$B_coeff


cos(X_new[1,1])
-sin(X_new[1,2])


################### gcv p=p n=300
n=500
r_n_hat_list=0
b1_hat_list=0
b2_hat_list=0

for(i in 1:500){
  X1=rnorm(n,mean=0,sd=2)
  X2=rnorm(n,mean=0,sd=2)
  X=cbind(X1,X2)
  e=rnorm(n,mean=0,sd=1)
  Y=sin(X1)+cos(X2)+e  
  
  
  Y_new=0
  X_new=rbind(c(1,1.3))
  h=c(1.5,1.5)
  
  

  req=req=local_ploy_linear_reg(X_new,k_f,h,X,Y)
  r_n_hat_list[i]=req$r_x
  b1_hat_list[i]=req$B_coeff[2,1]
  b2_hat_list[i]=req$B_coeff[3,1]
    
  print(i)
}

# b2
# true b2= -sin(X_new[2])

(mean(b2_hat_list+sin(X_new[2])))^2 #bias^2
var(b2_hat_list,na.rm = TRUE) #variance
mean((b2_hat_list+sin(X_new[2]))^2) #mse



# b1
# true b1= cos(X_new[1])

(mean(b1_hat_list-cos(X_new[1])))^2 #bias^2
var(b1_hat_list,na.rm = TRUE) #variance
mean((b1_hat_list-cos(X_new[1]))^2) #mse


