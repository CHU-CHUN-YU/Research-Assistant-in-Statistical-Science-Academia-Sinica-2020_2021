

##  simulation studies



# g_function #############################################

g_function <- function(x1,x2) {
diff=abs(x1-x2)
return(diff)
}
###############################################################
# L=c(0,log(log(n)),sqrt(log(n)),log(n))

###############################################################

simulation_data<-function(num,n,p){
r_star=t(matrix(c(1,1.5),1,p))
L=c(0,log(log(n)),sqrt(log(n)),log(n))
      
#############  X
X=matrix(0,n,p)
for (i_1 in 1:n){
  X[i_1,]=rbeta(p,2,2,ncp = 0)
}
############ a_star
a_star=0
for (i_2 in 1:n){
  a_star[i_2]<-(n-i_2)*L[num]/(n-1)
}
############# b_star
b_star=0
for (i_3 in 1:n-1){
  b_star[i_3]=a_star[i_3]
}
b_star[n]=0

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
#v1

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
#v2

##### W
w=cbind(zij,v1,v2)
      
##### theta
theta=matrix(c(r_star,t(a_star),t(b_star)),p+2*n,1)

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

adj=matrix(rbinom(n*n,1,as.vector(m)),n,n)

result=list(X=X,w=w,link_prob=link_prob,adj=adj,z=w[,1:2],theta=theta,v1=v1,v2=v2)
return(result)
}

#######################################################################
set.seed(100)
test1=simulation_data(num=2,n=3,p=2)
test1$adj
test1$X
test1$link_prob
test1$z
test1$theta
test1$v1

############################### parameter Estimation ####################################


#r_hat=c(1,1.5)
#a_hat=c(1,2,2)
#b_hat=c(1,2,2)

#sum(z_r*aij)
#sum(a_hat*d)
#sum(b_hat*b)
#sum(log(1+exp(w%*%theta)))

##########################################
#  log likelihood function()

data=test1
#theta=c(r_hat,a_hat,b_hat)

ln=function(theta,data){
  
  r_hat=theta[1:2]
  a_hat=theta[3:5]
  b_hat=theta[6:8]
  
#######one

    attach(data)
 
  aij=matrix(NaN,3,3)
  for (i in 1:3){
    for(j in 1:3){
      if (i==j)
      {aij[i,j]=NaN}
      else
      {aij[i,j]=adj[i,j]}
    }
  }
  aij<-na.omit(as.vector(t(aij)))
  z_r=0
  for (i in 1:6){
    z_r[i]=r_hat%*%z[i,] 
  }
### two & three
  d=0
  b=0
  for(i in 1:3){
    d[i]=sum(adj[i,])
    b[i]=sum(adj[,i])
  }
####four
  w=cbind(z,v1,v2)
  return(sum(z_r*aij)+sum(a_hat*d)+sum(b_hat*b)+sum(log(1+exp(theta %*% t(w)))))
}
#############################################################################
#r_hat=t(matrix(c(1,1.5),1,2))
#a_hat=c(1,2,2)
#b_hat=c(1,2,2)
ln(c(1,1.5,1,2,2,1,2,2),test1)







############################## Logistic Regression &
#############  Maximum Likelihood Estimation by Newton-Raphson Method ##############
############# ¤j¼Ë¥»©Ê½è #######################
## µLºI¶Z¶µ

n=6
data=simulation_data(num=2,n=n,p=2)

attach(data)

aij=matrix(NaN,n,n)
for (i in 1:n){
  for(j in 1:n){
    if (i==j)
    {aij[i,j]=NaN}
    else
    {aij[i,j]=adj[i,j]}
  }
}
aij<-na.omit(as.vector(t(aij)))


x1=data$z[,1]
x2=data$z[,2]
y=aij
data=data.frame(y,x1,x2)

attach(data)
beta=c(0,0)
beta_new=0
sum=0
e=1    #¥ýÀH«K³]¤@­Ó­È
delta=0.0001



while(e>0){
  #«ØH() u()-->u()«Ø¦¨function¤ñ¸û¤è«K
  ln=function(beta){  ##beta­n¤@°_¼g¦¨¦V¶q ÁÙ¬O¤À¶}b0 b1 ³£¥i
    #page9 ¬õÂI 
    for (i in 1:(n*(n-1))){
      sum[i]=y[i]*(beta%*%c(x1[i],x2[i]))-log(1+exp(beta%*%c(x1[i],x2[i]))) 
    }
    return(sum(sum))
  }
  
  ln0=function(beta){
    return((1/(2*delta))*(ln(beta+c(delta,0))-ln(beta-c(delta,0))))
  }
  
  ln1=function(beta){
    return((1/(2*delta))*(ln(beta+c(0,delta))-ln(beta-c(0,delta))))
  }
  
  ln00=function(beta){
    return((1/(2*delta))*(ln0(beta+c(delta,0))-ln0(beta-c(delta,0))))
  }
  
  ln01=function(beta){
    return((1/(2*delta))*(ln0(beta+c(0,delta))-ln0(beta-c(0,delta))))
  }
  
  
  ln10=function(beta){
    return((1/(2*delta))*(ln1(beta+c(delta,0))-ln1(beta-c(delta,0))))
  }
  ln11=function(beta){
    return((1/(2*delta))*(ln1(beta+c(0,delta))-ln1(beta-c(0,delta))))
  }
  
  u=function(beta){
    return(cbind(ln0(beta),ln1(beta)))
  }
  
  H=function(beta){
    return(cbind(c(ln00(beta),ln10(beta)),c(ln01(beta),ln11(beta))))
  }
  beta_new=beta-u(beta)%*%solve(H(beta)) #%*%¯x°}ªº¬Û­¼
  e=ln(beta_new)-ln(beta)
  beta=beta_new #²{¦b°_©l¦ì¸m´«¦¨·sªº¦ì¸m
  
}
beta_mle=beta_new




#  glm  model
model=glm(y~0+x1+x2, family="binomial", data=data ) #³Ì¤p¥­¤èªk

#### compare ####
#model$coefficients

print(n)
print(beta_mle)
print(model$coefficients)


###################################################################################3

################################# plot simulational data   ####################################

n=dim(test1$adj)[1]
X=cbind(matrix(1:n,n,1),test1$X)
r1<-graph.adjacency(test1$adj==1)

#Convert igraph object to a data frame in R
net1.edges <- as.data.frame(get.edgelist(r1))

# Converting the data to an igraph object:
net1 <- graph_from_data_frame(net1.edges,X, directed=T) 

# The labels are currently node IDs.
# Setting them to NA will render no labels:
V(net1)$label.color <- "black"

#change arrow size and edge color:
E(net1)$arrow.size <- .2
E(net1)$edge.color <- "gray80"


# We can also add a legend explaining the meaning of the colors we used:
#par(mfrow = c(2,1))


set.seed(12)
V(net1)$size=X[,2]*20 
plot(net1,vertex.label.cex=0.8) 



set.seed(12)
V(net1)$size <-X[,3]*20
plot(net1,vertex.label.cex=0.8) 



################################### plot    ##############################################
#


install.packages("igraph")
library(igraph)

set.seed(123)

traits <- read.table("c:\\Users\\user\\Desktop\\RA\\data\\ELattr.dat",header = F)
names(traits) <- c('seniority','status','gender','office','years','age','practice','law school')
traits[,1]<- paste('V',1:71,sep = '')

relation_1<-read.table("c:\\Users\\user\\Desktop\\RA\\data\\ELadv.dat",header = F)
#relation_2<-read.table("c:\\Users\\user\\Desktop\\RA\\data\\ELfriend.dat",header = F)
#relation_3<-read.table("c:\\Users\\user\\Desktop\\RA\\data\\ELwork.dat",header = F)


r1<-graph.adjacency(relation_1==1)

#class(r1)
#r2<-graph.adjacency(relation_2==1)
#r3<-graph.adjacency(relation_3==1)


#Convert igraph object to a data frame in R
net1.edges <- as.data.frame(get.edgelist(r1))


# Converting the data to an igraph object:
net1 <- graph_from_data_frame(net1.edges ,traits, directed=T) 


# Removing loops from the graph:#?ˆª?™¤?‡ªå·±è?Ÿè‡ªå·?
net1<- simplify(net1, remove.multiple = F, remove.loops = T) 

# Let's and reduce the arrow size and remove the labels:
#plot(net1, edge.arrow.size=.4,vertex.label=NA,edge.curved=.1)

# Generate colors base on media type:
colrs <- c("tomato", "dodgerblue")
V(net1)$color <- colrs[V(net1)$status]


# Compute node degree (#links) and use it to set node size:
# mode= all/in/out
#V(net1)$size <- strength(net1,mode='in')

# taking the log to improve it
V(net1)$size <- log(strength(net1,mode='in')) * 4 
# The labels are currently node IDs.
# Setting them to NA will render no labels:
V(net1)$label.color <- "black"

#V(net1)$label <- NA

#change arrow size and edge color:
E(net1)$arrow.size <- .2

E(net1)$edge.color <- "gray80"


V(net1)$label <- ifelse( strength(net1)>=40, V(net1)$name, NA )
# We can also add a legend explaining the meaning of the colors we used:
par(mar = c(0, 0, 0, 0))

set.seed(12)
plot(net1,vertex.label.cex=0.8) 
legend('right', c("partner","associate"), pch=21, pt.bg=colrs, pt.cex=2.5, bty="n", ncol=1)

#####################################################################
#seniority
#status (1=partner; 2=associate)
#gender (1=man; 2=woman)
#office (1=Boston; 2=Hartford; 3=Providence)
#years with the firm
#age
#practice (1=litigation; 2=corporate)
#law school (1: harvard, yale; 2: ucon; 3: other)













