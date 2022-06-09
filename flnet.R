library(flexclust)
library(grpreg)
library(MASS)
library(refund)
library(ROCR)
library(fda)
library(igraph)

set.seed(1) #set.seed
use_samople<-sample(1:1e5,6,replace = FALSE)
n1<-300 # training 
n2<-300 # testing 
n<-n1+n2 
N <- 200  # number of time points
t <- seq(0,1,by=1/N) 
J <- 2  # number of eigenfunctions
rangeval<-c(0,1) #range value
mfunctions<-1 # number of functional predictors
true_number_of_p<-1 # number of important functional predictors

#####################    x(t)       ####################

phi_list<-matrix(1,nrow = J+1,ncol = N+1)
for(j in 1:J){
  phi_list[j+1,] <- sqrt(2)*cos(pi*j*t)
}
set.seed(use_samople[1])
setsamplez1<-sample(1:1e5,mfunctions,replace = FALSE)
Z1<-list()

for(mset in 1:mfunctions){
  set.seed(setsamplez1[mset])
  Z1[[mset]]<-mvrnorm(n=n,mu=rep(0,(J+1)),Sigma = diag(rep(c(1:(J+1))^(-2))))
  
}


X<-matrix(0,nrow = n,ncol = mfunctions*(N+1))
X_list<-list()
for(mset in 1:mfunctions){
  for(i in 1:(n)){
    X[i,((mset-1)*(N+1)+1):(mset*(N+1))]<-t(phi_list)%*%(Z1[[mset]][i,])
  }
  X_list[[mset]]<-X[,((mset-1)*(N+1)+1):(mset*(N+1))]
}

#####################    beta(t)       ####################

beta<-matrix(0,N+1,mfunctions)
set.seed(2020)
beta_vector<-mvrnorm(1,rep(1,true_number_of_p*(J+1)),diag(1,true_number_of_p*(J+1)))
for(mfunctions_i in 1:true_number_of_p){
  beta[,mfunctions_i]<-t(phi_list)%*%(beta_vector[((J+1)*(mfunctions_i-1)+1):((J+1)*mfunctions_i)])
  
}
true_beta<-beta
#####################    y       ####################

yita<-rep(0,n)
for (mset in 1:mfunctions) {
  
  for(i in 1:(n)){
    ina<-rep(0,N)
    for(j in 1:N){
      
      ina[j]<-0.5*(X_list[[mset]][i,j]%*%beta[j,mset]+(X_list[[mset]][i,j+1]%*%beta[j+1,mset]))
    }
    yita[i]<-yita[i]+sum(ina)*1/N
  }
  
}
set.seed(use_samople[2])
e<-rnorm(n,0,0.1)
alpha_s<-0.1
set.seed(use_samople[3])
true_group_structure<-matrix(0,n,3)
true_group_structure[1:(n/3),1]<- 1
true_group_structure[(n/3+1):(n/3*2),2]<- 1
true_group_structure[(n/3*2+1):(n/3*3),3]<- 1
alpha_true<-c(rnorm(n/3,-1,alpha_s^2),rnorm(n/3,0,alpha_s^2),rnorm(n/3,1,alpha_s^2))
all_y<-alpha_true+yita+e
training_set_id<-c(1:100,201:300,401:500)
testing_set_id<-c(101:200,301:400,501:600)
y<-all_y[training_set_id]


#####################      A_matrix       ####################

p_b<-0.02
set.seed(use_samople[4])
A_vector<-c(rbinom(n/3*(n/3-1)/2*3,1,0.2),rbinom(n/3*n/3*3,1,p_b))
A_matrix11<-matrix(0,n/3,n/3)
A_matrix11[which(upper.tri(A_matrix11))]<-A_vector[1:(n/3*(n/3-1)/2)] 
A_matrix22<-matrix(0,n/3,n/3)
A_matrix22[which(upper.tri(A_matrix22))]<-A_vector[(1+n/3*(n/3-1)/2):(n/3*(n/3-1)/2*2)]
A_matrix33<-matrix(0,n/3,n/3)
A_matrix33[which(upper.tri(A_matrix22))]<-A_vector[(1+n/3*(n/3-1)/2*2):(n/3*(n/3-1)/2*3)]
A_matrix12<-matrix(A_vector[(1+n/3*(n/3-1)/2*3):(n/3*(n/3-1)/2*3+n/3*n/3)],n/3,n/3)
A_matrix13<-matrix(A_vector[(1+n/3*(n/3-1)/2*3+n/3*n/3):(n/3*(n/3-1)/2*3+n/3*n/3*2)],n/3,n/3)
A_matrix23<-matrix(A_vector[(1+n/3*(n/3-1)/2*3+n/3*n/3*2):(n/3*(n/3-1)/2*3+n/3*n/3*3)],n/3,n/3)

A_matrix11[lower.tri(A_matrix11)] <- t(A_matrix11)[lower.tri(A_matrix11)]
A_matrix22[lower.tri(A_matrix22)] <- t(A_matrix22)[lower.tri(A_matrix22)]
A_matrix33[lower.tri(A_matrix33)] <- t(A_matrix33)[lower.tri(A_matrix33)]
A_matrix21<-t(A_matrix12)
A_matrix31<-t(A_matrix13)
A_matrix32<-t(A_matrix23)

A_matrix<-rbind(cbind(A_matrix11,A_matrix12,A_matrix13),cbind(A_matrix21,A_matrix22,A_matrix23),cbind(A_matrix31,A_matrix32,A_matrix33))

 


#####################      fpca       ####################


sore_list<-list()
group<-c()
pca_list<-list()

for(mset in 1:mfunctions){
  nowdata<-X_list[[mset]][training_set_id,]
  #  Set up a B-spline basis  
  rangeval<-c(0,1)
  knots  <- t
  mynorder <- 4
  nbasis <- length(knots) + mynorder - 2
  hgtbasis <- create.bspline.basis(rangeval, nbasis, mynorder, knots)
  Lfdobj <- 2
  lambda <- 10^(-5)   
  growfdPar <- fdPar(hgtbasis, Lfdobj, lambda)
  datafd <- smooth.basis(t, t(nowdata), growfdPar)$fd
  daytemppcaobj <- pca.fd(datafd, nharm=10, growfdPar)
  varsum<-cumsum(daytemppcaobj$varprop)
  nhar<-min(which(varsum>0.95))
  group<-c(group,rep(mset,nhar))
  daytemppcaobj <- pca.fd(datafd, nharm=nhar, growfdPar)
  sore_list[[mset]]<-daytemppcaobj$scores
  pca_list[[mset]]<-predict(daytemppcaobj$harmonics,seq(0,1,by=1/200)) 
}


group<-factor(x=group)
p_group<-group
Sore<-matrix(0,nrow = n1,ncol = length(group))
for(mset in 1:mfunctions){
  Sore[,which(group==mset)]<-sore_list[[mset]]
}
 
####################    L matrix, D matrix ############
 


L_matrix<-matrix(0,n1,n1)
D_matrix<-matrix(0,n1,n1)
for(D_i in 1:length(training_set_id)){
  D_matrix[ D_i,D_i]<-sum(abs(A_matrix[training_set_id[D_i],training_set_id]))
}
L_matrix=D_matrix-A_matrix [training_set_id,training_set_id]

####################  fit model ########################

gamma_alpha<-c(0.01)
lambda_alpha<-c( 0.01 )
x_tuta<-cbind(diag(dim(Sore)[1]),Sore)
M_matrix<- rbind(cbind(L_matrix+gamma_alpha*diag(rep(1,dim(Sore)[1])),matrix(0,dim(Sore)[1],dim(Sore)[2])),matrix(0,dim(Sore)[2],dim(Sore)[2]+dim(Sore)[1]))
coef<-solve(t(x_tuta)%*%x_tuta+lambda_alpha*M_matrix)%*%(t(x_tuta)%*%y)


####################  fit high dimensional model ########################

mynetwork<-function(Sore,y,p_group,gamma_alpha,lambda_alpha,lambda_beta,L_matrix,mytol,max_iter){
  n<-dim(Sore)[1]
  p<- length(names(table(p_group)))
  
  beta<-rep(0,dim(Sore)[2])
  alpha<-rep(0,n)
  beta_new<-beta
  alpha_new<-alpha
  
  for(update_i in 1:max_iter){
    
    # beta-update
    new_y<-y- alpha_new 
    
    afit<-grpreg(Sore,  new_y, p_group, penalty="grSCAD",lambda = lambda_beta,gamma=3.7)
    beta_new<-afit$beta[-1]
    
    
    # alpha-update
    new_y<-y-Sore%*%beta_new 
    
    new_L_matrix<-L_matrix+gamma_alpha*diag(1,n,n)
    alpha_new<-solve(diag(1,n,n)+lambda_alpha*new_L_matrix)%*%new_y
    
    
    
    ########### r
    r<-c(alpha_new-alpha,beta_new-beta)
    nowtol<-norm(r,type="2")
    nowtol
    if(nowtol<mytol) break
    
    beta<- beta_new
    alpha<-alpha_new
    
  }
  return(list(alpha=alpha,beta=beta,update_i=update_i,nowtol=nowtol))
  
  
}
end_fit_network_oracle2<-mynetwork(Sore,y,p_group,gamma_alpha=0.01,lambda_alpha=0.01,lambda_beta=0.01, L_matrix,mytol=1e-5,max_iter=1e5)

