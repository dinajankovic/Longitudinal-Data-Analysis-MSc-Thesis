install.packages("rootSolve")
library(rootSolve)

install.packages("tibble")
library(tibble)


dataset = read.csv("ICHS.csv", header = TRUE)
attach(dataset)
head(dataset)

n = 250
m = 6
p = 4

cov2 = Gender
cov3 = VitaminA
cov4 = Age

y = Response

#We also need to introduce the intercept parameter beta_0
intercept = c(rep(1,len = n*m))
cov1 = intercept

dataset = add_column(dataset, cov1, .after = "Time")

set.seed(31417)   #this will always provide us with the same result
#in our case, the missingness percentage is 3.33% for seed = 31417
r = rbinom(n*m,1,0.95)
#Compute missingness percentage
r1 = r[r<1]
missing = length(r1)/(n*m)
paste("missingness percentage: ",round(missing*100, digits=2),"%")

#Introduce missing responses
for (i in (1:n*m)){
  if(r[i]<1) y[i] = NA
}

#Replace by 0 the missing responses
for (i in (1:n*m)){
  if(r[i]<1) y[i] = 0
}


mu = function(x){
  return(1/(1+exp(-x)))
}

mu_prime = function(x){
  return(exp(x)/((1+exp(x))^2))
}

#FINDING gamma_hat
logisticfunction_gamma = function(gamma){

  prod = cov1*gamma[1] + cov2*gamma[2] + cov3*gamma[3] + cov4*gamma[4];
  mean = 1/(1+exp(-prod));
  residual = r-mean;
  c(F1 = sum(cov1*residual), F2 = sum(cov2*residual), F3 = sum(cov3*residual), F4 = sum(cov4*residual)) }

gamma_hat = multiroot(f = logisticfunction_gamma, start = c(0,0,0,0))$root
gamma_hat


#compute estimated pi_{ij} using logit(pi_{ij})=gamma_hat^T X_{ij}
pi_hat = 1/(1+exp(-(gamma_hat[1]*cov1 + gamma_hat[2]*cov2 + gamma_hat[3]*cov3 + gamma_hat[4]*cov4)))

#compute y_{ij}^* (inverse probability weighted responses)
w = y*r/pi_hat

#compute estimated beta in the case of working indepedence
logisticfunction0 = function(b){

  prod = cov1*b[1] + cov2*b[2] + cov3*b[3] + cov4*b[4];
  mean = mu(prod);
  residual = w-mean;
  c(F1 = sum(cov1*residual), F2 = sum(cov2*residual), F3 = sum(cov3*residual), F4 = sum(cov4*residual)) }

beta_indep = multiroot(f = logisticfunction0, start = c(0,0,0,0))$root
beta_indep


#Create a list of 1500 matrices of dimension 6x3; only matrices List1[[k]] with k(mod6)=1 are defined
List1 = list()
for (k in 1:(n*m)) {
  if ((k %% m) == 1) {
    List1[[k]] = data.matrix(dataset[k:(k+m-1),4:7])
  }
}

#Extract matrices X_1,...,X_250 from List1 by mapping the indices from 1-250 to 1-1500
X = list()
for (i in 1:n){
  X[[i]] = List1[[m*i-m+1]]
}


#Create a list of 1500 6X1 matrices from y: only matrices y_list1[[k]] with k(mod6)=1 are defined
y_list <- list()
for (k in 1:(n*m)){
  if ((k %% m) == 1) {
    y_list[[k]] = data.matrix(y[k:(k+m-1)])
  }
}

#Create a list of 1500 6X1 matrices from r: only matrices r_list1[[k]] with k(mod6)=1 are defined
r_list = list()
for (k in 1:(n*m)) {
  if ((k %% m) == 1) {
    r_list[[k]] = data.matrix(r[k:(k+m-1)])
  }
}

#Create a list of 1500 6X1 matrices from w: only matrices w_list1[[k]] with k(mod6)=1 are defined
w_list = list()
for (k in 1:(n*m)) {
  if ((k %% m) == 1) {
    w_list[[k]] = data.matrix(w[k:(k+m-1)])
  }
}

#Extract matrices y_star1,...,y_star250 from w_list by mapping the indices from 1-250 to 1-1500
#y_star[[i]] is a 6X1 matrix with values y_{ij}^*,j=1,...,6
y_star = list()
for (i in 1:n) {
    y_star[[i]] = w_list[[m*i-m+1]]
}


#Create a list of 1500 6X1 matrices from pi_hat: only matrices pi_list1[[k]] with k(mod6)=1 are defined
pi_list = list()
for (k in 1:(n*m)) {
  if ((k %% m) == 1) {
    pi_list[[k]] = data.matrix(pi_hat[k:(k+m-1)])
  }
}


#Extract matrices pi1,...,pi250 from pi_list by mapping the indices from 1-250 to 1-1500
#pi[i]] is a 6X1 matrix with values pi_{ij}, j=1,...,6
pi = list()
for (i in 1:n) {
    pi[[i]] = pi_list[[6*i-5]]
}

#Compute inverse probabilities
#q[[i]] is a 6X1 matrix with values 1/pi_{ij}, j=1,...,6
q = list()
for (i in 1:n) {
   q[[i]] = 1/pi[[i]]-1
}


#Create a list of 250 6x1 matrices containing centered values when beta=beta_indep
#y_star_c[[i]] is a 6X1 matrix with values Y_{ij}^*-mu(X_{ij}^T beta_indep), j=1,...,6
y_star_c = list()
for (i in 1:n){
  y_star_c[[i]] = matrix(rep(0,len=m),nrow=m)
}


for (i in 1:n){
  for (j in 1:m){
          y_star_c[[i]][j,] = y_star[[i]][j,] - mu(t(X[[i]][j,]) %*% beta_indep)
  }
}

#Create a list of 250 6x1 matrices containing the variance of Y_{ij}^* when beta=beta_indep
#sigma[[i]]] is a 6X1 matrix with values sigma_{ij}^*(beta_indep), j=1,...,6
sigma = list()
for (i in 1:n){
  sigma[[i]] = matrix(rep(0,len = m),nrow = m)
}


for (i in 1:n){
  for (j in 1:m){
          sigma[[i]][j,] = mu_prime(t(X[[i]][j,]) %*% beta_indep)  +
                 (1/pi_list[[6*i-5]][j,]-1)*(mu_prime(t(X[[i]][j,]) %*% beta_indep) +
                 (mu(t(X[[i]][j,]) %*% beta_indep))^2  )
  }
}


#Create a list of 250 6x1 matrices containing standardized Y_{ij}^* when beta=beta_indep
#y_tilde[[i]]] is a 6X1 matrix with values (Y_{ij}^*-mu(X_{ij}^T beta_indep))/sqrt{sigma_{ij}^*(beta_indep)}, j=1,...,6
y_tilde = list()
for (i in 1:n){
  y_tilde[[i]] = matrix(rep(0,len=m),nrow=m)
}


for (i in 1:n){
  for (j in 1:m){
      y_tilde[[i]][j,] = y_star_c[[i]][j,] / sqrt( sigma[[i]][j,] )
  }
}

#Compute the estimated correlation matrix

#Case 1: 1-dependent
R_hat1 = matrix(0,m,m)
R_hat1[m,m] = 1
for (j in 1:(m-1)){
  R_hat1[j,j] = 1
  for (i in 1:n){
    R_hat1[j,j+1] = R_hat1[j,j+1] + (1/(n-p))*y_tilde[[i]][j]*y_tilde[[i]][j+1]
    R_hat1[j+1,j] =  R_hat1[j,j+1]
  }
}

#Find the inverse of R_hat1
R_hat1_inv = solve(R_hat1)
R_hat1_inv


#Case 2: exhangeable
N = n*m*(m-1)/2

alpha = 0
for (i in 1:n){
  for (k in 2:m){
    for (j in 1:(k-1)){
         alpha = alpha + (1/(N-p))*y_tilde[[i]][j]*y_tilde[[i]][k]
    }
  }
}


R_hat2 = matrix(0,m,m)
for (j in 1:m){
      R_hat2[j,j] = 1
}
for (j in 1:(m-1)){
      for (k in (j+1):m){
       R_hat2[j,k] = alpha
      }
}
for (j in 2:m){
      for (k in 1:(j-1)){
       R_hat2[j,k] = alpha
      }
}

#Find the inverse of R_hat2
R_hat2_inv = solve(R_hat2)
R_hat2_inv



#### OUR METHOD #####
#Case 1: 1-dependent (R = R_hat1)
logisticfunction1 = function(b){
  F = c(0,0,0,0)
  b = c(b[1],b[2],b[3],b[4])

  for (t in 1:p) {
    for (i in 1:n){
      for (k in 1:m){
        for (j in 1:m){


          F[t] = F[t] + X[[i]][j,t]*mu_prime(t(X[[i]][j,]) %*% b)   /  sqrt(   mu_prime(t(X[[i]][j,]) %*% b)  +
                 (1/pi_list[[m*i-m+1]][j,]-1)*(mu_prime(t(X[[i]][j,]) %*% b) +  (mu(t(X[[i]][j,]) %*% b))^2  ))*R_hat1_inv[j,k]*(y_star[[i]][k]-mu(t(X[[i]][j,]) %*% b))/sqrt(   mu_prime(t(X[[i]][j,]) %*% b)  +
                  (1/pi_list[[m*i-m+1]][k,]-1)*(mu_prime(t(X[[i]][j,]) %*% b) + (mu(t(X[[i]][j,]) %*% b))^2  ))

        }

      }
    }
  }

  c(F[1],F[2],F[3],F[4])


} #end of logisticfunction1

bhat1 = multiroot(f = logisticfunction1, start = c(beta_indep[1],beta_indep[2],beta_indep[3],beta_indep[4]))$root
bhat1


#exhangeable case: R = R_hat2
logisticfunction2 = function(b){
  F = c(0,0,0,0)
  b = c(b[1],b[2],b[3],b[4])

  for (t in 1:p) {
    for (i in 1:n){
      for (k in 1:m){
        for (j in 1:m){


          F[t] = F[t] + X[[i]][j,t]*mu_prime(t(X[[i]][j,]) %*% b)   /  sqrt(   mu_prime(t(X[[i]][j,]) %*% b)  +
                 (1/pi_list[[m*i-m+1]][j,]-1)*(mu_prime(t(X[[i]][j,]) %*% b) +  (mu(t(X[[i]][j,]) %*% b))^2  ))*R_hat2_inv[j,k]*(y_star[[i]][k]-mu(t(X[[i]][j,]) %*% b))/sqrt(   mu_prime(t(X[[i]][j,]) %*% b)  +
                  (1/pi_list[[m*i-m+1]][k,]-1)*(mu_prime(t(X[[i]][j,]) %*% b) + (mu(t(X[[i]][j,]) %*% b))^2  ))

        }

      }
    }
  }

  c(F[1],F[2],F[3],F[4])


} #end of logisticfunction2
bhat2 = multiroot(f = logisticfunction2, start = c(beta_indep[1],beta_indep[2],beta_indep[3],beta_indep[4]))$root
bhat2


A1_list = list()
for (i in 1:n){
  A1_list[[i]] = matrix(rep(0,len=m*m),nrow=m)
  for (j in 1:m){
    A1_list[[i]][j,j] = mu_prime(t(X[[i]][j,]) %*%  bhat1 )
  }
}

A2_list = list()
for (i in 1:n){
  A2_list[[i]] = matrix(rep(0,len=m*m),nrow=m)
  for (j in 1:m){
    A2_list[[i]][j,j] = mu_prime(t(X[[i]][j,]) %*%  bhat2 )
  }
}


D1_list = list()
for (i in 1:n){
  for (j in 1:m){
    D1_list[[i]] = A1_list[[i]]%*%X[[i]]
  }
}

D2_list = list()
for (i in 1:n){
  for (j in 1:m){
    D2_list[[i]] = A2_list[[i]]%*%X[[i]]
  }
}


prod1_list = list()
for (i in 1:n){
  prod1_list[[i]] = matrix(rep(0,len=m),nrow=m)

  for (j in 1:m){
    prod1_list[[i]][j,] = 0
    prod1_list[[i]][j,] = t(X[[i]][j,]) %*%  bhat1
    }
}

prod2_list = list()
for (i in 1:n){
  prod2_list[[i]] = matrix(rep(0,len=m),nrow=m)

  for (j in 1:m){
    prod2_list[[i]][j,] = 0
    prod2_list[[i]][j,] = t(X[[i]][j,]) %*%  bhat2
    }
}


mu1_list = list()
for (i in 1:n){
   mu1_list[[i]] =  mu(prod1_list[[i]])
}

mu2_list = list()
for (i in 1:n){
   mu2_list[[i]] =  mu(prod2_list[[i]])
}

mu1_prime_list = list()
for (i in 1:n){
    mu1_prime_list[[i]] = matrix(rep(0,len=m*m),nrow=m)
    mu1_prime_list[[i]] = mu_prime(prod1_list[[i]])
}

mu2_prime_list = list()
for (i in 1:n){
    mu2_prime_list[[i]] = matrix(rep(0,len=m*m),nrow=m)
    mu2_prime_list[[i]] = mu_prime(prod2_list[[i]])
}


A1_star_list = list()
for (i in 1:n){
  A1_star_list[[i]] = matrix(rep(0,len=m*m),nrow=m)
  for (j in 1:m){

  A1_star_list[[i]][j,j] =  mu1_prime_list[[i]][j,]  +
                             q[[i]][j,]*( mu1_prime_list[[i]][j,] +
                             mu1_list[[i]][j,]^2  )
  }
}

A2_star_list = list()
for (i in 1:n){
  A2_star_list[[i]] = matrix(rep(0,len=m*m),nrow=m)
  for (j in 1:m){

  A2_star_list[[i]][j,j] =  mu2_prime_list[[i]][j,]  +
                             q[[i]][j,]*( mu2_prime_list[[i]][j,] +
                             mu2_list[[i]][j,]^2  )
  }
}



V1_star_list = list()
for (i in 1:n){
   V1_star_list[[i]] = A1_star_list[[i]]^(1/2) %*% R_hat1 %*% A1_star_list[[i]]^(1/2)
}

V2_star_list = list()
for (i in 1:n){
   V2_star_list[[i]] = A2_star_list[[i]]^(1/2) %*% R_hat2 %*% A2_star_list[[i]]^(1/2)
}


Sigma1_star_list = list()
for (i in 1:n){
  for (j in 1:m){
    Sigma1_star_list[[i]] = (y_star[[i]]-mu1_list[[i]]) %*% t(y_star[[i]]-mu1_list[[i]])
  }
}

Sigma2_star_list = list()
for (i in 1:n){
  for (j in 1:m){
    Sigma2_star_list[[i]] = (y_star[[i]]-mu2_list[[i]]) %*% t(y_star[[i]]-mu2_list[[i]])
  }
}

M_hat1 = matrix(rep(0,len=p*p),nrow=p)
for (i in 1:n){
  M_hat1 = M_hat1 + t(D1_list[[i]]) %*% solve(V1_star_list[[i]]) %*% Sigma1_star_list[[i]] %*% solve(V1_star_list[[i]]) %*% D1_list[[i]]
}

M_hat2 = matrix(rep(0,len=p*p),nrow=p)
for (i in 1:n){
  M_hat2 = M_hat2 + t(D2_list[[i]]) %*% solve(V2_star_list[[i]]) %*% Sigma2_star_list[[i]] %*% solve(V2_star_list[[i]]) %*% D2_list[[i]]
}

H_hat1 = matrix(rep(0,len=p*p),nrow=p)
for (i in 1:n){
  H_hat1 = H_hat1 + t(D1_list[[i]]) %*% solve(V1_star_list[[i]]) %*% D1_list[[i]]
}

H_hat2 = matrix(rep(0,len=p*p),nrow=p)
for (i in 1:n){
  H_hat2 = H_hat2 + t(D2_list[[i]]) %*% solve(V2_star_list[[i]]) %*% D2_list[[i]]
}

B1 = solve(H_hat1) %*% M_hat1 %*% solve(H_hat1)
B2 = solve(H_hat2) %*% M_hat2 %*% solve(H_hat2)

#Compute standard error of estimator and p-value of two-sided test for beta=0
#1-dependent case

SE1 = c()
for (i in 1:p){
    SE1[i] = sqrt(B1[i,i])
}
SE1

test_stat1 = abs(bhat1/SE1)

p_value1 = c()
p_value1 = 2*(1-pnorm(test_stat1))
p_value1

#exchangeable case

SE2 = c()
for (i in 1:p){
    SE2[i] = sqrt(B2[i,i])
}
SE2

test_stat2 = abs(bhat2/SE2)

p_value2 = c()
p_value2 = 2*(1-pnorm(test_stat2))
p_value2



######## METHOD of Yi-Ma-Carroll ###############

#Create the 6x1 matrix D
# D_j is -(1/n)sum_{i=1}^{n} X_{ij}^T X_{ij} \mu'(X_{ij}^T beta_indep)

D = matrix(rep(0,len=m),nrow=1)

for (j in 1:m){
  for (i in 1:n){
    D[1,j]=D[1,j] + (-1/n)*t(X[[i]][j,])%*%t(t(X[[i]][j,]))%*%mu_prime(t(X[[i]][j,]) %*% beta_indep)
  }
}
D

#Create a list of 6x4 matrices Phi_1,...,Phi_250
#Phi[[i]]=Phi_i has j-th row X_{ij}^T (Y_{ij}^*- mu_{ij}(\beta_indep))
Phi = list()
PPhi = matrix(rep(0,len=m*p),nrow=m)

for (i in 1:n){
  for (j in 1:m){

    PPhi[j,] = y_star_c[[i]][j,]%*%t(X[[i]][j,])
  }
  Phi[[i]] = PPhi
}


#Compute the 6x6 matrix V
#V[j,k] is (1/n)sum_{i=1}^{250} \Phi_{ij}^T \Phi_{ik}
V = matrix(rep(0,len=m*m),nrow=m)

for (j in 1:m) {
  for (k in 1:m){
    for (i in 1:n){
      V[j,k] = V[j,k] + (1/n)*y_star_c[[i]][j,]%*%y_star_c[[i]][k,]%*%t(X[[i]][j,])%*%t(t(X[[i]][k,]))
    }
  }
}

#compute the inverse of V
V_inv = solve(V)

#Compute the 1xm vector K
#K=D V^{-1}
K = D%*%V_inv
K

#Defining new variables T
#T_{ij}^l=K_j X_{ij}^l, for l=1,2,3

ccov1 = c()
ccov2 = c()
ccov3 = c()
ccov4 = c()

for (i in 1:n){
  for (j in 1:m){
    ccov1[m*i-m+j] = K[j]*cov1[m*i-m+j]
    ccov2[m*i-m+j] = K[j]*cov2[m*i-m+j]
    ccov3[m*i-m+j] = K[j]*cov3[m*i-m+j]
    ccov4[m*i-m+j] = K[j]*cov4[m*i-m+j]
  }
}


#compute estimated beta using Yi-Ma-Carroll method
logisticfunction3 = function(b){

  prod = ccov1*b[1] + ccov2*b[2] + ccov3*b[3] + ccov4*b[4];
  mean = mu(prod);
  residual = w-mean;
  c(F1 = sum(ccov1*residual), F2 = sum(ccov2*residual), F3 = sum(ccov3*residual), F4 = sum(ccov4*residual)) }

bhat3 = multiroot(f = logisticfunction3, start = c(beta_indep[1],beta_indep[2],beta_indep[3],beta_indep[4]))$root
bhat3


#Start to compute matrices needed for the asymptotic variance of estimator

#Compute the pxp matrix Gamma_0
Gamma0 = matrix(rep(0,len=p*p),nrow=p)

for (i in 1:n){
  for (j in 1:m){

    Gamma0 = Gamma0 + (-1/n)*( K[,j] * (t(t(X[[i]][j,])) %*% t(X[[i]][j,]) ) ) * mu_prime(t(X[[i]][j,]) %*% beta_indep)[1,1]

  }
}
Gamma0

#Compute the inverse of Gamma0
Gamma0_inv = solve(Gamma0)
Gamma0_inv



#Compute the pxp matrix Sigma_psi
Sigma_psi = matrix(rep(0,len=p*p),nrow=p)
for (i in 1:n){
  Sigma_psi = Sigma_psi + (1/n)* ( t(Phi[[i]])%*%t(K)%*%K%*%Phi[[i]])

}


#Compute the pxp matrix Sigma_beta
Sigma_beta = Gamma0_inv%*%Sigma_psi%*%t(Gamma0_inv)
Sigma_beta

#Compute standard error of estimator and p-value of two-sided test for beta=0

SE3 = c()
for (i in 1:p){
    SE3[i] = sqrt(Sigma_beta[i,i]/n)
}
SE3

test_stat3 = abs(bhat3/SE3)
test_stat3


p_value3 = c()
p_value3 = 2*(1-pnorm(test_stat3))
p_value3
