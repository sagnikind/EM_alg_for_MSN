#______________________________________________________________________________________________________________

n = 1000 # number of samples
p = 3 # dimension

#______________________________________________________________________________________________________________

# functions

# W_phi
W = function(x){
  return(dnorm(x)/pnorm(x))
}

# inside Phi(skewing function)
s_f = function(x,a,B,l){
  return(as.numeric((t(l)%*%(x-a))/sqrt(1+t(x-a)%*%solve(B)%*%(x-a))))
}

#______________________________________________________________________________________________________________

# parameters
xi = c(-1,0,3)
Psi = 10*matrix(c(1,1.2,-.9,1.2,4,-.6,-.9,-.6,9),nrow = 3,byrow = TRUE)
eta = c(10,-5,0)
Omega = Psi + eta%*%t(eta)
lambda = as.numeric(1/sqrt(1+t(eta)%*%solve(Psi)%*%eta))*solve(Psi)%*%eta

#______________________________________________________________________________________________________________

# random vector generation 
library("MASS")
Y = mvrnorm(n,rep(0,p),Omega,empirical = FALSE)
X = matrix(nrow = n,ncol = p)
for(i in 1:n){
  if(rnorm(1)<s_f(Y[i,],rep(0,p),Omega,lambda)){
    X[i,] = Y[i,] + xi
  }else{
    X[i,] = -Y[i,] + xi
  }
}

#______________________________________________________________________________________________________________

# empirical mean and covariance

x_bar = as.numeric(t(X)%*%rep(1/n,n))
S = diag(rep(1/n,p))%*%t(X)%*%X - x_bar%*%t(x_bar)

#______________________________________________________________________________________________________________

# functions

# inside fun_q
q_f = function(x,a,B){
  return(as.numeric(as.numeric(1/sqrt(1+t(x-a)%*%solve(B)%*%(x-a)))*solve(B)%*%(x-a)))
}

# fun_q = E(Y|X=x)
fun_q = function(x,a,B,l){
  return(as.numeric(l+W(s_f(x,a,B,l))*q_f(x,a,B)))
}

# fun_P = E(YY^\top|X=x)
fun_P = function(x,a,B,l){
  return(l%*%t(l)+solve(B)+W(s_f(x,a,B,l))*
           (q_f(x,a,B)%*%t(l)+l%*%t(q_f(x,a,B))-s_f(x,a,B,l)*q_f(x,a,B)%*%t(q_f(x,a,B))))
}

# fun_t = E(W|X=x)
fun_t = function(x,a,B,l){
  return(as.numeric(t(x-a)%*%l + sqrt(1+t(x-a)%*%solve(B)%*%(x-a))*W(s_f(x,a,B,l))))
}

# fun_r = E(W^2|X=x)
fun_r = function(x,a,B,l){
  return(as.numeric((t(x-a)%*%l)^2 + (1+t(x-a)%*%solve(B)%*%(x-a)) + (t(x-a)%*%l)*sqrt(1+t(x-a)%*%solve(B)%*%(x-a))*W(s_f(x,a,B,l))))
}


# fun_s = E(WY|X=x)
fun_s = function(x,a,B,l){
  return(as.numeric(fun_t(x,a,B,l)*l+as.numeric(fun_r(x,a,B,l) - (t(x-a)%*%l)*fun_t(x,a,B,l))*q_f(x,a,B)*as.numeric(1/sqrt(1+t(x-a)%*%solve(B)%*%(x-a)))))
}

# fun_llh(log-likelihood)
fun_llh = function(a,B,l){
  h = 0
  for(i in 1:n){
    h = h + log(pnorm(s_f(X[i,],a,B,l)))
  }
  return(as.numeric(-n*0.5*(log(det(B))+sum(diag(S%*%solve(B)))+t(x_bar-a)%*%solve(B)%*%(x_bar-a)) + h))
}

#______________________________________________________________________________________________________________#

library("moments")
# initial parameters
ini_xi = x_bar
ini_Omega = S
ini_lambda = vector(length = p)
for(i in 1:p){
  ini_lambda[i] = skewness(X[,i])
}

#______________________________________________________________________________________________________________#

fin_lambda = rep(0,p)
P = matrix(rep(0,p^2),nrow = p)
PX = rep(0,p)
s = rep(0,p)
for(i in 1:n){
  fin_lambda = fin_lambda + fun_q(X[i,],ini_xi,ini_Omega,ini_lambda)
  P = P + fun_P(X[i,],ini_xi,ini_Omega,ini_lambda)
  PX = PX + t(X[i,])%*%fun_P(X[i,],ini_xi,ini_Omega,ini_lambda)
  s = s + fun_s(X[i,],ini_xi,ini_Omega,ini_lambda)
}
fin_lambda = fin_lambda/n
P = P/n
PX = as.vector(PX/n)
s = s/n

A = P - fin_lambda%*%t(fin_lambda)

library("expm")
fun_xi = function(a){
  H = S + (x_bar-a)%*%t(x_bar-a)
  B = solve(sqrtm(A))%*%sqrtm(sqrtm(A)%*%H%*%sqrtm(A))%*%solve(sqrtm(A))
  return(0.5*n*(sum(diag(solve(B)%*%(S + (x_bar - a)%*%t(x_bar -a)))) + sum(diag(B%*%P)) - t(fin_lambda)%*%B%*%fin_lambda
                + t(a)%*%P%*%a - 2*t(PX)%*%a + 2*t(s)%*%a))
}
t = nlm(fun_xi,ini_xi)
fin_xi = t$estimate
fin_Omega = solve(sqrtm(A))%*%sqrtm(sqrtm(A)%*%(S + (x_bar-fin_xi)%*%t(x_bar-fin_xi))%*%sqrtm(A))%*%solve(sqrtm(A))

while( abs((fun_llh(fin_xi,fin_Omega,fin_lambda) / fun_llh(ini_xi,ini_Omega,ini_lambda)) - 1) > 10^-9){
  ini_xi = fin_xi
  ini_Omega = fin_Omega
  ini_lambda = fin_lambda
  fin_lambda = rep(0,p)
  P = matrix(rep(0,p^2),nrow = p)
  PX = rep(0,p)
  s = rep(0,p)
  for(i in 1:n){
    fin_lambda = fin_lambda + fun_q(X[i,],ini_xi,ini_Omega,ini_lambda)
    P = P + fun_P(X[i,],ini_xi,ini_Omega,ini_lambda)
    PX = PX + t(X[i,])%*%fun_P(X[i,],ini_xi,ini_Omega,ini_lambda)
    s = s + fun_s(X[i,],ini_xi,ini_Omega,ini_lambda)
  }
  fin_lambda = fin_lambda/n
  P = P/n
  PX = as.vector(PX/n)
  s = s/n
  A = P - fin_lambda%*%t(fin_lambda)
  
  fun_xi = function(a){
    H = S + (x_bar-a)%*%t(x_bar-a)
    B = solve(sqrtm(A))%*%sqrtm(sqrtm(A)%*%H%*%sqrtm(A))%*%solve(sqrtm(A))
    return(0.5*n*(sum(diag(solve(B)%*%(S + (x_bar - a)%*%t(x_bar -a)))) + sum(diag(B%*%P)) - t(fin_lambda)%*%B%*%fin_lambda
                  + t(a)%*%P%*%a - 2*t(PX)%*%a + 2*t(s)%*%a))
  }
  t = nlm(fun_xi,ini_xi)
  fin_xi = t$estimate
  fin_Omega = solve(sqrtm(A))%*%sqrtm(sqrtm(A)%*%(S + (x_bar-fin_xi)%*%t(x_bar-fin_xi))%*%sqrtm(A))%*%solve(sqrtm(A))
}

fin_eta = as.numeric(as.numeric(1/sqrt(1+t(fin_lambda)%*%fin_Omega%*%fin_lambda))*fin_Omega%*%fin_lambda)
fin_Psi = fin_Omega - fin_eta%*%t(fin_eta)

fin_xi
fin_eta
fin_Psi




