#______________________________________________________________________________________________________________

n = 1000 # number of samples
p = 1 # dimension

#______________________________________________________________________________________________________________

# functions

# W_phi
W = function(x){
  return(dnorm(x)/pnorm(x))
}

# inside Phi(skewing function)
s_f = function(x,a,B,l){
  z = (x-a)/sqrt(B)
  return(l*sqrt(B)*z/sqrt(1+z^2))
}

#______________________________________________________________________________________________________________

# parameters
xi = -1
Psi = 1
eta = -0.5
Omega = Psi + eta^2
lambda = (1/sqrt(1+(eta)*(Psi^-1)*eta))*(Psi^-1)*eta

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

x_bar = mean(X)
S = sum(X^2)/n - x_bar^2

#______________________________________________________________________________________________________________

# functions

# fun_q = E(Y|X=x)
fun_q = function(x,a,B,l){
  z = (x-a)/sqrt(B)
  return(l + (z/sqrt(1+z^2))*W(s_f(x,a,B,l))/sqrt(B))
}

# fun_P = E(YY^\top|X=x)
fun_P = function(x,a,B,l){
  z = (x-a)/sqrt(B)
  return(l^2+1/B+W(s_f(x,a,B,l))*l*z*(2+z^2)/(sqrt(B)*(1+z^2)^1.5))
}

# fun_t = E(W|X=x)
fun_t = function(x,a,B,l){
  z = (x-a)/sqrt(B)
  return(l*z*sqrt(B) + sqrt(1+z^2)*W(s_f(x,a,B,l)))
}

# fun_r = E(W^2|X=x)
fun_r = function(x,a,B,l){
  z = (x-a)/sqrt(B)
  return(B*(l*z)^2 + 1 + z^2 + sqrt(B)*l*z*sqrt(1+z^2)*W(s_f(x,a,B,l)))
}

# fun_s = E(WY|X=x)
fun_s = function(x,a,B,l){
  z = (x-a)/sqrt(B)
  return((l*fun_t(x,a,B,l) + z*fun_r(x,a,B,l)/sqrt(B))/(1+z^2))
}

# fun_lllh(log-likelihood)
fun_llh = function(a,B,l){
  h = 0
  for(i in 1:n){
    h = h + log(pnorm(s_f(X[i],a,B,l)))
  }
  return(as.numeric(n*0.5*log(2/pi)-n*0.5*(log((B))+S/B+((x_bar-a)^2)/(B)) + h))
}

#______________________________________________________________________________________________________________#

library("moments")
# initial parameters
ini_xi = x_bar
ini_Omega = S
ini_lambda =  skewness(X[,1])


#______________________________________________________________________________________________________________#

fin_lambda = 0
for(i in 1:n){
  fin_lambda = fin_lambda + fun_q(X[i],ini_xi,ini_Omega,ini_lambda)
}
fin_lambda = fin_lambda/n

P = 0
PX = 0
PX_2 = 0
for(i in 1:n){
  P = P + fun_P(X[i],ini_xi,ini_Omega,ini_lambda)
  PX = PX + X[i]*fun_P(X[i],ini_xi,ini_Omega,ini_lambda)
  PX_2 = PX_2 + (X[i]^2)*fun_P(X[i],ini_xi,ini_Omega,ini_lambda)
}
P = P/n
PX = PX/n
PX_2 = PX_2/n

s = 0
sX = 0
for(i in 1:n){
  s = s + fun_s(X[i],ini_xi,ini_Omega,ini_lambda)
  sX = sX + X[i]*fun_s(X[i],ini_xi,ini_Omega,ini_lambda)
}
s = s/n
sX = sX/n

f_theta = function(a){
  b = sqrt(((x_bar-a)^2 + S)/(P - fin_lambda^2))
  return(0.5*n*((S + (x_bar - a)^2)/b 
                + (P - fin_lambda^2)*b
                + ( - 2*sX +2*a*s + PX_2 - 2*PX*a +P*a^2)))
}
t = nlm(f_theta,ini_xi)

fin_xi = t$estimate[1]
fin_Omega = sqrt(((x_bar-fin_xi)^2 + S)/(P - fin_lambda^2))

while( abs((fun_llh(fin_xi,fin_Omega,fin_lambda) / fun_llh(ini_xi,ini_Omega,ini_lambda)) - 1) > 10^-9){
  ini_xi = fin_xi
  ini_Omega = fin_Omega
  ini_lambda = fin_lambda
  fin_lambda = 0
  for(i in 1:n){
    fin_lambda = fin_lambda + fun_q(X[i],ini_xi,ini_Omega,ini_lambda)
  }
  fin_lambda = fin_lambda/n
  
  P = 0
  PX = 0
  PX_2 = 0
  for(i in 1:n){
    P = P + fun_P(X[i],ini_xi,ini_Omega,ini_lambda)
    PX = PX + X[i]*fun_P(X[i],ini_xi,ini_Omega,ini_lambda)
    PX_2 = PX_2 + (X[i]^2)*fun_P(X[i],ini_xi,ini_Omega,ini_lambda)
  }
  P = P/n
  PX = PX/n
  PX_2 = PX_2/n
  
  s = 0
  sX = 0
  for(i in 1:n){
    s = s + fun_s(X[i],ini_xi,ini_Omega,ini_lambda)
    sX = sX + X[i]*fun_s(X[i],ini_xi,ini_Omega,ini_lambda)
  }
  s = s/n
  sX = sX/n
  
  f_theta = function(a){
    b = sqrt(((x_bar-a)^2 + S)/(P - fin_lambda^2))
    return(0.5*n*((S + (x_bar - a)^2)/b 
                  + (P - fin_lambda^2)*b
                  + ( - 2*sX +2*a*s + PX_2 - 2*PX*a +P*a^2)))
  }
  t = nlm(f_theta,ini_xi)
  
  fin_xi = t$estimate[1]
  fin_Omega = sqrt(((x_bar-fin_xi)^2 + S)/(P - fin_lambda^2))
}

fin_eta = as.numeric(as.numeric(1/sqrt(1+(fin_lambda)^2*fin_Omega))*fin_Omega*fin_lambda)
fin_Psi = fin_Omega - fin_eta^2

fin_xi
fin_eta
fin_Psi
