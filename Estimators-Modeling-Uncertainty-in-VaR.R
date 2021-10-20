
#These codes perform the estimators used to obtain the GARCH (1,1) parameters as such proposed by the thesis.


library(rugarch)
library(tictoc)
library(parallel)
library(Metrics)
library(Rsolnp)
library(fBasics)
library(pracma)
library(PerformanceAnalytics)
library(zoo)
library(quantmod)
library(resample)
library(gdata)
library(NlcOptim)

### Auxiliary functions

 #f_ZeQ is the functions to obtain the standardized residuals (f_ZeQ) of GARCH models.
 
 #HZ calculates the conditional standard deviation.
 
 #.gl is the function to estimate either the degrees of freedom of Student-t distributions or the   shape parameter #of GED distribution.


#Functions to obtain the standardized residuals (f_ZeQ) and conditional sd (H).
f_ZeQ<-function(pars,X){
    
    a0<-pars[1]
    a1<-pars[2]
    b1<-pars[3]
    
    zez<-array(0,length(X))
    
    h<-array(0, length(X))
    
    h[1]<-a0/(1-(a1+b1))
    
    if(h[1]<=0){
     h[1]<-var(X)       #median(X^2) --> original Preminger (2017)
    }
    
    for (t in 2:length(X)) {
      
      h[t]= a0 + a1*(X[(t-1)]^2) + (b1*h[t-1]);
      
    }
    
    for (i in 1:length(h)) {
      zez[i]<-X[i]/(h[i]^0.5)
    }
    return(zez)
}

HZ<-function(pars,X){
  
  a0<-pars[1]
  a1<-pars[2]
  b1<-pars[3]
  
  H<-array(0, (length(X)+1))
  
  H[1]<-a0/(1-(a1+b1))
    
    if(H[1]<=0){
     H[1]<-var(X)       #median(X^2) --> original Preminger (2017)
    }
  
  for (t in 2:(length(X)+1)) {
    
  H[t]=a0+(a1*(X[(t-1)]^2))+(b1*H[t-1])
    
  }
  H<-H^0.5
  
  return(H)
}
########################################################################
########### Estimate the shape parameter of densities ##################
#This function allows to obtain the shape parameters that are required to
#forecast the conditional sd (ugarchforecast), and to estimate VaR.
#inputs:
# Z:standardized innovations
# dis:distribution : std or ged.
.gl<- function(Z, dis) {
  gl=fitdist(distribution = dis, Z,
             control = list(delta = 1e-10, tol = 1e-8, trace = 0))[[1]][["shape"]]
  
  return(gl)
}
###########################################################################
### Maximum likelihood estimator

MLE_f<- function(X) {
  
  
  spec1 = ugarchspec(variance.model = list(model = "sGARCH"),
                     mean.model=list(armaOrder=c(0,0), include.mean=FALSE),
                     distribution.model = dis, fixed.pars = list(shape=df))
  model1<-try(
    ugarchfit(data = X, spec = spec1, solver = "hybrid"),
    silent = TRUE
  )
  if (class(model1) == "try-error"){
    model1 <- ugarchfit(data = X, spec = spec1,
                        solver = "gosolnp",
                        solver.control = list(n.restarts = 10))
  }
  if (model1@fit$convergence != 0){
    cat("\n\nConvergence error!!\n\n")
  } 
  
  
  Param_ML<-model1@fit[["coef"]]
  
  return(Param_ML)
  
}

### Gaussian Quasi-Maximum likelihood estimator

GQML_f<- function(X) {
  
  spec1 = ugarchspec(variance.model = list(model = "sGARCH"),
                     mean.model=list(armaOrder=c(0,0), include.mean=FALSE),
                     distribution.model = "norm")   #fixed.pars = list(shape=df))
  model1<-try(
    ugarchfit(data = X, spec = spec1, solver = "hybrid"),
    silent = TRUE
  )
  if (class(model1) == "try-error"){
    model1 <- ugarchfit(data = X, spec = spec1,
                        solver = "gosolnp",
                        solver.control = list(n.restarts = 10))
  }
  if (model1@fit$convergence != 0){
      cat("\n\nConvergence error!!\n\n")
  } 
  
  
  Param_GQML<-model1@fit[["coef"]]
  
  return(Param_GQML)
  
}


### Student-t Estimator

stQML_f<-function(X){ 
  
  spec1 = ugarchspec(variance.model = list(model = "sGARCH"),
                     mean.model=list(armaOrder=c(0,0), include.mean=FALSE),
                     distribution.model = "std")   #fixed.pars = list(shape=df))
  model1<-try(
    ugarchfit(data = X, spec = spec1, solver = "hybrid"),
    silent = TRUE
  )
  if (class(model1) == "try-error"){
    model1 <- ugarchfit(data = X, spec = spec1,
                        solver = "gosolnp",
                        solver.control = list(n.restarts = 10))
  }
  if (model1@fit$convergence != 0){
      cat("\n\nConvergence error!!\n\n")
  } 
  
  
  
  Param_STD<-model1@fit[["coef"]]
  
  return(Param_STD)
}

### GED Estimator  

 GED_f<-function(X){ 
  
  spec1 = ugarchspec(variance.model = list(model = "sGARCH"),
                     mean.model=list(armaOrder=c(0,0), include.mean=FALSE),
                     distribution.model = "ged")   #fixed.pars = list(shape=df))
  model1<-try(
    ugarchfit(data = X, spec = spec1, solver = "hybrid"),
    silent = TRUE
  )
  if (class(model1) == "try-error"){
    model1 <- ugarchfit(data = X, spec = spec1,
                        solver = "gosolnp",
                        solver.control = list(n.restarts = 10))
  }
  if (model1@fit$convergence != 0){
      cat("\n\nConvergence error!!\n\n")
  } 
  
  
  
  Param_GED<-model1@fit[["coef"]]
 
  return(Param_GED)
}
  

### Non-Gaussian Quasi-Maximum likelihood estimator (Fan,Qi and Xiu, 2014)

########################################################################
  
  NGQML_f<-function(X){
  
  # First step, obtaining parameters of GQMLE    
    
  Param_GQML<-GQML_f(X)
  
  ##################################################################
  
  # Second step, estimating standirdized residuals (f_ZeQ), Fan,Qi and Xiu (2014)
  zez<-f_ZeQ(Param_GQML,X)
  
  
  # Obtaining the scale parameter eta
  
  ni<-function(par){
    eta=par
    nu=7    #7 original
    x=zez/eta
    f = (1+(x^2)/(nu-2))^-((nu+1)/2);
    ll=-mean(-log(eta)+log(f))
  }
  ans<- gosolnp(pars = 0.3, fun=ni,  
                LB = 0.00001, UB = 5,
                control = list(delta = 1e-10, tol = 1e-8, trace = 0),
                n.restarts = 5, n.sim=70)
  eta<-ans[["pars"]]
  
  # Third step, estimating the parameters through NGML estimator (Fan,Qi and Xiu, 2014)
  
  param_Fan<-function(pars){
    a0<-pars[1]
    a1<-pars[2]
    b1<-pars[3]
    sigma<-1
    nu<-7
    
    h<-array(0, length(X))
    h[1]=var(X)
    
    for (t in 2:length(X)) {
      
      h[t]=a0+(a1*X[(t-1)]^2)+(b1*h[t-1])
    }
    x=X/as.matrix(eta*sqrt(h*sigma))
    v<-as.matrix((eta*(sqrt(h*sigma))))
    
    x<-X/v
    f<-(1+(x^2)/(nu-2))^-((nu+1)/2)
    ll<--mean(-log(v)+log(f))
    
  }
  
  ans<- try(gosolnp(pars = c(Param_GQML[1], Param_GQML[2], Param_GQML[3]), fun=param_Fan,  
                LB = c(0.00001,0.0001,0.0001), UB = c(0.5,0.8,0.99),
                control = list(delta = 1e-10, tol =    1e-8, trace = 0),
                n.restarts = 5, n.sim=70)) #cluster = cl)
  
  if (ans[["convergence"]] != 0){
      cat("\n\nConvergence error!!\n\n")
  } 
  
  Param_NGQML<-ans$pars
  
  #z=f_ZeQ(Param_NGQML, X)
  #shape=.gl(z, "std")
  shape = 7
  
  Param_NGQML=c(Param_NGQML,shape=shape)
  
  return(Param_NGQML)
  
}



### Rank Estimator (Andrews, 2012)

  ##########################################################################################
  #Functions for estimating a GARCH(p,q) time series model R using ranks
  #Andrews (2012) 
  #see the original code at (https://faculty.wcas.northwestern.edu/~mea405/Rcode1.txt)
  ###########################################################################################
  Rank_f<-function(X){
    
    w <- length(X)
    
    p<-1
    q<-1
    n<-w
    
    for (i in 1:length(X)) {
    if (X[i]==0)(X[i]=mean(X))*1e-3
       }
    
    ################
    lambda <- function(x){
      if(x<1) {c<-(7*(sqrt(5/7)*qt((x+1)/2,7))^2-5)/((sqrt(5/7)*qt((x+1)/2,7))^2+5)}
      if(x==1) {c<-7}
      return(c)
    }
    ###################################################
    
    ###################################################
    #--------------------------------------------------
    # the weight function lambda squared.
    #-------------------------------------------------
    lambda2 <- function(x){
      c <- lambda(x)^2
      return(c)
    }
    ###################################################
    
    ###################################################
    #--------------------------------------------------
    # getD returns the value of D for a candidate parameter vector theta
    #-------------------------------------------------
    getD <- function(theta){
      #--------------------------------------------------
      # calculate sigma square
      sigmasq <- array(1, (n+max(0,q-p)))
      for (t in (p+1):n){
        sigmasq[t+max(0,q-p)]<-1+sum(theta[1:p]*((X[(t-1):(t-p)])^2))
        if(q>0){sigmasq[t+max(0,q-p)]<-sigmasq[t+max(0,q-p)]+
          sum(theta[(p+1):(p+q)]*sigmasq[(t+max(0,q-p)-1):(t+max(0,q-p)-q)])}
      }
      #--------------------------------------------------
      # calculate epsilon
      epsilon <- log(X^2)-log(abs(sigmasq[(1+max(0,q-p)):(n+max(0,q-p))]))
      #--------------------------------------------------
      # sort epsilon
      sortep <- sort(epsilon[p+1:n])
      #--------------------------------------------------
      # calculate D
      D<-0
      for (t in 1:(n-p)){
        D<-D+lambda(t/(n-p+1))*(sortep[t]-mean(sortep))
      }
      if(min(theta)<0){D=D+10000000}
      if((q>0) & (sum(theta[(p+1):(p+q)])>=1)){D=D+10000000}
      return(D)
    }
    ###################################################
    
    #---------------------------------------------------------------------------------
    #End of auxiliary functions
    #---------------------------------------------------------------------------------
    ###################################################
    #--------------------------------------------------
    # Main Body
    
    #--------------------------------------------------
    # generate 100 initial values for theta=(alpha1/alpha0, ...
    #, alphap/alpha0, beta1, ... ,betaq)
    # under the condition that the sum of alpha1, ... ,alphap, beta1, ... ,betaq is less than 1
    #--------------------------------------------------
    no_initial <- 100
    temp_alphabeta <-matrix(0,nr=no_initial,nc=p+q+1)
    temp_theta <-matrix(0,nr=no_initial,nc=p+q)
    count<-0
    repeat{
      alphabeta <- runif(p+q)
      if (sum(alphabeta)<1){
        count<-count + 1
        temp_alphabeta[count,] <- c(var(X)*(1-sum(alphabeta)),alphabeta)
      }
      if(count == no_initial)break 
    }
    
    for (i in 1:p){
      temp_theta[,i]<-temp_alphabeta[,(i+1)]/temp_alphabeta[,1]
    }
    if (q>0){
      for (i in (p+1):(p+q)){
        temp_theta[,i]<-temp_alphabeta[,(i+1)]
      }
    }
    #--------------------------------------------------
    # calculate D for the initial parameter vectors theta
    #--------------------------------------------------
    temp_D<- array(1000,no_initial)
    for (i in 1:no_initial){
      temp_D[i] <-getD(temp_theta[i,])
    }
    
    #--------------------------------------------------
    # find the 3 initial values with the smallest values for D
    #--------------------------------------------------
    store_D <- array(1000,3)
    store_theta <-matrix(0,nr=3,nc=p+q)
    for (j in 1:3){
      store_D[j] <- min(temp_D)
      store_theta[j,] <- temp_theta[which.min(temp_D),]
      temp_D[which.min(temp_D)] <- max(temp_D)+1
    }
    
    #--------------------------------------------------
    # Using the 3 initial values as starting points, find optimized values for theta
    # The Rank estimate has the smallest corresponding D value
    #--------------------------------------------------
    min_theta <- store_theta[1,]
    min_D <- store_D[1]
    for (j in 1:3){
      theta<-store_theta[j,]
      theta<-optim(par=theta,fn=getD)
      tempD<-theta$value
      if (tempD<min_D){
        min_D<-tempD
        min_theta<-theta$par
      }
    }
    
    sigma_sq <- array(1, (n+max(0,q-p)))
    for (t in (p+1):n){
      sigma_sq[t+max(0,q-p)]<-1+sum(min_theta[1:p]*((X[(t-1):(t-p)])^2))
      if(q>0){sigma_sq[t+max(0,q-p)]<-sigma_sq[t+max(0,q-p)]+
        sum(min_theta[(p+1):(p+q)]*sigma_sq[(t+max(0,q-p)-1):(t+max(0,q-p)-q)])}
    }
    ep<-log((X[(p+1):n])^2)-log(sigma_sq[(p+1+max(0,q-p)):(n+max(0,q-p))])
    
    omega<-mean(exp(ep))
    alpha1<-min_theta[1]*omega
    Param_Rank<-c(omega, alpha1,min_theta[2])
    names(Param_Rank)<-c("omega", "alpha1", "beta1")
    
    #z=f_ZeQ(Param_NGQML, X)
    #shape=.gl(z, "std")
    shape = 7
  
   Param_Rank=c(Param_Rank,shape=shape)
  
  
  return(Param_Rank)
  
  }


### Rank 2

 #This code uses the Rank approach (Andrews,2012) in the first step before obtaining the volatility scale parameter   #(eta) introduced by Fan et al. (2014). See (Andrews, 2014).


Rank_B<-function(X){
  
    w <- length(X)
  
    
    Param<-Rank_f(X) #First step <-Rank estimation (Andrews, 2012).
    
    #Henceforward, the Rank parameters are used to obtain eta and 
    #perform the NGQML(Fan et al. 2014)
    
    zez<-f_ZeQ(Param[1:3],X)
  
    
    ni<-function(par){
    eta=par
    nu=7    #7 original
    x=zez/eta
    f = (1+(x^2)/(nu-2))^-((nu+1)/2);
    ll=-mean(-log(eta)+log(f))
  }
  ans<- gosolnp(pars = 0.3, fun=ni,  
                LB = 0.00001, UB = 5,
                control = list(delta = 1e-10, tol = 1e-8, trace = 0),
                n.restarts = 5, n.sim=70)
  eta<-ans[["pars"]]
    
  
  param_Fan<-function(pars){
    a0<-pars[1]
    a1<-pars[2]
    b1<-pars[3]
    sigma<-1
    nu<-7
    
    h<-array(0, length(X))
    h[1]=var(X)
    
    for (t in 2:length(X)) {
      
      h[t]=a0+(a1*X[(t-1)]^2)+(b1*h[t-1])
    }
    x=X/as.matrix(eta*sqrt(h*sigma))
    v<-as.matrix((eta*(sqrt(h*sigma))))
    
    x<-X/v
    f<-(1+(x^2)/(nu-2))^-((nu+1)/2)
    ll<--mean(-log(v)+log(f))
    
  }
  
  ans<- gosolnp(pars = c(Param[1], Param[2], Param[3]), fun=param_Fan,  
                LB = c(0.00001,0.0001,0.0001), UB = c(0.5,0.8,0.99),
                control = list(delta = 1e-10, tol =1e-8, trace = 0),
                n.restarts = 5, n.sim=70) #cluster = cl)
  
  
  Param_NGQML<-ans$pars
  
  #z=f_ZeQ(Param_NGQML, X)
  #shape=.gl(z, "std")
  shape = 7
  
  Param_NGQML=c(Param_NGQML,shape=shape)
  
  
  return(Param_NGQML)
  
  }
    

### Least-squares non-Gaussian estimation (Preminger Storti, 2017)


   Prem_f1<-function(X){
  
  # First step, obtaining parameters of GQMLE
    
  Param_GQML<-GQML_f(X) 
  
   for (i in 1:length(X)) {
   if (X[i]==0)(X[i]=mean(X))#*1e-3) #*1e-6
    }
  
  ##################################################################
  # Second step, estimating standirdized residuals (f_ZeQ) (Preminger, Storti, 2017)
  
   zes<-f_ZeQ(Param_GQML,X)
  
  
    ces<-mean(log(zes^2)) # setting the trim function c 
  
  # Least-squares estimation (Preminger, Storti, 2017)
    
  LSGAR_f<-function(pars){
    
    a0<-pars[1]
    a1<-pars[2]
    b1<-pars[3]
    
    
    ldata=array(0,length(X))
    ldata=log(X^2)-ces
    
    h<-array(0, length(ldata))
    h[1]<-mean(X^2)
    
    for (t in 2:length(X)) {
      
      h[t]=a0+(a1*X[(t-1)]^2)+(b1*h[t-1])
    }
    
    LSGARobj<-array(0, length(ldata))
    LSGARobj<-(ldata-as.matrix(log(h)))
    #return(LSGARobj)
  }
  
 
  coef<-try(
    lsqnonlin(LSGAR_f, Param_GQML,  options = list(tau=1, tolx = 1e-20, tolg = 1e-20, 
                                                   maxeval = 10000)),
    silent = TRUE
  )
  if (class(coef) == "try-error"){
    coef <- fmincon(Param_GQML, fn=LSGAR_f, tol = 1e-06, maxfeval = 10000, maxiter = 5000)
  }
  
  Param_prem<-coef[[1]] 
  
  if(Param_prem[2]+Param_prem[3]>=1){
    Param_prem[1:3]<-Param_GQML[1:3]
  }
  
  for(i in 1:length(Param_prem)){
  
  if(Param_prem[i]<=0){
    Param_prem[i]=0.0001
  }
  }
  e<-(is.na(Param_prem))
  for(i in 1:length(Param_prem)){
  
  if(e[i] == "TRUE"){
    Param_prem[i]=Param_GQML[i]
  }
  }
  
  #z=f_ZeQ(Param_NGQML, X)
  #shape=.gl(z, "std")
  shape = 7
  
  Param_prem=c(Param_prem,shape=shape)
  names(Param_prem)<-c("omega", "alpha1", "beta1", "shape")
  
  return(Param_prem)
  
  }
