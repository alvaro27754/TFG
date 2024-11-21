################################################
#############Prueba SARAR Hetero################
################################################

rm(list=ls())

#cargar librerías
library(gamlss)
library(reshape)
library(ggplot2)
library(spatialreg)
library(sphet)
library(spdep)


#Se hará la simulación sobre la base de Columbus

data(oldcol)
W=COL.nb
matstand=nb2mat(W)
wt1<-matstand
n=length(COL.OLD$CRIME)

lw1 = mat2listw(wt1, style="W")
y = COL.OLD$CRIME
x = as.matrix(cbind(1, COL.OLD[,c(1,2,7, 8)]))
wt2 = wt1

hetero_sarar = function(x,z,y,wt1,wt2,threshold = 1e-4, max_iter = 200)
{
  diff_lambda = 10000;diff_rho = 10000 
  iter_count = 0
  lambda_fin=0
  rho_fin=0
  n <- length(y)
  B_f=B=(diag(length(y))-lambda_fin*wt2)
  A_f=A=(diag(length(y))-rho_fin*wt1)
  while(diff_lambda > threshold & diff_rho> threshold)
  {
    
    m1<-lm(B_f%*%A_f%*%y~B_f%*%x-1)
    vero_ini=function(rhos){
      beta_ini=solve(t(x)%*%t((diag(n)-rhos[2]*wt2))%*%(diag(n)-rhos[2]*wt2)%*%x)%*%
         t(x)%*%t((diag(n)-rhos[2]*wt2))%*%(diag(n)-rhos[2]*wt2)%*%(diag(n)-rhos[1]*wt1)%*%y #paso 1
      #res_ini=y-x%*%beta_ini #paso 2
      vv=t((diag(n)-rhos[1]*wt1)%*%y-x%*%beta_ini)%*%t((diag(n)-rhos[2]*wt2))%*%
         (diag(n)-rhos[2]*wt2)%*%((diag(n)-rhos[1]*wt1)%*%y-x%*%beta_ini)
      #Paso 3 maximizar la funci?n de log verosmilitud y obtener lambda
      lc=-n/2*log(pi)-1/2*n*log(sum(m1$residuals^2)/m1$df.residual)+log(det(diag(n)-rhos[2]*wt2))+log(det(diag(n)-rhos[1]*wt1))-1/2*vv
      return(-lc)
    }
    
    
    p_ini<-optim(c(rho_fin,lambda_fin),vero_ini,method="L-BFGS-B",lower = rep(-0.99999, 2), upper = rep(0.99999, 2))
    lambda_ini=p_ini$par[2]
    rho_ini=p_ini$par[1]
    
    
    B=(diag(n)-lambda_ini*wt2)
    A=(diag(n)-rho_ini*wt1)
    
    # paso 4 estimar los betas por minimos cuadrados generalizados
    
    m2<-lm(B%*%A%*%y~B%*%X-1)
    
    vero_fin=function(rhos){
      
      beta_ini=solve(t(x)%*%t((diag(length(y))-rhos[2]*wt2))%*%(diag(n)-rhos[2]*wt2)%*%x)%*%
         t(x)%*%t((diag(n)-rhos[2]*wt2))%*%(diag(n)-rhos[2]*wt2)%*%(diag(n)-rhos[1]*wt1)%*%y #paso 1
      #res_ini=y-x%*%beta_ini #paso 2
      vv=t((diag(n)-rhos[1]*wt1)%*%y-x%*%beta_ini)%*%t((diag(n)-rhos[2]*wt2))%*%
         (diag(n)-rhos[2]*wt2)%*%((diag(n)-rhos[1]*wt1)%*%y-x%*%beta_ini)
      #Paso 3 maximizar la funci?n de log verosmilitud y obtener lambda
      lc=-n/2*log(pi)-1/2*n*log(sum(m1$residuals^2)/m1$df.residual)+log(det(diag(n)-rhos[2]*wt2))+log(det(diag(n)-rhos[1]*wt1))-1/2*vv
      return(-lc)
    }
    
    p_fin<-optim(c(p_ini$par[1],p_ini$par[2]),vero_fin,method="L-BFGS-B",lower=rep(-0.99999, 2), upper = rep(0.99999, 2))
    
    
    lambda_fin=p_fin$par[2]
    rho_fin=p_fin$par[1]
    
    diff_lambda=lambda_fin-lambda_ini;diff_rho=rho_fin-rho_ini
    B_f=(diag(length(y))-lambda_fin*wt2)
    A_f=(diag(length(y))-rho_fin*wt1)
    
    m_f<-gamlss(B_f%*%A_f%*%y~B_f%*%X-1,~Z-1)
    
    iter_count = iter_count + 1
    if(iter_count > max_iter) {
      stop("No converge")
    }
    
  }
  
  coef = c("(Intercept)" = m_f$mu.coefficients[1], x1 = m_f$mu.coefficients[2], m_f$mu.coefficients[3], 
           "(rho)"=rho_fin[1],"(lambda)"=lambda_fin[1],
           "alpha_0"=m_f$sigma.coefficients[1], "alpha_1"=m_f$sigma.coefficients[2], "alpha_2"=m_f$sigma.coefficients[3])
  return(coef)
}

hetero_sarar = function(x,z,y,wt1,wt2,threshold = 1e-4, max_iter = 200)
{
   diff_lambda = 10000;diff_rho = 10000 
   iter_count = 0
   lambda_fin=0
   rho_fin=0
   B_f=B=(diag(length(y))-lambda_fin*wt2)
   A_f=A=(diag(length(y))-rho_fin*wt1)
   while(diff_lambda > threshold & diff_rho> threshold)
   {
      
      #Paso 1 estimar un GAMLSS para obtener una versión inicial de la varianza
      
      m1<-gamlss(B_f%*%A_f%*%y~B_f%*%x-1,~z)
      
      var2<-diag(c(predict(m1,what="sigma",type="response")^2))
      
      vero_ini=function(rhos){
         
         beta_ini=solve(t(x)%*%t((diag(length(y))-rhos[2]*wt2))%*%solve(var2)%*%(diag(length(y))-rhos[2]*wt2)%*%x)%*%t(x)%*%t((diag(length(y))-rhos[2]*wt2))%*%solve(var2)%*%(diag(length(y))-rhos[2]*wt2)%*%(diag(length(y))-rhos[1]*wt1)%*%y #paso 1
         #res_ini=y-x%*%beta_ini #paso 2
         vv=t((diag(length(y))-rhos[1]*wt1)%*%y-x%*%beta_ini)%*%t((diag(length(y))-rhos[2]*wt2))%*%solve(var2)%*%(diag(length(y))-rhos[2]*wt2)%*%((diag(length(y))-rhos[1]*wt1)%*%y-x%*%beta_ini)
         ve=var2;L=chol(ve)
         ldetW=2*sum(diag(log(L))) 
         #Paso 3 maximizar la funci?n de log verosmilitud y obtener lambda
         lc=-n/2*log(pi)-1/2*ldetW+log(det(diag(length(y))-rhos[2]*wt2))+log(det(diag(length(y))-rhos[1]*wt1))-1/2*vv
         return(-lc)
      }
      
      
      p_ini<-optim(c(0,0),vero_ini,method="L-BFGS-B",lower = rep(-0.99999, 2), upper = rep(0.99999, 2))
      lambda_ini=p_ini$par[2]
      rho_ini=p_ini$par[1]
      
      B=(diag(length(y))-lambda_ini*wt2)
      A=(diag(length(y))-rho_ini*wt1)
      
      # paso 4 estimar los betas por minimos cuadrados generalizados
      
      m2<-gamlss(B%*%A%*%y~B%*%x-1,~z)
      
      var2<-diag(c(predict(m2,what="sigma",type="response")^2))
      
      vero_fin=function(rhos){
         
         beta_ini=solve(t(x)%*%t((diag(length(y))-rhos[2]*wt2))%*%solve(var2)%*%(diag(length(y))-rhos[2]*wt2)%*%x)%*%t(x)%*%t((diag(length(y))-rhos[2]*wt2))%*%solve(var2)%*%(diag(length(y))-rhos[2]*wt2)%*%(diag(length(y))-rhos[1]*wt1)%*%y #paso 1
         #res_ini=y-x%*%beta_ini #paso 2
         vv=t((diag(length(y))-rhos[1]*wt1)%*%y-x%*%beta_ini)%*%t((diag(length(y))-rhos[2]*wt2))%*%solve(var2)%*%(diag(length(y))-rhos[2]*wt2)%*%((diag(length(y))-rhos[1]*wt1)%*%y-x%*%beta_ini)
         ve=var2;L=chol(ve)
         ldetW=2*sum(diag(log(L))) 
         #Paso 3 maximizar la funci?n de log verosmilitud y obtener lambda
         lc=-n/2*log(pi)-1/2*ldetW+log(det(diag(length(y))-rhos[2]*wt2))+log(det(diag(length(y))-rhos[1]*wt1))-1/2*vv
         return(-lc)
      }
      
      p_fin<-optim(c(p_ini$par[1],p_ini$par[2]),vero_fin,method="L-BFGS-B",lower=rep(-0.99999, 2), upper = rep(0.99999, 2))
      
      
      lambda_fin=p_fin$par[2]
      rho_fin=p_fin$par[1]
      
      diff_lambda=lambda_fin-lambda_ini;diff_rho=rho_fin-rho_ini
      B_f=(diag(length(y))-lambda_fin*wt2)
      A_f=(diag(length(y))-rho_fin*wt1)
      
      m_f<-gamlss(B_f%*%A_f%*%y~B_f%*%x-1,~z)
      
      iter_count = iter_count + 1
      if(iter_count > max_iter) {
         stop("No converge")
      }
      
   }
   
   coef = c(m_f$mu.coefficients,
            "(rho)"=rho_fin[1],"(lambda)"=lambda_fin[1],
            "alpha_0"=m_f$sigma.coefficients[1], "alpha_1"=m_f$sigma.coefficients[2], "alpha_2"=m_f$sigma.coefficients[3])
   return(coef)
}



m_3<-hetero_sarar(x=X,z=0*X,y=y,wt1=wt1,wt2=wt1)
m_3

sacsarlm(y~X-1, listw = lw1)
