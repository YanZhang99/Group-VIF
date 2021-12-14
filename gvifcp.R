####Detecting Multiple change point via group VIF method####
############################################################

gvifmcp<-function(y,X,l)
{
  #library("changepoint")
  call<-match.call()
  
  if(!is.numeric(l))stop("The length of the partition shoule be numeric")
  
  ###Input data###
  y <- as.vector(y)-mean(y)
  n<-length(y)
  a<-floor(n/l)-1
  q<-dim(X)[2]
  
  if (length(q)== 0) {
    X<-cbind(rep(1, n), X)
    q <- 1
  }
  
  ########### transform:cutting stage ###################
  
  
  X_temp<- NULL
  y_temp <- NULL
  X_temp[[1]] <- as.matrix(X[1:(n - a * l), ])
  y_temp[[1]] <- y[1:(n - a * l)]
  
  for (i in 2:(a+1)) {
    X_temp[[i]] <- as.matrix(X[(n - (a- i + 2) * l + 1):(n - (a -
                                                                i+1) * l), ])
    y_temp[[i]] <- y[(n - (a- i + 2) * l + 1):(n - (a- i+1) * l)]
  }
  
  K_temp <- matrix(0, nrow = a+1, ncol = 1, byrow = TRUE)
  
  X_temp1<-NULL
  for (i in 1:(a+1)) {X_temp1[[i]]<-kronecker(K_temp[i,, drop = FALSE], X_temp[[i]]) }
  
  X0_temp<-NULL
  X0_temp[[0]]<-NULL
  for(i in 1:a){X0_temp[[i]]<-rbind(X0_temp[[i-1]],X_temp1[[i]])}
  
  #initial parameter values###
  w = 0.05;dw = 0.05
  flag <- 0
  res <- NULL
  sigma<-NULL
  m<-0; gvif.cp<-NULL
  
  y.cp<-y_temp[[1]]
  X.cp<-rbind(X_temp[[1]],X_temp[[2]])
  X_new<-NULL
  gvif.cp<-cp<-NULL
  
  ###Detecting a Change-point in the segements### 
  ####Singe change point detection,library("changepoint")###
  ###for normal distribution data###
  
  for (i in 1:a) 
  { 
    alpha <- w/(1+ i- flag)
    y.cp<-c(y.cp,y_temp[[i+1]])
    
    X_new<-rbind(X0_temp[[i]],X_temp[[i+1]])
    
    H<-X.cp%*%solve(t(X.cp)%*%X.cp,tol = 1e-50)%*%t(X.cp)
    
    gamma2<-t(X_new)%*%(X_new)-t(X_new)%*%H%*%(X_new)
    
    res<-y.cp-H%*%(as.matrix(y.cp))
    
    sigma2<-t(res)%*%res/(n-(a-i)*l-(m+1)*q)
    T.cp<-t(t(X_new)%*%res)%*%solve(gamma2,tol=1e-50)%*%t(X_new)%*%res/sigma2
    
    qsq<-qchisq(1-alpha,q)
    
    if(is.na(qsq)){qsq<-0}
    
    if(T.cp >qsq)
    {
      #y1<-NULL
      y1<-c(y_temp[[i]],y_temp[[i+1]])
      cp<-cpts(cpt.meanvar(y1,penalty="Manual",pen.value=0.05,method="AMOC",test.stat = "Normal"))
     
       if(i==1)
      {cp<-cp}
     else
     {
       cp<-cp+(n-(a-i+3)*l+1)
      }
      gvif.cp<-c(gvif.cp,cp)
     
      
      
      ###updata X^(i+1)### 
      if(i!=a)
      { 
        X.cp<-cbind(rbind(X.cp,do.call(cbind,rep(list(X_temp[[i+2]]),m+1))),rbind(X_new,X_temp[[i+2]]))}
      else
      {X.cp<-X.cp}
      
      m<-m+1
      w <- w + dw 
      flag <- i 
      }
    
    else 
      
    {w <- w - alpha / (1 - alpha) 
    
    if(i!=a)
      
    { X.cp<-rbind(X.cp,do.call(cbind,rep(list(X_temp[[i+2]]),m+1))) }  ####if no cp in this segment, update X.cp
    else
    {X.cp<-X.cp}
    
    }
    
    if(w <= 0)break }
  
###deal with redundant change point###
 gcp<-NULL
 #print(gvif.cp)
if (m>1)
{
   
  raw_data<-gvif.cp
  flag_data<-rep(1,m)
 for (k in 1:(m-1))
 { 
   if (flag_data[k] < 1) 
     next
  for (j in (k+1):(m))
  {
    ab<-abs(raw_data[j]-raw_data[k])
    if(ab<=3*l|ab<=10) 
      { flag_data[j] = 0}
    else
      break
  }
 }
  for (k in 1:(m))
  {
    if (flag_data[k] > 0)
    {
      gcp<-append(gcp, raw_data[k])
    }
  }
  
 }
  else{gcp<-gvif.cp}
  
 cpm<-length(gcp)
  return(list(cpset=gcp,cpnumber=cpm))
}


