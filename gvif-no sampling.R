gvif <- function(y, X, group, w0 = 0.05, dw = 0.05, trace = TRUE, mode = c("dense", "sparse")){
  
  mf <- match.call()
  m <- match("mode", names(mf), 0L)
  if(m == 0)mode <- "sparse"
  m<- n
  
  y <- as.vector(y)- mean(y);
  n <- length(y)
  X<-as.matrix(X)
  X<- t(t(X) - apply(X,2,mean))
  if (nrow(X) != length(y)) 
    stop("X and y do not have the same number of observations")
  
  # initials
  sel <- rep(NA, 0)
  J<-length(group)
  res <- y - mean(y)
  sigma <-sd(res)
  w <- w0
  flag <- 0
  X.sel <- cbind(rep(1, n), X[, sel])
  
  
  for(j in 1:J){
    if(mode == "dense")alpha <- w/(1 + j - flag) else alpha <- w/2/j
    
    X.can<-as.matrix(X[,group[[j]]])
    res.can<-t(X.can)%*%res
   
    
    Hm.sam<-X.sel%*%solve(t(X.sel)%*%X.sel,tol=1e-50)%*%t(X.sel)
    
    gamma2.sam<-as.matrix(t(X.can)%*%X.can-t(X.can)%*%Hm.sam%*%X.can)
    
    Tm<-t(res.can)%*%solve(gamma2.sam,tol=1e-50)%*%res.can/sigma^2
    
    dfX.can<-dim(X.can)[2]
    if(Tm > qchisq(1-alpha,dfX.can))
    {
      X.sel <- cbind(X.sel, X.can)
      
      lm.sel <- lm(y ~ X.sel - 1)
      
      if(!all(!is.na(lm.sel$coef)))next
      
      sel <- c(sel, j)
      
      res <- summary(lm.sel)$residuals
      
      sigma<-sqrt(t(res)%*%res)/sqrt(n - dim(X.sel)[2])
      
      if(mode == "dense")w <- w + dw else w <- w + dw - alpha;
      
      flag <- j; }
    
    else if(mode == "dense")w <- w - alpha / (1 - alpha) else  w <- w - alpha;
    
    if(w <= 0)break;
  }
  
  return(list(select = sel,selectmatrix =X.sel[,-1]));
}

