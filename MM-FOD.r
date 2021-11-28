###The proposed MM-FOD method
library(fdapace)
#################################################################################
### function to compute the test statistics
fun_statistic=function(n,p,y,y_sample) 
{

T_stat <- rep(0,n)
R <- FPCA(y_sample$Ly,y_sample$Lt,list(nRegGrid=p,FVEthreshold=0.9))   
mu <- R$mu
lambda <- R$lambda
phi <- R$phi 
D <- dim(phi)[2]
xi=matrix(0,n,D)
for (i in 1:n) 
{   
      for (k in 1:D) {   
        xi[i,k]<- t(matrix((y$Ly[[i]]-mu)))%*%matrix(phi[,k])}
      T_stat[i] <- sum(xi[i,]^2/lambda)
}

  list(T_stat=T_stat,D=D)
}
#################################################################################
#################################################################################
### function to choose a clean set
fun_h_clean <- function(n,p,y,iter,k_sub)
{

T_matrix <- matrix(0,iter,n)
for (m in 1:iter) 
{
      ind <- sample(c(1:n),floor(n/2), replace = FALSE)
      y_sample <- list()
      for (k in 1:length(ind)) {
        y_sample$Lt[[k]] <- y$Lt[[ind[k]]]
        y_sample$Ly[[k]] <- y$Ly[[ind[k]]]
      }
   T_matrix[m,]=fun_statistic(n,p,y,y_sample)$T_stat 
}
  T_max <- apply(T_matrix, 2, max)
  T_min <- apply(T_matrix, 2, min)
  clean_max=order(T_max,decreasing=FALSE)[1:floor(n*k_sub)]
  clean_min=order(T_min,decreasing=FALSE)[1:floor(n*k_sub)]
  cleanset <- intersect(clean_max,clean_min)
list(cleanset=cleanset)
}


######################################################################################
### function to test outliers through FDR control procedure
fun_test=function(n,p,y,y_clean, alpha0)
{  
  T.test=rep(0,n)
  Q_T=rep(0,n)  

  AA=fun_statistic(n,p,y,y_clean)
  T.test=AA$T_stat
  D=AA$D
  th_median=median(T.test)/qchisq(0.5,D)  

  T.test=T.test/th_median
  for (i in 1:n)
  {Q_T[i]=qnorm((pchisq(T.test[i],D)+1)/2)}
  
  ############  find the threshold t0 by iteration 
  len.tt=100
  U_thr=sqrt(2*log(n))
  t_it=seq(0,U_thr,length=len.tt)
  t_fdr=rep(0,len.tt)
  for (it in 1:len.tt)
  {
    t_fdr[it]=2*n*(1-pnorm(t_it[it]))/max(length(which(abs(Q_T)>=t_it[it])),1)
  }
  print(T.test)
  print(Q_T)
  
  if (length(which(t_fdr<=alpha0))>0){
    t0=t_it[min(which((t_fdr-alpha0)<=0))]}
  else
  {t0=sqrt(2*log(n))}
  print(t_fdr)
  print(t0)  
  ###############  perform outlier detection
  out_set=which(Q_T>t0)
  list(out_set=out_set)
  }
  
 
 
 
 
 
 
 
 
 
 
###################################################### 
##### an illustrate example
########################################################

 
##############  given parameters for the MM-FOD method
iter=100  ## number of sampling 
k_sub=0.6   ## ratio of remaining clean set
#y ##input: functional data with n*p matrix
## preprocessing of functional data
yy=Sparsify(y, tt, p)  ## translate the matrix to list
### perform MM-FOD algorithm
#### (1) find clean set
clean_hat=fun_h_clean(n,p,yy,iter,k_sub)$cleanset
y_clean <- list()
for (i in 1:length(clean_hat)) 
{
  y_clean$Lt[[i]] <- yy$Lt[[clean_hat[i]]]
  y_clean$Ly[[i]] <- yy$Ly[[clean_hat[i]]]
}
####### (2) testing with FDR control
 out_hat=fun_test(n,p,yy,y_clean, alpha0)$out_set
### Output: the set of outliers
##########################################################