#----------------------------------------------------------
# Function for the ABC algorithm
# And matrices needed to generate data from the JR-NMM
#----------------------------------------------------------

#'@rdname exp_matJR
#'@title exp_matJR
#'@description exponential matrix eAt for the JR-NMM
#'@return exponential matrix eAt for the JR-NMM
#'@export
exp_matJR<-function(t,a,b)
{
  eat<-exp(-a*t)
  eatt<-exp(-a*t)*t
  ebt<-exp(-b*t)
  ebtt<-exp(-b*t)*t
  ret<-diag(c(eat+a*eatt,eat+a*eatt,ebt+b*ebtt,eat-a*eatt,eat-a*eatt,ebt-b*ebtt))
  ret[1,4]<-  ret[2,5]<-eatt
  ret[3,6]<-ebtt
  ret[4,1]<-ret[5,2]<--a^2*eatt
  ret[6,3]<--b^2*ebtt
  return (ret)
}

#'@rdname cov_matJR
#'@title cov_matJR
#'@description covariance matrix C for the JR-NMM
#'@return covariance matrix C for the JR-NMM
#'@export
cov_matJR<-function(t,sigma,a,b)
{
  em2at<-exp(-2*a*t)
  em2bt<-exp(-2*b*t)
  e2at<-exp(2*a*t)
  e2bt<-exp(2*b*t)
  sigma<-sigma^2
  ret<-diag(c(em2at*(e2at-1-2*a*t*(1+a*t))*sigma[4]/(4*a^3),
              em2at*(e2at-1-2*a*t*(1+a*t))*sigma[5]/(4*a^3),
              em2bt*(e2bt-1-2*b*t*(1+b*t))*sigma[6]/(4*b^3),
              em2at*(e2at-1-2*a*t*(a*t-1))*sigma[4]/(4*a),
              em2at*(e2at-1-2*a*t*(a*t-1))*sigma[5]/(4*a),
              em2bt*(e2bt-1-2*b*t*(b*t-1))*sigma[6]/(4*b)))
  ret[1,4]<-em2at*t^2*sigma[4]/2
  ret[2,5]<-em2at*t^2*sigma[5]/2
  ret[3,6]<-em2bt*t^2*sigma[6]/2
  ret[4,1]<-em2at*t^2*sigma[4]/2
  ret[5,2]<-em2at*t^2*sigma[5]/2
  ret[6,3]<-em2bt*t^2*sigma[6]/2
  return (ret)
}

#----------------------------------------------------------
# ABC for the JR-NMM:
# Estimation of theta=(sig,mu,C) based on simulated reference data
#----------------------------------------------------------
#'@rdname ABC_Parallel_JRNMM_sigmuC
#'@title ABC_Parallel_JRNMM_sigmuC
#'@description Run 1 iteration from the sdbmsABC after sampling tilde.theta=(sigma,mu,C). Possibly faster code
#'@return d, tilde.theta
#'@export
ABC_Parallel_JRNMM_sigmuC<-function(refY,spec_refY,spx,T,h,M,w,startv,A,B,a,b,v0,r,vmax){
  #sample a theta value from the prior
  sig<-runif(1,1300,2700)
  mu<-runif(1,160,280)
  C<-runif(1,129,141)
  #setting
  grid<-seq(from=0,to=T,by=h)

  #hyperparameters
  Lsupport<-1000
  span_val<-5*T

  #simulated reference data
  #spectral density derived from the reference data
  idx_P<-2:length(spx)

  #conditionally on theta, simulate a new artificial dataset Y
  dm<-exp_matJR(h,a,b)
  cm<-t(chol(cov_matJR(h,c(0,0,0,0.01,sig,1),a,b)))
  Y<-Splitting_JRNMM_output_Cpp(h,startv,grid,dm,cm,mu,C,A,B,a,b,v0,r,vmax)
  #Y<-sol[2,]-sol[3,] #output process
  spec_Y<-spectrum(Y,log="no",span=span_val,plot=FALSE)$spec

  dist_vec<-rep(0,M)
  for (j in 1:M){
    #Integrate absolute error (IAE): invariant spectral density
    spec_refYj<-spec_refY[j,]
    idy_P<-abs(spec_refYj-spec_Y)
    IAE1<-((spx[idx_P]-spx[idx_P-1]) %*% (idy_P[idx_P]+idy_P[idx_P-1]) / 2)

    #Integrate absolute error (IAE): invariant density
    #Find the common support
    startSupp<-min(refY[j,],Y)
    endSupp<-max(refY[j,],Y)
    Dens_refY<-density(refY[j,],kernel="gaussian",from=startSupp,to=endSupp,n=Lsupport,adjust=1)
    invDens_refY<-Dens_refY$y
    dsx<-Dens_refY$x
    idx_D<-2:length(dsx)
    invDens_Y<-density(Y,kernel="gaussian",from=startSupp,to=endSupp,n=Lsupport,adjust=1)$y
    idy_D<-abs(invDens_refY-invDens_Y)
    IAE2<-((dsx[idx_D]-dsx[idx_D-1]) %*% (idy_D[idx_D]+idy_D[idx_D-1]) / 2)

    #combined distance using the weight
    dist_vec[j]<-IAE1+w*IAE2
  }

  #median of the calculated distances
  Dist<-median(dist_vec)

  #return the distance together with the sampled theta value
  ret<-array(0,dim=c(1,4,1))
  ret[,,1]<-c(Dist,sig,mu,C)
  return(ret)
}
