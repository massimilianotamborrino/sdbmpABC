#----------------------------------------------------------
# Run this file to obtain the results (kept ABC posterior samples)
# of the proposed ABC algorithm.
# The kept marginal posterior samples are stored
# in the text file "Marginal_posteriors".
#----------------------------------------------------------

#----------------------------------------------------------
# The user can specify the setting: epsilon, N, M, w, T, h
#----------------------------------------------------------

cut<-10^3 #number of kept samples (relates to the threshold epsilon)
N<-2.5*10^6 #number of iterations
M<-30  #number of observed datasets (corresponds to paths of the output process)
w<-1930.17 #weight for the distance calculation
T<-200 #time interval for the datasets
h<-0.002 #time step (corresponds to Delta)

#----------------------------------------------------------
numDist<-1
grid<-seq(from=0,to=T,by=h) #time grid

#----------------------------------------------------------

#theta_true: parameters used to simulate the reference data
sig<-2000.0
mu<-220.0
C<-135.0

#fixed model coefficients
A<-3.25
B<-22.0
a<-100.0
b<-50.0
v0<-6.0
vmax<-5.0
r<-0.56

#starting value X(0)
#The process X converges exponentially fast to its invariant regime.
#Thus, the choice of X(0) is not crucial.
#Its impact on the distribution of X diminishes exponentially fast.
startv<-c(0.08,18,15,-0.5,0,0)

#----------------------------------------------------------
#generation of the M reference datasets
refData<-matrix(0, nrow = M, ncol = length(grid))
for(j in 1:M)  refData[j,]<-Splitting_JRNMM_output_Cpp(h,startv,grid,exp_matJR(h,a,b),t(chol(cov_matJR(h,c(0,0,0,0.01,sig,1),a,b))),mu,C,A, B, a, b, v0, r, vmax)

#----------------------------------------------------------

#calculation of the M corresponding invariant spectral densities
span_val<-5*T
specDens<-spectrum(refData[1,],log="no",span=span_val,plot=FALSE)
spx<-specDens$freq*(1/h)
stepP<-diff(spx)[1]
Lspec<-length(spx)
specRef<-matrix(0,nrow=M,ncol=Lspec)
for (i in 1:M)   specRef[i,]<-spectrum(refData[i,],log="no",span=span_val,plot=FALSE)$spec


#----------------------------------------------------------
ncl<-detectCores()
cl<-makeCluster(ncl)
registerDoSNOW(cl)
start_time<-Sys.time()
merge_d<-foreach(i=1:N,.combine='rbind',.packages = c('sdbmsABC')) %dopar% {
  ABC_Parallel_JRNMM_sigmuC(refData,specRef,spx,T,h,M,w,startv,A, B, a, b, v0, r, vmax)
}
end_time<-Sys.time()
stopCluster(cl)
merge_d1<-merge_d[,1:4]
#----------------------------------------------------------

#sort the theta values sampled from the prior
#with respect to the derived corresponding distances
sort_d1<-merge_d1[order(merge_d1[,1]),]
#keep all the sorted sig values
sort_sig=sort_d1[,2]
#keep all the sorted mu values
sort_mu=sort_d1[,3]
#keep all the sorted C values
sort_C=sort_d1[,4]

#----------------------------------------------------------

#store the kept (cut/percentile) (marginal) sig samples
#into the txt file "Marginal_posteriors"
sigMat<-matrix(0, nrow = numDist, ncol = cut)
sigMat[1,]<-sort_sig[1:cut]
write(t(sigMat), file = "sig.txt",ncolumns = cut,sep = " ")

#----------------------------------------------------------

#store the kept (cut/percentile) (marginal) mu samples
#into the txt file "Marginal_posteriors"
muMat<-matrix(0, nrow = numDist, ncol = cut)
muMat[1,]<-sort_mu[1:cut]
write(t(muMat), file = "mu.txt",ncolumns = cut,sep = " ")

#----------------------------------------------------------

#store the kept (cut/percentile) (marginal) C samples
#into the txt file "Marginal_posteriors"
CMat<-matrix(0, nrow = numDist, ncol = cut)
CMat[1,]<-sort_C[1:cut]
write(t(CMat), file = "C.txt",ncolumns = cut,sep = " ")

#----------------------------------------------------------

#runtime of the algorithm
end_time-start_time





