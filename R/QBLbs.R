#' Quantitative Bayesian Lasso for Detecting Rare (or Common) Haplotype Association Correcting for Population Stratification
#'
#' \code{QBLbs} is designed to identify rare (or common) haplotype effects associated with a quantitative trait.
#' The software is able to adjust for population stratification using principal components (PCs). 
#' The input does not allow missing observations and subjects with missing data are removed. 
#' The function returns an object containing posterior inferences after the burn-in period.
#' 
#' @param dat The data set should consist of n rows and 1+n.cov+2p (allelic format) or 1+n.cov+p (genotypic format) columns, 
#' where n is the number of individuals, n.cov is the number of covariates and PCs, 
#' and p is the number of SNPs. The first column is the trait, followed by covariates and PCs.
#' If in allelic format, the other 2*p columns are allels with one column for each allele of the single-locus genotypes.
#' If in genotypic format, the other p columns are genotypes with one column for each SNP.
#' Missing data are not allowed and will be removed.
#' 
#' @param numSNPs The number of SNPs which form the haplotypes.
#' @param allelic TRUE if in allelic format and FALSE if in genotypic format; default is TRUE.
#' @param baseline The name of the baseline haplotype (e.g., "h00000"); default is the one with the maximum frequency.
#' @param interaction TRUE if consider the interaction between covariates and SNPs; only considers interaction with a maximum of the first 2 covariates; FALSE if else.
#' @param cov The number of covariates and PCs.
#' @param pooling.tol Haplotypes whose frequency is below pooling tolerance will be pooled together; default is 0
#' @param zero.tol Tolerance for haplotype frequencies below which haplotypes are assumed not to exist; default is 1/(20*n) where n is the number of individuals.
#' @param a First hyperparameter of the Gamma(\eqn{a},\eqn{b}) prior for regression coefficients,
#'    \strong{\eqn{\beta}}. The prior variance of \eqn{\beta} is 2/\eqn{\lambda^2}. The Gamma prior parameters a and b 
#'    are formulated such that the mean and variance of
#'    the Gamma distribution are \eqn{a/b} and \eqn{a/b^2}. The default value of a is 20.
#'
#' @param b Second hyperparameter of the Gamma(a,b) distribution described above; default
#'    is 20.
#' @param start.beta Starting value of all regression coefficients, \strong{\eqn{\beta}};
#'    default is 0.01.
#' @param lambda Starting value of the \eqn{\lambda} parameter described above;
#'    default is 1.
#' @param D Starting value of the D parameter, which is the within-population
#'    inbreeding coefficient; default is 0.
#' @param seed Seed to be used for the MCMC in Bayesian Lasso; default is a
#'    random seed. If exact same results need to be reproduced, seed should be
#'    fixed to the same number.
#' @param e A (small) number \eqn{\epsilon} in the null hypothesis of no association,
#'    \eqn{H_0: |\beta| \le \epsilon}. The default is 0.1. Changing e from the default of 0.1 may necessitate choosing a
#'    different threshold for Bayes Factor (one of the outputs) to infer
#'    association. 
#' @param burn.in Burn-in period of the MCMC sampling scheme; default is 20000.
#' @param num.it Total number of MCMC iterations including burn-in; default is
#'    50000.
#' @param sigmas Starting value of the \eqn{\sigma^2} parameter, which is the random error; default is 1.
#' @param asig First hyperparameter of the Inverse-Gamma(\eqn{asig}, \eqn{bsig}) prior for the random error,
#'    \eqn{\sigma^2}. The prior mean and variance of \eqn{\sigma^2} is \eqn{bsig/(asig-1)} and 
#'    \eqn{bsig^2/(asig-1)^2(asig-2)}, respectively. The default value of \eqn{asig} is 2.
#'
#' @param bsig Second hyperparameter of Inverse-Gamma(\eqn{asig}, \eqn{bsig}) described above; default
#'    is 1.
#' @param CC The constant value set for the proposal distribution for updating haplotype
#' frequencies, \strong{\eqn{f}}. The proposal distribution of the MH algorithm is 
#' Dirichlet\eqn{(a_1,...,a_p)} with \eqn{a_1+...+a_p=CC}; default is 1000.
#' @param dumpConvergence TRUE if check convergence; FALSE if else; default is FALSE.
#' @param n.chains The number of chains for checking convergence; default if 1.
#' @param plot.interval The interval between two posterior samples that are extracted to draw trace plots.
#'
#' @return Return a list with the following components:
#' \describe{
#' \item{BF}{Bayes Factors for all regressors (haplotypes, covariates and PCs).}
#'
#' \item{beta}{The coefficient estimates for all regressors.}
#'
#' \item{CI.beta}{The 95% credible intervals for all regressors.}
#'
#' \item{CI.lambda}{The 95% credible intervals for \eqn{\lambda}.}
#'
#' \item{CI.D}{The 95% credible intervals for D.}
#' 
#' \item{freq}{The estimated haplotype frequencies.}
#'
#'}
#' 
#'  @examples
#'  data(QBLbsData)
#'  out<-QBLbs(QBLbsData, numSNPs=5, cov=5, burn.in=600, num.it=1000)
#'  
#' @export
#'
#' @useDynLib LBL QBLmcmc
#'
#' 
QBLbs<-function(dat, numSNPs=5, allelic=TRUE, baseline = "missing", interaction=F,cov=1, pooling.tol=0,zero.tol="missing",
                            a = 20, b = 20, start.beta = 0.01, lambda = 1, D = 0, seed = NULL, e = 0.1, burn.in = 20000, 
                            num.it = 50000,sigmas=1,asig=2,bsig=1, CC=1000,
                            dumpConvergence=F,n.chains=1,plot.interval=NULL)
{
  
  if(!require("hapassoc")){
    install.packages("hapassoc")
    library(hapassoc)
  }
  if(!require("coda")){
    install.packages("coda")
    library(coda)
  }
  if(!require("smoothmest")){
    install.packages("smoothmest")
    library(smoothmest)
  }
  
  if (missing(dat)) {
    stop(paste("Must provide data!\n\n"))
  }
  
  if(allelic){
    if(ncol(dat)!=(1+cov+2*numSNPs)){
      stop(paste("The data should have 1+cov+2*numSNPs columns.\n\n"))
    }
  }else{
    if(ncol(dat)!=(1+cov+numSNPs)){
      stop(paste("The data should have 1+cov+numSNPs columns.\n\n"))
    }
  }
  
  if(zero.tol=="missing"){
    haplos.new.list<-pre.hapassoc(dat, numSNPs=numSNPs, pooling.tol=pooling.tol, allelic=allelic,verbose=F)
  }else{
    haplos.new.list<-pre.hapassoc(dat, numSNPs=numSNPs, pooling.tol=pooling.tol,zero.tol=zero.tol, allelic=allelic,verbose=F)
  }
  
  if (is.null(seed)==FALSE) set.seed(seed)
  haplos.names <- names(haplos.new.list$initFreq)
  freq <- haplos.new.list$initFreq
  if (baseline=="missing")
  {
    baseline <- haplos.names[which.max(freq)]
  }
  
  column.subset <- colnames(haplos.new.list$haploDM) != baseline
  n.cov<-dim(haplos.new.list$nonHaploDM)[2] - 1
  
  if (cov!=0){
    cov.data<-haplos.new.list$nonHaploDM[,-1]    
    names.cov<-names(haplos.new.list$nonHaploDM)[-1]  
    
    if(interaction==T)
    {
      if(cov>1){
        # only considers interaction with the first 2 cov
        cov.data.int <- cov.data ## create and add dummy interaction columns
        cov.data.int<-cbind(haplos.new.list$haploDM[, column.subset]*cov.data[,1],
                            haplos.new.list$haploDM[, column.subset]*cov.data[,2],
                            cov.data.int)
  
        t2<-dim(cbind(haplos.new.list$haploDM[, column.subset]*cov.data[,1],haplos.new.list$haploDM[, column.subset]*cov.data[,2]))[2]
        colnames(cov.data.int)[1:t2]<-c(paste(names(haplos.new.list$haploDM[, column.subset,drop=F]), names.cov[1], sep=""),paste(names(haplos.new.list$haploDM[, column.subset,drop=F]), names.cov[2], sep=""))
        colnames(cov.data.int)[c(t2+1,t2+2)]<-names.cov
        
        hdat <- cbind(haplos.new.list$nonHaploDM[,1], haplos.new.list$haploDM[, column.subset], cov.data.int)     
        colnames(hdat) <- c(colnames(haplos.new.list$nonHaploDM)[1], colnames(haplos.new.list$haploDM[, column.subset,drop=F]), colnames(cov.data.int))

      }else{ #cov=1
        
        cov.data.int <- cov.data
        cov.data.int<-cbind(haplos.new.list$haploDM[, column.subset]*cov.data,cov.data.int)
        t2<-dim(haplos.new.list$haploDM[, column.subset,drop=F]*cov.data)[2]
        
        colnames(cov.data.int)[1:t2]<-paste(names(haplos.new.list$haploDM[, column.subset,drop=F]), names.cov, sep="")
        colnames(cov.data.int)[t2+1]<-names.cov
        
        hdat <- cbind(haplos.new.list$nonHaploDM[,1], haplos.new.list$haploDM[, column.subset], cov.data.int)     
        colnames(hdat) <- c(colnames(haplos.new.list$nonHaploDM)[1], colnames(haplos.new.list$haploDM[, column.subset,drop=F]), colnames(cov.data.int))
      }
      

    }else{ #interaction=F

      hdat <- cbind(haplos.new.list$nonHaploDM[,1], haplos.new.list$haploDM[, column.subset], cov.data)
      colnames(hdat) <- c(colnames(haplos.new.list$nonHaploDM)[1], colnames(haplos.new.list$haploDM[, column.subset,drop=F]), 
                          names(haplos.new.list$nonHaploDM)[-1])
    }
    
  }else{ #cov=0
    
    hdat <- cbind(haplos.new.list$nonHaploDM[,1],haplos.new.list$haploDM[, column.subset])
    colnames(hdat) <- c(colnames(haplos.new.list$nonHaploDM)[1],colnames(haplos.new.list$haploDM[, column.subset,drop=F]))
    
  }
  
  hap.name=colnames(hdat[,2:dim(haplos.new.list$haploDM)[2]])
  
  ####### prepare data for MCMC ######
  ID <- haplos.new.list$ID
  N <- sum(haplos.new.list$wt)
  y <- as.numeric(hdat[, 1])
  x <- data.matrix(hdat[, -1,drop=F])  
  colnames(x)<-NULL
  
  freq.new<-freq[names(haplos.new.list$haploDM[, column.subset,drop=F])]
  freq.new[length(freq.new)+1]<-freq[baseline]
  freq.new<-as.vector(freq.new)
  num.haplo.id<-as.vector(table(ID))
  x.length<-as.integer(dim(x)[2]) 
  
  h.length<-dim(haplos.new.list$haploDM)[2]
  freq.data<-haplos.new.list$haploDM[, column.subset,drop=F] 
  freq.num<-matrix(rep(NA,2*(dim(freq.data)[1])),ncol=2)
  for(i in 1:dim(freq.data)[1])
  {
    for(j in 1:dim(freq.data)[2])
    {
      if(freq.data[i,j]==2)
      {
        freq.num[i,1]=j
        freq.num[i,2]=j
        break
      }
      if(freq.data[i,j]==1)
      {
        if(is.na(freq.num[i,1])==T) freq.num[i,1]=j
        else
        {
          freq.num[i,2]=j
          break
        }
      }
    }
    if(is.na(freq.num[i,1])==T)
    {
      freq.num[i,1]=dim(freq.data)[2]+1
      freq.num[i,2]=dim(freq.data)[2]+1
    }
    if(is.na(freq.num[i,2])==T) freq.num[i,2]=dim(freq.data)[2]+1
  }
  
  
  beta=rep(start.beta, x.length+1) #initial all beta_j to a small number, x.length excludes baseline hap
  beta.out<-numeric((num.it-burn.in)*(x.length+1)) 
  lambda.out<-numeric(num.it-burn.in)
  sigmas.out<-numeric(num.it-burn.in)
  freq.out<-numeric((num.it-burn.in)*(h.length)) 
  D.out<-numeric(num.it-burn.in) 
  up.xz<-numeric(x.length*N)
  BDIC.out<-numeric(num.it-burn.in) 
  
  tmp<-apply(freq.num, 1, paste, collapse="")
  uniq.freq.num<-freq.num[match(unique(tmp),tmp),] 
  num.uniq.hap.pair<-nrow(uniq.freq.num)
  
  #dyn.load("QBLbs.so")
  
  if(dumpConvergence==T){

    if(is.null(plot.interval)==T){
      plot.interval<-(num.it-burn.in)%/%300
    }
    
    my.draws<-vector("list",n.chains) 
    beta.matrix<-vector("list",n.chains)
    lambda.vector<-matrix(numeric(n.chains*(num.it-burn.in)),ncol=n.chains)
    
    name1<-colnames(hdat[, -1,drop=F])
    names<-c("Intercept",name1,"lambda")          
    
    for(i in 1:n.chains){
      start.lambda<-rgamma(1,a,b)
      start.beta<-rdoublex(1,mu=0,lambda=1/start.lambda)
      beta<-rep(start.beta, x.length+1)
      out<-.C("QBLmcmc", PACKAGE="LBL", x=as.double(x), 
              n=as.integer(dim(x)[1]), as.double(y), 
              as.integer(N), as.integer(num.haplo.id), 
              as.integer(x.length), as.integer(h.length), as.integer(freq.num), 
              as.integer(uniq.freq.num), as.integer(num.uniq.hap.pair), #082421
              as.double(freq.new), as.double(D), as.double(sigmas),
              as.double(beta), as.double(a), as.double(b), as.double(asig),
              as.double(bsig),  as.double(start.lambda), as.integer(num.it), 
              as.integer(burn.in), beta.out=as.double(beta.out), 
              lambda.out=as.double(lambda.out),  sigmas.out=as.double(sigmas.out),
              freq.out=as.double(freq.out),D.out=as.double(D.out),as.integer(n.cov),
              up.xz=as.double(up.xz),as.integer(CC),BDIC.out=as.double(BDIC.out))                       
      
      beta.outn<-matrix(out$beta.out,nrow=num.it-burn.in, byrow=TRUE)
      lambda.outn<-matrix(out$lambda.out,nrow=num.it-burn.in, byrow=TRUE)
      
      beta_lambda_result<-cbind(beta.outn,lambda.outn)
      colnames(beta_lambda_result)<-names
      
      my.draws[[i]]<-mcmc(beta_lambda_result)
      beta.matrix[[i]]<-beta.outn
      lambda.vector[,i]<-lambda.outn  
      
    }

    mh.list <- mcmc.list(my.draws)
    
    if(n.chains==1){
      D_Raf<-raftery.diag(my.draws[[1]],q=0.025,r=0.005,s=0.95)
      num.beta<-dim(beta.matrix[[1]])[2]
      acceptance.rate<-1-rejectionRate(my.draws[[1]])
      
      #Convergence Numeric Summary
      sink("QBLDiagnosticNumericalSummary.txt")
      cat("Gelman Diagnostic needs at least 2 chains!\n\n\n")
      cat("Raftery Diagnostic \n")
      print(D_Raf)
      cat("Acceptance Rate \n")
      cat(acceptance.rate,"\n")
      sink()
      
      #Convergence Graphical Summary
      pdf(file = "QBLDiagnosticPlot.pdf")
      
      par(mfrow=c(3,2))
      
      #draw trace plot for beta
      seq.get<-seq(1,length(lambda.vector[,1]),plot.interval)
      
      for (i in 1: num.beta){
        plot(seq.get,beta.matrix[[1]][seq.get,i],type="l",col=1,ylim=c(-2,2),
             ylab=paste("beta",i,sep=""),xlab="iteration")
        abline(h=0,col=(n.chains+1))
      }
      
      #draw trace plot for lambda
      plot(seq.get,lambda.vector[seq.get,1],type="l",col=1,ylim=c(-2,2),
           ylab="lambda",xlab="iteration")
      abline(h=0,col=(n.chains+1))
      dev.off()
      
    }else if(n.chains>1){
      D_Gel<-gelman.diag(mh.list)
      
      D_Raf<-vector("list",n.chains)
      acceptance.rate<-vector("list",n.chains)
      
      num.beta<-dim(beta.matrix[[1]])[2]
      
      for (i in 1:n.chains){
        D_Raf[[i]]<-raftery.diag(my.draws[[i]],q=0.025,r=0.005,s=0.95)
        acceptance.rate[[i]]<-1-rejectionRate(my.draws[[i]]) 
      }
      
      #Convergence Numeric Summary
      sink("QBLDiagnosticNumericalSummary.txt")
      cat("Gelman Diagnostic \n")
      print(D_Gel)
      cat("\n\n\n")
      
      for (i in 1:n.chains){
        cat("Chain ",i,"'s Result \n",sep="")
        cat("Raftery Diagnostic \n")
        print(D_Raf[[i]])
        cat("Acceptance Rate \n")
        print(acceptance.rate[[i]])
        cat("\n\n\n")
      }
      sink()
      
      #Convergence Graphical Summary
      pdf(file = "QBLDiagnosticPlot.pdf")
      
      # draw Gelman plot
      gelman.plot(mh.list)
      
      par(mfrow=c(3,2))
      
      #draw trace plot for beta
      seq.get<-seq(1,length(lambda.vector[,1]),plot.interval)
      
      
      for (i in 1: num.beta){
        plot(seq.get,beta.matrix[[1]][seq.get,i],type="l",col=1,ylim=c(min(beta.matrix[[1]][seq.get,i])-0.5,max(beta.matrix[[1]][seq.get,i])+0.5),
             ylab=names[i],xlab="iteration")
        for (j in 2:n.chains){ 
          lines(seq.get,beta.matrix[[j]][seq.get,i],col=j)
        }
        abline(h=0,col=(n.chains+1))
      }
      
      #draw trace plot for lambda
      plot(seq.get,lambda.vector[seq.get,1],type="l",col=1,ylim=c(min(lambda.vector[seq.get,1])-0.5,max(lambda.vector[seq.get,1])+0.5),ylab="lambda"
           ,xlab="iteration")
      for (i in 2: n.chains){
        lines(seq.get,lambda.vector[seq.get,i],col=i,ylab="lambda")
      }
      abline(h=0,col=(n.chains+1))
      
      dev.off()
      
    }
    
  }else{ #dumpConvergence=F
    
    out<-.C("QBLmcmc", x=as.double(x), 
            n=as.integer(dim(x)[1]), as.double(y), 
            as.integer(N), as.integer(num.haplo.id), 
            as.integer(x.length), as.integer(h.length), as.integer(freq.num), 
            as.integer(uniq.freq.num), as.integer(num.uniq.hap.pair), #082421
            as.double(freq.new), as.double(D), as.double(sigmas),
            as.double(beta), as.double(a), as.double(b), as.double(asig),
            as.double(bsig),  as.double(lambda), as.integer(num.it), 
            as.integer(burn.in), beta.out=as.double(beta.out), 
            lambda.out=as.double(lambda.out),  sigmas.out=as.double(sigmas.out),
            freq.out=as.double(freq.out),D.out=as.double(D.out),as.integer(n.cov),
            up.xz=as.double(up.xz),as.integer(CC),BDIC.out=as.double(BDIC.out))                       
    
    BDIC.out<-out$BDIC.out
    BDIC_result=-2*mean(BDIC.out)+2*var(BDIC.out)

    ####### output ########
    beta.out<-matrix(out$beta.out,nrow=num.it-burn.in, byrow=TRUE)
    freq.out<-matrix(out$freq.out,nrow=num.it-burn.in, byrow=TRUE)
    sigmas.out<-out$sigmas.out
    lambda.out<-out$lambda.out
    D.out<-out$D.out
    
    ci.beta<-numeric((x.length+1)*2)
    ci.lambda<-numeric(2)
    ci.D<-numeric(2)
    post.mean.beta<-numeric(x.length+1)
    ci.freq<-numeric((h.length)*2)
    post.mean.freq<-numeric(h.length)
    k<-1
    for (i in 1:(x.length+1))
    {
      ci.beta[k]<-quantile(beta.out[,i], probs=0.025)
      ci.beta[k+1]<-quantile(beta.out[,i], probs=0.975)
      k<-k+2
      post.mean.beta[i]<- mean(beta.out[,i])
    }
    
    ci.lambda<-c(quantile(out$lambda.out, probs=0.025),quantile(out$lambda.out, probs=0.975))
    ci.D<-c(quantile(out$D.out, probs=0.025),quantile(out$D.out, probs=0.975))
    
    k<-1
    for (i in 1:(h.length))
    {
      ci.freq[k]<-quantile(freq.out[,i], probs=0.025)
      ci.freq[k+1]<-quantile(freq.out[,i], probs=0.975)
      k<-k+2
      post.mean.freq[i]<- mean(freq.out[,i])
    }
    
    prob.alt<-numeric(x.length+1)
    BF<-numeric(x.length+1)
    
    for (i in 1:(x.length+1))
    {
      prob.alt[i]<-length(beta.out[,i][abs(beta.out[,i]) > e])/(num.it-burn.in)
    }
    
    
    for (i in 1:(x.length+1)) 
    {
      prior.prob<-(b/(e+b))^a
      prior.odds<-prior.prob/(1-prior.prob)
      if (prob.alt[i]<=(100*prior.odds)/(100*prior.odds+1))
      {
        BF[i]<-round((prob.alt[i]/(1-prob.alt[i]))/prior.odds,4)
      } else
        BF[i]<-">100"
    }
    
    ci.beta<-matrix(ci.beta,nrow=x.length+1, ncol=2, byrow=TRUE)
    ci.OR<-data.frame(exp(ci.beta))
    OR<-exp(post.mean.beta)
    OR<-round(OR,4)
    ci.OR<-round(ci.OR,4)
    post.mean.freq<-round(post.mean.freq,4)
    ci.freq<-round(ci.freq,4)
    ci.lambda<-round(ci.lambda,4)
    ci.D<-round(ci.D,4)
    
    name1<-colnames(hdat[, -1,drop=F])
    pname<-c(baseline,name1)         
    names(BF)<-pname
    names(OR)<-pname
    names(post.mean.beta)<-pname
    
    ci.beta<-data.frame(pname,ci.beta)
    colnames(ci.beta)<-c("Name", "Lower", "Upper")
    
    ci.freq<-matrix(ci.freq,nrow=(h.length), ncol=2, byrow=TRUE)
    ci.freq<-data.frame(c(colnames(hdat[, 2:h.length]),baseline), ci.freq)
    colnames(ci.freq)<-c("Hap", "Lower", "Upper")
    names(post.mean.freq)<-c(colnames(hdat[, 2:h.length,drop=F]), baseline)   
    ci.OR<-data.frame(pname,ci.OR) 
    colnames(ci.OR)<-c("Name","Lower","Upper")
    names(ci.lambda)<-c("Lower","Upper")
    names(ci.D)<-c("Lower","Upper")
    
    ans <- list(BF = BF, beta = post.mean.beta, CI.beta = ci.beta, CI.lambda = ci.lambda, CI.D = ci.D,freq=post.mean.freq)
    return(ans)
    rm(out) 
    rm(BDIC.out) 
    rm(freq.out) 
    rm(beta.out)
    rm(D.out)
    rm(lambda.out)
    rm(sigmas.out)
    
  } 
  
}
