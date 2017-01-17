#' Inbreeding coefficient, initial estimator
#' @param s.a.f is the within-sample allele frequency, length M
#' @param a.f is the population-level allele frequency, length M
#' @param num.bins is the number of bins used in the estimate, deaults to 10.
#' @return an estimate of the sample's inbreeding coefficient
#' @examples calc.f.reg(sample.af,pop.af);
#' @export
#' @author Jack O'Brien, \email{jobrien@@bowdoin.edu}

calc.f.ini = function(s.a.f,a.f,num.bins =10)
{
	seq.bins = seq(0.05,0.5,length.out=num.bins+1);
	
	point.mat = matrix(rep(0,20),ncol=2);

	exp.het = a.f*(1-a.f)*2;
	obs.het = s.a.f*(1-s.a.f)*2;

	for(j in 2:length(seq.bins))
	{
		this.bin = which((exp.het <= seq.bins[j])&(exp.het > seq.bins[j-1]))
		mean.exp = mean(exp.het[this.bin]);
		mean.obs = mean(obs.het[this.bin]);
		point.mat[j-1,]=c(mean.exp,mean.obs);
	}

	f.val = 1-lm(point.mat[,2]~point.mat[,1]+0)[[1]][1]

	return(f.val)
}

#' Inbreeding coefficient, regressed estimator
#' @param s.a.f is the within-sample allele frequency, length M
#' @param a.f is the population-level allele frequency, length M
#' @return an estimate of the sample's inbreeding coefficient
#' @examples calc.f.reg(sample.af,pop.af);
#' @export
#' @author Jack O'Brien, \email{jobrien@@bowdoin.edu}

calc.f.reg = function(s.a.f,a.f)
{
	exp.het = a.f*(1-a.f)*2;
	obs.het = s.a.f*(1-s.a.f)*2;
	corr.2=lm(obs.het~exp.het+0)[[1]][1]
	f.val = 1-corr.2;

	return(f.val)
}

#' Inbreeding coefficient, direct estimator
#' @param s.a.f is the within-sample minor allele frequency, length M
#' @param a.f is the population-level minor allele frequency, length M
#' @return an estimate of the sample's inbreeding coefficient
#' @examples calc.f.fst(sample.maf,pop.maf);
#' @export
#' @author Jack O'Brien, \email{jobrien@@bowdoin.edu}

calc.f.fst = function(s.a.f,a.f)
{
	exp.het = a.f*(1-a.f)*2;
	obs.het = s.a.f*(1-s.a.f)*2;
	corr.het = 1-mean(obs.het)/mean(exp.het);
	return(corr.het)
}

#' MLE allele frequencies
#' @param data.set - number of ref and non read counts in array
#' @return A vector of allele frequencies
#' @examples calc.mle.p(abind(ref,non,along=3));
#' @author Jack O'Brien, \email{jobrien@@bowdoin.edu}
#' @export

calc.mle.p = function(data.set)
{
	
	M = dim(data.set)[2];
	N = dim(data.set)[3];

	ref = rowSums(data.set[,,1]);
	non = rowSums(data.set[,,2]);

	mle.vals = non/(ref+non);

	return(mle.vals)
}

#' Density of the Beta-binomial distribution using the p/f parameterization
#' @param K - number of non-reference read counts
#' @param N - total number of read counts
#' @param p - non-reference allele frequency
#' @param f - inbreeding coefficient
#' @return A log-probability of the density
#' @author Jack O'Brien, \email{jobrien@@bowdoin.edu}
#' @examples dbetabinom(10,100,0.1,0.3);
#' @export

dbetabinom.pf = function(K,N,p,f)
{
	 if (!requireNamespace("VGAM", quietly = TRUE)) {
	    stop("Pkg needed for this function to work. Please install it.",
	      call. = FALSE)
  	}
	alpha = (1-f)/f*p
	beta  = (1-f)/f*(1-p)
	prob =  dbetabinom.ab(K,N,alpha,beta,log=TRUE)
	return(prob)	
}

#' Internal MCMC function, initializes Markov chain 
#' @param data.set : an array of M x N x 2 of SNP read counts
#' @return An initialized state of the chain
#' @export
#' @author Jack O'Brien, \email{jobrien@@bowdoin.edu}

init.model = function(data.set)
{
	
	M = dim(data.set)[1];	
	N = dim(data.set)[2];

	f = rep(1/2,N);
	p = runif(M)
		
	mle.vals = calc.mle.p(data.set)
	p        = mle.vals;

	llk.mat  = array(NA,dim=c(M,N));
	model = list(M,N,f,p,mle.vals,llk.mat,0,0,0);
	names(model) = c("M","N","f","p","mles","llk.mat","llk","lpr","lprop")
	model$llk.mat = calc.llk.mat(model,data.set)

	model$llk = sum(model$llk.mat)
	model$lpr = calc.lpr(model);

	return(model)
	
}

#' Internal MCMC functions, calculates log-prior probability
#' @param Model iteration

calc.lpr = function(model)
{
	return(0);
}

#' Internal MCMC function, proposes new allele frequency value
#' @param model iteration and data set

propose.p=function(model,data.set,nu=20)
{
	new.model = model;
	snp = sample(1:model$M,1);

	mle.val = model$mles[snp];

	alpha = mle.val*nu;
	beta  = (1-mle.val)*(nu)

	new.model$p[snp] = rbeta(1,alpha,beta);

	new.model$lprop   	= dbeta(model$p[snp],alpha,beta,log=TRUE)-dbeta(new.model$p[snp],alpha,beta,log=TRUE);
	new.model$lpr     	= calc.lpr(new.model);
	new.model$llk.mat[snp,] = calc.llk.p.j(snp,new.model,data.set);
	new.model$llk     	= sum(new.model$llk.mat)

	return(new.model)
}

#' Internal MCMC function, proposes new f value
#' @param model iteration and data set

propose.f=function(model,data.set)
{
	new.model = model;
	sam = sample(1:model$N,1);

	alpha = beta = 1;

	new.model$f[sam] = rbeta(1,alpha,beta);
	
	new.model$lprop = 0;
	new.model$lpr   = calc.lpr(new.model);
	new.model$llk.mat[,sam] = calc.llk.f.i(sam,new.model,data.set)
	new.model$llk     	= sum(new.model$llk.mat)

	return(new.model)
}

#' Likelihood function
#' @param Model iteration
#' @examples
#' 	llk = calc.llk(model)
#' @export

calc.llk = function(model)
{
	return(sum(llk.mat))
}

#' Internal likelihood function
#' @param model, data set

calc.llk.mat=function(model,data.set)
{
	tmp = lapply(1:model$N,calc.llk.f.i,model,data.set)
	mat = matrix(unlist(tmp),ncol=model$N,byrow=FALSE)
	return(mat)
}

#' Internal likelihood function
#' @param element, model, data set

calc.llk.f.i = function(i,model,data.set)
{
	return(dbetabinom.pf(data.set[,i,2],data.set[,i,1]+data.set[,i,2],model$p,model$f[i]))
}

#' Internal likelihood function
#' @param element, model, data set

calc.llk.p.j = function(j,model,data.set)
{
	return(dbetabinom.pf(data.set[j,,2],data.set[j,,1]+data.set[j,,2],model$p[j],model$f))
}

#' Internal MCMC function
#' @param model iteration, data.set, p

propose.model = function(model,data.set,p=TRUE)
{
	test = runif(1);
	#model$M/(model$M+model$N)
	if (p==TRUE)
	{
		if (test < 0.2)
		{	new.model = propose.p(model,data.set)
		}else
		{	new.model = propose.f(model,data.set)	}
	}else
	{
		new.model = propose.f(model,data.set)
	}
	return(new.model)
}

#' MCMC function for simultaneous estimation of allele frequency and inbreeding coefficients from SNP read count data
#'
#' @param data.set : data set array of size M x N x 2 of read counts
#' @param num.iter : number of iterations, defaults to 1000
#' @param thin : the number of iterations to thin the MCMC by, defaults to 10
#' @param p : Boolean indicating if allele frequency is to be sampled, defaults to TRUE   
#' @return A list of states of the model for each (thinned) iteration of the MCMC
#' @keywords MCMC, f-statistic, inbreeding coefficient
#' @author Jack O'Brien, \email{jobrien@@bowdoin.edu}
#' @export
#' @examples 
#' chain = mcmc.fstat(data.set,1e4,100);
#' @export

mcmc.fstat = function(data.set,num.iter=1000,thin=10,p=TRUE)
{

	curr = init.model(data.set)

	chain = rep(list(),num.iter/thin)

	for(i in 1:num.iter)
	{
		prop = propose.model(curr,data.set,FALSE)
		llk.ratio = prop$llk - curr$llk;
		lprop.ratio     = prop$lprop;
		
		alpha = llk.ratio+lprop.ratio;

		if (log(runif(1))<alpha)
		{
			curr = prop;
		}

		if (i %% thin == 0)
		{
			print(i/thin)
			chain[[i/thin]] = curr;
		}
	}		
	return(chain)
}


