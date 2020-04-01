# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

#' Calibration of pvalues in a slow multi-hypothesis setting
#'
#' Compute the null distribution of test statistics on one gaussian variable.
#' Useful if testing a large number of variables at once since it allows
#'   running permutations only once behorehand rather than for every variable.
#'   Used in conjunction with "\code{get_calibrated_pvalues}"
#' @importFrom stats rnorm
#' @inheritParams aziz.test
#' @param doall Default=TRUE. All permutations are performed
#' @return A vector that can be used in \code{get_calibrated_pvalues()}
#' @seealso \code{\link{get_calibrated_pvalues}}, \code{\link{aziz.test}}
#' @examples
#' y = c(rep(1,200),rep(0,200))
#' x = rnorm(400)
#' calibration = calibrate_test(y,rep=100)
#' es = aziz.test(y,x,rep=0)$es #No need for permutations, pvalues computed from calibration
#' get_calibrated_pvalues(calibration,es)
#'
#' @export
calibrate_test = function(y,w=NULL,rep=10000000,doall=TRUE,unidirectional=0,flatten=0.5,ignoremax=0){
  e = rnorm(length(y));
  res_m <- aziz.test(y,e,w,rep=rep,doall=doall,unidirectional=unidirectional,
                     flatten=flatten,ignoremax=ignoremax);
  return(sort(res_m$perm));
}

#' Use the calibration of pvalues in a slow multi-hypothesis setting
#'
#' Compute the pvalues from a single set of permutations obtained from \code{calibrate_test}.
#' Useful if testing a large number of variables at once since it allows
#'   running permutations only once behorehand rather than for every variable.
#'   Used in conjunction with "\code{calibrate_test}"
#'
#' @param calibration Output of function \code{calibrate_test()}
#' @param es1 Max Enrichment score given by function \code{aziz.test()} $es.
#'   A vector containing the max enrichment scores from many variables is acceptable
#' @return calibrated pvalue(s) corresponding to the max enrichment score(s) given
#' @examples
#' y = c(rep(1,200),rep(0,200))
#' x = rnorm(400)
#' calibration = calibrate_test(y,rep=100)
#' es = aziz.test(y,x,rep=0)$es #No need for permutations, pvalues computed from calibration
#' get_calibrated_pvalues(calibration,es)
#' @seealso \code{\link{calibrate_test}}, \code{\link{aziz.test}}
#' @export
get_calibrated_pvalues = function(calibration,es1){
  return(1- findInterval(es1,calibration)/length(calibration));
  #return(1- findInterval(es1,calibration,left.open=TRUE)/length(calibration)); #Might be better but does not work on hpf
}


#' Statistical test for heterogeneous effects
#'
#' Main function running the statistical test looking for heterogeneous effects/
#'   aberration enrichment. Takes a vector of case/control labels (\code{y}) and a vector
#'   of numeric measurements (\code{x}) to be tested for association with case/control status.
#'   For example, in a clinical trial setting \code{y} can indicate individuals
#'   on a drug vs placebo and \code{x} can be a change in disease severity measurement
#'   from baseline. This test will return a pvalue indicating drug efficacy and is more
#'   powerful than other test in a heterogeneous effects setting.
#'   Another usage example is in -omics data where \code{y} would indicate
#'   disease vs healthy control and \code{x} could be a gene's expression vector accross samples.
#'
#' @param y A binary vector of sample labels (cases=1, controls=0).
#' @param x A numerical vector. Variable tested for association. Preferably continuous
#' @param w Default = NULL. Optional numerical vector of weights.
#'   1 means all weights are equal to 1 and only the ordering is considered.
#'   If NULL (default), a standardisation of x is used to calculate the weights
#'   giving larger weights to aberrations of larger magnitude.
#' @param rep Default=100000. Number of permutations to be used to calculate p-values.
#' @param doall Default=FALSE. Logical. If TRUE all \code{rep} permutations are performed.
#'   If FALSE only enough permutations are performed to get accurate pvalues.
#'   Variable that are clearly not associated need only a 100 permutations.
#' @param eps Default = 0.000000001. Small numeric value. Standard deviation of
#'   the gaussian node added to x before ordering samples. In the case of equalities,
#'   this ensures the ordering is not biased. Adjust lower if x has low variability.
#' @param unidirectional Default = 0. Can be 0, 1 or -1. 0 is for  testing both
#'   directions of effect. 1 is for testing cases<controls and -1 is for
#'   testing cases>controls.
#' @param flatten Default = 0.5. Numeric value recommended between 0 and 1.
#'   If weights are not given, we take the max of flatten and the absolute
#'   value of the Z-score of \code{x} as the weights (Default behavior).
#' @param ignoremax Default=0. Optional value indicating if we should ignore
#'   the first few values when selecting the maximal enrichment score.
#' @return A result object with the following fields: (for clarity use \code{\link{print_summary}})
#' \describe{
#' \item{es}{Max enrichment score.}
#' \item{pval}{Permutation pvalue, if permutations were performed.}
#' \item{oddcas}{Proportion of cases in the aberrant interval driving the max enrichment score.
#'   This is described as the proportion r in the main paper.}
#' \item{direction}{direction of the effect. 1: cases<controls, 2: cases>controls.}
#' \item{oddratio}{Odds ratio of being in the aberrant interval for cases/controls.
#'   Equal to \code{oddcas} divided by the same calculation on controls.}
#' }
#' Other info fields (Can be useful ):
#' \describe{
#' \item{esm}{Max enrichment score in both directions.}
#' \item{esind}{Index of the Max enrichment score in both directions.
#'   can also be interpreted the number of samples in the aberrant interval.}
#' \item{ncas}{Number of cases in the aberrant interval.}
#' \item{escurve}{A vector of the computed standardized enrichment scores at all positions.}
#' \item{perm}{A vector of all max enrichment scores obtained in permutations.}
#' }
#' @examples
#' y = c(rep(1,200),rep(0,200))
#' x = rnorm(400)
#'
#' res = aziz.test(y,x,rep=100) #run 100 permutations to calculate pvalue
#' print_summary(res)
#'
#' #Inducing an aberration enrichment signal by perturbing some of the cases
#' x[1:20]=x[1:20]-3;
#' res2 = aziz.test(y,x,rep=100)
#' print_summary(res2)
#' @export
aziz.test=function(y,x,w=NULL,rep=100000,doall=FALSE,eps=0.000000001,
                   unidirectional=0,flatten=0.5,ignoremax=0){
  checked=T;
  if (length(y)<20 | length(x)<20) {print("Error : Sample size too smalll"); checked=F; }
  if (checked & length(y)!=length(x)) {print("Error : x and y should be of the same size"); checked=F; }
  if (checked & !is.numeric(x)) {print("Error : x should be numeric"); checked=F; }
  if (checked & !is.numeric(y)) {print("Error : y should be a numeric vector of 0s and 1s"); checked=F; }
  if (checked & !all(y %in% c(0,1))) {print("Error : y should only contain 0 and 1"); checked=F; }
  if (checked & ((sum(y==0)<2) | (sum(y==1)<2) )) {print("Error : y should contain multiple 0 and multiple 1"); checked=F; }
  tab=table(x);
  if (checked & length(tab)<2 ) {print("Error : x should not be constant"); checked=F; }
  if (checked & max(tab)>1 ) {print("Warning : Duplicates are present in x. x should preferably be continuous."); }

  res=NULL;
  if (checked) {
    res=aziz.test.main(y,x,w=w,rep=rep,doall=doall,eps=eps,
                       unidirectional=unidirectional,flatten=flatten,ignoremax=ignoremax)
  }
  return(res)
}

# -Full aziz.test function after checking inputs validity
#' Statistical test for heterogeneous effects
#'
#' Main function running the statistical test looking for heterogeneous effects/
#'   aberration enrichment. Takes a vector of case/control labels and a vector
#'   of measurements (numeric) to be tested for association.
#'
#' @noRd
aziz.test.main=function(pheno,pred,w=NULL,rep=100000,doall=FALSE,eps=0.000000001,unidirectional=0,
                   flatten=0.5,ignoremax=0,approx=1,weighting=TRUE,trimmed=0.75){# negative pred are used as weights by default
  prep = preprocess_weights(pheno,pred,w=w,eps=eps,weighting=weighting,flatten=flatten)
  ord  = prep$ord; signmult = prep$signmult; w = prep$w;
  esnorm    = generate_es_normalizer_generic(w[ord],pheno,approx,weighting=weighting);
  esnormrev = generate_es_normalizer_generic(w[rev(ord)],pheno,approx,weighting=weighting);
  esa0 = function(s,o,norm)estransformer(cumsum( (s*w)[o]),norm,trimmed=trimmed); #Enrichment score is computed and standardized here.
  esa_pos = function(s)(as.numeric(esa0(s,ord,esnorm))); #Positive direction (underexpression)
  esa_neg = function(s)(as.numeric(esa0(s,rev(ord),esnormrev))); #Negative direction
  esa_pos_m = function(s)max(esa_pos(s));esa_neg_m = function(s)max(esa_neg(s));
  if (ignoremax>0){ esa_pos_m=function(s)max(esa_pos(s)[-(1:ignoremax)]); esa_neg_m=function(s)max(esa_neg(s)[-(1:ignoremax)]); }
  esa_2d = function(s)max(c(esa_pos_m(s),esa_neg_m(s))); #Max of both directions
  esa_applied = switch(as.character(unidirectional),'1'=esa_pos_m,'-1'=esa_neg_m,esa_2d)#both directions by default
  #esa_applied = esa_2d; if (unidirectional==1)  esa_applied = esa_pos_m; if (unidirectional==-1) esa_applied = esa_neg_m;
  esperm = permute(rep,signmult,esa_applied,doall);

  dir1 = esa_pos(signmult);dir2 = esa_neg(signmult);
  esm = c(max(dir1),max(dir2)); esind = c(which.max(dir1),which.max(dir2));
  direction = 1; escurve = dir1; if(esm[2]>esm[1]) {direction=2; escurve=dir2; } # WARNING The ignoremax parameter could cause discrepancies here
  vip = 1:esind[direction]; if (direction==2) vip = (length(pheno)+1-esind[direction]):length(pheno);
  ncas = sum(pheno[ord][vip]); nctr = esind[direction]-ncas;
  ncastotal = sum(pheno); nctrtotal = length(pheno)-ncastotal;
  oddcas = ncas/ncastotal; oddctr = nctr/nctrtotal; oddratio = oddcas/oddctr;
  #summary=paste("Score:",format(esperm$real,digits=3),"direction:",direction,"r:",format(oddcas,digits=3),"pval:",esperm$pval)
  #details=paste("indmax:",esind[direction],"cases:",ncas,"odd ratio:",format( oddratio ,digits=3))
  return(list(es=esperm$real,esm=esm,esind=esind,pval=esperm$pval,direction=direction,escurve=escurve,
              vip=vip,perm=esperm$perm,ncas=ncas,oddratio=oddratio,oddcas=oddcas))
}

#' Preprocessing to calculate weights
#'
#' If weights are not given (default), the x vector is standardized
#'It is also flattened so that no value is below \code{flatten}
#' @importFrom stats rnorm
#' @noRd
preprocess_weights=function(pheno,pred,w=NULL,eps=0.000000001,weighting=TRUE,flatten=0){
  if (eps>0) pred = scale(pred+rnorm(length(pred),0,eps));#Add noise to avoid equal values biasing the ranking in case of equality
  if (length(w)==0) { w = abs(pred); if(flatten) w[which(w<flatten)] = flatten;}
  if (length(w)==1) w = rep(1,length(pred));
  ord = order(pred,decreasing=FALSE);
  ind1 = which(pheno==1);ind0 = which(pheno==0);
  signmult = pheno; signmult[ind0] = -1;
  if(weighting) { signmult[ind1] = 1/length(ind1); signmult[ind0] = -1/length(ind0);}#Weighting based on proportion of cases and controls. Makes the ES mean 0 under null
  return(list(ord=ord,w=w,signmult=signmult))
}


#' Run the permutations
#'
#' It will do all permutations if doall=T,
#'   or it will gradually increase permutations by a factor of 10 until
#'   reaching an acceptable accuracy at estimating pval or max permutations
#' @noRd
permute=function(rep, phen,f,doall=FALSE){#rep assumed >100 or 0
  real = f(phen)
  if (rep){
    lrep = 100; perm = sapply(1:lrep,function(x){f(sample(phen,length(phen)))});
    pval = sum(perm>=real)/length(perm);
    lrepall = c(1000,10000,100000,1000000,10000000); lrepall = lrepall[which(lrepall<rep)]
    for (lrep in lrepall) {
      if (!doall & pval<(100/lrep)){#pval < 0.1 (100/1000)  -> Do lrep=1000
        perm = sapply(1:lrep,function(x){f(sample(phen,length(phen)))});
        pval = sum(perm>=real)/length(perm);
      }
    }
    if ((pval< 100/rep & rep>100)| doall){
      perm = sapply(1:rep,function(x){f(sample(phen,length(phen)))})
      pval = sum(perm>=real)/length(perm);
    }
  }else {pval = 1; perm=c(); }
  return(list(real=real,perm=perm,pval=pval))
}

