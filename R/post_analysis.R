
#' Reformats multiple results into one table
#'
#' If running multiple variables (and storing in a list),
#' this transform the list of results to one coherent data.frame
#'
#' @param res_esa listed outputs of multiple calls to \code{aziz.test()} on multiple variables
#' @return A data.frame containing all results in an accessible presentation
#' @export
reformat_results = function(res_esa){
  es1 = rep(0,length(res_esa));
  es_sign = rep(-1,length(res_esa));
  es_indmax = rep(-1,length(res_esa));
  r = rep(1,length(res_esa)); or = r; ncas = r-1; pval = r;#initializing
  for (i in 1:length(res_esa)){
    es1[i] = (res_esa[[i]])$es;
    es_sign[i] = (res_esa[[i]])$direction;
    es_indmax[i] = (res_esa[[i]])$esind[es_sign[i]];
    r[i]    = (res_esa[[i]])$oddcas; or[i] = (res_esa[[i]])$oddratio;
    ncas[i] = (res_esa[[i]])$ncas; pval[i] = (res_esa[[i]])$pval;
  }#WARNING the direction attribute was not previously returned
  return(data.frame(es=es1,es_sign=es_sign,es_indmax=es_indmax,r=r,or=or,ncases=ncas,pval=pval))
}

#' Print a formatted version of the result details
#'
#' Print a formatted version of the result details
#' @param x output of the \code{aziz.test()} function
#' @seealso \code{\link{print_summary}}
#' @export
print_details=function(x){
  print(paste("indmax:",x$esind[x$direction],"; cases:",x$ncas,"; odds ratio:",format(x$oddratio,digits=3)))
}

#' Print a formatted version of the results
#'
#' Print a formatted version of the results for clarity
#'
#' @param x output of the \code{aziz.test()} function
#' @seealso \code{\link{print_details}} for more info
#' @export
print_summary=function(x){
  pvalprint=as.character(x$pval); if (x$pval< 1/length(x$perm))pvalprint="< 1/#permutations"
  print(paste("Score:",format(x$es,digits=3),"; direction:",x$direction,
              "; r:",format(x$oddcas,digits=3),"; pval:",pvalprint))
}

# -Return original value corresponding to max ES
#'
#' Return original value corresponding to max ES
#'
#' @noRd
aberrant_expression=function(finalexpri,di,esmi){
  if (di==1) return(finalexpri[order(finalexpri)[esmi]]);
  return(finalexpri[order(finalexpri,decreasing = T)[esmi]])
}

# -Returns which samples are in the aberrant interval
#'
#' Returns index of samples that are in the aberrant interval
#'
#' @noRd
which_aberrant_internal=function(finalexprj,finalexpri,di,esmi){
  ipoint = aberrant_expression(finalexpri,di,esmi)
  if (di==1) return(which(finalexprj<= ipoint)); #direct sense
  if (di==2) return(which(finalexprj>= ipoint)); #antisense
  return(NULL)
}

#' Returns which samples are in the aberrant interval defined by \code{test.aziz()}
#'
#' Returns index of samples that are in the aberrant interval defined by \code{test.aziz()}
#' Can take the same samples or a new set of previously unseen samples
#' @param xi Numerical vector. Can be the same as the tested variable x
#'   or it can be a new set of unseen samples.
#' @param x Numerical vector of the variable tested by \code{test.aziz()}
#' @param res Result of running \code{test.aziz()}
#' @return indexes of samples in \code{xi} that are within the aberrant interval
#' @examples
#' y = c(rep(1,200),rep(0,200))
#' x = rnorm(400)
#' #Inducing an aberration enrichment signal by perturbing some of the cases
#' x[1:20]=x[1:20]-3;
#' res2 = aziz.test(y,x,rep=20000)
#' print_summary(res2)
#' which_aberrant(x,x,res2)
#' which_aberrant(c(-5,1.5,-2.5,-0.5,2),x,res2)#testing if new values are within the aberrant interval
#' @export
which_aberrant=function(xi,x,res){
  return(which_aberrant_internal(xi,x,res$direction,res$esind[res$direction]))
}

#' Polynomial prediction of the mapping from test statistic to p-value
#'
#' Used to predict p-values for associations that received a zero p-value from permutations.
#' Works on pvalues under the inverse of the number of permutations.
#' @importFrom stats lm
#' @importFrom stats predict
#' @noRd
correct_zero_pvalues = function(trainx,trainy,es1,esa_m,accuracythresh,maxpvaltrain=0.05){
  regdata = data.frame(x=trainx,y=trainy);
  ind = which(regdata[,2]>= log(accuracythresh) & regdata[,2]< log(maxpvaltrain));
  res_lm = lm(y~ poly(x,3),data=regdata[ind,]);
  predict_pval = function(esx)exp(predict(res_lm,data.frame(x=esx)));
  esa = predict_pval(es1);
  esa_cor = esa_m;
  toreplace = which(esa_m<accuracythresh);
  if (length(toreplace)) esa_cor[toreplace] = sapply(esa[toreplace],function(x)min(x,accuracythresh))
  return(esa_cor);
}

