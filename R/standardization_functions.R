
#' Standardizing the enrichment score. Main function
#'
#' Divide by sd under the null and optionally trim the score vector
#' Trimming is 3/4 by default. Just to avoid potential instability at the end
#' @noRd
estransformer=function(x0,esnorm,trimmed=0.75){#only return the 3/4 by default
  x = (x0/esnorm$sd)
  x = x[1:round(trimmed*length(x))];
  return(x)
}

#' Standardizing the enrichment score. Estimating null distribution
#'
#' \code{weighting} corrects for case-control imbalance within the steps
#' @noRd
generate_es_normalizer_generic=function(w,pheno,approx=1,weighting=TRUE){
  if (approx==1 & weighting) return(generate_es_normalizer(w,pheno))
  if (approx==1) return(generate_es_normalizer_old(w,pheno))
  return(generate_es_normalizer2(w,pheno))
}

#' Standardizing the enrichment score. Estimating null distribution
#'
#' \code{weighting} corrects for case-control imbalance within the steps
#' @noRd
generate_es_normalizer=function(w,pheno){ #w is the ranked expression
  p = length(which(pheno==1))/length(pheno)
  cumsumw = cumsum(w)
  cumsumwsq = cumsum(w^2);
  m = 0
  s = ( (1/length(which(pheno==1))) + (1/length(which(pheno==0))) )*sqrt(p*(1-p)/(length(pheno)-1))*sapply(1:length(w),function(n)sqrt( length(pheno)*cumsumwsq[n] - (cumsumw[n])^2 ))
  return(list(mean=m,sd=s))
}

#' Standardizing the enrichment score. old function. Obsolete now
#'
#' This is used if \code{weighting} is False. Mean under the null is non-zero
#' @noRd
generate_es_normalizer_old=function(w,pheno){ #w is the ranked expression
  p = length(which(pheno==1))/length(pheno)
  cumsumw = cumsum(w)
  cumsumwsq = cumsum(w^2);
  m = cumsumw*(2*p -1)
  s = 2*sqrt(p*(1-p)/(length(pheno)-1))*sapply(1:length(w),function(n)sqrt( length(pheno)*cumsumwsq[n] - (cumsumw[n])^2 ))
  return(list(mean=m,sd=s))
}

#' Standardizing the enrichment score. obsolete function
#'
#' Estimates variance under the null by permutations. Obsolete
#' @importFrom stats sd
#' @noRd
generate_es_normalizer2=function(w,pheno,nperm=2000){
  p = length(which(pheno==1))/length(pheno);
  cummeanw = cumsum(w)/(1:length(w));
  m = cummeanw*(2*(1:length(w))*p -(1:length(w)))
  resp = mapply(function(k)cumsum(2*(sample(pheno,length(pheno))-0.5)*w),1:nperm)
  s = rep(1,length(w));
  s[1:(length(w)/2)] = apply(resp[1:(length(w)/2),],1,sd)
  m = apply(resp,1,mean)
  return(list(mean=m,sd=s))
}

#' Unused. What if we try chisq at every position and take the best
#'
#' Unused. What if we try chisq at every position and take the best
#' @importFrom stats chisq.test
#' @noRd
best_chisq=function(pheno,e){#Unused. What if we try chisq at every position and take the best
  tp = 1; js = 1;
  for (j in 10:(length(pheno)-10)){
    ca = sum(pheno[order(e)][1:j])
    t = chisq.test(cbind(c(ca,length(which(pheno==1))-ca),c(j-ca,length(which(pheno==0))-j+ca)))$p.value
    if (t<tp){tp=t; js=j;}
  }
  print(paste(js,tp))
}

