#' aziztest: A package for finding associations in heterogeneous setting (aberration enrichment)
#'
#' This package contains the statistical test presented in Mezlini et al. (2020)
#'   "Finding associations in a heterogeneous setting: Statistical test for aberration enrichment"
#'   <https://www.biorxiv.org/content/10.1101/2020.03.23.002972v2>.
#'   It is used to detect associations that are beyond the broad pattern of comparing
#'   averages between all cases and all controls. Instead it looks for a heterogeneous association
#'   where only some of the cases present the signal of interest while the majority
#'   are indistinguishable from controls.
#'   For example, in a clinical trial setting our test can be used to assess
#'   treatment efficacy in a context of heterogeneous treatment effect, where
#'   the drug works well on only some of the patients.
#'   Another usage example is in -Omics data where a relevant gene's dysregulation is
#'   present in only some of the disease cases.
#'
#' The main function is the \code{aziz.test()} function used to test for
#' heterogeneous associations/ aberration enrichment.
#'
#' @section aziztest extra functions:
#' If you are testing multiple variables at once (such as all genes in a gene expression dataset),
#' you can store the results in a list and the reformat it into an easy to use data.frame using
#' the function \code{reformat_results()}.
#'
#' In the context of a large number of variables, calibration can be used to speed up p-value calculation
#' with functions \code{calibrate_test()} and \code{get_calibrated_pvalues()}.
#'
#' @docType package
#' @name aziztest
NULL
