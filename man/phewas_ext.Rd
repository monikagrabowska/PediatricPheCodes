\name{phewas_ext}
\alias{phewas_ext}
\title{PheWAS analysis using multiple methods}

\description{
%%This function performs a PheWAS analysis using multiple methods, with optional adjustment for user-specified covariates. By default, uses a generalized linear model.

\arguments{
\item{phenotypes}{
%%The names of the outcome variables (i.e., the phecodes).}

\item{genotypes}{
%%The names of the prediction variables (e.g., PLINK SNP genotype).}

\item{data}{
%%A data frame containing all of the variables to be used in the PheWAS analysis.}

\item{covariates}{
%%The names of the covariates.}

\item{cores}{
%%The number of cores to be used in the parallel socket cluster implementation.}

\item{additive.genotypes}{
%%Identifies if additive genotypes are being supplied (default is TRUE).}

\item{method}{
%%Determines the statistical method used to check associations (generalized linear model = “glm”, conditional logistic regression = “clogit”, likelihood ratio test = “lrt”, or Firth’s logistic regression = “logistf”).}

\item{strata}{
%%Name of the strata column (only necessary for conditional logistic regression).}

\item{factor.contrasts}{
%%Contrasts used to generate names for conditional logistic regression (by default, uses 	the contrasts returned by the contr.phewas() function).}

\item{return.models}{
%%Returns a list of the complete models (default is FALSE).}

\item{min.records}{
%%The minimum number of records necessary to perform a test. For logistic regression, there must be at least this number of cases and at least this number of controls. The default value is 20.}

\item{MASS.confint.level}{
%%Calculates a confidence interval at the specified level. Logistic models will report OR CIs. Confidence intervals are not calculated by default.}
}

\value{
%%Returns a data frame containing the results of the PheWAS analysis (e.g., a data frame with the 	following columns: phenotype, snp, covariates, beta, SE, OR, p, type, n_total, n_cases, 		n_controls, HWE_p, allele_freq, n_no_snp, formula, expanded_formula, note).}

\references{
\url{
%%https://github.com/PheWAS/PheWAS/blob/master/R/phewas_ext.R}
\url{
%%https://github.com/PheWAS/PheWAS/blob/master/man/phewas_ext.Rd}
}

\note{
%%Code was not modified from the phewas_ext() function in the original PheWAS R package.
}
