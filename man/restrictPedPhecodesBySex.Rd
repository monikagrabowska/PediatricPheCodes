\name{restrictPedPhecodesBySex}
\alias{restrictPedPhecodesBySex}
\title{Applies sex-specific phecode restrictions}

\description{
%%  Called within the createPedPhenotypes() function. This function applies the sex restrictions necessary to properly evaluate the sex-specific pediatric phecodes.}

\usage{
restrictPedPhecodesBySex(x)}

\arguments{
\item{phenotypes}{
%% A data frame with boolean phenotypes for individuals in the population. The first column should contain the individual IDs; the remaining columns should indicate case/control/exclude status for each of the pediatric phecodes.}

\item{id.sex}{
%% A data frame with two columns: id and sex (specified as “M” or “F”).}
}

\value{
%% Returns a data frame with all the same columns as the input data frame except for Sex. The data frame is also filtered to contain only elements of the selected sex in the function’s argument.}

\references{
\url{
%% https://github.com/PheWAS/PheWAS/blob/master/R/restrictPhecodesBySex.R}
\url{
%% https://github.com/PheWAS/PheWAS/blob/master/man/restrictPhecodesBySex.Rd}
}

\note{
%%  We modified the restrictPhecodesBySex() function from the original PheWAS R package to use 	the sex-specific pediatric phecodes by default.}
