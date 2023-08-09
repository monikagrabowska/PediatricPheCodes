\name{createPedPhenotypes}
\alias{createPedPhenotypes}
%- Also NEED an '\alias' for EACH other topic documented here.
\description{
%%  This function translates ICD-9-CM and ICD-10-CM codes into pediatric phecodes. It takes in a data frame with four columns: id, vocabulary_id (specified as “ICD9CM” and/or “ICD10CM”), code, and index. The function returns a data frame with boolean phenotypes for the individuals in the input file – the first column contains the individual IDs and the remaining columns indicate case/control/exclude status for each of the pediatric phecodes (case = TRUE, control = FALSE, exclude = NA).
}
\usage{
createPedPhenotypes(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{id.vocab.code.index}{

%% A data frame with four columns: id, vocabulary_id, code, and index. The vocabulary_id must be specified as “ICD9CM” or “ICD10CM”. The code column must be a character column to ensure that all ICD codes are properly mapped to phecodes. Leading and trailing zeros should be preserved for ICD codes in the code column. The index can be a code count or a diagnosis date. If dates are supplied, only distinct ICD code/date pairs are counted for an individual.

}

\item{id.sex}{

%% A data frame with two columns: id and sex (specified as “M” or “F”). Needed for proper evaluation of the sex-specific pediatric phecodes.

}

\item{min.code.count){

%% The minimum ICD code count for an individual to be considered a case 				(default = 2).

}

\item{full.population.ids} {

%% A list of IDs for the complete population. This is needed to ascertain control status for individuals with no ICD codes in the input data frame. By default, uses the number of unique person IDs in the id.vocab.code.index data frame.

}
}
\details{

}
\value{
%%  The function returns a data frame with boolean phenotypes for the individuals in the input file – the first column contains the individual IDs and the remaining columns indicate case/control/exclude status for each of the pediatric phecodes (case = TRUE, control = FALSE, exclude = NA).
}
\references{
 \url{

  %%https://github.com/PheWAS/PheWAS/blob/master/R/createPhenotypes.R

  }

  \url{

  %%https://github.com/PheWAS/PheWAS/blob/master/man/createPhenotypes.Rd

  }
}

