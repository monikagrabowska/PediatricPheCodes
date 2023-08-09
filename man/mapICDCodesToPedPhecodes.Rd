\name{mapICDCodesToPedPhecodes}
\alias{mapICDCodesToPedPhecodes}
%- Also NEED an '\alias' for EACH other topic documented here.

\description{
%%  Called within the createPedPhenotypes() function. This function takes the ICD-9-CM and ICD-10-CM codes in the id.vocab.code.index data frame that was passed into createPedPhenotypes() and maps them to pediatric phecodes. It returns a data frame with three columns: the id and index columns from the original id.vocab.code.index data frame, along with a new phecode column.
}
\usage{
mapICDCodesToPedPhecodes(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{

  %%The id.vocab.code.index data frame originally passed into the createPedPhenotypes() function.

}
}
\value{
%%  This function returns a dataframe with three columns, 1. the original ICD codes, 2. The original index column, 3. A new column representing the corresponding phecode for each icd code.
}
\references{
  \url{

  %%https://github.com/PheWAS/PheWAS/blob/master/R/mapCodesToPhecodes.R

  }

  \url{

  %%https://github.com/PheWAS/PheWAS/blob/master/man/mapCodesToPhecodes.Rd

}
}
\note{
%%  Code not modified from the mapCodesToPhecodes() function in the original PheWAS R package.
}


