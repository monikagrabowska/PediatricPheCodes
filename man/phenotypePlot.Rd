\name{phenotypePlot}

\alias{phenotypePlot}

\description{

%%Called within the plotManhattan() function that generates a PheWAS Manhattan plot with phenotypes on the x-axis and p-values on a -log10 scale on the y-axis. This function enables some additional plotting options compared to the plotManhattan() function.

}

\arguments{

\item{d}{

%%The data frame passed into the plotManhattan() function, which is typically the data 	frame of PheWAS results returned by the phewas_ext() function.

}

\item{suggestive.line}{

%%The p-value at which to draw a blue line to show a significance threshold. The default value is 0.05. An NA value results in no line.

}

\item{significant.line}{

%%The p-value at which to draw a red line to show an adjusted significance threshold. The default applies a Bonferroni correction. An NA value results in no line.

}

\item{annotate.phenotype}{

%%Determines whether phenotype descriptions (corresponding to phecodes with p < annotate.level or phecodes specified in annotate.list) are annotated on the Manhattan plot. The default is TRUE.

}

\item{annotate.angle}{

%%The angle used for the annotation text (default = 0).

}

\item{annotate.size}{

%%The size of the annotation text (default = 3).

}

\item{annotate.level}{

%%The p-value threshold for annotating phenotypes. The default value is that provided in the significant.line argument.

}

\item{annotate.list}{

%%A character vector of phecodes to be annotated on the plot regardless of significance.

}

\item{lc.labels}{

%%Determines if the labels should be displayed in lower case (default is FALSE).

}

\item{y.axis.interval}{

%%The interval used for the y-axis of the Manhattan plot (default value = 5).

}

\item{size.x.labels}{

%%The size of the x-axis labels (default = 9).

}

\item{size.y.labels}{

%%The size of the y-axis labels (default = 9).

}

}

\value{

%% This function returns a phenotype-significance Manhattan plot to its parent function

}

\references{

\url{

%%https://github.com/PheWAS/PheWAS/blob/master/R/phenotypePlot.R

}

\url{

%%https://github.com/PheWAS/PheWAS/blob/master/man/PheWAS_Plotting.Rd

}

}

\note{

%% Code was only minimally modified from the phenotypePlot() function in the original PheWAS R 	package to use the pediatric phecode descriptions and categorization by default.

}
