% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Coxmos_multiblock_functions.R
\name{deleteNearZeroCoefficientOfVariation.mb}
\alias{deleteNearZeroCoefficientOfVariation.mb}
\title{deleteNearZeroCoefficientOfVariation.mb}
\usage{
deleteNearZeroCoefficientOfVariation.mb(X, LIMIT = 0.1)
}
\arguments{
\item{X}{List of numeric matrices or data.frames. Explanatory variables. Qualitative variables must
be transform into binary variables.}

\item{LIMIT}{Numeric. Cutoff for minimum variation. If coefficient is lesser than the limit, the
variables are removed because not vary enough (default: 0.1).}
}
\value{
A list of three objects.
\code{X}: A list with as many blocks as X input, but with the variables filtered.
\code{variablesDeleted}: A list with as many blocks as X input, with the name of the variables that have been removed.
\code{coeff_variation}: A list with as many blocks as X input, with the coefficient of variation per variable.
}
\description{
Filters out variables from a dataset that exhibit a coefficient of variation below a
specified threshold, ensuring the retention of variables with meaningful variability.
}
\details{
The \code{deleteNearZeroCoefficientOfVariation} function is a pivotal tool in data preprocessing,
especially when dealing with high-dimensional datasets. The coefficient of variation (CoV) is a
normalized measure of data dispersion, calculated as the ratio of the standard deviation to the mean.
In many scientific investigations, variables with a low CoV might be considered as offering limited
discriminative information, potentially leading to noise in subsequent statistical analyses. By
setting a threshold through the \code{LIMIT} parameter, this function provides a systematic approach to
identify and exclude variables that do not meet the desired variability criteria. The underlying
rationale is that variables with a CoV below the set threshold might not contribute significantly
to the variability of the dataset and could be redundant or even detrimental for certain analyses.
The function returns a modified dataset, a list of deleted variables, and the computed coefficients
of variation for each variable. This comprehensive output ensures that researchers are well-informed
about the preprocessing steps and can make subsequent analytical decisions with confidence.
}
\examples{
data("X_multiomic")
X <- X_multiomic
filter <- deleteNearZeroCoefficientOfVariation.mb(X, LIMIT = 0.1)

}
\author{
Pedro Salguero Garcia. Maintainer: pedsalga@upv.edu.es
}
