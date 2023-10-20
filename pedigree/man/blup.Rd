\name{blup}
\alias{blup}
\title{Function to calculate breeding values using an animal model}
\description{
  Fit an animal model to data, use a given variance ratio
  (\eqn{\alpha = \frac{\sigma^2_e}{\sigma^ 2_a}}). Calculate inverse of the
  additive genetic relationship matrix using function \code{makeInv()}
  of this package.
}
\usage{
blup(formula, ped, alpha,trim = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{formula of the model, do not include the random effect
    due to animal (generally ID).}
  \item{ped}{\code{data.frame} with columns corresponding to ID, SIRE,
    DAM and the columns in the formula.}
  \item{alpha}{Variance ratio (\eqn{\frac{\sigma^2_e}{\sigma^ 2_a}}).}
  \item{trim}{If \code{TRUE}, trims the pedigree using the available
    phenotype data using function \code{\link[pedigree]{trimPed}}.}
}
\value{
  Vector of solutions to the model, including random animal effects.}
\seealso{
  \code{\link[HaploSim]{SamplePedigree}}, \code{\link[pedigree]{gblup}},
   \code{\link[pedigree]{makeAinv}},\code{\link[pedigree]{blup}}
}
\examples{
example(gblup)
sol <- blup(P~1,ped = ped,alpha = 1/h2 - 1)
}
\keyword{utilities}
