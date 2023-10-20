\name{calcG}
\alias{calcG}
\title{Function to calculate a relationship matrix from marker data
  (usually allele count data), G matrix.}
\description{
  Function to calculate a relationship matrix from marker data. Option
  to return the inverse of matrix. Inverse calculated using
  \code{Matrix} package.
}
\usage{
calcG(M, data = NULL,solve = FALSE)
}
\arguments{
  \item{M}{Matrix of marker genotypes, usually the count of one of the
    two SNP alleles at each markers (0, 1, or 2).}
  \item{data}{Optional logical vector which can tell of which
    individuals we have phenotypes.}
  \item{solve}{Logic, if TRUE then function returns the inverse of the
    relationship matrix.}
}
\value{
  Matrix of class \code{dgeMatrix}.}
\seealso{
  \code{\link[HaploSim]{SamplePedigree}}, \code{\link[pedigree]{gblup}},
   \code{\link[pedigree]{makeAinv}},\code{\link[pedigree]{blup}}
}
\examples{
example(gblup)
G <- calcG(M)
Ginv <- calcG(M,solve = TRUE)
}
\keyword{utilities}
