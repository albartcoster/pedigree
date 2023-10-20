\name{calcInbreeding}
\alias{calcInbreeding}
\title{Calculates inbreeding coefficients for individuals in a pedigree.}
\description{
  Calculates inbreeding coefficients of individuals in a pedigree.}
\usage{
calcInbreeding(ped)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ped}{\code{data.frame} with three columns: id,id parent1,id parent2}
}
\value{
  Logical.
}
\examples{
id <- 1:6
dam <- c(0,0,1,1,4,4)
sire <- c(0,0,2,2,3,5)
ped <- data.frame(id,dam,sire)
(F <- calcInbreeding(ped))
}
\keyword{utilities}

