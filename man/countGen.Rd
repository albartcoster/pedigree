\name{countGen}
\alias{countGen}
\title{Count generation number for each individual in a pedigree.}
\description{Counts generation number for individuals in a pedigreee.
  }
\usage{
countGen(ped)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ped}{\code{data.frame} with three columns: id,id parent1,id parent2}
}
\value{
  Numeric vector
}
\examples{
id <- 1:5
dam <- c(0,0,1,1,4)
sire <- c(0,0,2,2,3)
ped <- data.frame(id,dam,sire)
(gens <- countGen(ped))
}
\keyword{utilities}

