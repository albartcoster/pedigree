\name{makeA}
\alias{makeA}
\title{Makes the A matrix for a part of a pedigree}
\description{
  Makes the A matrix for a part of a pedigree and stores it in a file called
  \code{A.txt}.}
  \usage{
makeA(ped,which)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ped}{\code{data.frame} with three columns: id,id parent1,id parent2 }
  \item{which}{Logical vector specifying between which indiduals
    additive genetic relationship is required. Goes back through the
    whole pedigree but only for subset of individuals.
  }
}
\value{
  Logical.
}
\examples{
id <- 1:6
dam <- c(0,0,1,1,4,4)
sire <- c(0,0,2,2,3,5)
ped <- data.frame(id,dam,sire)
makeA(ped,which = c(rep(FALSE,4),rep(TRUE,2)))
A <- read.table("A.txt")
}
\keyword{utilities}

