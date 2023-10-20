\name{add.Inds}
\alias{add.Inds}
\title{Function to add missing individuals to a pedigree}
\description{
  Function add.Inds() adds missing individuals to a pedigree and returns the
  complete pedigree as a data.frame with the same headers as the
  original pedigree. Remeber to check for errors beforehand with
  function \code{errors.ped}. Unknown parents should be coded as NA.
}
\usage{
add.Inds(ped)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ped}{\code{data.frame} with three columns: id,id parent1,id parent2 }
}
\value{
  data.frame of three columns with identical header as input.
}
\author{Albart Coster, Albart.Coster@wur.nl}
\seealso{ \code{\link{orderPed}}}
\examples{
ID <- 3:5
DAM <- c(1,1,3)
SIRE <- c(2,2,4)
pedigree <- data.frame(ID,DAM,SIRE)
pedigree <- add.Inds(pedigree)
}
\keyword{utilities}
