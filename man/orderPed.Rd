\name{orderPed}
\alias{orderPed}
\title{Orders a pedigree}
\description{
  Orders a pedigree so that offspring follow parents.
}
\usage{
orderPed(ped)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ped}{\code{data.frame} with three columns: id,id parent1,id parent2 }
}
\value{
  numerical vector
}
\examples{
id <- 1:6
dam <- c(0,0,1,1,4,4)
sire <- c(0,0,2,2,3,5)
pedigree <- data.frame(id,dam,sire)
(ord <- orderPed(pedigree))
pedigree <- pedigree[6:1,]
(ord <- orderPed(pedigree))
pedigree <- pedigree[order(ord),]
pwrong <- pedigree
pwrong[1,2] <- pwrong[6,1]
}
\keyword{utilities}

