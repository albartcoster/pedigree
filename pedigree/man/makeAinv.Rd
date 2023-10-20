\name{makeAinv}
\alias{makeAinv}
\title{Makes inverted A matrix for a pedigree}
\description{
  Makes inverted A matrix for a pedigree and stores it in a file called
  \code{Ainv.txt}.}
  \usage{
makeAinv(ped)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ped}{\code{data.frame} with three columns: id,id parent1,id parent2 }
}
\value{
  Logical.
}
\examples{
id <- 1:6
dam <- c(0,0,1,1,4,4)
sire <- c(0,0,2,2,3,5)
ped <- data.frame(id,dam,sire)
makeAinv(ped)
Ai <- read.table("Ainv.txt")
nInd <- nrow(ped)
Ainv <- matrix(0,nrow = nInd,ncol = nInd)
Ainv[as.matrix(Ai[,1:2])] <- Ai[,3]
dd <- diag(Ainv)
Ainv <- Ainv + t(Ainv)
diag(Ainv) <- dd
}
\keyword{utilities}

