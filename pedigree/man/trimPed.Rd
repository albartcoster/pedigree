\name{trimPed}
\alias{trimPed}
\title{Function to trim a pedigree based on available data}
\description{
  Trims a pedigree given a vector of data. Branches without data are
  trimmed off the pedigree.
}
\usage{
trimPed(ped, data,ngenback = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ped}{\code{data.frame} with three columns: id,id parent1,id parent2}
  \item{data}{TRUE-FALSE vector. Specifies if data for an individual is
    available.}
  \item{ngenback}{Number of generations back. Specifies the number of
    generations to keep before the individuals with data.}
}
\value{
  Logical vector specifying if an individual should stay in the
  pedigree.}
\examples{
id <- 1:5
dam <- c(0,0,1,1,4)
sire <- c(0,0,2,2,3)
data <- c(FALSE,FALSE,TRUE,FALSE,FALSE)
ped <- data.frame(id,dam,sire)
yn <- trimPed(ped,data)
ped <- ped[yn,]
}
\keyword{utilities}
