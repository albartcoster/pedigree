\name{countOff}
\Rdversion{1.1}
\alias{countOff}
\title{Function that counts the number of offspring (and following
  generations for each individual in a pedigree.}
\description{
  Function to count the number of offspring for each individual in a
  pedigree. With loops, offspring of later generations will be counted
  several times.}
\usage{
countOff(ped)
}
\arguments{
  \item{ped}{\code{data.frame} with three columns: id,id parent1,id parent2}
}
\value{Numeric vector with number of offspring for each individual in
  the pedigree.}
\author{
  Albart Coster
}
\examples{
example(countGen)
countOff(ped)
}
\keyword{utilities}
