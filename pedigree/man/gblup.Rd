\name{gblup}
\alias{gblup}
\title{Function to calculate breeding values using an animal model and a
  relationship matrix calculated from the markers (G matrix)}
\description{
  Fit an animal model to data, use a given variance ratio
  (\eqn{\alpha = \frac{\sigma^2_e}{\sigma^ 2_a}}). Calculate genetic
  relationship matrix using the function \code{calcG} of this package.
}
\usage{
gblup(formula, data, M, lambda)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{formula of the model, do not include the random effect
    due to animal (generally ID).}
  \item{data}{\code{data.frame} with columns corresponding to ID and the
    columns mentioned in the formula.}
  \item{M}{Matrix of marker genotypes, usually the count of one of the
    two SNP alleles at each markers (0, 1, or 2).}
  \item{lambda}{Variance ratio (\eqn{\frac{\sigma^2_e}{\sigma^ 2_a}})}
}
\value{
  Vector of solutions to the model, including random animal effects.}
\seealso{
  \code{\link[HaploSim]{SamplePedigree}}, \code{\link[pedigree]{gblup}},
   \code{\link[pedigree]{makeAinv}},\code{\link[pedigree]{blup}}
}
\examples{
h2 <- 0.5
example(SamplePedigree)
ped <- phList$ped
hList <- phList$hList
qtlList <- ListQTL(hList = hList,frqtl = 0.1,sigma2qtl = 1)
qtl <- tapply(unlist(qtlList),list(rep(names(qtlList),times = unlist(lapply(qtlList,length))),
                   unlist(lapply(qtlList,function(x)seq(1,length(x))))),mean,na.rm = TRUE)
qtl <- melt(qtl)
names(qtl) <- c("POS","TRAIT","a")
HH <- getAll(hList,translatePos = FALSE)
rownames(HH) <- sapply(hList,function(x)x@hID)
QQ <- HH[,match(qtl$POS,colnames(HH))]
g <- QQ%*%matrix(qtl$a,ncol = 1)
ped$G <- with(ped,g[match(hID0,rownames(g))]+g[match(hID1,rownames(g))])
sigmae <- sqrt(var(ped$G)/h2 - var(ped$G))
ped$P <- ped$G + rnorm(nrow(ped),0,sigmae)
M <- with(ped,HH[match(hID0,rownames(HH)),] + HH[match(hID1,rownames(HH)),]) 
rownames(M) <- ped$ID
sol <- gblup(P~1,data = ped[,c('ID','P')],M = M,lambda = 1/h2 - 1)
}
\keyword{utilities}
