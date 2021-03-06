\name{listDD}
\alias{listDD}
\title{List direct descendants for all nodes of a tree}
\description{
  The function \code{listDD} lists the direct descendants from each node
  of a tree. The tree can be of class \code{\link[ape:read.tree]{phylo}},
  \linkS4class{phylo4} or \linkS4class{phylo4d}.
}
\usage{
listDD(x, nameBy=c("label","number"))
}
\arguments{
  \item{x}{A tree of  class \code{\link[ape:read.tree]{phylo}},
    \linkS4class{phylo4} or \linkS4class{phylo4d}.}
  \item{nameBy}{a character string indicating whether the returned list
    must be named by node labels ("label") or by node numbers ("number").}
}
\value{
  A list whose components are vectors of named nodes (or tips) for a
  given internal node.
}
\author{ Thibaut Jombart \email{tjombart@imperial.ac.uk} }
\seealso{\code{\link{listTips}} which lists the tips descending from
  each node. \cr
  
  \code{\link{treePart}} which defines partitions of tips
  according to the tree topology.}
\examples{
if(require(ape) & require(phylobase)){
## make a tree
x <- as(rtree(20),"phylo4")
plot(x,show.node=TRUE)
listDD(x)
}
}
\keyword{manip}
