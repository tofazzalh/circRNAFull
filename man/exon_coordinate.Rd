\name{exon_coordinate}
\alias{exon_coordinate}
\docType{data}
\title{
  A sample text file containing the exon cordinates obtained from reference genome annotation file.
}
\description{
  This data is a sample data containing the exon coordinates of chromosome 21.
}
\usage{data("exon_coordinate")}
\format{
  A data frame with 8939 observations on the following 7 variables.
  \describe{
    \item{\code{V1}}{chromosome}
    \item{\code{V2}}{start of exon}
    \item{\code{V3}}{end of exon}
    \item{\code{V4}}{transcript name}
    \item{\code{V5}}{transcript with exon index}
    \item{\code{V6}}{strand}
    \item{\code{V7}}{transcript with chromosome}
  }
}
\details{
  This is an example exon coordinate file.
}
\source{
  This example data is obtained from human genome annotation file.
}
\references{
  R package circRNAFull.
}
\examples{
  data(exon_coordinate)
}
\keyword{datasets}
