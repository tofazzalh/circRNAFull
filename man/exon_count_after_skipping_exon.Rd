\name{exon_count_after_skipping_exon}
\alias{exon_count_after_skipping_exon}
\docType{data}
\title{
  A sample exon count file after deleting skip exon.
}
\description{
  This data is a sample exon count data deleting skip exon.
}
\usage{data("exon_count_after_skipping_exon")}
\format{
  A data frame with 65 observations on the following 8 variables.
  \describe{
    \item{\code{V1}}{circle ID}
    \item{\code{V2}}{transcript name}
    \item{\code{V3}}{exon number}
    \item{\code{V4}}{chromosome}
    \item{\code{V5}}{start of exon}
    \item{\code{V6}}{end of exon}
    \item{\code{V7}}{strand}
    \item{\code{V8}}{number of supporting reads}
  }
}
\details{
  This is an example exon count file after deleting skip exon.
}
\source{
  This example data is exon count data after deleting skip exon.
}
\references{
  R package circRNAFull.
}
\examples{
  data(exon_count_after_skipping_exon)
  ## maybe str(exon_count_after_skipping_exon) ; plot(exon_count_after_skipping_exon) ...
}
\keyword{datasets}
