\name{skip_exon}
\alias{skip_exon}
\docType{data}
\title{
  A sample skip exon file.
}
\description{
  This data is a sample data containing the skipping exons of circRNAs.
}
\usage{data("skip_exon")}
\format{
  A data frame with 2 observations on the following 7 variables.
  \describe{
    \item{\code{circle_id}}{circRNA id}
    \item{\code{transcript}}{transcript name}
    \item{\code{skip_exon_start}}{start of skip exon}
    \item{\code{skip_exon_end}}{end of skip exon}
    \item{\code{intron_start}}{start of intron}
    \item{\code{intron_end}}{end of intron}
    \item{\code{unique_reads}}{name of supporting reads}
  }
}
\details{
  This is an example skip exon file.
}
\source{
  This example data is obtained from \code{\link[circRNAFull]{skip_exon_detection}}.
}
\references{
  R package circRNAFull.
}
\examples{
  data(skip_exon)
}
\keyword{datasets}
