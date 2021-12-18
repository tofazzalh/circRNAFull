\name{circle_ID}
\alias{circle_ID}
\docType{data}
\title{
  A sample file containing the circle ids, transcript name and spanning reads of circRNAs.
}
\description{
  This data is a sample data containing the circle ids, transcript name and spanning reads of circRNAs.
}
\usage{data("circle_ID")}
\format{
  A data frame with 90 observations on the following 3 variables.
  \describe{
    \item{\code{V1}}{circRNA id}
    \item{\code{V2}}{transcript name}
    \item{\code{V3}}{spanning reads of the circRNAs}
  }
}
\details{
  This is an example circle id file.
}
\source{
  This example data is obtained from \code{\link[circRNAFull]{extract_circle_ids}}.
}
\references{
  R package circRNAFull.
}
\examples{
  data(circle_ID)
}
\keyword{datasets}
