\name{CIRCexplorer_output}
\alias{CIRCexplorer_output}
\docType{data}
\title{
  Sample output of the circRNA prediction tool CIRCexplorer
}
\description{
  This data is a sample data produced by CIRCexplorer.
}
\usage{data("CIRCexplorer_output")}
\format{
  A data frame with 91 observations on the following 18 variables.
  \describe{
    \item{\code{V1}}{chromosome}
    \item{\code{V2}}{start of junction}
    \item{\code{V3}}{end of junction}
    \item{\code{V4}}{circula RNA/junction reads}
    \item{\code{V5}}{flag to indicate realignment of fusion junctions}
    \item{\code{V6}}{strand}
    \item{\code{V7}}{thickStart: no meaning}
    \item{\code{V8}}{thickEnd: no meaning}
    \item{\code{V9}}{itemRgb}
    \item{\code{V10}}{number of exons}
    \item{\code{V11}}{exon sizes}
    \item{\code{V12}}{exon offsets}
    \item{\code{V13}}{Number of junction reads}
    \item{\code{V14}}{circular RNA type}
    \item{\code{V15}}{name of gene}
    \item{\code{V16}}{name of isoform}
    \item{\code{V17}}{exon or intron index}
    \item{\code{V18}}{left intron/right intron}
  }
}
\details{
  This is an example output from the circRNA prediction tool CIRCexplorer.
}
\source{
  This example data is produced by CIRCexplorer.
}
\references{
  R package circRNAFull.
}
\examples{
  data(CIRCexplorer_output)
  ## maybe str(CIRCexplorer_output) ; plot(CIRCexplorer_output) ...
}
\keyword{datasets}
