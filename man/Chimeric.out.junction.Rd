\name{Chimeric.out.junction}
\alias{Chimeric.out.junction}
\docType{data}
\title{
A simple junction file obtained from the chimeric alignment produced by STAR
}
\description{
This data is a sample data produced by STAR aligner.
}
\usage{data("Chimeric.out.junction")}
\format{
  A data frame with 7092 observations on the following 14 variables.
  \describe{
    \item{\code{V1}}{chromosome of the donor}
    \item{\code{V2}}{first base of the intron of the donor}
    \item{\code{V3}}{strand of the donor}
    \item{\code{V4}}{chromosome of the acceptor}
    \item{\code{V5}}{first base of the intron of the acceptor}
    \item{\code{V6}}{strand of the acceptor}
    \item{\code{V7}}{ junction type}
    \item{\code{V8}}{repeat length to the left of the junction}
  	\item{\code{V9}}{repeat length to the right of the junction}
	  \item{\code{V10}}{read name}
	  \item{\code{V11}}{first base of the first segment}
    \item{\code{V12}}{CIGAR of the first segment}
    \item{\code{V13}}{first base of the second segment}
    \item{\code{V14}}{CIGAR of the second segment}
  }
}
\details{
This is an example junction file obtained from the chimeric alignment.
}
\source{
This example data is produced by STAR aligner.
}
\references{
R package circRNAFull.
}
\examples{
data(Chimeric.out.junction)
## maybe str(Chimeric.out.junction) ; plot(Chimeric.out.junction) ...
}
\keyword{datasets}
