% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/circRNAFull.R
\name{extract_reads_from_bam}
\alias{extract_reads_from_bam}
\title{Extracting individual alignment file for each circRNA from chimeric alignment bam file}
\usage{
extract_reads_from_bam(circle_id, bamfile, outfolder)
}
\arguments{
\item{circle_id}{A data frame containing circle ids, transcript name and spanning reads of circRNAs (obtained from function \code{\link[circRNAFull]{extract_circle_ids}})}

\item{bamfile}{A chimeric alignment bam file read as a data frame}

\item{outfolder}{The name of the output directory}
}
\value{
Individual alignment files for each circRNAs
}
\description{
This function generates alignment file for each circRNA from the chimeric alignment produced by STAR
}
\examples{
#loading an example circle_id file
data(circle_ID)
circle_id<-circle_ID

\donttest{
\dontrun{
#Please upload your chimeric alignment bam file.Suppose you have the chimeric
#alignemnt bam file 'Chimeric.out.bam' in you working directory.
#Then run:
bam <- scanBam("Chimeric.out.bam")
bamfile <- as.data.frame(bam)
}
}

#Creating an output directory
outfolder<-tempdir()

\donttest{
\dontrun{
#Extracting individual alignment file for each circRNAs. The individual alignment file
#will be generated in the 'outfolder' directory.
extract_reads_from_bam<-function(circle_id, bamfile, outfolder)
}
}
}
