% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{correctIntegerCN}
\alias{correctIntegerCN}
\title{FUNCTION_TITLE}
\usage{
correctIntegerCN(
  cn,
  segs,
  callColName = "event",
  purity,
  ploidy,
  cellPrev,
  maxCNtoCorrect.autosomes = NULL,
  maxCNtoCorrect.X = NULL,
  correctHOMD = TRUE,
  minPurityToCorrect = 0.2,
  gender = "male",
  chrs = c(1:22, "X")
)
}
\arguments{
\item{cn}{PARAM_DESCRIPTION}

\item{segs}{PARAM_DESCRIPTION}

\item{callColName}{PARAM_DESCRIPTION, Default: 'event'}

\item{purity}{PARAM_DESCRIPTION}

\item{ploidy}{PARAM_DESCRIPTION}

\item{cellPrev}{PARAM_DESCRIPTION}

\item{maxCNtoCorrect.autosomes}{PARAM_DESCRIPTION, Default: NULL}

\item{maxCNtoCorrect.X}{PARAM_DESCRIPTION, Default: NULL}

\item{correctHOMD}{PARAM_DESCRIPTION, Default: TRUE}

\item{minPurityToCorrect}{PARAM_DESCRIPTION, Default: 0.2}

\item{gender}{PARAM_DESCRIPTION, Default: 'male'}

\item{chrs}{PARAM_DESCRIPTION, Default: c(1:22, "X")}
}
\value{
OUTPUT_DESCRIPTION
}
\description{
FUNCTION_DESCRIPTION
}
\details{
DETAILS
}
\examples{
\dontrun{
if(interactive()){
 #EXAMPLE1
 }
}
}
