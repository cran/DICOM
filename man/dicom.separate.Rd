\name{dicom.separate}
\alias{dicom.separate}
\title{Read All DICOM Files in a Directory}
\description{
  All DICOM files are imported and a text file summarizing their content
  recorded.
}
\usage{
dicom.separate(path, debug=TRUE)
}
\arguments{
  \item{path}{Path name to the DICOM directory.}
  \item{debug}{Flag to provide limited progress information.}
}
%\value{}
%\details{}
\references{
  Digital Imaging and Communications in Medicine (DICOM)\cr
  \url{http://medical.nema.org}
}
\seealso{
  \code{\link{dicom.info}}
}
%\examples{}
\author{Brandon Whitcher}
\keyword{file}