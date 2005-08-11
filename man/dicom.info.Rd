\name{dicom.info}
\alias{dicom.info}
\title{Read Single DICOM File}
\description{
  All information, both header and image, is read into a list structure
  from a DICOM file.
}
\usage{
dicom.info(fname, endian="little", flipud=TRUE)
}
\arguments{
  \item{fname}{File name}
  \item{endian}{Binary encoding; default is \code{"little"}.}
  \item{flipud}{Flip image vertically; default is \code{TRUE}.}
}
\value{
  
}
\details{
  
}
\references{
  Digital Imaging and Communications in Medicine (DICOM)
  \url{http://medical.nema.org}
  \url{http://en.wikipedia.org/wiki/Digital_Imaging_and_Communications_in_Medicine}
}
\seealso{
  \code{\link{dicom.separate}}
}
%\examples{}
\author{Brandon Whitcher}
\keyword{file}
