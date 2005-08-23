\name{extract.hdr}
\alias{extract.hdr}
\title{Extract Single Field from DICOM Headers}
\description{
  A particular DICOM field is extracted for a collection of DICOM
  headers.
}
\usage{
extract.hdr(hdrs, string, numeric=TRUE, names=FALSE)
}
\arguments{
  \item{hdrs}{List object of DICOM headers.}
  \item{string}{DICOM field name.}
  \item{numeric}{Logical; values are converted to numbers when
    \code{TRUE}.}
  \item{names}{Logical; file names are kept with elements of the
    vector.}
}
\value{
  Vector of values from the requested DICOM field.
}
\details{
  The DICOM field is extracted from each DICOM header and placed into a
  vector.
}
\references{
  Digital Imaging and Communications in Medicine (DICOM)
  \url{http://medical.nema.org}
}
\seealso{
  \code{\link{dicom.info}}
}
%\examples{}
\author{B. Whitcher}
\keyword{misc}