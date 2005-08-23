\name{write.dicom.hdr}
\alias{write.dicom.hdr}
\title{Write DICOM Table to ASCII File}
\description{
  A wrapper to \code{write.table} specifically for DICOM tables.
}
\usage{
write.dicom.hdr(dtable, filename, ...)
}
\arguments{
  \item{dtable}{The DICOM table.}
  \item{filename}{Name of the file to be created.}
  \item{...}{Additional parameters to be passed to \code{write.table}.}
}
\value{
  None.
}
\details{
  This function is a straightforward wrapper to \code{write.table}.
}
\references{
  Digital Imaging and Communications in Medicine (DICOM)
  \url{http://medical.nema.org}
}
\seealso{
  \code{\link{write.table}}
}
%\examples{}
\author{B. Whitcher}
\keyword{file}