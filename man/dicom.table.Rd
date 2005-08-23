\name{dicom.table}
\alias{dicom.table}
\title{Construct Data Frame from DICOM Headers}
\description{
  A data frame is created given the valid DICOM fields provided by the
  user.
}
\usage{
dicom.table(hdrs, fields, numeric=rep(TRUE,length(fields)))
}
\arguments{
  \item{hdrs}{List object of DICOM headers.}
  \item{fields}{Vector of character strings; i.e., valid DICOM field
    names.}
  \item{numeric}{Logical vector; each field is converted to numeric when
    \code{TRUE}.}    
}
\value{
  Data frame where the rows correspond to images and the columns to
  DICOM fields.
}
\details{
  The function \code{extract.hdr} is called for each DICOM field and the
  resulting vectors are rearranged into a convenient matrix.
}
\references{
  Digital Imaging and Communications in Medicine (DICOM)
  \url{http://medical.nema.org}
}
\seealso{
  \code{\link{extract.hdr}}
}
%\examples{}
\author{B. Whitcher}
\keyword{misc}