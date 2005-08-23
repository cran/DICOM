\name{str2time}
\alias{str2time}
\alias{str2date}
\title{Convert DICOM Time/Date Entry}
\description{
  The DICOM time entry (TM) is converted into two alternative formats: a
  text version of the original format and a number in seconds.  The
  DICOM date entry (DA) is converted into a simple alternative format.
}
\usage{
str2time(tt)
str2date(dd)
}
\arguments{
  \item{tt}{TM field from a DICOM header.}
  \item{dd}{DA field from a DICOM header.}
}
\value{
  For "TM", a list structure containing two fields
  \item{txt}{A text version of the time where colons have been inserted
    for readability.}
  \item{time}{Time in seconds from midnight.}
  for "DA", a simple character string.
}
\details{
  DICOM "TM" format consists of a string of characters of the format
  \code{hhmmss.frac}; where \code{hh} contains hours (range "00" -
  "23"), \code{mm} contains minutes (range "00" - "59"), \code{ss}
  contains seconds (range "00" - "59"), and \code{frac} contains a
  fractional part of a second as small as 1 millionth of a second (range
  000000 - 999999).  A 24 hour clock is assumed.  Midnight can be
  represented by only 0000 since 2400 would violate the hour
  range.  The string may be padded with trailing spaces.  Leading and
  embedded spaces are not allowed.  One or more of the components
  \code{mm}, \code{ss}, or \code{frac} may be unspecified as long as
  every component to the right of an unspecified component is also
  unspecified.  If \code{frac} is unspecified the preceding "." may not
  be included.  \code{Frac} shall be held to six decimal places or less
  to ensure its format conforms to the ANSI HISPP MSDS Time common data
  type.  Examples:
  \item{1.}{070907.0705 represents a time of 7 hours, 9 minutes and
    7.0705 seconds.}
  \item{2.}{1010 represents a time of 10 hours, and 10 minutes.}
  \item{3.}{021 is an invalid value.}
  Notes: For reasons of backward compatibility with versions of this
  standard prior to V3.0, it is recommended that implementations also
  support a string of characters of the format \code{hh:mm:ss.frac} for
  this VR.

 DICOM "DA" format A string of characters of the format yyyymmdd; where
 yyyy shall contain year, mm shall contain the month, and dd shall
 contain the day. This conforms to the ANSI HISPP MSDS Date common data
 type. Example:
 \item{1.}{19930822  would represent August 22, 1993.}
 Notes: For reasons of backward compatibility with versions of this
 standard prior to V3.0, it is recommended that implementations also
 support a string of characters of the format yyyy.mm.dd for this VR.
}
\references{
  Digital Imaging and Communications in Medicine (DICOM)
  \url{http://medical.nema.org}
  \url{http://en.wikipedia.org/wiki/Digital_Imaging_and_Communications_in_Medicine}
}
\seealso{
  \code{\link{dicom.info}}
}
%\examples{}
\author{B. Whitcher}
\keyword{misc}