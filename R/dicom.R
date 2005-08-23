dicom.info <- function(fname, endian="little", flipud=TRUE) {
  ##
  ## "The default DICOM Transfer Syntax, which shall be supported by
  ## all AEs, uses Little Endian encoding and is specified in Annex
  ## A.1." (PS 3.5-2004, page 38)
  ##
  ## PS 3.5-2004, Sect 7.1.2: Data Element Structure with Explicit VR
  ## Explicit VRs store VR as text chars in 2 bytes.
  ## VRs of OB, OW, SQ, UN, UT have VR chars, then 0x0000, then 32 bit VL:
  ##
  ## +-----------------------------------------------------------+
  ## |  0 |  1 |  2 |  3 |  4 |  5 |  6 |  7 |  8 |  9 | 10 | 11 |
  ## +----+----+----+----+----+----+----+----+----+----+----+----+
  ## |<Group-->|<Element>|<VR----->|<0x0000->|<Length----------->|<Value->
  ##
  ## Other Explicit VRs have VR chars, then 16 bit VL:
  ##
  ## +---------------------------------------+
  ## |  0 |  1 |  2 |  3 |  4 |  5 |  6 |  7 |
  ## +----+----+----+----+----+----+----+----+
  ## |<Group-->|<Element>|<VR----->|<Length->|<Value->
  ##
  ## Implicit VRs have no VR field, then 32 bit VL:
  ##
  ## +---------------------------------------+
  ## |  0 |  1 |  2 |  3 |  4 |  5 |  6 |  7 |
  ## +----+----+----+----+----+----+----+----+
  ## |<Group-->|<Element>|<Length----------->|<Value->

  data(dicom.VRfields)
  data(dicom.fields)
  ## Some shortcuts...
  DICM.gr <- dicom.fields$group
  DICM.el <- dicom.fields$element
  VR.code <- dicom.VRfields$code
  ## Open connection
  fid <- file(fname, "rb")
  ## First 128 bytes are not used
  skip <- readBin(fid, integer(), n=128, size=1)
  ## Next four bytes spell "DICM"
  DICM <- readChar(fid, n=4)
  hdr <- NULL
  pixel.data <- FALSE
  while(!pixel.data) {
    group <- dec2hex(readBin(fid, integer(), size=2, endian=endian), 4)
    element <- dec2hex(readBin(fid, integer(), size=2, endian=endian), 4)
    if(group %in% "7FE0" && element %in% "0010")
      pixel.data <- TRUE
    vrstr <- readChar(fid, n=2)
    if(any(index <- VR.code %in% vrstr)) {
      VR <- dicom.VRfields[index, ]
      if(VR$bytes == 0) {
        ## This is an OB, OW, SQ, UN or UT, with 32-bit VL
        skip <- readBin(fid, integer(), size=2)
        length <- readBin(fid, integer(), size=4, endian=endian)
        if(pixel.data) {
          bytes <- as.numeric(hdr[hdr[,3] %in% "BitsAllocated", 6]) / 8
          img <- readBin(fid, integer(), length, size=bytes, endian=endian)
          nr <- as.numeric(hdr[hdr[,3] %in% "Rows", 6])
          nc <- as.numeric(hdr[hdr[,3] %in% "Columns", 6])
        } else {
          skip <- readBin(fid, integer(), length, size=1)
        }
        value <- NA
      } else {
        ## This is an explicit VR with 16-bit VL
        length <- readBin(fid, integer(), size=2, endian=endian)
        if(any(c("UL","US") %in% VR$code)) {
          ## "Unsigned Long" and "Unsigned Short"
          value <- readBin(fid, integer(), n=length/VR$bytes,
                           size=VR$bytes, signed=FALSE, endian=endian)
          value <- paste(value, collapse=" ")
        } else {
          if(any(c("SL","SS") %in% VR$code)) {
            ## "Signed Long" and "Signed Short"
            value <- readBin(fid, integer(), n=length/VR$bytes,
                             size=VR$bytes, signed=TRUE, endian=endian)
            value <- paste(value, collapse=" ")
          } else {
            ## Trim trailing white space
            value <- sub(" +$", "", readChar(fid, length))
            ## Replace all "\\"s with " "
            value <- gsub("[\\]", " ", value)
            ## Remove all non {a-zA-Z} characters with white space
            ## value <- gsub("[^{a-zA-Z}]", " ", value)
            if(VR$code %in% "UI")
              value <- sub("\\0", "", value) # Remove trailing \0
          }
        }
      }
    } else {
      ## Implicit VR with 32-bit VL
      stop("\"dicom.info\" cannot handle implicit VR.")
    }
    if(any((index1 <- DICM.gr %in% group) &
           (index2 <- DICM.el %in% element))) {
      hdr <- rbind(hdr,
                   c(group, element, dicom.fields[(index1 & index2),"name"],
                     VR$code, length, value))
    } else {
      hdr <- rbind(hdr, I(c(group, element, "Unknown", VR$code,
                            length, value)))
    }
  }
  close(fid)

  hdr <- data.frame(group = I(hdr[,1]), element = I(hdr[,2]),
                    name = I(hdr[,3]), code = I(hdr[,4]),
                    length = I(hdr[,5]), value = I(hdr[,6]))
  hdr$group <- as.character(hdr$group)
  hdr$element <- as.character(hdr$element)
  hdr$name <- as.character(hdr$name)
  hdr$code <- as.character(hdr$code)
  hdr$length <- as.numeric(hdr$length)
  hdr$value <- as.character(hdr$value)

  img <- matrix(img[1:(nc*nr)], nc, nr)
  if(flipud)
    img <- img[,nr:1]

  list(hdr = hdr, img = img)
}

dicom.separate <- function(path, debug=TRUE) {
  filenames <- system(paste("ls", path), intern=TRUE)
  nfiles <- length(filenames)
  headers <- images <- vector("list", nfiles)
  names(images) <- names(headers) <- sub("\\.", "", filenames)
  for(i in 1:nfiles) {
    if(debug && (i %% 100 == 0))
      cat("  ", i, "files processed...", fill=TRUE)
    dcm <- dicom.info(paste(path, filenames[i], sep="/"))
    images[[i]] <- dcm$img
    headers[[i]] <- dcm$hdr
  }
  list(hdr = headers, img = images)
}

dicom2analyze <- function(hdr, img) {
  A <- vector("list")
  A$size.of.hdr <- as.integer(348)
  A$endian <- "little"
  A$extents <- integer(1)
  A$session.error <- integer(1)
  ## dim
  dim <- c(4, nrow(img), ncol(img),
           ifelse(!is.na(dim(img)[3]), dim(img)[3], 1),
           ifelse(!is.na(dim(img)[4]), dim(img)[4], 1), 0, 0, 0)
  A$dim <- as.integer(dim)
  ## datatype
  datatype <- switch(hdr$value[hdr$name %in% "BitsAllocated"],
                     "8" = 4, "16" = 8, "32" = 16, "64" = 64)
  A$datatype <- as.integer(datatype)
  ##
  A$bitpix <- as.integer(hdr$value[hdr$name %in% "BitsAllocated"])
  A$dim.un0 <- integer(1)
  ## pixdim
  pixdim23 <- unlist(strsplit(hdr$value[hdr$name %in% "PixelSpacing"], " "))
  pixdim4 <- hdr$value[hdr$name %in% "SliceSpacing"]
  pixdim5 <- hdr$value[hdr$name %in% "RepetitionTime"]
  A$pixdim <- as.numeric(c(0, pixdim23,
                           ifelse(!length(pixdim4), 0, pixdim4),
                           ifelse(!length(pixdim5), 0, pixdim5),
                           0, 0, 0))
  ##
  A$vox.offset <- 0
  A$cal.max <- 0
  A$cal.min <- 0
  A$compressed <- 0
  A$verified <- 0
  A$glmax <- A$glmin <- integer(1)
  A$views <- integer(1)
  A$vols.added <- integer(1)
  A$start.field <- integer(1)
  A$field.skip <- integer(1)
  A$omax <- A$omin <- A$smax <- A$smin <- integer(1)
  ## descrip
  descrip <- hdr$value[hdr$name %in% "SeriesDescription"]
  if((n <- nchar(descrip)) < 80)
    A$descrip <- paste(descrip, paste(rep(" ", 80-n-1), collapse=""))
  else
    A$descrip <- substring(descrip, 1, 80)
  ##
  A$db.type <- paste(rep(" ", 10), collapse="")
  A$db.name <- paste(rep(" ", 18), collapse="")
  A$regular <- " "
  A$hkey.un0 <- " "
  A$vox.units <- paste(rep(" ", 4), collapse="")
  A$cal.units <- paste(rep(" ", 8), collapse="")
  A$aux.file <- paste(rep(" ", 24), collapse="")

  A$orient <- " "
  ## originator
  if(any(index <- hdr$name %in% "RequestingPhysician")) {
    originator <- hdr$value[index]
    if((n <- nchar(originator)) < 10)
      A$originator <- paste(originator, paste(rep(" ", 10-n-1), collapse=""))
    else
      A$originator <- substring(originator, 1, 10)
  }
  else
    A$originator <- paste(rep(" ", 10), collapse="")
  ##
  A$generated <- paste(rep(" ", 10), collapse="")
  ## scannum
  scannum <- hdr$value[hdr$name %in% "StudyID"]
  if((n <- nchar(scannum)) < 10)
    A$scannum <- paste(scannum, paste(rep(" ", 10-n-1), collapse=""))
  else
    A$scannum <- substring(scannum, 1, 10)
  ## patient.id
  patient.id <- hdr$value[hdr$name %in% "PatientName"]
  if((n <- nchar(patient.id)) < 10)
    A$patient.id <- paste(patient.id, paste(rep(" ", 10-n-1), collapse=""))
  else
    A$patient.id <- substring(patient.id, 1, 10)
  ## exp.date
  exp.date <- hdr$value[hdr$name %in% "StudyDate"]
  if((n <- nchar(exp.date)) < 10)
    A$exp.date <- paste(exp.date, paste(rep(" ", 10-n-1), collapse=""))
  else
    A$exp.date <- substring(exp.date, 1, 10)
  ## exp.time
  exp.time <- hdr$value[hdr$name %in% "StudyTime"]
  if((n <- nchar(exp.time)) < 10)
    A$exp.time <- paste(exp.time, paste(rep(" ", 10-n-1), collapse=""))
  else
    A$exp.time <- substring(exp.time, 1, 10)
  ##
  A$hist.un0 <- paste(rep(" ", 3), collapse="")
  return(A)
}

extract.hdr <- function(hdrs, string, numeric=TRUE, names=FALSE) {
  out.list <- lapply(hdrs,
                     function(hdr) {
                       if(sum(index <- hdr$name %in% string) > 0)
                         hdr$value[index]
                       else
                         NA
                       })
  out.names <- names(out.list)
  out.vec <- unlist(out.list)
  if(numeric) {
    out.vec <- as.numeric(out.vec)
    names(out.vec) <- out.names
  }
  return(out.vec)
}

dicom.table <- function(hdrs, fields, numeric=rep(TRUE,length(fields))) {
  N <- length(fields)
  mat <- NULL
  for(i in 1:N)
    mat <- cbind(mat, extract.hdr(hdrs, fields[i], numeric[i]))
  dimnames(mat)[[2]] <- fields
  as.data.frame(mat)
}

write.dicom.hdr <- function(dtable, filename, ...)
  write.table(dtable, filename, quote=FALSE, sep="\t", ...)

str2time <- function(tt) {
  tt <- as.numeric(tt)
  hh <- as.integer(trunc(tt / 10000))
  tt <- tt %% 10000
  mm <- as.integer(trunc(tt / 100))
  ss <- tt %% 100
  list(txt = sprintf("%02i:%02i:%08.5f", hh, mm, ss),
       time = 3600*hh + 60*mm + ss)
}
  
str2date <- function(dd)
  format(as.Date(dd, "%Y%m%d"), "%d %b %Y")
