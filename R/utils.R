
#####################
# Utils             #
#####################

#'@title GetDType
#'@description Get HDF5 data types.
#'@author mojaveazure
#'@param x An R object or string describing HDF5 datatype.
#'@return The corresponding HDF5 data type.
#'@import hdf5r
#'@seealso \link{hdf5r::h5types}
#'@export
GetDType<-function(x) {
  return(switch(
    EXPR = class(x = x),
    'numeric' = h5types$double,
    'integer' = h5types$int,
    'character' = H5T_STRING$new(size = Inf),
    'logical' = H5T_LOGICAL$new(),
    stop(paste("Unknown data type:", class(x = x)))
  ))
}

#'@title get_dspace
#'@description Get the HDF5 dataspace interface object given x.
#'@param x The dataspace interface name. Either scalar or simple.
#'@return The corresponding HDF5 dataspace interface object.
GetDSpace<-function(x) {
  dspaces<-c("scalar", "simple")
  if(!("scalar" %in% dspaces)) {
    stop(paste("Wrong dspace. Choose either scalar or simple."))
  }
  return (H5S$new(type = x, dims = NULL, maxdims = NULL))
}

AsGenericType <- function(object) {
  # Set type
  if(is.character(x = object)) {
    object <- as.character(x = object)
  } else if(is.numeric(x = object)) {
    if(object@axis == 0) {
      if(ncol(x = object) < 2) {
        object <- as.numeric(x = object)
      }
    } else if(object@axis == 1) {
      if(nrow(x = object) < 2) {
        object <- as.numeric(x = object)
      }
    }
  } else {
    stop("Error: type not recognized.")
  }
  return (object)
}

#' @title ChunkPoints
#' @description  Generate chunk points.
#' @author mojaveazure
#' @param data.size How big is the data being chunked.
#' @param chunk.size How big should each chunk be.
#' @return A matrix where each column is a chunk, row 1 is start points, row 2 is end points.
#' @export
ChunkPoints<-function(data.size, chunk.size) {
  return(vapply(
    X = 1L:ceiling(data.size / chunk.size),
    FUN = function(i) {
      return(c(
        start = (chunk.size * (i - 1L)) + 1L,
        end = min(chunk.size * i, data.size)
      ))
    },
    FUN.VALUE = numeric(length = 2L)
  ))
}

#'@title CompressGZIPBase64Encode
#'@description Gzip compress and Base 64 encode a character
#'@param c The character to compress
CompressGZIPBase64Encode<-function(c) {
  return (base64enc::base64encode(what = memCompress(from = c, type = "gzip")))
}

#'@title DecompressGZIPBase64Decode
#'@description Base 64 decode and decompress a character
#'@param gzb64c The character to decompress
DecompressGZIPBase64Decode<-function(gzb64c) {
  return (rawToChar(memDecompress(from = base64enc::base64decode(what = gzb64c), type = "gzip", asChar = F), multiple = F))
}
