
#####################
# Writers           #
#####################

setGeneric("write", 
           function(object, verbose = FALSE, ...) standardGeneric("write"), signature = c("object", "verbose", "...")
)

#' @import hdf5r
setMethod("write", "SCopeLoom", function(object, verbose = FALSE) {
  
  tryCatch({
    h5 <- H5File$new(filename = object@file.path, mode = "w")

    # Writing the matrix
    write(object = object@matrix
          , file = h5
          , chunk.size = object@options$chunk.size
          , display.progress = object@options$display.progress
          , verbose = verbose)
  
    # Writing the global attributes
    # attrs
    for(global.attr in object@attrs) {
      write(object = global.attr, file = h5, verbose = verbose)
    }
    
    # Writing the row attributes
    # row_attrs
    if(verbose) {
      print("Adding row attributes...")
    }
    h5$create_group("row_attrs")
    for(row.attr in object@row.attrs) {
      write(object = row.attr, file = h5, verbose = verbose)
    }
    
    # Writing the column attributes
    # col_attrs
    if(verbose) {
      print("Adding column attributes...")
    }
    h5$create_group("col_attrs")
    for(column.attr in object@col.attrs) {
      write(object = column.attr, file = h5, verbose = verbose)
    }
    # col_edges
    print("Adding columns edges...")
    h5$create_group("col_edges")
    # row_edges
    print("Adding row edges...")
    h5$create_group("row_edges")
    # layers
    print("Adding layers...")
    h5$create_group("layers")
    # finalize
    gc()
    h5$flush()
    h5$close_all()
  }, warning = function(w) {
    warning(w)
  }, error = function(e) {
    if (file.exists(object@file.path))
      file.remove(object@file.path)
    stop(e)
  }, finally = {

  })
})

setMethod("write", "MainMatrix", function(object, ...) {
  args <- list(...)
  loom <- args$file
  chunk.size <- args$chunk.size
  display.progress <- args$display.progress
  
  if(verbose) print("write MainMatrix")
  
  if(!object@flush) {
    if(loom$mode == "r") stop("File open as read-only.")
    
    # Check if any NA values
    if(any(is.na(x = object))) {
      stop("Please make sure that the expression matrix (dgem) does not contain any NA values.")
    }
    
    row.names(x = object)<-NULL
    colnames(x = object)<-NULL
    dtype<-GetDType(x = object[1, 1])
    loom$create_dataset(
      name = 'matrix',
      dtype = dtype,
      dims = rev(x = dim(x = object))
    )
    chunk.points<-ChunkPoints(
      data.size = dim(x = object)[2],
      chunk.size = chunk.size
    )
    if (display.progress) {
      pb<-utils::txtProgressBar(char = '=', style = 3)
    }
    
    for (col in 1:ncol(x = chunk.points)) {
      row.start <- chunk.points[1, col]
      row.end <- chunk.points[2, col]
      loom[['matrix']][row.start:row.end, ] <- t(x = as.matrix(x = object[, row.start:row.end]))
      if(display.progress) {
        utils::setTxtProgressBar(pb = pb, value = col / ncol(x = chunk.points))
      }
    }
    flush(loom = loom)
    object@flush <- TRUE
  }
  
})

setMethod("write", "SCopeTree", function(object, ...) {
  # Level 1
  write(object = object@L1, ...)
  # Level 2
  write(object = object@L2, ...)
  # Level 3
  write(object = object@L3, ...)
})


setMethod("write", "MetaData", function(object, ...) {
  meta.data<-list(object@annotations, 
                  object@metrics, 
                  object@embeddings, 
                  object@clusterings, 
                  object@regulonThresholds)
  # Convert list to JSON
  meta.data.json<-rjson::toJSON(meta.data)
  # Compress (gzip) JSON 
  compressed.meta.data<-CompressGZIPBase64Encode(c = as.character(meta.data.json))
  ga.meta.data <- new(Class = "GAMetaData", compressed.meta.data)
  write(object = ga.meta.data, ...)
})

setMethod("write", "LoomGlobalAttribute", function(object, ...) {
  args <- list(...)
  loom = args$file
  
  if(verbose) {
    print(paste0(object@key, " - ( ", class(object), " )"))
  }
  
  # Get Dtype
  dtype<-guess_dtype(x = object)
  
  if(!object@flush) {
    loom$create_attr(attr_name = object@key, robj = object, dtype = dtype, space = GetDSpace(x = "scalar"))
    gc()
    flush(loom = loom)
    object@flush <- TRUE
  }
})

setMethod("write", "LoomRowAttribute", function(object, ...) {
  args <- list(...)
  loom = args$file
  
  if(verbose) {
    print(paste0(object@key, " - ( ", class(object), " )"))
  }
  
  # Get Dtype
  if(is.character(x = object)) {
    dtype <- guess_dtype(x = object, string_len = "estimate")
  } else {
    dtype <- guess_dtype(x = object)
  }
  
  # Transpose
  object <- t(x = object)
  
  # if(verbose) print(dtype)
    
  # Create the row dataset
  if(!object@flush) {
    loom$create_dataset(name = paste0("row_attrs/", object@key), robj = AsGenericType(object = object), dtype = dtype)
    gc()
    flush(loom = loom)
    object@flush <- TRUE
  }
})

setMethod("write", "Embeddings", function(object, ...) {
  # Default Embedding
  e.dflt <- object[[sapply(X = object, FUN = function(e) return (e@default))]]
  e.dflt <- as.matrix(x = e.dflt[,])
  colnames(e.dflt) <- NULL
  ca.e.dflt <- new(Class = "CADefaultEmbedding", t(x = e.dflt))
  write(object = ca.e.dflt, ...)
  # Embeddings_X
  e.x <- do.call(what = "cbind", args = lapply(X = object, FUN = function(e) {
    return (e[,1])
  }))
  ca.e.x <- new(Class = "CAEmbeddingsX", t(x = e.x))
  write(object = ca.e.x, ...)
  # Embeddings_Y
  e.y <- do.call(what = "cbind", args = lapply(X = object, FUN = function(e) {
    return (e[,2])
  }))
  ca.e.y <- new(Class = "CAEmbeddingsY", t(x = e.y))
  write(object = ca.e.y, ...)
})

setMethod("write", "LoomColumnAnnotation", function(object, ...) {
  args <- list(...)
  if(verbose) {
    print(paste0(object@key, " - ( ", class(object), " 2)"))
  }
  # Write the column attribute
  callNextMethod()
})

setMethod("write", "LoomColumnMetric", function(object, ...) {
  args <- list(...)
  if(verbose) {
    print(paste0(object@key, " - ( ", class(object), " 2)"))
  }
  # Write the column attribute
  callNextMethod()
})

setMethod("write", "LoomColumnAttribute", function(object, ...) {
  args <- list(...)
  loom = args$file

  if(verbose) {
    print(paste0(object@key, " - ( ", class(object), " )"))
  }
  # Get the Dtype
  if(is.character(x = object)) {
    dtype <- guess_dtype(x = object, string_len = "estimate")
  } else {
    dtype <- guess_dtype(x = object)
  }
  
  # if(verbose) print(dtype)
  
  if(!object@flush) {
    # Create the dataset
    loom$create_dataset(name = paste0("col_attrs/", object@key), robj = AsGenericType(object = object), dtype = dtype)
    # Flush
    gc()
    flush(loom = loom)
    object@flush <- TRUE
  }
})

