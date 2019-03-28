
#####################
# Writers           #
#####################

setGeneric("write", 
           function(object, verbose = FALSE, ...) standardGeneric("write"), signature = c("object", "verbose", "...")
)

#' @import hdf5r
setMethod("write", "SCopeLoom", function(object, verbose = FALSE) {
  
  o <- tryCatch({
    if(file.exists(object@file.path)) {
      h5 <- H5File$new(filename = object@file.path, mode = "r+")
    } else {
      h5 <- H5File$new(filename = object@file.path, mode = "w")
    }

    # Writing the matrix
    matrix <- write(object = object@matrix
        , file = h5
        , chunk.size = object@options$chunk.size
        , display.progress = object@options$display.progress
        , verbose = verbose)
    object@matrix <- matrix
    
    # Writing the global attributes
    # attrs
    new.attrs <- list()
    for(global.attr.idx in seq_along(along.with = object@attrs)) {
      global.attr <- object@attrs[[global.attr.idx]]
      ga <- write(object = global.attr, file = h5, verbose = verbose)
      new.attrs[[names(x = object@attrs)[[global.attr.idx]]]] <- ga
    }
    object@attrs <- new.attrs
    
    # Writing the row attributes
    # row_attrs
    if(verbose) {
      print("write row attributes...")
    }
    # create group
    if(!h5$link_exists(name = "row_attrs")) {
      h5$create_group("row_attrs")
    }
    new.row.attrs <- list()
    for(row.attr.idx in seq_along(along.with = object@row.attrs)) {
      row.attr <- object@row.attrs[[row.attr.idx]]
      ra <- write(object = row.attr, file = h5, verbose = verbose)
      new.row.attrs[[names(x = object@row.attrs)[[row.attr.idx]]]] <- ra
    }
    object@row.attrs <- new.row.attrs
    
    # Writing the column attributes
    # col_attrs
    if(verbose) {
      print("write column attributes...")
    }
    # create group
    if(!h5$link_exists(name = "col_attrs")) {
      h5$create_group("col_attrs")
    }
    new.col.attrs <- list()
    for(col.attr.idx in seq_along(along.with = object@col.attrs)) {
      col.attr <- object@col.attrs[[col.attr.idx]]
      ca <- write(object = col.attr, file = h5, verbose = verbose)
      new.col.attrs[[names(x = object@col.attrs)[[col.attr.idx]]]] <- ca
    }
    object@col.attrs <- new.col.attrs
    # col_edges
    print("write columns edges...")
    # create group
    if(!h5$link_exists(name = "col_edges")) {
      h5$create_group("col_edges")
    }
    # row_edges
    print("write row edges...")
    # create group
    if(!h5$link_exists(name = "row_edges")) {
      h5$create_group("row_edges")
    }
    # layers
    print("write layers...")
    # create group
    if(!h5$link_exists(name = "layers")) {
      h5$create_group("layers")
    }
    # finalize
    gc()
    h5$flush()
    h5$close_all()
    return (object)
  }, warning = function(w) {
    warning(w)
  }, error = function(e) {
    if (file.exists(object@file.path))
      file.remove(object@file.path)
    stop(e)
  }, finally = {
    
  })
  return (o)
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
    
    row.names(x = object) <- NULL
    colnames(x = object) <- NULL
    dtype <- GetDType(x = object[1, 1])
    loom$create_dataset(
      name = 'matrix',
      dtype = dtype,
      dims = rev(x = dim(x = object))
    )
    chunk.points <- ChunkPoints(
      data.size = dim(x = object)[2],
      chunk.size = chunk.size
    )
    if (display.progress) {
      pb <- utils::txtProgressBar(char = '=', style = 3)
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
  } else {
    print("Skip.")
  }
  
  return (object)
})

#####################
# Global Attributes

setMethod("write", "SCopeTree", function(object, ...) {
  if(!object@flush) {
    # Level 1
    l1 <- write(object = object@L1, ...)
    # Level 2
    l2 <- write(object = object@L2, ...)
    # Level 3
    l3 <- write(object = object@L3, ...)
    object@flush <- TRUE
  } else {
    print("Skip.")
  }
  return (object)
})


setMethod("write", "MetaData", function(object, ...) {
  if(verbose) print("write MetaData")
  
  if(!object@flush) {
    meta.data<-list(object@annotations, 
                    object@metrics, 
                    object@embeddings, 
                    object@clusterings, 
                    object@regulonThresholds)
    # Convert list to JSON
    meta.data.json<-rjson::toJSON(x = meta.data)
    # Compress (gzip) JSON 
    compressed.meta.data <- CompressGZIPBase64Encode(c = as.character(x = meta.data.json))
    loom.meta.data <- new(Class = "LoomMetaData", compressed.meta.data)
    write(object = loom.meta.data, ...)
    object@flush <- TRUE
  } else {
    print("Skip.")
  }
  return (object)
})

setMethod("write", "LoomGlobalAttribute", function(object, ...) {
  args <- list(...)
  loom = args$file
  
  if(verbose) {
    print(paste0(object@key))
  }
  
  # Get Dtype
  dtype <- guess_dtype(x = object)
  
  if(!object@flush) {
    loom$create_attr(attr_name = object@key, robj = object, dtype = dtype, space = GetDSpace(x = "scalar"))
    gc()
    flush(loom = loom)
    object@flush <- TRUE
  }
  return (object)
})

#####################
# Row Attributes

setMethod("write", "Regulons", function(object, ...) {
  if(verbose) print("write Regulons")
  
  ca.e.x <- new(Class = "LoomEmbeddingsX", t(x = e.x))
  write(object = ca.e.x, ...)
})

setMethod("write", "LoomRowAttribute", function(object, ...) {
  args <- list(...)
  loom = args$file
  
  if(verbose) {
    print(paste0(object@key))
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
  loom$create_dataset(name = paste0("row_attrs/", object@key), robj = AsGenericType(object = object), dtype = dtype)
  gc()
  flush(loom = loom)
})

setMethod("write", "RowAttribute", function(object, ...) {
  args <- list(...)
  loom = args$file
  
  if(!object@flush) {
    if(ncol(x = object) > 1) {
      data <- t(x = as.matrix(x = data.frame(object, stringsAsFactors = F)))
    } else {
      if(is.character(x = object[[1]])) {
        data <- as.character(x = object[[1]])
      } else if(is.numeric(x = object[[1]])) {
        data <- as.numeric(x = object[[1]])
      } else {
        stop("Not implemented")
      }
    }
    lra <- new(Class = "LoomRowAttribute", key = object@key, data)
    write(object = lra, verbose = verbose, ...)
    object@flush <- TRUE
  } else {
    print("Skip.")
  }
  return (object)
})

#####################
# Column Attributes

setMethod("write", "Embeddings", function(object, ...) {
  if(verbose) print("write Embeddings")
  
  # Check if any embedding hasn't been flushed yet
  if(any(!sapply(X = object, FUN = function(e) return (!e@flush)))) {
    print("Skip.")
  } else {
    # Default Embedding
    e.dflt <- object[[sapply(X = object, FUN = function(e) return (e@default))]]
    e.dflt <- as.matrix(x = e.dflt[,])
    colnames(e.dflt) <- NULL
    ca.e.dflt <- new(Class = "LoomDefaultEmbedding", t(x = e.dflt))
    write(object = ca.e.dflt, ...)
    # Embeddings_X
    e.x <- do.call(what = "cbind", args = lapply(X = object, FUN = function(e) {
      return (e[,1])
    }))
    ca.e.x <- new(Class = "LoomEmbeddingsX", t(x = e.x))
    write(object = ca.e.x, ...)
    # Embeddings_Y
    e.y <- do.call(what = "cbind", args = lapply(X = object, FUN = function(e) {
      return (e[,2])
    }))
    ca.e.y <- new(Class = "LoomEmbeddingsY", t(x = e.y))
    write(object = ca.e.y, ...)
    # Set flush
    for(e.idx in seq_along(along.with = list(x = object))) {
      object[[e.idx]]@flush <- TRUE
    }
  }
  return (object)
})

setMethod("write", "LoomColumnAnnotation", function(object, ...) {
  args <- list(...)
  if(verbose) {
    print(paste0(object@key, " - ( ", class(object), " 2)"))
  }
  # Write the column attribute
  callNextMethod()
})

setMethod("write", "ColumnAnnotation", function(object, ...) {
  
  if(verbose) {
    print(paste0(object@key))
  }
  
  if(!object@flush) {
    if(ncol(x = object) > 1) {
      data <- t(x = as.matrix(x = data.frame(object, stringsAsFactors = F)))
    } else {
      data <- as.character(x = object[[1]])
    }
    lca <- new(Class = "LoomColumnAnnotation", key = object@key, data)
    write(object = lca, ...)
    object@flush <- TRUE
  } else {
    print("Skip")
  }
  return (object)
})

setMethod("write", "LoomColumnMetric", function(object, ...) {
  args <- list(...)
  if(verbose) {
    print(paste0(object@key, " - ( ", class(object), " 2)"))
  }
  # Write the column attribute
  callNextMethod()
})

setMethod("write", "ColumnMetric", function(object, ...) {
  
  if(verbose) {
    print(paste0(object@key))
  }
  
  if(!object@flush) {
    if(ncol(x = object) > 1) {
      data <- t(x = as.matrix(x = data.frame(object, stringsAsFactors = F)))
    } else {
      data <- as.character(x = object[[1]])
    }
    lcm <- new(Class = "LoomColumnMetric", key = object@key, data)
    write(object = lcm, ...)
    object@flush <- TRUE
  } else {
    print("Skip.")
  }
  return (object)
})

setMethod("write", "LoomColumnAttribute", function(object, ...) {
  args <- list(...)
  loom = args$file
  
  # Get the Dtype
  if(is.character(x = object)) {
    dtype <- guess_dtype(x = object, string_len = "estimate")
  } else {
    dtype <- guess_dtype(x = object)
  }
  
  # if(verbose) print(dtype)
  
  # Create the dataset
  loom$create_dataset(name = paste0("col_attrs/", object@key), robj = AsGenericType(object = object), dtype = dtype)
  # Flush
  gc()
  flush(loom = loom)
})

setMethod("write", "ColumnAttribute", function(object, ...) {
  args <- list(...)

  if(verbose) {
    print(paste0(object@key))
  }
  
  if(!object@flush) {
    if(ncol(x = object) > 1) {
      data <- t(x = as.matrix(x = data.frame(object, stringsAsFactors = F)))
    } else {
      if(is.character(x = object[[1]])) {
        data <- as.character(x = object[[1]])
      } else if(is.numeric(x = object[[1]])) {
        data <- as.numeric(x = object[[1]])
      } else {
        stop("Not implemented")
      }
    }
    lca <- new(Class = "LoomColumnAttribute", key = object@key, data)
    write(object = lca, verbose = verbose, ...)
    object@flush <- TRUE
  } else {
    print("Skip.")
  }
  return (object)
})

