#######################
# Class Definitions
#######################
source(file = "R/utils.R")
source(file = "R/writers.R")
source(file = "R/validators.R")
source(file = "R/creators.R")

ROW_AXIS = 0
COLUMN_AXIS = 1

#' The SCopeLoom Class
#'
#' The SCopeLoom object is a R Loom representation of single-cell expression data for R.
#' 
#' @slot attrs A list of global attributes for this Loom
#' @slot col.attrs A list of column attributes for this Loom
#' @slot row.attrs A list of row attributes for this Loom
#' @slot col.edges A list of column edges for this Loom
#' @slot row.edges A list of row edges for this Loom
#'
#' @name SCopeLoom-class
#' @rdname SCopeLoom-class
#' @exportClass SCopeLoom
#' @import hdf5r
#'
SCopeLoom <- setClass(
  Class = 'SCopeLoom',
  slots = c(
    matrix = 'MainMatrix',
    attrs = 'list',
    col.attrs = 'list',
    row.attrs = 'list',
    col.edges = 'list',
    row.edges = 'list',
    layers = 'list',
    file.path = 'character',
    h5 = 'H5File',
    options = 'list'
  )
)

MainMatrix <- setClass(
  Class = 'MainMatrix',
  slots = c(
    flush = 'logical'
  ),
  prototype = prototype(
    flush = FALSE
  ),
  contains = 'dgCMatrix'
)

LoomGlobalAttribute <- setClass(
  Class = 'LoomGlobalAttribute',
  slots = c(
    key = 'character',
    flush = 'logical'
  ),
  prototype = prototype(
    flush = FALSE
  ),
  contains = 'character'
)

GlobalAttribute <- setClass(
  Class = 'GlobalAttribute',
  slots = c(
    key = 'character',
    flush = 'logical'
  ),
  prototype = prototype(
    flush = FALSE
  )
)

LoomRowAttribute <- setClass(
  Class = 'LoomRowAttribute',
  slots = c(
    axis = 'numeric',
    key = 'character'
  ),
  prototype = prototype(
    axis = ROW_AXIS
  ),
  contains = 'matrix'
)

RowAttribute <- setClass(
  Class = 'RowAttribute',
  slots = c(
    axis = 'numeric',
    key = 'character',
    flush = 'logical'
  ),
  prototype = prototype(
    axis = ROW_AXIS,
    flush = FALSE
  ),
  contains = 'data.frame'
)

LoomColumnAnnotation <- setClass(
  Class = 'LoomColumnAnnotation',
  contains = 'LoomColumnAttribute'
)

ColumnAnnotation <- setClass(
  Class = 'ColumnAnnotation',
  contains = 'ColumnAttribute'
)

LoomColumnMetric <- setClass(
  Class = 'LoomColumnMetric',
  contains = 'LoomColumnAttribute'
)

ColumnMetric <- setClass(
  Class = 'ColumnMetric',
  contains = 'ColumnAttribute'
)

LoomColumnAttribute <- setClass(
  Class = 'LoomColumnAttribute',
  slots = c(
    axis = 'numeric',
    key = 'character'
  ),
  prototype = prototype(
    axis = COLUMN_AXIS
  ),
  contains = 'matrix'
)

ColumnAttributeDataFrame <- setClass(
  Class = 'ColumnAttributeDataFrame',
  slots = c(
    axis = 'numeric',
    flush = 'logical'
  ),
  prototype = prototype(
    axis = COLUMN_AXIS,
    flush = FALSE
  ),
  contains = 'data.frame'
)

ColumnAttributeList <- setClass(
  Class = 'ColumnAttributeList',
  slots = c(
    axis = 'numeric',
    flush = 'logical'
  ),
  prototype = prototype(
    axis = COLUMN_AXIS,
    flush = FALSE
  ),
  contains = 'list'
)

### MetaData

LoomMetaData <- setClass(
  Class = "LoomMetaData",
  prototype = prototype(
    key = "MetaData"
  ),
  contains = "LoomGlobalAttribute"
)

MetaData <- setClass(
  Class = "MetaData",
  slots = c(
    annotations = 'list',
    metrics = 'list',
    embeddings = 'list',
    clusterings = 'list',
    regulonThresholds = 'list'
  ),
  prototype = prototype(
    annotations = list(),
    metrics = list(),
    embeddings = list(),
    clusterings = list(),
    regulonThresholds = list()
  ),
  contains = 'GlobalAttribute'
)

### Tree

LoomSCopeTreeLevel <- setClass(
  Class = "LoomSCopeTreeLevel",
  slots = c(
    key = 'character'
  ),
  contains = 'LoomGlobalAttribute'
)

SCopeTree <- setClass(
  Class = "SCopeTree",
  slots = c(
    L1 = 'LoomSCopeTreeLevel',
    L2 = 'LoomSCopeTreeLevel',
    L3 = 'LoomSCopeTreeLevel'
  ),
  contains = 'GlobalAttribute'
)

### Embeddings

LoomDefaultEmbedding <- setClass(
  Class = 'LoomDefaultEmbedding',
  prototype = prototype(
    key = 'Embedding'
  ),
  contains = 'LoomColumnAttribute'
)

LoomEmbeddingsX <- setClass(
  Class = 'LoomEmbeddingsX',
  prototype = prototype(
    key = 'Embeddings_X'
  ),
  contains = 'LoomColumnAttribute'
)

LoomEmbeddingsY <- setClass(
  Class = 'LoomEmbeddingsY',
  prototype = prototype(
    key = 'Embeddings_Y'
  ),
  contains = 'LoomColumnAttribute'
)

#' The Embedding Class (extends list)
#'
#' The Embeddings object is a representatial class for the embeddings of the observations that will stored as 
#' column attributes in the loom file through:
#' - CADefaultEmbedding (default embedding)
#' - CAEmbeddingsX (X coordinates of the additional embedding)
#' - CAEmbeddingsY (Y coordinates of the additional embeddings)
#' 
#' @name Embeddings-class
#' @rdname Embeddings-class
#' @exportClass Embeddings
#'
Embeddings <- setClass(
  Class = 'Embeddings',
  contains = 'list'
)

setMethod("initialize","Embedding", function(.Object, .Data, ...) {
  callNextMethod(.Object, .Data, ...)
})

#' The Embedding Class (extends data.frame)
#'
#' The Embedding object is a representatial class for a given embedding contained in the Embeddings Class.
#' 
#' @slot name A character to set the name of this embedding
#' @slot default A boolean to set this embedding as default 
#'
#' @name Embedding-class
#' @rdname Embedding-class
#' @exportClass Embedding
#'
Embedding <- setClass(
  Class = 'Embedding',
  slots = c(
    name = 'character',
    default = 'logical'
  ),
  prototype = prototype(
    default = FALSE
  ),
  contains = 'ColumnAttributeDataFrame'
)

### Regulons

#' The LoomRegulons Class (extends list)
#'
LoomRegulons <- setClass(
  Class = 'LoomRegulons',
  contains = 'LoomRowAttribute'
)

#' The Regulons Class (extends list)
#'
#' The Regulons object is a representatial class for the regulons that will stored as 
#' row attributes in the loom file through:
#' - LoomRegulons
#' 
#' @name Regulons-class
#' @rdname Regulons-class
#' @exportClass Regulons
#'
Regulons <- setClass(
  Class = 'Regulons',
  contains = 'ColumnAttributeList'
)
