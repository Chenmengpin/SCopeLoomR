#######################
# Class Definitions
#######################

ROW_AXIS = 0
COLUMN_AXIS = 1

setGeneric("write", 
           function(object) standardGeneric("write")
)

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
#'
SCopeLoom <- setClass(
  Class = 'SCopeLoom',
  slots = c(
    matrix = 'matrix',
    attrs = 'list',
    col.attrs = 'list',
    row.attrs = 'list',
    col.edges = 'list',
    row.edges = 'list',
    layers = 'list',
    options = 'list'
  )
)

CreateSCopeLoom <- function(
  title = NULL
  , genome = NULL
  , dgem
  , embedding = NULL
  , embedding.name = NULL
  , tree = NULL
  # , fbgn.gn.mapping.file.path = NULL
  , chunk.size = 1000
  , display.progress = T
) {
  loom <- new(Class = "SCopeLoom")
  cn<-colnames(dgem)
  rn<-row.names(dgem)
  ##########
  # Matrix #
  ##########
  is.dgem.sparse<-F
  # Check the type of the sparse matrix
  # convert to dgCMatrix if necessary to speedup populating the matrix slot
  if(class(dgem) == "dgTMatrix") {
    print("Converting to dgCMatrix...")
    dgem<-methods::as(object = dgem, Class = "dgCMatrix")
  }
  # Check if sparse matrix
  if(sum(class(x = dgem) %in% c("dgCMatrix", "dgTMatrix")) > 0) {
    is.dgem.sparse<-T
  }
  loom@matrix <- dgem
  
  #####################
  # Global Attributes #
  #####################
  loom@attrs[["title"]] <- new(Class = "LoomGlobalAttribute", key = "title", name = title)
  loom@attrs[["genome"]] <- new(Class = "LoomGlobalAttribute", key = "Genome", name = genome)
  loom@attrs[["tree"]] <- tree
  loom@attrs[["R.version"]] <- new(Class = "LoomGlobalAttribute", key = "R.version", name = as.character(R.version.string))
  loom@attrs[["CreationDate"]] <- new(Class = "LoomGlobalAttribute", key = "CreationDate", name = as.character(Sys.time()))
  loom@attrs[["MetaData"]] <- new(Class = "MetaData", key = "MetaData")
  
  #####################
  # Column Attributes #
  #####################
  
  ## CellID
  loom@col.attrs[["CellID"]] <- t(x = as.matrix(x = cn))
  ## Embeddings (Add the default embedding)
  if(!is.null(x = default.embedding)) {
    print("Adding default embedding...")
    if(is.null(x = default.embedding.name)) {
      default.embedding.name<-"Embedding"
    }
    embeddings <- new(Class = "Embeddings")
    embeddings[[default.embedding.name]] <- new(
      Class = "Embedding"
      , name = default.embedding.name
      , default = TRUE)
    loom@col.attrs[["Embeddings"]] <- embeddings
  } else {
    warning("No default embedding set for the loom. You'll not be able to visualize it in SCope.")
  }
  
  ## nUMI
  # Check if matrix is raw counts
  if(any(dgem%%1!=0)) {
    print("Adding default metrics nUMI...")
    if(is.dgem.sparse) {
      nUMI<-Matrix::colSums(dgem)
    } else {
      nUMI<-colSums(dgem)
    }
  } else {
    warning("Default metric nUMI was not added because the input matrix does not seem to be UMI counts.")
  }
  loom@col.attrs[["nUMI"]] <- nUMI
  
  ## nGene
  print("Adding default metrics nGene...")
  if(is.dgem.sparse) {
    nGene<-Matrix::colSums(dgem > 0)
  } else {
    nGene<-colSums(dgem > 0)
  }
  loom@col.attrs[["nGene"]] <- nGene
  
  #####################
  # Column Attributes #
  #####################
  
  ## Gene
  loom@row.attrs[["Gene"]] <- t(x = as.matrix(x = rn))
  
  # Check if Flybase gene mapping is not empty
  if(!is.null(fbgn.gn.mapping.file.path)) {
    AddFBgn(loom = loom, dgem = dgem, fbgn.gn.mapping.file.path = fbgn.gn.mapping.file.path)
  }
  
  # Options
  loom@options[["chunk.size"]] <- chunk.size
  loom@options[["display.progress"]] <- display.progress
}

setMethod("write", "CreateSCopeLoom", function(object) {
  loom<-H5File$new(filename = file.name, mode = "w")
})

LoomGlobalAttribute <- setClass(
  Class = 'LoomGlobalAttribute',
  slots = c(
    key = 'character'
  ),
  contains = 'character'
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

setMethod("write", "LoomColumnAttribute", function(object) {
  print(object@key)
})

### MetaData

GAMetaData <- setClass(
  Class = "MetaData",
  prototype = c(
    key = "MetaData"
  ),
  contains = "LoomGlobalAttribute"
)

setMethod("write", "MetaData", function(object) {
  meta.data<-list(object@annotations, 
                  object@metrics, 
                  object@embeddings, 
                  object@clusterings, 
                  object@regulonThresholds)
  # Convert list to JSON
  meta.data.json<-rjson::toJSON(meta.data)
  # Compress (gzip) JSON 
  compressed.meta.data<-compress_gzb64(c = as.character(meta.data.json))
  ga.meta.data <- new(Class = "GAMetaData", compressed.meta.data)
  write(object = ga.meta.data)
})

MetaData <- setClass(
  Class = "MetaData",
  slots = c(
    annotations = 'list',
    metrics = 'list',
    embeddings = 'list',
    clusterings = 'list',
    regulonThresholds = 'list'
  )
)

### Tree

SCopeTreeLevel <- setClass(
  Class = "TreeLevel",
  contains = 'LoomGlobalAttribute'
)

SCopeTree <- setClass(
  Class = "SCopeTree",
  slots = c(
    L1 = 'SCopeTreeLevel',
    L2 = 'SCopeTreeLevel',
    L3 = 'SCopeTreeLevel'
  )
)

### Embeddings

CADefaultEmbedding <- setClass(
  Class = 'CADefaultEmbedding',
  prototype = prototype(
    key = 'Embedding'
  ),
  contains = 'LoomColumnAttribute'
)

CAEmbeddingsX <- setClass(
  Class = 'CAEmbeddingsX',
  prototype = prototype(
    key = 'Embeddings_X'
  ),
  contains = 'LoomColumnAttribute'
)

CAEmbeddingsY <- setClass(
  Class = 'CAEmbeddingsY',
  prototype = prototype(
    key = 'Embeddings_Y'
  ),
  contains = 'LoomColumnAttribute'
)

setMethod("write", "Embeddings", function(object) {
  # Default Embedding
  e.dflt <- object@Embeddings[lapply(X = object@Embeddings, FUN = function(e) return (e@default))]
  e.dflt <- as.matrix(x = e.dflt[,])
  colnames(e.dflt) <- NULL
  ca.e.dflt <- new(Class = "CADefaultEmbedding", e.dflt)
  write(object = ca.e.dflt)
  # Embeddings_X
  e.x <- do.call(what = "cbind", args = lapply(X = object@Embeddings, FUN = function(e) {
    return (e[,1])
  }))
  ca.e.x <- new(Class = "CAEmbeddingsX", e.x)
  write(object = ca.e.x)
  # Embeddings_Y
  e.y <- do.call(what = "cbind", args = lapply(X = object@Embeddings, FUN = function(e) {
    return (e[,2])
  }))
  ca.e.y <- new(Class = "CAEmbeddingsX", e.y)
})

Embeddings <- setClass(
  Class = 'Embeddings',
  contains = 'list'
)

setValidity("Embedding", function(object) {
  y_dims <- dim(x = object)[2]
  # Check if 2 Y-coordinates
  if (y_dims != 2) {
    paste("Embedding has Y-dimension length of ", y_dims, ". Should be 2", sep = "")
  } else if(any(is.na(x = object))) { # Check if any NA values
    "Embedding contains NA values. No NAs are allowed"
  } else {
    TRUE
  }
})

Embedding <- setClass(
  Class = 'Embedding',
  slots = c(
    name = 'character',
    default = 'logical'
  ),
  prototype = prototype(
    default = FALSE
  ),
  contains = 'data.frame'
)

#' #'@title AddFBgn
#' #'@description Add the Flybase gene as a row attribute to the given .loom file handler.
#' #'@param loom                       The loom file handler.
#' #'@param dgem                       A matrix of the gene expression with M genes as rows and N cells as columns.
#' #'@param fbgn.gn.mapping.file.path  A N-by-2 data.frame containing the mapping between the Flybase gene and the gene symbol.
#' #'@export
#' AddFBgn<-function(loom
#'                    , dgem
#'                    , fbgn.gn.mapping.file.path) {
#'   if(loom$mode=="r") stop("File open as read-only.")
#'   fbgn.gn.mapping<-utils::read.table(file = fbgn.gn.mapping.file.path, header = F, sep = "\t", quote = '', stringsAsFactors = F)
#'   colnames(fbgn.gn.mapping)<-c("FBgn",RA_GENE_NAME)
#'   tmp<-data.frame(row.names(dgem))
#'   colnames(tmp)<-RA_GENE_NAME
#'   genes<-merge(x = tmp, y = fbgn.gn.mapping, by = RA_GENE_NAME)
#'   add_row_attr(loom = loom, key = "FBgn", value = genes$FBgn)
#' }

df<-data.frame("a"=c(2,2,3,4), "b"=c(1,2,3,4))
row.names(x = df)<-c("r1","r2","r3","r4")
e <- new("Embedding", name = "test", df)
e.df<-as.data.frame(x = as.matrix(x = e))