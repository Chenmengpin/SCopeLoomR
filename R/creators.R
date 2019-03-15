
#####################
# Creators          #
#####################

CreateSCopeLoom <- function(
  file.path = NULL
  , title = NULL
  , genome = NULL
  , dgem
  , embedding = NULL
  , embedding.name = NULL
  , tree = NULL
  , chunk.size = 1000
  , display.progress = T
) {
  loom <- new(Class = "SCopeLoom")
  # Set the file path
  loom@file.path <- file.path
  
  cn <- colnames(x = dgem)
  rn <- row.names(x = dgem)
  
  ##########
  # Matrix #
  ##########
  # Check the type matrix
  # Convert to dgCMatrix if necessary to speedup populating the matrix slot
  if(class(x = dgem) %in% c("matrix", "dgTMatrix")) {
    print("Converting to dgCMatrix...")
    dgem <- methods::as(object = dgem, Class = "dgCMatrix")
  }
  print("Adding main matrix...")
  loom@matrix <- new(Class = "MainMatrix", dgem)
  
  #####################
  # Global Attributes #
  #####################
  loom <- AddGlobalAttribute(object = loom, name = "Title", data = title)
  loom <- AddGlobalAttribute(object = loom, name = "Genome", data = genome)
  loom <- AddGlobalAttribute(object = loom, name = "SCopeTree", data = tree)
  loom <- AddGlobalAttribute(object = loom, name = "SCopeLoomR.version", data = as.character(x = packageVersion("SCopeLoomR")))
  loom <- AddGlobalAttribute(object = loom, name = "R.version", data = R.version.string)
  loom <- AddGlobalAttribute(object = loom, name = "CreationDate", data = as.character(x = Sys.time()))
  loom <- AddGlobalAttribute(object = loom, name = "MetaData", data = new(Class = "MetaData"))

  #####################
  # Column Attributes #
  #####################
  ## CellID
  print("Adding CellID column attribute...")
  loom <- AddColumnAttribute(object = loom, name = "CellID", data = cn)
  
  ## Embeddings (Add the default embedding)
  if(!is.null(x = embedding)) {
    loom <- AddEmbedding(object = loom, name = embedding.name, data = embedding, is.default = TRUE)
  } else {
    warning("No default embedding set for the loom. You'll not be able to visualize it in SCope.")
  }
  
  ## nUMI
  loom <- AddNbUMI(object = loom)
  ## nGenes
  loom <- AddNbGene(object = loom)
  
  #####################
  # Row Attributes    #
  #####################
  
  ## Gene
  loom <- AddRowAttribute(object = loom, name = "Gene", data = rn)

  # Options
  loom@options[["chunk.size"]] <- chunk.size
  loom@options[["display.progress"]] <- display.progress
  
  return (loom)
}

CreateSCopeTree <- function(level.1
                            , level.2
                            , level.3) {
  tree <- new(Class = "SCopeTree"
              , L1 = new(Class = "LoomSCopeTreeLevel", key = "SCopeTreeL1", level.1)
              , L2 = new(Class = "LoomSCopeTreeLevel", key = "SCopeTreeL2", level.2)
              , L3 = new(Class = "LoomSCopeTreeLevel", key = "SCopeTreeL3", level.3))
  return (tree)
}

CreateLoomGlobalAttribute <- function(object
                               , name
                               , data) {
  if(is.character(x = data)) {
    return (new(Class = "LoomGlobalAttribute", key = name, data))
  } else {
    return (data)
  }
}

AddGlobalAttribute <- function(object
                               , name
                               , data) {
  object@attrs[[name]] <- CreateLoomGlobalAttribute(object = object, name = name, data = data)
  return (object)
}

AddRowAttribute <- function(object, name, data) {
  object@row.attrs[[name]] <- CreateRowAttribute(name = name, data = data)
  return (object)
}

CreateRowAttribute <- function(name, data) {
  ra <- data.frame(data, stringsAsFactors = F)
  colnames(x = ra) <- name
  return (new(Class = "RowAttribute", key = name, ra))
}

AddAnnotation <- function(object, name, data) {
  # Add column attribute
  object@col.attrs[[name]] <- new(Class = "ColumnAnnotation", key = name, data)
  # Add associated global attribute
  gmd.annotations <- slot(object = object@attrs$MetaData, "annotations")
  gmd.annotations[[length(x = gmd.annotations)+1]] <- as.list(as.character(x = unique(x = values)))
  slot(object = object@attrs$MetaData, "annotations") <- gmd.annotations
  return (object)
}

AddNbGene <- function(object) {
  # Check if matrix is raw counts
  print("Adding default metrics nGene...")
  n.gene <- Matrix::colSums(x = object@matrix > 0)
  n.gene <- data.frame("n.Gene" = n.gene)
  object <- AddMetric(object = object, name = "nGene", data = n.gene)
  return (object)
}

AddNbUMI <- function(object) {
  # Check if matrix is raw counts
  if(!any(object@matrix%%1!=0)) {
    print("Adding default metrics nUMI...")
    n.umi <- Matrix::colSums(x = object@matrix)
  } else {
    warning("Default metric nUMI was not added because the input matrix does not seem to be UMI counts.")
  }
  n.umi <- data.frame("n.UMI" = n.umi)
  object <- AddMetric(object = object, name = "nUMI", data = n.umi)
  return (object)
}

AddMetric <- function(object, name, data) {
  # Add column attribute
  ca <- new(Class = "ColumnAttribute", data)
  object@col.attrs[[name]] <- new(Class = "ColumnMetric", key = name, ca)
  # Add associated global attribute
  gmd.metrics <- slot(object = object@attrs$MetaData, "metrics")
  gmd.metrics[[length(x = gmd.metrics)+1]]<-list(name = name)
  slot(object = object@attrs$MetaData, "metrics") <- gmd.metrics
  return (object)
}

AddEmbedding <- function(object, name, data, is.default = FALSE) {
  if(is.default) {
    print("Adding default embedding...")
  }
  if(is.null(x = name)) {
    name<-"Embedding"
  }
  if(!("Embeddings" %in% names(x = object@col.attrs))) {
    embeddings <- new(Class = "Embeddings")
  } else {
    embeddings <- object@col.attrs$Embeddings
  }
  embeddings[[name]] <- new(
    Class = "Embedding"
    , name = name
    , default = is.default
    , data)
  object@col.attrs[["Embeddings"]] <- embeddings
  return (object)
}

AddColumnAttribute <- function(object, name, data) {
  object@col.attrs[[name]] <- CreateColumnAttribute(name = name, data = data)
  return (object)
}

CreateColumnAttribute <- function(name, data) {
  ca <- data.frame(data, stringsAsFactors = F)
  colnames(x = ca) <- name
  return (new(Class = "ColumnAttribute", key = name, ca))
}


