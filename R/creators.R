
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
  
  cn<-colnames(dgem)
  rn<-row.names(dgem)
  
  ##########
  # Matrix #
  ##########
  # Check the type matrix
  # Convert to dgCMatrix if necessary to speedup populating the matrix slot
  if(class(x = dgem) %in% c("matrix", "dgTMatrix")) {
    print("Converting to dgCMatrix...")
    dgem <- methods::as(object = dgem, Class = "dgCMatrix")
  }
  loom@matrix <- new(Class = "MainMatrix", dgem)
  
  #####################
  # Global Attributes #
  #####################
  loom <- AddLoomGlobalAttribute(object = loom, name = "title", data = title)
  loom <- AddLoomGlobalAttribute(object = loom, name = "Genome", data = genome)
  loom <- AddGlobalAttribute(object = loom, name = "tree", data = tree)
  loom <- AddLoomGlobalAttribute(object = loom, name = "SCopeLoomR.version", data = as.character(x = packageVersion("SCopeLoomR")))
  loom <- AddLoomGlobalAttribute(object = loom, name = "R.version", data = as.character(R.version.string))
  loom <- AddLoomGlobalAttribute(object = loom, name = "CreationDate", data = as.character(Sys.time()))
  loom <- AddGlobalAttribute(object = loom, name = "MetaData", data = new(Class = "MetaData"))

  #####################
  # Column Attributes #
  #####################
  
  ## CellID
  loom@col.attrs[["CellID"]] <- new(Class = "LoomColumnAttribute", key = "CellID", t(x = as.matrix(x = cn)))
  
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
  loom@row.attrs[["Gene"]] <- new(Class = "LoomRowAttribute", key = "Gene", as.matrix(x = rn))
  
  # Options
  loom@options[["chunk.size"]] <- chunk.size
  loom@options[["display.progress"]] <- display.progress
  
  return (loom)
}

CreateSCopeTree <- function(level.1
                            , level.2
                            , level.3) {
  tree <- new(Class = "SCopeTree"
              , L1 = new(Class = "SCopeTreeLevel", key = "SCopeTreeL1", level.1)
              , L2 = new(Class = "SCopeTreeLevel", key = "SCopeTreeL2", level.2)
              , L3 = new(Class = "SCopeTreeLevel", key = "SCopeTreeL3", level.3))
  return (tree)
}

AddLoomGlobalAttribute <- function(object
                               , name
                               , data) {
  object@attrs[[name]] <- new(Class = "LoomGlobalAttribute", key = name, data)
  return (object)
}

AddGlobalAttribute <- function(object
                               , name
                               , data) {
  object@attrs[[name]] <- data
  return (object)
}

AddAnnotation <- function(object, name, data) {
  # Add column attribute
  object@col.attrs[[name]] <- new(Class = "LoomColumnAnnotation", key = name, as.character(x = t(x = as.matrix(x = data))))
  # Add associated global attribute
  gmd.annotations <- slot(object = object@attrs$MetaData, "annotations")
  gmd.annotations[[length(x = gmd.annotations)+1]]<-as.list(as.character(unique(values)))
  slot(object = object@attrs$MetaData, "annotations") <- gmd.annotations
  return (object)
}

AddNbGene <- function(object) {
  # Check if matrix is raw counts
  print("Adding default metrics nGene...")
  n.gene <- Matrix::colSums(x = object@matrix > 0)
  object <- AddMetric(object = object, name = "nGene", data = t(x = n.gene))
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
  object <- AddMetric(object = object, name = "nUMI", data = t(x = n.umi))
  return (object)
}

AddMetric <- function(object, name, data) {
  # Add column attribute
  object@col.attrs[[name]] <- new(Class = "LoomColumnMetric", key = name, data)
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


