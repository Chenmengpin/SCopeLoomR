
#####################
# Validators        #
#####################

setValidity(Class = "LoomColumnAnnotation", function(object) {
  if(length(x = unique(x = value)) > 245) {
    "Cannot add column annotation with more than 245 unique values."
  } else {
    TRUE
  }
})

setValidity("Embedding", function(object) {
  y.dims <- dim(x = object)[2]
  # Check if 2 Y-coordinates
  if (y.dims != 2) {
    paste("Embedding has Y-dimension length of ", y.dims, ". Should be 2", sep = "")
  } else if(any(is.na(x = object))) { # Check if any NA values
    "Embedding contains NA values. No NAs are allowed"
  } else {
    TRUE
  }
})