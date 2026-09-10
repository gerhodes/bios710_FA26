## parser function for extracting CDS locations
parse_cds <- function(x) {
  ## strings we need to parse will be something like:
  ### "[17594092-17594331, 17598616-17598823, 17604859-17604996, 17608971-17609143]" 
  ## removes brackets
  x <- gsub("\\[|\\]", "", x)
  ## split on commas
  pieces <- strsplit(x, ",\\s*")[[1]]
  ## loop through each interval
  do.call(
    rbind,
    lapply(pieces, function(y) {
      ## split on dash
      coords <- strsplit(y, "-")[[1]]
      ## separate into start and end columns
      data.frame(
        start = as.integer(coords[1]),
        end   = as.integer(coords[2])
      )
    })
  )
}
