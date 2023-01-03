slideFunct <- function(data, window, step){
  total <- length(data)
  spots <- seq(from=1, to=total, by=step)
  as.data.frame(pbsapply((1+window):(length(spots)-window),function(w) {mean(as.numeric(na.omit(data[spots[w-window]:(spots[w]+window)]))) }),stringsAsFactors=FALSE)
  #return(result)
}

