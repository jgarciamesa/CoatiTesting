library(seqinr)
library(stringr)

find_stop = function(filename) {
  fas = read.fasta(filename,set.attributes = FALSE,forceDNAtolower = FALSE)
  for(s in fas[]) {
    for(i in 1:(length(s)%/%3)) {
      if(str_c(s[(3*i-2):(3*i)],collapse = "") %in% c("TAG","TAA","TGA") & 3*i != length(s)) {
        print(paste("Found stop codon in",filename,",position",(3*i-2),". Length",length(s)))
      }
    }
  }
}

if(!interactive()) {
  ARGS = commandArgs(trailing=TRUE)
  find_stop(ARGS)
}
