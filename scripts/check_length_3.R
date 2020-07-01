library(seqinr)

ARGS=commandArgs(trailing = TRUE)
f = read.fasta(ARGS, set.attributes = FALSE)
if(length(f[[1]])%%3!=0 || length(f[[2]])%%3!=0){print(ARGS)}
