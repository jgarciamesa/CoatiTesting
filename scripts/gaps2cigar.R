library(seqinr)
library(stringr)


gaps_cigar = function(gap_table, line, model) {
	fasta = read.table(gap_table, header=FALSE, sep="\n", stringsAsFactors=FALSE, skip=as.integer(line)-1, nrows=1)
	
	if(file.exists(paste("aln/",model,"/",fasta,sep=""))) {
		seqs = read.fasta(paste("aln/",model,"/",fasta,sep=""), set.attributes = TRUE, forceDNAtolower = FALSE)
	} else if(file.exists(paste("aln/coati/",fasta,sep=""))) {
		seqs = read.fasta(paste("aln/coati/",fasta,sep=""), set.attributes = TRUE, forceDNAtolower = FALSE)
	} else {
		f = list.files(".",fasta,recursive=TRUE,include.dirs=TRUE)[1]
		seqs = read.fasta(fasta,set.attributes = TRUE, forceDNAtolower = FALSE)
	}
	s1 = getSequence(seqs)[[1]]
	s2 = getSequence(seqs)[[2]]

	cigar = "" 	# cigar string
	current = ""
	nxt = ""
	count = 0

	stopifnot(length(s1) == length(s2))

	for(i in 1:length(s1)) {
		 if(s1[i] == '-'){
			if(current == "I") {
				count = count + 1
				next
			} else {
			  nxt = "I"
			  }
		} else if(s2[i] == '-') {
			if(current == "D") {
				count = count + 1
				next
			} else {
			  nxt = "D"
			  }
		} else if((s1[i] != '-') & (s2[i] != '-')) {
			if(current == "M") {
				count = count + 1
				next
			} else {
			  nxt = "M"
			  }
		}

	  	if(count > 0) {
			cigar = paste(cigar,paste(count,current, sep = ""), sep = "")
	  	}
		current = nxt
		count = 1
	  
	}
	
	cigar = paste(cigar,paste(count,current,sep=""),sep="")
	print(cigar)

}

if(!interactive()) {
	ARGS = commandArgs(trailing = TRUE)
	gaps_cigar(ARGS[1], ARGS[2], ARGS[3])
}
