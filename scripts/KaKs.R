library(seqinr)
library(stringr)

#################### code modified from Reed A Cartwright ######################
ka = function(a) { if(class(a) != "alignment") {return(NA)}
	ka = kaks(a)$ka[1]
	ka
}

ks = function(a) {
	if(class(a) != "alignment") {return(NA)}
	ks = kaks(a)$ks[1]
	ks
}

omega = function(a) { #dN/dS
	k = kaks(a)
	ka = k$ka[1]
	ks = k$ks[1]
	if(ka < 0 || ks < 0) {
		return(0)
	} else if(ka == 0 && ks == 0) {
		return(0)
	} else if(ks == 0) {
		return(10)
	} else {
		return(ka/ks)
	}
}

ps_accuracy = function(ref,w) {
	tp = length(intersect(which(w>1),which(ref>1)))
	fp = length(intersect(which(w>1),which(ref<1)))
	fn = length(intersect(which(w<1),which(ref>1)))
	return(2*tp/(2*tp+fp+fn))
}

ns_accuracy = function(ref,w) {
	tp = length(intersect(which(w<1),which(ref<1)))
	fp = length(intersect(which(w<1),which(ref>1)))
	fn = length(intersect(which(w>1),which(ref<1)))
	return(2*tp/(2*tp+fp+fn))
}

process_aln = function(d,fasta) {
	ret = list()

	for(f in fasta) {
		aln = try(read.alignment(str_c(d,"/",f),"fasta"),silent=TRUE)
		if(class(aln) == "try-error") {
			ret[[f]] = list()
		} else {
			seq = aln$seq %>% str_split(pattern="")
			ung = (seq[[1]] %in% c("a","c","g","t","n"))
			seq = lapply(seq,subset,ung) # subset = "["
			aln$seq = lapply(seq,str_c,collapse="")
			ret[[f]] = aln
		}
	}
	ret
}

################################################################################

k_metric = function() {#filename) {
	# read alignments
	#ref = list.files("data/ref_alignments",pattern = "*.fasta")
	#marg = list.files("coati/aln/ref/coati_marg",pattern="*.fasta")
	#toy = list.files("coati/aln/ref/coati_toy",pattern="*.fasta")
	#dna = list.files("coati/aln/ref/coati_dna",pattern="*.fasta")
	#mafft = list.files("coati/aln/ref/mafft",pattern="*.fasta")
	prank = list.files("coati/aln/ref/prank",pattern="*.fasta")
	#clustalo = list.files("coati/aln/ref/clustalo",pattern="*.fasta")

	# remove insertions w.r.t. model organism (human)
	ref = process_aln("data/ref_alignments",prank)
	marg = process_aln("coati/aln/ref/marg",prank)
	toy = process_aln("coati/aln/ref/toy",prank)
	dna = process_aln("coati/aln/ref/dna",prank)
	mafft = process_aln("coati/aln/ref/mafft",prank)
	clustalo = process_aln("coati/aln/ref/clustalo",prank)
	prank = process_aln("coati/aln/ref/prank",prank)

	# calculate ka
	ref_ka = sapply(ref,ka)
	marg_ka = sapply(marg,ka)
	toy_ka = sapply(toy,ka)
	dna_ka = sapply(dna,ka)
	mafft_ka = sapply(mafft,ka)
	prank_ka = sapply(prank,ka)
	clustalo_ka = sapply(clustalo,ka)

	# calculate ks
	ref_ks = sapply(ref,ks)
	marg_ks = sapply(marg,ks)
	toy_ks = sapply(toy,ks)
	dna_ks = sapply(dna,ks)
	mafft_ks = sapply(mafft,ks)
	prank_ks = sapply(prank,ks)
	clustalo_ks = sapply(clustalo,ks)

	# calculate omega
	ref_omega = sapply(ref,omega)
	marg_omega = sapply(marg,omega)
	toy_omega = sapply(toy,omega)
	dna_omega = sapply(dna,omega)
	mafft_omega = sapply(mafft,omega)
	prank_omega = sapply(prank,omega)
	clustalo_omega = sapply(clustalo,omega)

	# calculate root mean-square error of ka
	marg_rmse_ka = sqrt(sum((ref_ka-marg_ka)^2)/length(ref_ka))
	toy_rmse_ka = sqrt(sum((ref_ka-toy_ka)^2)/length(ref_ka))
	dna_rmse_ka = sqrt(sum((ref_ka-dna_ka)^2)/length(ref_ka))
	mafft_rmse_ka = sqrt(sum((ref_ka-mafft_ka)^2)/length(ref_ka))
	prank_rmse_ka = sqrt(sum((ref_ka-prank_ka)^2)/length(ref_ka))
	clustalo_rmse_ka = sqrt(sum((ref_ka-clustalo_ka)^2)/length(ref_ka))

	# calculate root mean-square error of ks
	marg_rmse_ks = sqrt(sum((ref_ks-marg_ks)^2)/length(ref_ks))
	toy_rmse_ks = sqrt(sum((ref_ks-toy_ks)^2)/length(ref_ks))
	dna_rmse_ks = sqrt(sum((ref_ks-dna_ks)^2)/length(ref_ks))
	mafft_rmse_ks = sqrt(sum((ref_ks-mafft_ks)^2)/length(ref_ks))
	prank_rmse_ks = sqrt(sum((ref_ks-prank_ks)^2)/length(ref_ks))
	clustalo_rmse_ks = sqrt(sum((ref_ks-clustalo_ks)^2)/length(ref_ks))

	# calculate accuracy of positive selection
	marg_pos = ps_accuracy(ref_omega,marg_omega)
	toy_pos = ps_accuracy(ref_omega,toy_omega)
	dna_pos = ps_accuracy(ref_omega,dna_omega)
	mafft_pos = ps_accuracy(ref_omega,mafft_omega) 
	prank_pos = ps_accuracy(ref_omega,prank_omega)
	clustalo_pos = ps_accuracy(ref_omega,clustalo_omega)

	# calculate accuracy of negative selection
	marg_neg = ns_accuracy(ref_omega,marg_omega)
	toy_neg = ns_accuracy(ref_omega,toy_omega)
	dna_neg = ns_accuracy(ref_omega,dna_omega)
	mafft_neg = ns_accuracy(ref_omega,mafft_omega) 
	prank_neg = ns_accuracy(ref_omega,prank_omega)
	clustalo_neg = ns_accuracy(ref_omega,clustalo_omega)

	print("dna, toy, marg, prank, mafft, clustalo")
	print(paste("ka root mean-squared error:",dna_rmse_ka,toy_rmse_ka,marg_rmse_ka,prank_rmse_ka,mafft_rmse_ka,clustalo_rmse_ka))
	print(paste("ks root mean-squared error:",dna_rmse_ks,toy_rmse_ks,marg_rmse_ks,prank_rmse_ks,mafft_rmse_ks,clustalo_rmse_ks))
	print(paste("Accuracy of + selection:",dna_pos,toy_pos,marg_pos,prank_pos,mafft_pos,clustalo_pos))
	print(paste("Accuracy of - selection:",dna_neg,toy_neg,marg_neg,prank_neg,mafft_neg,clustalo_neg))

}

if(!interactive()) {
	#ARGS = commandArgs(trailing=TRUE)
	k_metric()#ARGS)
}


