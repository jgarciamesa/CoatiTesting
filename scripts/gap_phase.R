# gap_phase.R
# Analize gap phases from pairwise alignments

library(seqinr)

gap_analysis = function(file_list,path) {
  
  total_gaps = 0       # total number of gaps
  total_gaps3 = 0      # total number of gaps w/ length multiple of 3
  phase = c(0,0,0)     # number of gaps on phase 2,0,1
  phase3 = c(0,0,0)    # number of gaps w/ length multiple of the on phase 2,0,1
  
  alignments = read.csv(file = file_list, header = FALSE)
  
  for(i in 1:nrow(alignments)) {
    
    aln = read.fasta(file = paste0(path,alignments[i,]), forceDNAtolower = FALSE)
    
    seqs = rbind(getSequence(aln)[[1]],getSequence(aln)[[2]])
    
    for(i in 1:2) {
      if(any(seqs[i,] == "-")) {
        # find indexes of gaps and group them by continuous positions (# of groups = # of gaps)
        pos_gaps = split(which(seqs[i,]=="-"),cumsum(c(1,diff(which(seqs[i,]=="-")) != 1)))
        
        # for each gap found
        for(j in 1:length(pos_gaps)) {
          if(length(pos_gaps[[j]]) %% 3 == 0){
            total_gaps3 = total_gaps3 + 1
            phase3[(pos_gaps[[j]][1]%%3)+1] = phase3[(pos_gaps[[j]][1]%%3)+1] + 1
          } else {
            total_gaps = total_gaps + 1
            phase[(pos_gaps[[j]][1]%%3)+1] = phase[(pos_gaps[[j]][1]%%3)+1] + 1
          }
        }
      }
    }
  }
  
  
  
  print(paste0("Total gaps w/ length multiple of 3: ", total_gaps3))
  print(paste0("Phase Zero: ", phase3[2], "    Phase One: ", phase3[3], "    Phase Two: ", phase3[1]))
  print(paste0("Total gaps w/ length NOT multiple of 3: ", total_gaps))
  print(paste0("Phase Zero: ", phase[2], "    Phase One: ", phase[3], "    Phase Two: ", phase[1]))
}

analysis_per_alignment = function(file_list,path1,path2) {
  
  total_gaps_1 = 0       # total number of gaps
  total_gaps3_1 = 0      # total number of gaps w/ length multiple of 3
  phase_1 = c(0,0,0)     # number of gaps on phase 2,0,1
  phase3_1 = c(0,0,0)    # number of gaps w/ length multiple of the on phase 2,0,1
  
  total_gaps_2 = 0       # total number of gaps
  total_gaps3_2 = 0      # total number of gaps w/ length multiple of 3
  phase_2 = c(0,0,0)     # number of gaps on phase 2,0,1
  phase3_2 = c(0,0,0)    # number of gaps w/ length multiple of the on phase 2,0,1

  alignments = read.csv(file = file_list, header = TRUE, sep = "\t")
  
  more_gaps = 0
  less_gaps = 0
  mcoati = 0
  ecm = 0

  for(i in 1:nrow(alignments)) {
    
    aln1 = read.fasta(file = paste0(path1,alignments[i,1]), forceDNAtolower = FALSE)
    aln2 = read.fasta(file = paste0(path2,alignments[i,1]), forceDNAtolower = FALSE)
    
    seqs1 = rbind(getSequence(aln1)[[1]],getSequence(aln1)[[2]])
    seqs2 = rbind(getSequence(aln2)[[1]],getSequence(aln2)[[2]])

	g1 = 0
	g2 = 0

    g1 = length(split(which(seqs1[1,]=="-"),cumsum(c(1,diff(which(seqs1[1,]=="-")) != 1))))
    g1 = g1 + length(split(which(seqs1[2,]=="-"),cumsum(c(1,diff(which(seqs1[2,]=="-")) != 1))))
    g2 = length(split(which(seqs2[1,]=="-"),cumsum(c(1,diff(which(seqs2[1,]=="-")) != 1))))
    g2 = g2 + length(split(which(seqs2[2,]=="-"),cumsum(c(1,diff(which(seqs2[2,]=="-")) != 1))))

    
#    for(k in 1:2) {
#      if(any(seqs1[k,] == "-", seqs2[k,] == "-")) {
#        # find indexes of gaps and group them by continuous positions (# of groups = # of gaps)
#        pos_gaps1 = split(which(seqs1[k,]=="-"),cumsum(c(1,diff(which(seqs1[k,]=="-")) != 1)))
#        pos_gaps2 = split(which(seqs2[k,]=="-"),cumsum(c(1,diff(which(seqs2[k,]=="-")) != 1)))
#
#		#TODO: this currently compares seq1 of aln1 with seq1 of aln2, then seq2 of aln1 and seq2 of aln2.
#		#      instead it should compare number of gaps per alignment, not seq to seq comparison
#        if(length(pos_gaps1) != length(pos_gaps2)) {
#
#		  g1 = g1 + length(pos_gaps1)
#		  g2 = g2 + length(pos_gaps2)
#        
#          # for each gap found in first alignment
#          for(j in 1:length(pos_gaps1)) {
#            if(length(pos_gaps1[[j]]) %% 3 == 0){
#              total_gaps3_1 = total_gaps3_1 + 1
#              phase3_1[(pos_gaps1[[j]][1]%%3)+1] = phase3_1[(pos_gaps1[[j]][1]%%3)+1] + 1
#            } else {
#              total_gaps_1 = total_gaps_1 + 1
#              phase_1[(pos_gaps1[[j]][1]%%3)+1] = phase_1[(pos_gaps1[[j]][1]%%3)+1] + 1
#            }
#          } 
#          
#          # for each gap found in second alignment
#          for(j in 1:length(pos_gaps2)) {
#            if(length(pos_gaps2[[j]]) %% 3 == 0){
#              total_gaps3_2 = total_gaps3_2 + 1
#              phase3_2[(pos_gaps2[[j]][1]%%3)+1] = phase3_2[(pos_gaps2[[j]][1]%%3)+1] + 1
#            } else {
#              total_gaps_2 = total_gaps_2 + 1
#              phase_2[(pos_gaps2[[j]][1]%%3)+1] = phase_2[(pos_gaps2[[j]][1]%%3)+1] + 1
#            }
#          } 
#        }
#      }
#    }
	if(g1 > g2) {
		#print(alignments[i,1])
		if(alignments[i,2] < alignments[i,3]){
			more_gaps = more_gaps + 1
			mcoati = mcoati + 1
			print("mcoati better with more gaps")
		} else {
			less_gaps = less_gaps + 1
			ecm = ecm + 1
			print("ecm better with less gaps")
		}
	} else if(g1 < g2) {
		#print(alignments[i,1])
		if(alignments[i,2] > alignments[i,3]){
			more_gaps = more_gaps + 1
			ecm = ecm + 1
			print("ecm better with more gaps")
		} else {
			less_gaps = less_gaps + 1
			mcoati = mcoati + 1
			print("mcoati better with less gaps")
		}
	}
	#if(g1 == g2) {
		#print(alignments[i,])
		#print(paste0("same length: ",g1))
	#}
  }
   print(paste0("more gaps: ",more_gaps,"; less gaps: ",less_gaps))
   print(paste0("mcoati was better: ", mcoati, "; ecm was better: ",ecm))
  
  
  print("--Analysis of alignments with different number of gaps--") 
  print(path1)
  print(paste0("Total gaps w/ length multiple of 3: ", total_gaps3_1))
  print(paste0("Phase Zero: ", phase3_1[2], "    Phase One: ", phase3_1[3], "    Phase Two: ", phase3_1[1]))
  print(paste0("Total gaps w/ length NOT multiple of 3: ", total_gaps_1))
  print(paste0("Phase Zero: ", phase_1[2], "    Phase One: ", phase_1[3], "    Phase Two: ", phase_1[1]))
  print("------------------------------------------------------------------------------------")
  print(path2)
  print(paste0("Total gaps w/ length multiple of 3: ", total_gaps3_2))
  print(paste0("Phase Zero: ", phase3_2[2], "    Phase One: ", phase3_2[3], "    Phase Two: ", phase3_2[1]))
  print(paste0("Total gaps w/ length NOT multiple of 3: ", total_gaps_2))
  print(paste0("Phase Zero: ", phase_2[2], "    Phase One: ", phase_2[3], "    Phase Two: ", phase_2[1]))
}


if(!interactive()) {
  ARGS = commandArgs(trailingOnly = TRUE)
  #gap_analysis(ARGS[1],ARGS[2])
  analysis_per_alignment(ARGS[1], ARGS[2], ARGS[3])
}
