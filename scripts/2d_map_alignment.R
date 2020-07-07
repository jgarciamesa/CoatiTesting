# 2d_map_alignment.R
# Create a 2d map from two alignments
# Usage: Rscript --vanilla 2d_map_alignment.R alignment1.fasta alignment2.fasta

library(seqinr)

twoDmap = function(filename1, filename2) {
	# Read both alignments
	aln_temp1 = read.fasta(file = filename1, forceDNAtolower = FALSE)
	aln_temp2 = read.fasta(file = filename2, forceDNAtolower = FALSE)

	# Organize alignmnets into reference sequences (aln1) and second sequences (aln2)
	aln1 = aln2 = list()
	aln1[[1]] = getSequence(aln_temp1[[1]])
	aln1[[2]] = getSequence(aln_temp2[[1]])
	aln2[[1]] = getSequence(aln_temp1[[2]])
	aln2[[2]] = getSequence(aln_temp2[[2]])

	# Remove gaps for labeling the rows and columns on the map
	aln1_no_gaps = aln1[[1]][aln1[[1]] != "-"]
	aln2_no_gaps = aln2[[1]][aln2[[1]] != "-"]

	# Length of sequences without gaps
	L = length(aln1_no_gaps)
	l = length(aln2_no_gaps)

	# Initialize tex file which will contain map
	map = paste0(filename1, "_2d.tex")
	if(file.exists(map)) {
		file.remove(map)
	}

	f = file(map, open = "a" )

	# Write tex file config lines
	writeLines("\\documentclass[tikz]{standalone}", f)
	write("\\usepackage{tikz}", file=f,append=TRUE)
	write("\\usetikzlibrary{arrows}", file=f,append=TRUE)
	write("\\begin{document}", file=f,append=TRUE)
	write("\\begin{tikzpicture}", file=f,append=TRUE)
	write(paste0("\\draw[very thin, color=gray!50, step=0.5] (0,0) grid (", l/2+1, ", ",
	 L/2+1, ");"),file=f,append=TRUE)

	# Draw ref sequence
	for(i in seq(from = L/2+0.5, to = 1, by = -0.5)) {
		write(paste0("\\node at (",0.25,",",i-0.75,") {",
			aln1_no_gaps[L+1-(i-0.5)*2],"};"), file = f, append = TRUE)
	}

	# Draw second sequence
	for(i in seq(from = l/2+0.5, to = 1, by = -0.5)) {
		write(paste0("\\node at (",1.75+l/2-i,",",L/2+0.75,") {",
			aln2_no_gaps[l+1-(i-0.5)*2],"};"), file = f, append = TRUE)
	}


	# Set colors for different alignments
	colors = c("blue","orange")

	for(j in 1:2) {
		# Starting position to draw aligment
		x = 0.7 + (j-1) * 0.1
		y = L/2 + 0.2 + (j-1) * 0.1
		for(i in 1:length(aln1[[j]])) {
			if(aln1[[j]][i] == "-") {
				# insertion (-->)
				write(paste0("\\draw[->,",colors[j],", thick] (",x,",",y,") -- (",
					x+0.5,",",y,");"), file = f, append = TRUE)
				x = x + 0.5
			} else if(aln2[[j]][i] == "-") {
				# deletion (down)
				write(paste0("\\draw[->,",colors[j],", thick] (",x,",",y,") -- (",
					x,",",y-0.5,");"), file = f, append = TRUE)
				y = y - 0.5
			} else {
				# match/mismatch
				write(paste0("\\draw[->,",colors[j],", thick] (",x,",",y,") -- (",
					x+0.5,",",y-0.5,");"), file = f, append = TRUE)
				x = x + 0.5
				y = y - 0.5
			}

		}
	}


	# Draw legend and end latex document
	write("\\node at (1.3,-0.5) {\\textcolor{blue}{Marginal Coati}};", file=f,append=TRUE)
	write("\\node at (0.55,-1) {\\textcolor{orange}{ECM}};", file=f,append=TRUE)
	write("\\end{tikzpicture}", file=f,append=TRUE)
	write("\\end{document}", file=f,append=TRUE)

	# Close connection to file
	close(f)
}






if(!interactive()) {
	ARGS = commandArgs(trailingOnly = TRUE)
	twoDmap(ARGS[1], ARGS[2])
}
