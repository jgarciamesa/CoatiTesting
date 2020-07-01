# Get sequences from Ensembl

library(httr)
library(jsonlite)
library(xml2)
library(seqinr)
library(stringr)

get_canonical_transcript = function(request) {
	max_l = 2001/3 # max length in aa
	r = fromJSON(toJSON(content(request)))
	if(r$Transcript$Translation$length[r$Transcript$is_canonical==1] > max_l) {
		return("MAX_LENGTH_EXCEPTION")
	} else{
		return(r$Transcript[["id"]][r$Transcript$is_canonical==1])
	}
}

getseq_main = function(gene_table, out_folder, line) {
	server = "https://rest.ensembl.org/"

	geneID = read.table(gene_table, header=TRUE, sep="\t", stringsAsFactors=FALSE, skip=as.integer(line)-1, nrows=1)

	sequences = c()
	id_seq = array(dim=ncol(geneID))
	# once the new get_geneId is run, remove `-1`
	for(i in 1:ncol(geneID)) {
		req = GET(paste(server,"lookup/id/",geneID[i],"?expand=1",sep = ""))
		stop_for_status(req)

		id = get_canonical_transcript(req)
		if(id == "MAX_LENGTH_EXCEPTION") {
			break
		}
		# improve appending for efficiency
		id_seq[i] = id
		req = GET(paste(server,"sequence/id/",id,"?type=cds",sep = ""),
			content_type("text/x-fasta"))
		stop_for_status(req)
		sequences = append(sequences, content(req))
	}
	if(id != "MAX_LENGTH_EXCEPTION") {
		write.table(sequences, file = paste(out_folder,"/",geneID[1],".fasta",sep=""),
			quote = FALSE, row.names = FALSE, col.names = FALSE)
	}
}



if(!interactive()) {
	ARGS = commandArgs(trailing=TRUE)
	getseq_main(gene_table=ARGS[1], out_folder=ARGS[2], line=ARGS[3])
}
