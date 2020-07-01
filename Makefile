
.PHONY: download_geneid download_fasta gap_cigar reference d_seq no_gaps_reference

RSCRIPT = Rscript --vanilla

SHELL:=/bin/bash

MODEL = coati
SPECIES = gorilla

################################################################################
# Download gene ID table from ensembl                                          #
################################################################################

download_geneid: scripts/get_geneId.R
	$(RSCRIPT) $<

################################################################################
# Download first N sequences from ENSEMBL (if small size)                      #
################################################################################
N = 20
GENE_LIST = ./raw_data/$(SPECIES)_geneId.tsv
LIST_ALL := $(shell echo {1..$(N)})
FASTA = $(addsuffix .fasta,$(addprefix raw_data/fasta_files/,$(LIST_ALL)))

download_fasta: $(FASTA)
	@mv raw_data/fasta_files/* raw_data/$(SPECIES)/

# download fasta files from ENSEMBL
raw_data/fasta_files/%.fasta: scripts/get_sequences.R
	@echo "Downloading "$*" out of $(N)."
	@$(shell $(RSCRIPT) $< $(GENE_LIST) $(shell dirname $@) $*)

################################################################################
# Create cigar strings with indel types for simulation                         #
################################################################################

GAPS_FILE = ./data/$(SPECIES)/gaps_$(SPECIES).csv
GAPS_N := $(shell cat $(GAPS_FILE) | wc -l | cut -f1)
GAPS_LIST := $(shell echo {1..$(GAPS_N)})
GAPS := $(addsuffix .gap,$(GAPS_LIST))

gap_cigar: scripts/gaps2cigar.R $(GAPS)

# extract gap information and encode it using CIGAR strings
%.gap: scripts/gaps2cigar.R
	@$(shell $(RSCRIPT) $< $(GAPS_FILE) $* ${MODEL} | cut -d '"' -f 2 >> data/$(SPECIES)/gaps_$(SPECIES)_cigar.csv)

################################################################################
# Generate alignments with biologically inferred indel information (cigar str) #
################################################################################

CIGAR = ./data/$(SPECIES)/gaps_$(SPECIES)_cigar.csv
SEQS := $(shell cat ./data/$(SPECIES)/nogaps_$(SPECIES).csv)
REF_ALIG := $(addprefix data/ref_alignments/,$(addsuffix .ref,$(SEQS)))

reference: $(CIGAR) $(REF_ALIG)
	@mkdir -p data/$(SPECIES)/ref_alignments
	@mv data/ref_alignments/* data/$(SPECIES)/ref_alignments

data/ref_alignments/%.ref: scripts/simulate2.R scripts/write_fasta.R
	@echo $*
	$(shell $(RSCRIPT) $< $(SPECIES) raw_data/$(SPECIES)/$* > data/ref_alignments/$* )

################################################################################
# Create reference alignments with no gaps for testing                         #
################################################################################

no_gaps_reference:
	$(shell bash scripts/rem_gaps_ref.sh $(SPECIES))

################################################################################
# Clean pipeline results except gene id list and raw fasta downloads		   #
################################################################################

.PHONY: clean_pipeline

clean_pipeline:
	rm -f data/$(SPECIES)/*.csv
	rm -f data/$(SPECIES)/no_gaps_ref/*
	rm -f data/$(SPECIES)/ref_alignments/*
	rm -f aln/{coati,mcoati,dna,ecm,mecm,prank,mafft,clustalo}/*
	rm -f aln/ref/{coati,mcoati,dna,ecm,mecm,prank,mafft,clustalo}/*
	rm -f data/dseq.csv
	rm -f data/dseq_summary.csv

