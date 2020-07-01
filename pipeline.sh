## TODO: change species and model here and on all 3 Makefiles
species=gorilla
model=mafft #coati mcoati dna ecm mecm prank mafft clustalo
models=(coati mcoati dpmcoati dna ecm mecm prank mafft clustalo)

# download geneid list
#  modify scripts/get_geneID.R with desired species
 @tput setaf 11; echo Downloading ${species} gene id list
make download_geneid

# download sequences
#  modify Makefile with desired species
 @tput setaf 11; echo Downloading ${species} fasta sequences 
make download_fasta

# align to identify which ones have gaps
tput setaf 11; echo Initial alignment of ${species} fasta sequences
tput setaf 15
for m in ${models[*]}
do
	echo Initial alignment with ${m}
	make -f Makefile_aln -j8 align_${m}
done

# identify which alignments have gaps
tput setaf 11; echo Identify which alignments have gaps
tput setaf 15
bash scripts/gaps.sh ${species} ${models[*]}

# create gap CIGAR strings
tput setaf 11; echo Create biologically-like gap maps
tput setaf 15
make gap_cigar

# create list of no gap sequences to be converted into reference alignments
tput setaf 11; echo Create list of sequences with no gaps to become reference alignments
tput setaf 15
bash scripts/nogaps.sh ${species}

# create reference alignments
tput setaf 11; echo Introduce gaps and make reference alignments
tput setaf 15
make reference

# remove gaps from reference
tput setaf 11; echo Remove gaps from true alignments generating benchmark data set
tput setaf 15
make no_gaps_reference
#for file in data/${species}/no_gaps_ref/*; do mv $file data/${species}/no_gaps_ref/$(basename $file .ref); done

# align reference alignments
tput setaf 11; echo Align benchmark data set sequences
for m in ${models[*]}
do
	echo Aligning with ${m} model
	make -f Makefile_aln -j8 reference_${m}
done
tput setaf 15

# compute statistics
tput setaf 11; echo TODO: compute statistics
tput setaf 15
make -f Makefile_stats d_seq
R --vanilla --slave -e 'data=read.csv("data/dseq.csv");cat(sum(data[,2])/length(data[,2]))' >> data/dseq_summary.csv
Rscript --vanilla scripts/number_alignments.R
Rscript --vanilla scripts/kaks.R ${species} ${models[*]}
