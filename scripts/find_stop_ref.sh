mv -f ../data/stop_ref.txt ../data/stop_ref.txt.bak
for filename in $(ls ../data/ref_alignments/)
do
	Rscript --vanilla find_stop_codon.R ../data/ref_alignments/${filename} >> ../data/stop_ref.txt
done

