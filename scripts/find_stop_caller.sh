mv -f ../data/stop.txt ../data/stop.txt.bak
for dir in $(ls ../coati/aln/ref/)
do
	for filename in $(ls ../coati/aln/ref/${dir}/)
	do
		Rscript --vanilla find_stop_codon.R ../coati/aln/ref/${dir}/${filename} >> ../data/stop.txt
	done
done

