# create a file with the name of the pairwise aligned ENSEMBL sequences
#  that contain gaps

if [ $# -lt 2 ]
then
	echo "At least two arguments are required"
fi

species="$1"

shift 1

mkdir -p "data/${species}/"
rm -f "data/${species}/gaps_${species}_tmp.csv"

for arg in $@
do
	for file in aln/${arg}/*.fasta
	do
		c=$(grep -c '-' ${file})
		if [ ${c} -gt 0 ]
		then
			echo $(basename ${file}) >> "data/${species}/gaps_${species}_tmp.csv"
		fi
	done
done

cat "data/${species}/gaps_${species}_tmp.csv"| sort -n | uniq > "data/${species}/gaps_${species}.csv"
rm "data/${species}/gaps_${species}_tmp.csv"

