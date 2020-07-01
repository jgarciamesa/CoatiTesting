# create a file with the name of the pairwise aligned ENSEMBL sequences
#  that DON'T contain gaps

if [ $# -lt 1 ]
then
	echo "At least one argument is required"
fi

species=$1

rm -f data/${species}/nogaps_${species}.csv

for file in raw_data/${species}/*
do
	c=$(grep -c $(basename ${file}) data/${species}/gaps_${species}.csv)
	if [ ${c} -eq 0 ]
	then
		echo $(basename ${file}) >> data/${species}/nogaps_${species}.csv
	fi
done

