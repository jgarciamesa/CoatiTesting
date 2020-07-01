# gap3.sh
# Calculate proportion of indel phases of length 3 gaps
# Usage: bash gap3.sh directory

rm -f "data/gap3_phase_$(basename $1).txt"
for line in $(grep -r -h "[A,C,G,T]\-\-\-[A,C,G,T]\|[A,C,G,T]\-\-[\-]$" "$1")
do
	echo `expr index $line ---` >> "data/gap3_phase_$(basename $1).txt"
done
