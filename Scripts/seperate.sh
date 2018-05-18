#/bin/bash

# Run where mpi output files (.out) are
# data_R.awk file called with size of message and algorithm
# give <file> as $1 if filename=file.out

#SCRIPT=/media/storage/jagpreet/Dropbox/Personal/Research/ETH/scripts_js
SCRIPT=/home/mohit/Desktop/Thesis/plottingcollectiveproject

alg=$1

for f in $alg/*; do
	for s in 1 2 3; do
		echo $f		
		awk -v S=$s -f $SCRIPT/sep.awk $f >> $f\_$s
		#mv $f\_$s $alg/
	done
done
