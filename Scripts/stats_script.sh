#!/bin/bash 
SCRIPT=/home/mohit/Desktop/Thesis/plottingcollectiveproject/Testing

alg=$1

for f in $alg/*; do
	dest=$(echo $f| cut -d "/" -f 2)
#	echo $dest
#	echo plot_data_$f\.txt
	Rscript $SCRIPT/calculate_stats.R $f plot_data_$dest\.txt
	mv plot_data_$dest\.txt stats_$alg/
done
	
	
