#/bin/bash

# Run where mpi output files (.out) are
# data_R.awk file called with size of message and algorithm
# give <file> as $1 if filename=file.out

#SCRIPT=/media/storage/jagpreet/Dropbox/Personal/Research/ETH/scripts_js
SCRIPT=/home/mohit/Desktop/Thesis/plottingcollectiveproject/Testing

alg=$1

for s in 1048576; do
	if [[ "$alg" == "MPI_BCAST" ]]
		then
			echo "nP runid wtime" > $alg\_$s
    
    		else 
			echo "nP chunks runid wtime" > $alg\_$s 
    	fi
    for f in sampleresult.out; do
	awk -v S=$s -v ALG=$alg -f $SCRIPT/data_R.awk $f >> $alg\_$s
	mv $alg\_$s $alg/
    done
done
