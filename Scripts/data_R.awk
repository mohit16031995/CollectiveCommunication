# -v ALG is the name of the algorithm for which data need to be extracted
# -v S is size of message for the plot and grouping of data
BEGIN {
    temp = ALG;
    s = S;
    flag = 0;
    nP = 0;
    chunks = 0;
}
{   
    if($1==temp && $3==s) {
	if (temp=="MPI_BCAST") {
		flag=1; nP=$2;
	}
	else {
		flag=1;
		nP=$2;
		chunks=$4;
	}	
	
    }
    else if($1=="Run" && flag) {
	if (temp=="MPI_BCAST") {
		print nP" "$2" "$4
	}
	else {
		print nP" "chunks" "$2" "$4
	}
    }
    else if ($1=="Tree") {
		
    }
    else if ($1=="Logs,") {
	
    }
   else {
	flag=0;
   }
}
END {
}
