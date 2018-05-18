# -v ALG is the name of the algorithm for which data need to be extracted
# -v S is size of message for the plot and grouping of data
BEGIN {
    flag = 0;
    s = S;
    print "nP chunks runid wtime";		
}
{   
    if($3==1) {	
        if(flag==3)	
            flag = 0;
	flag=flag+1;
    }
    if(flag==s)
        print $0;
}
END {
}
