#/bin/bash

FILE=/home/mohit/Desktop/simulator/goal_2Tree_reduce.c
IFS=$'\n' read -d '' -r -a lines < daintpy.out
IFS=$'\n' read -d '' -r -a lines3 < readValues.txt

L=1623300
#L=1082200
#L=1059947

O=100
#O=66
G=100
#G=66
#G=103
#echo "Binary Pipeline" >> sim_binary.txt
# 65536 131072 262144 524288 1048576 3145728 5242880 7340032 9437184
# 32 64 128 188 256 366 512 
#for s in 1048576; do
	#for p in 32 64 128 188 256 366 512; do	
	for k in "${lines3[@]}"; do
		s=$(echo $k | cut -d " " -f 2)
		p=$(echo $k | cut -d " " -f 1)
		ch=$(echo $k | cut -d " " -f 3)
		mintime=1000000.0
		optimal=0
		prev=999999
		echo "s = $s p = $p"
		for ch in {1..70}; do
			chunks=$ch
	#	for i in "${lines[@]}"; do
	#		cs=$(echo $i | cut -d " " -f 1)
	#		ovalue=$(echo $i | cut -d " " -f 2)
	#		y=10000
	#		o=$(echo $ovalue*$y | bc)
	#		o=${o%.*}
			o=4262100
			#o=2841400
			g=$o
	#		chunks=$(($s / $cs))
	#		if [ $chunks == 0 ]; then
	#			break
	#		fi
	#		if [ $chunks == $prev ]; then
	#			continue			
	#		fi
	#		echo "chunks = $chunks"
	#		prev=$chunks
	#		if (("$chunks" == "$ch")); then
				echo "$s $p"
				gcc $FILE -lm && ./a.out $p $s $chunks > output.goal
				LogGOPSim-master/tests/testsim/txt2bin -i output.goal -o out.bin
				LogGOPSim-master/tests/testsim/LogGOPSim -L $L -g $g -o $o -O $O -G $G -f out.bin > res.txt 
				IFS=$'\n' read -d '' -r -a lines2 < res.txt
				for j in "${lines2[@]}"; do			
					time=$(echo $j)
				done
				echo "time = $time"
				echo "mintime = $mintime"
				if (( $(echo "$mintime > $time" | bc -l) )); then
					mintime=$time
					optimal=$chunks
				fi
				#echo "p = $p s = $s Chunks = $chunks time = $time" >> binaryreduceall.txt
#			fi
		done
				echo "p = $p s = $s optimalChunks = $optimal time = $mintime" >> 2Treereduceoptimalchunks.txt		
	done
#done

