#!/bin/bash -l
#SBATCH --ntasks=256
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:50:00
#SBATCH --constraint=mc
#SBATCH --switches=1

#export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

IFS=$'\n' read -d '' -r -a lines < cut2.txt
len=${#lines[@]}
#echo $len

for i in "${lines[@]}"; do
#	line="${line[$i]}"
#	echo $i	
	np=$(echo $i | cut -d " " -f 1)
	msize=$(echo $i | cut -d " " -f 2)
	c1=$(echo $i | cut -d " " -f 3)
	c2=$(echo $i | cut -d " " -f 4)
	c3=$(echo $i | cut -d " " -f 5)
#	c4=$(($c1/2))	
	c4=$(((c1*7)/4))
	c5=$(((c2*7)/4))
	c6=$(((c3*7)/4))

	
			
	echo "MPI_reduce $np $msize"	
	srun -N $np -n $np -c $SLURM_CPUS_PER_TASK --cpu_bind=rank --ntasks-per-node $SLURM_NTASKS_PER_NODE ./ex_mpi_reduce $msize
	echo "binomial_reduce $np $msize"	
	srun -N $np -n $np -c $SLURM_CPUS_PER_TASK --cpu_bind=rank --ntasks-per-node $SLURM_NTASKS_PER_NODE ./ex_binomial_reduce $msize
	echo "reduceScatter_reduce $np $msize"	
	srun -N $np -n $np -c $SLURM_CPUS_PER_TASK --cpu_bind=rank --ntasks-per-node $SLURM_NTASKS_PER_NODE ./ex_reduceScatter_reduce $msize
	echo "2TreeComplete_reduce $np $msize $c1"	
	srun -N $np -n $np -c $SLURM_CPUS_PER_TASK --cpu_bind=rank --ntasks-per-node $SLURM_NTASKS_PER_NODE ./ex_2treecomplete_reduce $msize $c1
	echo "2TreeComplete_reduce2 $np $msize $c4"	
	srun -N $np -n $np -c $SLURM_CPUS_PER_TASK --cpu_bind=rank --ntasks-per-node $SLURM_NTASKS_PER_NODE ./ex_2treecomplete_reduce $msize $c4
	echo "2TreeSanders_reduce $np $msize $c2"	
	srun -N $np -n $np -c $SLURM_CPUS_PER_TASK --cpu_bind=rank --ntasks-per-node $SLURM_NTASKS_PER_NODE ./ex_2treesanders_reduce $msize $c2
	echo "2TreeSanders_reduce2 $np $msize $c5"	
	srun -N $np -n $np -c $SLURM_CPUS_PER_TASK --cpu_bind=rank --ntasks-per-node $SLURM_NTASKS_PER_NODE ./ex_2treesanders_reduce $msize $c5
	echo "binary_reduce $np $msize $c3"	
	srun -N $np -n $np -c $SLURM_CPUS_PER_TASK --cpu_bind=rank --ntasks-per-node $SLURM_NTASKS_PER_NODE ./ex_binary_reduce $msize $c3
	echo "binary_reduce2 $np $msize $c6"	
	srun -N $np -n $np -c $SLURM_CPUS_PER_TASK --cpu_bind=rank --ntasks-per-node $SLURM_NTASKS_PER_NODE ./ex_binary_reduce $msize $c6

#	echo "2Treereduce_false2 $np $msize $c4"	
#	srun -N $np -n $np -c $SLURM_CPUS_PER_TASK --cpu_bind=rank --ntasks-per-node $SLURM_NTASKS_PER_NODE ./2treereduce_false $msize $c4
	#echo "2treecomplete_reduce $np $msize $c1"	
	#srun -N $np -n $np -c $SLURM_CPUS_PER_TASK --cpu_bind=rank --ntasks-per-node $SLURM_NTASKS_PER_NODE ./ex_2treecomplete_reduce $msize $c1
#	echo "2treecomplete_reduce2 $np $msize $c4"	
#	srun -N $np -n $np -c $SLURM_CPUS_PER_TASK --cpu_bind=rank --ntasks-per-node $SLURM_NTASKS_PER_NODE ./2treereduce $msize $c4
	#echo "binary_reduce $np $msize $c2"	
	#srun -N $np -n $np -c $SLURM_CPUS_PER_TASK --cpu_bind=rank --ntasks-per-node $SLURM_NTASKS_PER_NODE ./ex_binary_reduce $msize $c2
#	echo "linear_reduce $np $msize $c3"	
#	srun -N $np -n $np -c $SLURM_CPUS_PER_TASK --cpu_bind=rank --ntasks-per-node $SLURM_NTASKS_PER_NODE ./linearreduce $msize $c3
	#echo "MPI_reduce $np $msize"	
	#srun -N $np -n $np -c $SLURM_CPUS_PER_TASK --cpu_bind=rank --ntasks-per-node $SLURM_NTASKS_PER_NODE ./ex_mpi_reduce $msize
	#echo "MPIuser_reduce $np $msize"	
	#srun -N $np -n $np -c $SLURM_CPUS_PER_TASK --cpu_bind=rank --ntasks-per-node $SLURM_NTASKS_PER_NODE ./ex_mpiuser_reduce $msize
	#echo "2treecomplete_allreduce $np $msize $c1"	
	#srun -N $np -n $np -c $SLURM_CPUS_PER_TASK --cpu_bind=rank --ntasks-per-node $SLURM_NTASKS_PER_NODE ./ex_2treecomplete_allreduce $msize $c1
	#echo "2treecomplete_allreduceoptimal $np $msize $c1"	
	#srun -N $np -n $np -c $SLURM_CPUS_PER_TASK --cpu_bind=rank --ntasks-per-node $SLURM_NTASKS_PER_NODE ./ex_2treecomplete_allreduceoptimal $msize $c1
	#echo "2treecomplete_allreducebinary $np $msize $c1"	
	#srun -N $np -n $np -c $SLURM_CPUS_PER_TASK --cpu_bind=rank --ntasks-per-node $SLURM_NTASKS_PER_NODE ./ex_2treecomplete_allreducebinary $msize $c1
#	echo "binary_allreduce $np $msize $c2"	
#	srun -N $np -n $np -c $SLURM_CPUS_PER_TASK --cpu_bind=rank --ntasks-per-node $SLURM_NTASKS_PER_NODE ./ex_binary_allreduce $msize $c2
	#echo "linear_allreduce $np $msize $c3"	
	#srun -N $np -n $np -c $SLURM_CPUS_PER_TASK --cpu_bind=rank --ntasks-per-node $SLURM_NTASKS_PER_NODE ./linearallreduce $msize $c3
	#echo "MPI_allreduce $np $msize"	
	#srun -N $np -n $np -c $SLURM_CPUS_PER_TASK --cpu_bind=rank --ntasks-per-node $SLURM_NTASKS_PER_NODE ./ex_mpi_allreduce $msize
	#echo "MPIuser_allreduce $np $msize"	
	#srun -N $np -n $np -c $SLURM_CPUS_PER_TASK --cpu_bind=rank --ntasks-per-node $SLURM_NTASKS_PER_NODE ./ex_mpiuser_allreduce $msize
#	cc3=$(($c3/10))
#	echo "cc3 = $cc3"
#	for j in -1 0 1 2 3; do
#		echo "j = $j"
#		pp=$(($cc3*$j))
#		echo "linearPipeline $np $msize $(($c3-$pp))"
#		srun -N $np -n $np -c $SLURM_CPUS_PER_TASK --cpu_bind=rank --ntasks-per-node $SLURM_NTASKS_PER_NODE ./linearPipeline $msize $(($c3-$pp))
#	done
#	for add in 1 0 -1 -2 -3; do
		#cc1=$(($c1/10))
		#for l in -1 1 2 3 4; do
		#	pp=$(($cc1*$l))
#		echo "bintree $np $msize $(($c1+$add))"
#		srun -N $np -n $np -c $SLURM_CPUS_PER_TASK --cpu_bind=rank --ntasks-per-node $SLURM_NTASKS_PER_NODE ./bintree $msize $(($c1+$add))
		#done
#		echo "bin2tree $np $msize $(($c2+$add))"
#		srun -N $np -n $np -c $SLURM_CPUS_PER_TASK --cpu_bind=rank --ntasks-per-node $SLURM_NTASKS_PER_NODE ./bin2tree $msize $(($c2+$add))
#		echo "linearPipeline $np $msize $(($c2+$add))"
#		srun -N $np -n $np -c $SLURM_CPUS_PER_TASK --cpu_bind=rank --ntasks-per-node $SLURM_NTASKS_PER_NODE ./2TreeBdcast $np $msize $(($c2+$add))
#	cc2=$(($c2/10))
#	for j in -1 1 2 3 4; do
#		pp=$(($cc2*$j))
#		echo "2TreeComplete $np $msize $(($c2+($add*2)))"
#		srun -N $np -n $np -c $SLURM_CPUS_PER_TASK --cpu_bind=rank --ntasks-per-node $SLURM_NTASKS_PER_NODE ./2TreeComplete $msize $(($c2+($add*2)))
#	done
#	done
#	cc4=$(($c4/10))
#	for k in -1 1 2 3 4; do
#		pp=$(($cc4*$k))
#		echo "2TreeBdcast $np $msize $(($c4+($add*2)))"
#		srun -N $np -n $np -c $SLURM_CPUS_PER_TASK --cpu_bind=rank --ntasks-per-node $SLURM_NTASKS_PER_NODE ./2TreeBdcast $np $msize $(($c4+($add*2)))
#	done
done
