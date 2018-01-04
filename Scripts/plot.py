#!/usr/bin/python

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt



# open the data


bin2tree_1048576_1_data = np.genfromtxt("STATS/stats_bin2tree/plot_data_bin2tree_1048576_1.txt", dtype=None, delimiter=' ', names=True)
bin2tree_1048576_2_data = np.genfromtxt("STATS/stats_bin2tree/plot_data_bin2tree_1048576_2.txt", dtype=None, delimiter=' ', names=True)
bin2tree_1048576_3_data = np.genfromtxt("STATS/stats_bin2tree/plot_data_bin2tree_1048576_3.txt", dtype=None, delimiter=' ', names=True)

TreeBdcast_1048576_1_data = np.genfromtxt("STATS/stats_2TreeBdcast/plot_data_2TreeBdcast_1048576_1.txt", dtype=None, delimiter=' ', names=True)
TreeBdcast_1048576_2_data = np.genfromtxt("STATS/stats_2TreeBdcast/plot_data_2TreeBdcast_1048576_2.txt", dtype=None, delimiter=' ', names=True)
TreeBdcast_1048576_3_data = np.genfromtxt("STATS/stats_2TreeBdcast/plot_data_2TreeBdcast_1048576_3.txt", dtype=None, delimiter=' ', names=True)


bintree_1048576_1_data = np.genfromtxt("STATS/stats_bintree/plot_data_bintree_1048576_1.txt", dtype=None, delimiter=' ', names=True)
bintree_1048576_2_data = np.genfromtxt("STATS/stats_bintree/plot_data_bintree_1048576_2.txt", dtype=None, delimiter=' ', names=True)
bintree_1048576_3_data = np.genfromtxt("STATS/stats_bintree/plot_data_bintree_1048576_3.txt", dtype=None, delimiter=' ', names=True)

MPI_BCAST_1048576_data = np.genfromtxt("STATS/stats_MPI_BCAST/plot_data_MPI_BCAST_1048576.txt", dtype=None, delimiter=' ', names=True)


#optimal segments
#pip_data = np.genfromtxt("pip/plot_data_1048576-2.txt", dtype=None, delimiter=' ', names=True)

#seg2_data = np.genfromtxt("seg21/plot_data_1048576.txt", dtype=None, delimiter=' ', names=True)
#simple_data = np.genfromtxt("simple/plot_data_1048576.txt", dtype=None, delimiter=' ', names=True)
#btree_data = np.genfromtxt("btree/plot_data_1048576.txt", dtype=None, delimiter=' ', names=True)
#bin2tree_data = np.genfromtxt("plot_data_1048576.txt", dtype=None, delimiter=' ', names=True)
#bintree_data = np.genfromtxt("bintree/plot_data_1048576-3.txt", dtype=None, delimiter=' ', names=True)


scale=1e3

# convert from seconds to milliseconds
TreeBdcast_1048576_1_data["median"] = TreeBdcast_1048576_1_data["median"] * scale 
TreeBdcast_1048576_1_data["CINNormhigh"] = TreeBdcast_1048576_1_data["CINNormhigh"] * scale
TreeBdcast_1048576_1_data["CINNormlow"] = TreeBdcast_1048576_1_data["CINNormlow"] * scale

TreeBdcast_1048576_2_data["median"] = TreeBdcast_1048576_2_data["median"] * scale 
TreeBdcast_1048576_2_data["CINNormhigh"] = TreeBdcast_1048576_2_data["CINNormhigh"] * scale
TreeBdcast_1048576_2_data["CINNormlow"] = TreeBdcast_1048576_2_data["CINNormlow"] * scale

TreeBdcast_1048576_3_data["median"] = TreeBdcast_1048576_3_data["median"] * scale 
TreeBdcast_1048576_3_data["CINNormhigh"] = TreeBdcast_1048576_3_data["CINNormhigh"] * scale
TreeBdcast_1048576_3_data["CINNormlow"] = TreeBdcast_1048576_3_data["CINNormlow"] * scale


bin2tree_1048576_1_data["median"] = bin2tree_1048576_1_data["median"] * scale 
bin2tree_1048576_1_data["CINNormhigh"] = bin2tree_1048576_1_data["CINNormhigh"] * scale
bin2tree_1048576_1_data["CINNormlow"] = bin2tree_1048576_1_data["CINNormlow"] * scale

bin2tree_1048576_2_data["median"] = bin2tree_1048576_2_data["median"] * scale 
bin2tree_1048576_2_data["CINNormhigh"] = bin2tree_1048576_2_data["CINNormhigh"] * scale
bin2tree_1048576_2_data["CINNormlow"] = bin2tree_1048576_2_data["CINNormlow"] * scale

bin2tree_1048576_3_data["median"] = bin2tree_1048576_3_data["median"] * scale 
bin2tree_1048576_3_data["CINNormhigh"] = bin2tree_1048576_3_data["CINNormhigh"] * scale
bin2tree_1048576_3_data["CINNormlow"] = bin2tree_1048576_3_data["CINNormlow"] * scale

bintree_1048576_1_data["median"] = bintree_1048576_1_data["median"] * scale 
bintree_1048576_1_data["CINNormhigh"] = bintree_1048576_1_data["CINNormhigh"] * scale
bintree_1048576_1_data["CINNormlow"] = bintree_1048576_1_data["CINNormlow"] * scale

bintree_1048576_2_data["median"] = bintree_1048576_2_data["median"] * scale 
bintree_1048576_2_data["CINNormhigh"] = bintree_1048576_2_data["CINNormhigh"] * scale
bintree_1048576_2_data["CINNormlow"] = bintree_1048576_2_data["CINNormlow"] * scale

bintree_1048576_3_data["median"] = bintree_1048576_3_data["median"] * scale 
bintree_1048576_3_data["CINNormhigh"] = bintree_1048576_3_data["CINNormhigh"] * scale
bintree_1048576_3_data["CINNormlow"] = bintree_1048576_3_data["CINNormlow"] * scale


MPI_BCAST_1048576_data["median"] = MPI_BCAST_1048576_data["median"] * scale 



#pip_data["median"] = pip_data["median"] * scale
#pip_data["CINNormhigh"] = pip_data["CINNormhigh"] * scale
#pip_data["CINNormlow"] = pip_data["CINNormlow"] * scale
#print pip_data["median"]

#seg2_data["median"] = seg2_data["median"] * scale
#seg2_data["CINNormhigh"] = seg2_data["CINNormhigh"] * scale
#seg2_data["CINNormlow"] = seg2_data["CINNormlow"] * scale
#print seg2_data["median"]

#simple_data["median"] = simple_data["median"] * scale
#simple_data["CINNormhigh"] = simple_data["CINNormhigh"] * scale
#simple_data["CINNormlow"] = simple_data["CINNormlow"] * scale
#print simple_data["median"]

#btree_data["median"] = btree_data["median"] * scale
#btree_data["CINNormhigh"] = btree_data["CINNormhigh"] * scale
#btree_data["CINNormlow"] = btree_data["CINNormlow"] * scale
#print btree_data["median"]

#bin2tree_data["median"] = bin2tree_data["median"] * scale
#bin2tree_data["CINNormhigh"] = bin2tree_data["CINNormhigh"] * scale
#bin2tree_data["CINNormlow"] = bin2tree_data["CINNormlow"] * scale
#print bin2tree_data["median"]

#bintree_data["median"] = bintree_data["median"] * scale
#bintree_data["CINNormhigh"] = bintree_data["CINNormhigh"] * scale
#bintree_data["CINNormlow"] = bintree_data["CINNormlow"] * scale
#print bintree_data["median"]


min_x = 8
max_x = 540


#bcast-only timings
TreeBdcast_1048576_1_CI = np.vstack((TreeBdcast_1048576_1_data["median"] - TreeBdcast_1048576_1_data["CINNormhigh"], TreeBdcast_1048576_1_data["CINNormlow"] - TreeBdcast_1048576_1_data["median"]))
plt.errorbar(x=TreeBdcast_1048576_1_data["nP"], y=TreeBdcast_1048576_1_data["median"], yerr=TreeBdcast_1048576_1_CI, fmt="r-", label="2TreeBdcast_1")

TreeBdcast_1048576_2_CI = np.vstack((TreeBdcast_1048576_2_data["median"] - TreeBdcast_1048576_2_data["CINNormhigh"], TreeBdcast_1048576_2_data["CINNormlow"] - TreeBdcast_1048576_2_data["median"]))
plt.errorbar(x=TreeBdcast_1048576_2_data["nP"], y=TreeBdcast_1048576_2_data["median"], yerr=TreeBdcast_1048576_2_CI, fmt="r--", label="2TreeBdcast_2")

TreeBdcast_1048576_3_CI = np.vstack((TreeBdcast_1048576_3_data["median"] - TreeBdcast_1048576_3_data["CINNormhigh"], TreeBdcast_1048576_3_data["CINNormlow"] - TreeBdcast_1048576_3_data["median"]))
plt.errorbar(x=TreeBdcast_1048576_3_data["nP"], y=TreeBdcast_1048576_3_data["median"], yerr=TreeBdcast_1048576_3_CI, fmt="r:", label="2TreeBdcast_3")

bin2tree_1048576_1_CI = np.vstack((bin2tree_1048576_1_data["median"] - bin2tree_1048576_1_data["CINNormhigh"], bin2tree_1048576_1_data["CINNormlow"] - bin2tree_1048576_1_data["median"]))
plt.errorbar(x=bin2tree_1048576_1_data["nP"], y=bin2tree_1048576_1_data["median"], yerr=bin2tree_1048576_1_CI, fmt="g-", label="bin2tree_1")

bin2tree_1048576_2_CI = np.vstack((bin2tree_1048576_2_data["median"] - bin2tree_1048576_2_data["CINNormhigh"], bin2tree_1048576_2_data["CINNormlow"] - bin2tree_1048576_2_data["median"]))
plt.errorbar(x=bin2tree_1048576_2_data["nP"], y=bin2tree_1048576_2_data["median"], yerr=bin2tree_1048576_2_CI, fmt="g--", label="bin2tree_2")

bin2tree_1048576_3_CI = np.vstack((bin2tree_1048576_3_data["median"] - bin2tree_1048576_3_data["CINNormhigh"], bin2tree_1048576_3_data["CINNormlow"] - bin2tree_1048576_3_data["median"]))
plt.errorbar(x=bin2tree_1048576_3_data["nP"], y=bin2tree_1048576_3_data["median"], yerr=bin2tree_1048576_3_CI, fmt="g:", label="bin2tree_3")

bintree_1048576_1_CI = np.vstack((bintree_1048576_1_data["median"] - bintree_1048576_1_data["CINNormhigh"], bintree_1048576_1_data["CINNormlow"] - bintree_1048576_1_data["median"]))
plt.errorbar(x=bintree_1048576_1_data["nP"], y=bintree_1048576_1_data["median"], yerr=bintree_1048576_1_CI, fmt="m-", label="bintree_1")

bintree_1048576_2_CI = np.vstack((bintree_1048576_2_data["median"] - bintree_1048576_2_data["CINNormhigh"], bintree_1048576_2_data["CINNormlow"] - bintree_1048576_2_data["median"]))
plt.errorbar(x=bintree_1048576_2_data["nP"], y=bintree_1048576_2_data["median"], yerr=bintree_1048576_2_CI, fmt="m--", label="bintree_2")

bintree_1048576_3_CI = np.vstack((bintree_1048576_3_data["median"] - bintree_1048576_3_data["CINNormhigh"], bintree_1048576_3_data["CINNormlow"] - bintree_1048576_3_data["median"]))
plt.errorbar(x=bintree_1048576_3_data["nP"], y=bintree_1048576_3_data["median"], yerr=bintree_1048576_3_CI, fmt="m:", label="bintree_3")

MPI_BCAST_1048576_CI = np.vstack((MPI_BCAST_1048576_data["median"] - MPI_BCAST_1048576_data["CINNormhigh"], MPI_BCAST_1048576_data["CINNormlow"] - MPI_BCAST_1048576_data["median"]))
plt.errorbar(x=MPI_BCAST_1048576_data["nP"], y=MPI_BCAST_1048576_data["median"], yerr=MPI_BCAST_1048576_CI, fmt="b-", label="MPI_BCAST")



#pip-only timings
#pip_CI = np.vstack((pip_data["median"] - pip_data["CINNormhigh"],
#                     pip_data["CINNormlow"] - pip_data["median"]))
#plt.errorbar(x=pip_data["nP"], y=pip_data["median"], yerr=pip_CI, fmt="go-", label="PIP")

#seg2-only timings
#seg2_CI = np.vstack((seg2_data["median"] - seg2_data["CINNormhigh"],
 #                    seg2_data["CINNormlow"] - seg2_data["median"]))
#plt.errorbar(x=seg2_data["nP"], y=seg2_data["median"], yerr=seg2_CI, fmt="m+--", label="Seg2")

#Simple-only timings
#simple_CI = np.vstack((simple_data["median"] - simple_data["CINNormhigh"],
#                     simple_data["CINNormlow"] - simple_data["median"]))
#plt.errorbar(x=simple_data["nP"], y=simple_data["median"], yerr=simple_CI, fmt="b+-", label="Simple")

#bin2tree-only timings
#bin2tree_CI = np.vstack((bin2tree_data["median"] - bin2tree_data["CINNormhigh"],
 #                    bin2tree_data["CINNormlow"] - bin2tree_data["median"]))
#plt.errorbar(x=bin2tree_data["nP"], y=bin2tree_data["median"], yerr=bin2tree_CI, fmt="b+-", label="bin2tree-ran")


#Binomial tree-only timings
#btree_CI = np.vstack((btree_data["median"] - btree_data["CINNormhigh"],
#                     btree_data["CINNormlow"] - btree_data["median"]))
#plt.errorbar(x=btree_data["nP"], y=btree_data["median"], yerr=btree_CI, fmt="y+-", label="BMLTree")

#bintree_CI = np.vstack((bintree_data["median"] - bintree_data["CINNormhigh"],
   #                  bintree_data["CINNormlow"] - bintree_data["median"]))
#plt.errorbar(x=bintree_data["nP"], y=bintree_data["median"], yerr=bintree_CI, fmt="y+-", label="bintree")


# labels, etc.
plt.legend(loc="upper left")
plt.xlabel("# of processors")
plt.ylabel("execution time [ms]")
#plt.xlim(min_x, max_x)
plt.ylim(ymin=0)
x_values = [8, 16, 32, 64, 128, 150, 256, 512]


plt.xscale('log', basex=2)

figtitle = 's=1Mbyte'
plt.suptitle(figtitle, fontsize=16)

# write
plt.savefig("newgraph.png")
