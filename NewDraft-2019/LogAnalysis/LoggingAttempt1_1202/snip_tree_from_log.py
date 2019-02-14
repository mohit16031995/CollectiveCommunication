sourceFile = "reduce_sandersOld_8procs5mb9chunks.out"
f = open(sourceFile)
treeFile = "tree_"+sourceFile
fout = open(treeFile,mode='w')
for line in f:
	if line.startswith("rank"):
		fout.write(line)
fout.close()
