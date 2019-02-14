sourceFiles = ["reduce_complete_8procs5mb9chunks.out"]
for sourceFile in sourceFiles:
	f = open(sourceFile)
	treeFile = "recvsSortedNormalized_"+sourceFile
	fout = open(treeFile,mode='w')
	lines = f.readlines()
	lines = [line for line in lines if line.startswith("Log")]
	lines = sorted(lines, key=lambda log: float(log.split()[2]))
	tsprev = 1
	tsprevNorm = 0 
	diff = 0
	for line in lines:
		line = line.split()
		ts = float(line[2])*10000
		if (ts-tsprev >= 2):
			diff = diff + (ts-tsprev)/1 - 2
			print(diff)
		tsNorm = ts - diff
		line[2] = str(tsNorm)
		line[2] = line[2][:line[2].find('.')+6]
		interval = str(tsNorm-tsprevNorm)
		interval = interval[:interval.find('.')+6]
		line[2] = line[2] + " " + interval
		tsprev = ts
		tsprevNorm = tsNorm
		fout.write(' '.join(line)+'\n')
	fout.close()
	f.close()