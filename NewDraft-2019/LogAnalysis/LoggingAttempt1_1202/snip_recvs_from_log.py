sourceFiles = ["run6_recvsOnly_reduce_sandersOld_8procs5mb9chunks_new.out","run11_recvsOnly_reduce_sandersOld_8procs5mb9chunks_new.out","run16_recvsOnly_reduce_sandersOld_8procs5mb9chunks_new.out"]
for sourceFile in sourceFiles:
	f = open(sourceFile)
	treeFile = "recvsSorted_"+sourceFile
	fout = open(treeFile,mode='w')
	lines = f.readlines()
	lines = [line[:line.find("run")] for line in lines if line.startswith("Log")]
	lines = sorted(lines, key=lambda log: float(log.split()[2]))
	tsprev = 0
	for line in lines:
		line = line.split()
		#line[2] = str(float(line[2])*10000)
		ts = float(line[2])*10000
		line[2] = str(ts)
		line[2] = line[2][:line[2].find('.')+6]
		interval = str(ts-tsprev)
		interval = interval[:interval.find('.')+6]
		line[2] = line[2] + " " + interval
		tsprev = ts
		#line[2] = line[2][4:]                      #WARNING#####################
		#line[2] = line[2][:2] + '.' + line[2][2:]
		fout.write(' '.join(line)+'\n')
	fout.close()
	f.close()
