sourceFiles = ["reduce_sandersOld_8procs5mb9chunks_new.out"]
for sourceFile in sourceFiles:
	f = open(sourceFile)
	treeFile = "recvsOnly_"+sourceFile
	treeFiles = ["run6_"+treeFile,"run11_"+treeFile,"run16_"+treeFile]
	fouts = [open(treeFile,mode='w') for treeFile in treeFiles]
	lines = f.readlines()
	lines = [line for line in lines if line.startswith("Log")]
	#lines = sorted(lines, key=lambda log: float(log.split()[2]))
	#tsprev = 0
	for line in lines:
		line = line.split()
		run = int(line[11])
		if run==6:
			fout = fouts[0]
		if run==11:
			fout = fouts[1]
		if run==16:
			fout = fouts[2]
		#line[2] = str(float(line[2])*10000)
		#ts = float(line[2])*10000
		#line[2] = str(ts)
		#line[2] = line[2][:line[2].find('.')+6]
		#interval = str(ts-tsprev)
		#interval = interval[:interval.find('.')+6]
		#line[2] = line[2] + " " + interval
		#tsprev = ts
		#line[2] = line[2][4:]                      #WARNING#####################
		#line[2] = line[2][:2] + '.' + line[2][2:]
		fout.write(' '.join(line)+"\n")
	for fout in fouts:
		fout.close()
	f.close()
