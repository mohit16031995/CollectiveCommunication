sourceFiles = ["/home/blackcat/Desktop/TwoTreeAlgorithm/NewDraft-2019/LogAnalysis/allruns_SandersOld_16_5mb_9",
"/home/blackcat/Desktop/TwoTreeAlgorithm/NewDraft-2019/LogAnalysis/allruns_Complete_16_5mb_9",
"/home/blackcat/Desktop/TwoTreeAlgorithm/NewDraft-2019/LogAnalysis/allruns_SandersFlipped_16_5mb_9"]
prefix = "/home/blackcat/Desktop/TwoTreeAlgorithm/NewDraft-2019/LogAnalysis/allruns_"
for sourceFile in sourceFiles:
	f = open(sourceFile)
	out = sourceFile[len(prefix):]
	outFiles = ["allLogs_run6_"+out,"allLogs_run11_"+out,"allLogs_run16_"+out]
	fouts = [open(outFile,mode='w') for outFile in outFiles]
	lines = f.readlines()
	lines = [line for line in lines if line.startswith("Log")]
	#lines = sorted(lines, key=lambda log: float(log.split()[2]))
	for line in lines:
		#line = line.split()
		#print(line)
		lineAsList = line.split()
		#print(lineAsList)
		run = lineAsList.index("run")
		run = int(lineAsList[run+1])
		#print(run)
		if run==6:
			fout = fouts[0]
		if run==11:
			fout = fouts[1]
		if run==16:
			fout = fouts[2]
		fout.write(line)
	for fout in fouts:
		fout.close()
	f.close()
