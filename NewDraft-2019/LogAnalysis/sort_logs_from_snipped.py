sourceFiles = ["allLogs_run6_SandersFlipped_16_5mb_9",
"allLogs_run11_SandersFlipped_16_5mb_9",
"allLogs_run16_SandersFlipped_16_5mb_9",
"allLogs_run6_SandersOld_16_5mb_9",
"allLogs_run11_SandersOld_16_5mb_9",
"allLogs_run16_SandersOld_16_5mb_9",
"allLogs_run6_Complete_16_5mb_9",
"allLogs_run11_Complete_16_5mb_9",
"allLogs_run16_Complete_16_5mb_9"]
for sourceFile in sourceFiles:
	f = open(sourceFile)
	out = "sorted_"+sourceFile
	fout = open(out,mode='w')
	lines = f.readlines()
	lines = sorted(lines, key=lambda log: float(log.split()[1]))
	for line in lines:
		fout.write(line)
	fout.close()
	f.close()
