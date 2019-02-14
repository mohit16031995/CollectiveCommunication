sourceFiles = ['/home/blackcat/Desktop/TwoTreeAlgorithm/NewDraft-2019/LogAnalysis/sorted_allLogs_run6_Complete_8_5mb_9',
'/home/blackcat/Desktop/TwoTreeAlgorithm/NewDraft-2019/LogAnalysis/sorted_allLogs_run6_Complete_16_5mb_9',
'/home/blackcat/Desktop/TwoTreeAlgorithm/NewDraft-2019/LogAnalysis/sorted_allLogs_run11_Complete_8_5mb_9',
'/home/blackcat/Desktop/TwoTreeAlgorithm/NewDraft-2019/LogAnalysis/sorted_allLogs_run16_Complete_8_5mb_9',
'/home/blackcat/Desktop/TwoTreeAlgorithm/NewDraft-2019/LogAnalysis/sorted_allLogs_run6_SandersOld_8_5mb_9',
'/home/blackcat/Desktop/TwoTreeAlgorithm/NewDraft-2019/LogAnalysis/sorted_allLogs_run11_Complete_16_5mb_9',
'/home/blackcat/Desktop/TwoTreeAlgorithm/NewDraft-2019/LogAnalysis/sorted_allLogs_run16_Complete_16_5mb_9',
'/home/blackcat/Desktop/TwoTreeAlgorithm/NewDraft-2019/LogAnalysis/sorted_allLogs_run6_SandersOld_16_5mb_9',
'/home/blackcat/Desktop/TwoTreeAlgorithm/NewDraft-2019/LogAnalysis/sorted_allLogs_run11_SandersOld_8_5mb_9',
'/home/blackcat/Desktop/TwoTreeAlgorithm/NewDraft-2019/LogAnalysis/sorted_allLogs_run16_SandersOld_8_5mb_9',
'/home/blackcat/Desktop/TwoTreeAlgorithm/NewDraft-2019/LogAnalysis/sorted_allLogs_run11_SandersOld_16_5mb_9',
'/home/blackcat/Desktop/TwoTreeAlgorithm/NewDraft-2019/LogAnalysis/sorted_allLogs_run16_SandersOld_16_5mb_9',
'/home/blackcat/Desktop/TwoTreeAlgorithm/NewDraft-2019/LogAnalysis/sorted_allLogs_run6_SandersFlipped_8_5mb_9',
'/home/blackcat/Desktop/TwoTreeAlgorithm/NewDraft-2019/LogAnalysis/sorted_allLogs_run6_SandersFlipped_16_5mb_9',
'/home/blackcat/Desktop/TwoTreeAlgorithm/NewDraft-2019/LogAnalysis/sorted_allLogs_run11_SandersFlipped_8_5mb_9',
'/home/blackcat/Desktop/TwoTreeAlgorithm/NewDraft-2019/LogAnalysis/sorted_allLogs_run16_SandersFlipped_8_5mb_9',
'/home/blackcat/Desktop/TwoTreeAlgorithm/NewDraft-2019/LogAnalysis/sorted_allLogs_run11_SandersFlipped_16_5mb_9',
'/home/blackcat/Desktop/TwoTreeAlgorithm/NewDraft-2019/LogAnalysis/sorted_allLogs_run16_SandersFlipped_16_5mb_9',]
prefix = '/home/blackcat/Desktop/TwoTreeAlgorithm/NewDraft-2019/LogAnalysis/'
for sourceFile in sourceFiles:
	f = open(sourceFile)
	treeFile = "processed_"+sourceFile[len(prefix):]
	fout = open(treeFile,mode='w')
	lines = f.readlines()
	lines = [line[:line.find("in run")] for line in lines]
	tsprev = 0
	#shift=""
	for line in lines:
		line = line.split()
		#line[2] = str(float(line[2])*10000)
		ts = float(line[1])*10000
		line[1] = str(ts)
		line[1] = line[1][:line[1].find('.')+6]
		gap = len(line[1])-line[1].find('.')
		if gap!=6:
			line[1] = line[1]+'0'*(6-gap)
		interval = ts-tsprev
		if interval>=1:
			fout.write("\n\n\n")
		#	shift = shift+"\t"
		interval = str(interval)
		interval = interval[:interval.find('.')+6]
		gap = len(interval)-interval.find('.')
		if gap!=6:
			interval = interval+'0'*(6-gap)
		line[1] = line[1] + " " + interval
		tsprev = ts
		if line[0].startswith("LogISend"):
			to = line.index("to")
			line[to] = "to  "
			line[4] = "     ISEND       "
			line[5] = line[6] = ""
		if line[0].startswith("LogIRecv"):
			line[4] = "IRECV            "
			line[5] = line[6] = ""
		if line[0].startswith("LogReceived"):
			line[4] = "          RECVD  "
		line[0] = line[0]+'\t'#+shift
		line[2] = "   "+line[2]
		fout.write(' '.join(line)+'\n')
	fout.close()
	f.close()
