sourceFiles = ["/home/blackcat/Desktop/TwoTreeAlgorithm/NewDraft-2019/LogAnalysis/sorted_allLogs_run6_Complete_8_5mb_9",
"/home/blackcat/Desktop/TwoTreeAlgorithm/NewDraft-2019/LogAnalysis/sorted_allLogs_run11_Complete_8_5mb_9",
"/home/blackcat/Desktop/TwoTreeAlgorithm/NewDraft-2019/LogAnalysis/sorted_allLogs_run16_Complete_8_5mb_9",
"/home/blackcat/Desktop/TwoTreeAlgorithm/NewDraft-2019/LogAnalysis/sorted_allLogs_run6_SandersOld_8_5mb_9",
"/home/blackcat/Desktop/TwoTreeAlgorithm/NewDraft-2019/LogAnalysis/sorted_allLogs_run11_SandersOld_8_5mb_9",
"/home/blackcat/Desktop/TwoTreeAlgorithm/NewDraft-2019/LogAnalysis/sorted_allLogs_run16_SandersOld_8_5mb_9",
"/home/blackcat/Desktop/TwoTreeAlgorithm/NewDraft-2019/LogAnalysis/sorted_allLogs_run6_SandersFlipped_8_5mb_9",
"/home/blackcat/Desktop/TwoTreeAlgorithm/NewDraft-2019/LogAnalysis/sorted_allLogs_run11_SandersFlipped_8_5mb_9",
"/home/blackcat/Desktop/TwoTreeAlgorithm/NewDraft-2019/LogAnalysis/sorted_allLogs_run16_SandersFlipped_8_5mb_9",
"/home/blackcat/Desktop/TwoTreeAlgorithm/NewDraft-2019/LogAnalysis/sorted_allLogs_run6_Complete_16_5mb_9",
"/home/blackcat/Desktop/TwoTreeAlgorithm/NewDraft-2019/LogAnalysis/sorted_allLogs_run11_Complete_16_5mb_9",
"/home/blackcat/Desktop/TwoTreeAlgorithm/NewDraft-2019/LogAnalysis/sorted_allLogs_run16_Complete_16_5mb_9",
"/home/blackcat/Desktop/TwoTreeAlgorithm/NewDraft-2019/LogAnalysis/sorted_allLogs_run6_SandersOld_16_5mb_9",
"/home/blackcat/Desktop/TwoTreeAlgorithm/NewDraft-2019/LogAnalysis/sorted_allLogs_run11_SandersOld_16_5mb_9",
"/home/blackcat/Desktop/TwoTreeAlgorithm/NewDraft-2019/LogAnalysis/sorted_allLogs_run16_SandersOld_16_5mb_9",
"/home/blackcat/Desktop/TwoTreeAlgorithm/NewDraft-2019/LogAnalysis/sorted_allLogs_run6_SandersFlipped_16_5mb_9",
"/home/blackcat/Desktop/TwoTreeAlgorithm/NewDraft-2019/LogAnalysis/sorted_allLogs_run11_SandersFlipped_16_5mb_9",
"/home/blackcat/Desktop/TwoTreeAlgorithm/NewDraft-2019/LogAnalysis/sorted_allLogs_run16_SandersFlipped_16_5mb_9"]
prefix = "/home/blackcat/Desktop/TwoTreeAlgorithm/NewDraft-2019/LogAnalysis/sorted_allLogs_"
filters = ["LogEndTime:",
"LogReceived:"]
for sourceFile in sourceFiles:
	f = open(sourceFile)
	outSuffix = sourceFile[len(prefix):]
	lines = f.readlines()
	for filt in filters:
		filtered_lines = [line for line in lines if line.startswith(filt)]
		outFile = filt+"-"+outSuffix
		fout = open(outFile,mode='w')
		for line in filtered_lines:
			fout.write(line)
		fout.close()
	f.close()
