sourceFiles = ['/home/blackcat/Desktop/TwoTreeAlgorithm/NewDraft-2019/LogAnalysis/8_5mb_9/SandersFlipped/run16/allLogs_run16_SandersFlipped_8_5mb_9']
#prefix="/home/blackcat/Desktop/TwoTreeAlgorithm/NewDraft-2019/LogAnalysis/"
for sourceFile in sourceFiles:
	f = open(sourceFile)
	splitPoint = sourceFile.rfind('/')
	treeFile = sourceFile[:splitPoint+1]+"byrank_withSentTimes_processed_"+sourceFile[splitPoint+1:]
	fout = open(treeFile,mode='w')
	lines = f.readlines()
	lines = [line[:line.find("in run")].split() for line in lines]	
	lines = sorted(lines, key=lambda log: float(log[1]))
	lines = sorted(lines, key=lambda log: int(log[3]))
	tsprev = 0
	rankprev = lines[1][3]
	splitFileName = treeFile.split('_')
	chunks = int(splitFileName[-1])
	nP = int(splitFileName[-3])
	timings = [[0 for k in xrange(chunks)] for j in xrange(nP)]		#Map <source_rank,chunk_no> to time
	#print(timings)
	#shift=""
	linesToPrint = []
	for line in lines:
		#line = line.split()
		#line[2] = str(float(line[2])*10000)
		rank = int(line[3])
		if rank!=rankprev:
			tsprev = 0
			#fout.write("\n\n"+"-"*50+str(rank)+"-"*40+"\n\n")
			linesToPrint.append(["\n\n"+"-"*50+str(rank)+"-"*40+"\n"])
		ts = float(line[1])*10000
		line[1] = str(ts)
		line[1] = line[1][:line[1].find('.')+6]
		gap = len(line[1])-line[1].find('.')
		if gap!=6:
			line[1] = line[1]+'0'*(6-gap)
		interval = ts-tsprev
		if interval>=1:
			#fout.write("\n")
			linesToPrint.append([])
			#shift = shift+"\t"
		interval = str(interval)
		interval = interval[:interval.find('.')+6]
		gap = len(interval)-interval.find('.')
		if gap!=6:
			interval = interval+'0'*(6-gap)
		line[1] = line[1] + " " + interval
		tsprev = ts
		rankprev = rank
		if line[0].startswith("LogISend"):
			to = line.index("to")
			line[to] = "to  "
			line[4] = "     ISEND       "
			line[5] = line[6] = ""
			timings[rank][int(line[-3])] = ts
		if line[0].startswith("LogIRecv"):
			line[4] = "IRECV            "
			line[5] = line[6] = ""
		if line[0].startswith("LogReceived"):
			line[4] = "          RECVD  "
		line[0] = line[0]+'\t'#+shift
		line[2] = "   "+line[2]
		linesToPrint.append(line)
		#print(line)
	for line in linesToPrint:
		if line:
			if line[0].startswith("LogReceived"):
				recvTime = float(line[1][:line[1].find(' ')])
				source = int(line[-1])
				chunk = int(line[-4])
				sentTime = timings[source][chunk]
				suffix = "sent at "+str(sentTime)+", "+str(recvTime-sentTime)+" ago"
				line.append(suffix)
				print(line)
	for line in linesToPrint:
		fout.write(' '.join(line)+"\n")
	fout.close()
	f.close()
