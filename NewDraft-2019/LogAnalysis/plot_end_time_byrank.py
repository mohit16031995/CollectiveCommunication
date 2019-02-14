import matplotlib.pyplot as plt
import numpy as np

data_in = ["LogEndTime:-run16_Complete_16_5mb_9",
"LogEndTime:-run11_SandersOld_16_5mb_9"]
p=16

data1 = open(data_in[1])
lines = data1.readlines()
x1 = []
y1 = []
for line in lines:
	splitLine = line.split()
	x1.append(int(splitLine[splitLine.index("rank")+1]))
	y1.append(float(splitLine[1]))
X1 = np.arange(p)
x1 = [a+0.25 for a in x1]
plt.bar(x1, y1, width = 0.25, label = 'TwoTreeSanders')
plt.xticks(X1+0.25,X1)

data1 = open(data_in[0])
lines = data1.readlines()
x1 = []
y1 = []
for line in lines:
	splitLine = line.split()
	x1.append(int(splitLine[splitLine.index("rank")+1]))
	y1.append(float(splitLine[1]))
X1 = np.arange(p)
plt.bar(x1, y1, width = 0.25, label = 'TwoTreeComplete')
plt.xticks(X1+0.125,X1)

plt.ylabel("Time of Completion")
plt.xlabel("Process Rank")
plt.title('5MB')
plt.legend()
plt.show()
