import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

data_in = ["LogEndTime:-run16_Complete_16_5mb_9",
"LogEndTime:-run11_SandersOld_16_5mb_9"]
p=16
f,ax = plt.subplots()

data1 = open(data_in[1])
lines = data1.readlines()
x2 = []
y2 = []
for line in lines:
	splitLine = line.split()
	x2.append(splitLine[splitLine.index("rank")+1])
	y2.append(float(splitLine[1]))
X2 = np.arange(p)
plt.bar(X2, y2, width = 0.25, label = 'TwoTreeSanders')

data1 = open(data_in[0])
lines = data1.readlines()
x1 = []
y1 = []
for line in lines:
	splitLine = line.split()
	x1.append(splitLine[splitLine.index("rank")+1])
	y1.append(float(splitLine[1]))
X1 = np.arange(p)
plt.bar(X1, y1, width = 0.25, bottom = 0.0, label = 'TwoTreeComplete')


index=X1
ticks=[]
for i in range(p):
	c=x1[i]
	d=x2[i]
	ticks.append(d+","+c)

plt.xticks(index, ticks)

plt.ylabel("Time of Completion")
plt.xlabel("Process Rank")
plt.title('5MB')
plt.legend()
plt.show()
