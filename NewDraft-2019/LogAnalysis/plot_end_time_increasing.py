import matplotlib.pyplot as plt
import numpy as np

data_in = ["LogEndTime:-run16_Complete_16_5mb_9",
"LogEndTime:-run11_SandersOld_16_5mb_9"]
p=16
width = 0.3

data1 = open(data_in[1])
lines = data1.readlines()
x2 = []
y2 = []
for line in lines:
	splitLine = line.split()
	x2.append(splitLine[splitLine.index("rank")+1])
	y2.append(float(splitLine[1]))
X2 = np.arange(p) 
plt.bar(X2, y2, width = width, label = 'TwoTreeSanders')

data1 = open(data_in[0])
lines = data1.readlines()
x1 = []
y1 = []
for line in lines:
	splitLine = line.split()
	x1.append(splitLine[splitLine.index("rank")+1])
	y1.append(float(splitLine[1]))
X1 = np.arange(p) + width
plt.bar(X1, y1, width = width, label = 'TwoTreeComplete')


index=[]
ticks=[]
for i in range(p):
	a=X1[i]
	b=X2[i]
	c=x1[i]
	d=x2[i]
	index.append(a)
	index.append(b)
	ticks.append(c)
	ticks.append(d)

plt.xticks(index, ticks)

plt.ylabel("Time of Completion")
plt.xlabel("Process Rank")
plt.title('5MB')
plt.legend()
plt.show()
