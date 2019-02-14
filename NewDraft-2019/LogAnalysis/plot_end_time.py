import matplotlib.pyplot as plt
import numpy as np

data_in = ["LogEndTime:-run16_Complete_8_5mb_9",
"LogEndTime:-run16_SandersOld_8_5mb_9"]

data1 = open(data_in[0])
lines = data1.readlines()
x1 = []
y1 = []
for line in lines:
	splitLine = line.split()
	x1.append(splitLine[splitLine.index("rank")+1])
	y1.append(float(splitLine[1]))
X1 = np.arange(8)
plt.bar(X1, y1, width = 0.25, label = 'TwoTreeComplete')

data1 = open(data_in[1])
lines = data1.readlines()
x2 = []
y2 = []
for line in lines:
	splitLine = line.split()
	x2.append(splitLine[splitLine.index("rank")+1])
	y2.append(float(splitLine[1]))
X2 = np.arange(8) + 0.25
plt.bar(X2, y2, width = 0.25, label = 'TwoTreeSanders')
plt.legend()

index=[]
ticks=[]
for i in range(8):
	a=X1[i]
	b=X2[i]
	c=x1[i]
	d=x2[i]
	index.append(a)
	index.append(b)
	ticks.append(c)
	ticks.append(d)

plt.xticks(index, ticks)

plt.show()
