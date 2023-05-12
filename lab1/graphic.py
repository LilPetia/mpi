import matplotlib.pyplot as plt
import csv
import numpy as np
x = []
y = []
time_sum = 0.0
with open('output.csv','r') as csvfile:
    lines = csv.reader(csvfile, delimiter=',')
    for  row in lines:
        for i in range(len(row)):
          try:
            if (float(row[i]) > 0.00002): continue
            #print(row[i])
            y.append(np.float64(row[i])* 1000000000)
            time_sum += y[-1]
            x.append(i)
          except ValueError:
            pass
print("Среднее время:", time_sum / len(x) / 2)
plt.plot(x, y, color = 'g',
         marker = 'o',label = "time in ns")
  
plt.xticks(rotation = 25)
plt.xlabel('count')
plt.ylabel('time in ns')
# plt.title('Weather Report', fontsize = 20)
plt.grid()
plt.legend()
plt.show()