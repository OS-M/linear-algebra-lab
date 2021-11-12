import matplotlib.pyplot as plt
import numpy as np

lines = open("task3.txt", "r").readlines()

x = np.array(lines[0].strip().split(' '))
y = [np.array(x.strip().split(' '), dtype=float) for x in lines[2:]]
names = np.array(lines[1].split('/'))

fig = plt.figure(figsize=(10, 5))

ax = fig.add_subplot(1, 1, 1)
ax.set_xticklabels([int(float(i) * 100) for i in x])
ax.title.set_text('Duration')
ax.plot(x, y[0], label=names[0])
ax.plot(x, y[1], label=names[1])
ax.set_xlabel("Matrix sparsity, %")
ax.set_ylabel("Duration, sec.")
plt.legend()

# plt.show()
plt.savefig('task3.png')
