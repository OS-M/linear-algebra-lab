import matplotlib.pyplot as plt
import numpy as np

lines = open("task1.txt", "r").readlines()

x = np.array(lines[0].strip().split(' '))
y = [np.array(x.strip().split(' '), dtype=float) for x in lines[2:]]
names = np.array(lines[1].split('/'))

fig = plt.figure(figsize=(10, 5))

ax = fig.add_subplot(1, 2, 1)
ax.title.set_text('Prepare')
ax.plot(x, y[0], label=names[0])
ax.plot(x, y[1], label=names[1])
ax.set_xlabel("Matrix size")
ax.set_ylabel("Duration, sec.")
plt.legend()

ax = fig.add_subplot(1, 2, 2)
ax.title.set_text('Solve')
ax.plot(x, y[2], label=names[2])
ax.plot(x, y[3], label=names[3])
ax.set_xlabel("Matrix size")
ax.set_ylabel("Duration, sec.")
plt.legend()

# plt.show()
plt.savefig('task1.png')
