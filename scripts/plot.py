import sys
import re
from matplotlib import pyplot as plt
import numpy as np


def linspace(start, stop, n):
    if n == 1:
        yield stop
        return
    h = (stop - start) / (n - 1)
    for i in range(n):
        yield start + h * i
# from import


if len(sys.argv) > 1:
    filepath = sys.argv[1]
    filepath_spline = sys.argv[2]
else:
    filepath = "/home/manuel/hrg/spline_interpolation/output/trajectory.txt"
    filepath_spline = "/home/manuel/hrg/spline_interpolation/output/spline.txt"


with open(filepath) as fin:
    lines = [l.strip() for l in fin.readlines()]


with open(filepath_spline) as fin:
    lines_spline = [l.strip() for l in fin.readlines()]


print('nr of samples; ', len(lines))
print('nr of samples; ', len(lines_spline))


data = [[float(e) for e in line.split(' ')] for line in lines]
data_spline = [[float(e) for e in line.split(' ')] for line in lines_spline]

entries = []
for i in range(3):
    entries.append([d[i] for d in data])

entries_spline = []
for i in range(3):
    entries_spline.append([d[i] for d in data_spline])


# 2d
fig, axs = plt.subplots(3)
fig.suptitle('Cartesian')
for j in range(3):
    axs[j].plot(list(linspace(0, len(entries_spline[0]),
                              len(entries[0]))), entries[j])
    axs[j].plot(entries_spline[j])
# plt.show()


# 3d
# mpl.rcParams['legend.fontsize'] = 10
fig3d = plt.figure()
ax = fig3d.gca(projection='3d')
ax.plot(entries[0], entries[1], entries[2], ':',
        color='gray')
ax.plot(entries_spline[0], entries_spline[1],
        entries_spline[2])
ax.scatter(entries[0], entries[1], entries[2], color='red')

plt.show()
