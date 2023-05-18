#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from pandas import *

fig = plt.figure()
ax = plt.axes(projection = '3d')

with open('random_on_hemisphere.csv', mode='r') as f:
    data = read_csv(f)

    ax.scatter(data['X'].tolist(), data['Y'].tolist(), data['Z'].tolist())


plt.show()
