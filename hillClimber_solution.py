from numba import jit
import numpy as np
from bokeh.plotting import figure, show
import matplotlib.pyplot as plt
import time


a = np.loadtxt('tsp.txt')
x = np.array([])
y = np.array([])
for item in a:
    x = np.append(x, item[0])
    y = np.append(y, item[1])
sizedata = 1000
times = 100000

idx1 = np.random.randint(sizedata, size=1)
idx2 = np.random.randint(sizedata, size=1)

min_a = np.array([])

@jit
def cal_distance(a):
    distance = 0
    for j in range(1, len(a)):
        distance += np.linalg.norm(a[j] - a[j - 1])
    distance += np.linalg.norm(a[len(a)] - a[0])
    return distance

distance_original = cal_distance(a)

@jit
def cal1(distance_original):
    dismin1 = np.array([])
    for looptimes in range(times):
        idx1 = np.random.randint(sizedata, size=1)
        idx2 = np.random.randint(sizedata, size=1)

        a[[idx1, idx2], :] = a[[idx2, idx1], :]
        distance1 = cal_distance(a)

        if distance1 < distance_original:
            a[[idx1, idx2], :] = a[[idx2, idx1], :]
        else:
            distance_original = distance1

        print(distance_original, looptimes)
        dismin1 = np.append(dismin1, np.array([distance_original]))
        if looptimes in range(300, 500):

    return distance_original, a, dismin1
def movie(x, y):
    plt.ion()
    plt.figure(3)
    plt.scatter(x, y, s=0.8)
    plt.plot(x, y, linewidth=0.2)
    plt.pause(1e-3)

start = time.time()

distance3, a3, dismin3 = cal1(distance_original)


x=[]
y=[]
for i in range(len(a3)):
    x.append([a3[i][0]])
    y.append([a3[i][1]])
print(x)
print(y)
def movie():
    plt.ion()
    plt.figure(3)
    plt.scatter(x, y, s=0.8)
    plt.plot(x, y, linewidth=0.4)
    plt.pause(1e-3)
movie()
