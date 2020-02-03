from numba import jit
import numpy as np
from bokeh.plotting import figure, show
import matplotlib.pyplot as plt


a = np.loadtxt('tsp.txt')
x = np.array([])
y = np.array([])
for item in a:
    x = np.append(x, item[0])
    y = np.append(y, item[1])
times = 100000



min_a = np.array([])

def cal1():
    dismin1 = np.array([])
    min_a = np.array([])
    for looptimes in range(times):
        np.random.shuffle(a)
        distance = cal_distance(a, looptimes)
        if looptimes == 0:
            dismin = distance
        if distance > dismin:
            dismin = distance
            min_a = a
        print('looptimes',looptimes)
        print('dismin', dismin)
        if looptimes % 1 == 0:
            dismin1 = np.append(dismin1, np.array([dismin]))
            #print(dismin1)

    return (dismin1, dismin, min_a)

@jit
def cal_distance(a, looptimes):
    distance = 0
    for j in range(1, len(a)):
        distance += np.linalg.norm(a[j] - a[j - 1])
    return distance


# start = time.time()
# print(cal1())
# stop = time.time()
# print("time=", stop - start)
dismin1, min1, min_a1 = cal1()
dismin2, min2, min_a2 = cal1()
dismin3, min3, min_a3 = cal1()
dismin4, min4, min_a4 = cal1()
dismin5, min5, min_a5 = cal1()

def varN(n):
    return np.var(np.array([dismin1[n], dismin2[n],dismin3[n],dismin4[n],dismin5[n]]))

# var1 = varN(20000)
# var2 = varN(40000)
# var3 = varN(60000)
# var4 = varN(80000)
# y_error = [var1, var2, var3, var4]
# x1 = [20000, 40000, 60000, 80000]
# y1 = [dismin3[20000],dismin3[40000],dismin3[60000],dismin3[80000], ]


def creatx1(x):
    x1 = np.array([])
    for i in range(x):
        x1 = np.append(x1, [i])
    return x1

m = np.array([])
n = np.array([])


for item in min_a3:
    m = np.append(m, item[0])
    n = np.append(n, item[1])

p1 = figure(title="Random search minimum distance = {}".format(min3))
p1.circle(x, y)
p1.line(m, n, line_color = "steelblue", line_width = 0.5)
show(p1)

# x11 = creatx1(10000)
# plt.plot(x11, dismin3, '-g')
# plt.errorbar(x1, y1, y_error, fmt='o', capsize=5)
# plt.title('Random search minimum distance = {}'.format(min3))
# plt.show()
np.savetxt('randomsearch_long_orderpoint_x.txt',m)
np.savetxt('randomsearch_long_orderpoint_y.txt',n)
np.savetxt('randomsearch_long_distance3.txt', dismin3)
np.savetxt('randomsearch_long_distance1.txt', dismin1)
np.savetxt('randomsearch_long_distance2.txt', dismin2)
np.savetxt('randomsearch_long_distance4.txt', dismin4)
np.savetxt('randomsearch_long_distance5.txt', dismin5)


