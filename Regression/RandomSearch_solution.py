import numpy as np
import numexpr as ne
from random import random
from numba import jit
from bokeh.plotting import figure, show
import math
import time
import matplotlib.pyplot as plt
import functools

time1 = 500

def cos(q):
    return math.cos(q)
def sin(q):
    return math.sin((q))

operator = ['+', '-', '*', '/', 'cos', 'sin']
a = np.loadtxt('data.txt')
x0 = np.array([])
y0 = np.array([])
for i, item in enumerate(a):
    if i % 3 == 0:
        x0 = np.append(x0, item[0])
        y0 = np.append(y0, item[1])

def constant1():
    b = random.random()
    a = random.random()
    while a == 0:
        a = random.random()
    if b < 0.5:
        return str(10 * a)
    else:
        return '(' + '{}'.format(-10 * a) + ')'

#@jit
def generateHeap(heap):
    single = False
    Double = False
    for i in range(1, len(heap)):
        constant = str(constant1())
        c = ['', operator[np.random.randint(0, len(operator))], constant, 'x']
        parent = heap[i//2]
        if i == 1:
            heap[1] = c[1]
            continue
        if parent.isdigit() or parent == 'x':
            continue
        if parent != '':
            r1 = c[np.random.randint(1, len(c))]
            r2 = c[np.random.randint(2, len(c))]
            if Double:
                if i < 128:
                    heap[i] = r1
                else: heap[i] = r2
                Double = False
                continue
            if single:
                single = False
                continue
            if parent in ['sin', 'cos']:
                if i < 128:
                    heap[i] = r1
                else:
                    heap[i] = r2
                single = True
            if parent in ['+', '-', '*', '/']:
                if i < 128:
                    heap[i] = r1
                else:
                    heap[i] = r2
                Double = True
    #print('hh', heap)
    return heap

#@jit
def makevaluate(heap):
    for i in range(len(heap)-1, 0, -1):
        if heap[i] != '':
            if heap[i//2] in ['sin', 'cos']:
                heap[i//2] = heap[i//2] + '(' + heap[i] + ')'

            else:
                heap[i//2] = heap[i-1] + heap[i//2] + heap[i]
                heap[i-1] = ''
    return heap[1]

#@jit
def fitness(x0, func):
    sums = 0
    for i in range(len(x0)):
        func1 = func.replace('x', str(x0[i]))
        try:
            fit_single = abs(y0[i] - eval(func1))
        except:
            fit_single = 1.5
        sums += fit_single
    return sums/len(x0)

#@jit
def run_Random(times, x0):
    funcgood = ''
    goodfit = 100000
    shortest_fitness = []
    for i in range(times):
        heap = [''] * 256
        heap_new = generateHeap(heap)
        func = makevaluate(heap_new)
        #print(func)
        fitn = fitness(x0, func)
        if fitn < goodfit:
            goodfit = fitn
            print(goodfit)
            funcgood= func
        shortest_fitness.append(goodfit)
        print('fg', funcgood)
        print('fitnnnnnnnnn', goodfit)
        print('times=', i)
    np.savetxt('rs_goodfitness.txt', shortest_fitness)
    return funcgood, fitn

start = time.time()
func_good = []
func_good.append(str(run_Random(time1, x0)[0]))
print('000',type(func_good))
np.savetxt('rs_bestfunction.txt', func_good, fmt='%s')
end = time.time()
print('usetime=', end-start)
#
# print('fgggg', func_good)
# f = exec(func_good)
# print('f', f)

# y = eval(func_good.replace('x', str(0.3333)))
# # print('yy', y)
# list_new_y = np.array([])
# for i in x0:
#     y = eval(func_good.replace('x', str(i)))
#     list_new_y = np.append(list_new_y, y)
#     # print('i',i)
# # print('ly', list_new_y)
# # print(type(list_new_y))
# m = np.array([])
# for item in list_new_y:
#     m = np.append(m, item)


# TOOLS = "pan,wheel_zoom,box_zoom,reset,save,box_select"
# p1 = figure(title="Random Search", tools=TOOLS)
# p1.circle(x0, y0)
# p1.circle(x0, list_new_y, color = 'green')
# show(p1)



# fig_short = plt.figure()
# ax1 = fig_short.add_subplot(1, 1, 1)
# ax1.plot(x0, y0, lw=0.5, marker='o', markersize=1, mfc='w')
# ax1.set_title('dididi_distance~~')