import numpy as np
import numexpr as ne
from random import random
from numba import jit
from bokeh.plotting import figure, show
import math
import time
import matplotlib.pyplot as plt
# from Assignment2.random_search.RandomSearch_hw2 import generateHeap, makevaluate, fitness, constant1

def cos(q):
    return math.cos(q)
def sin(q):
    return math.sin((q))
operator = ['+', '-', '*', '/', 'cos', 'sin']
time1 = 10000
a = np.loadtxt('data.txt')
x0 = np.array([])
y0 = np.array([])
for i, item in enumerate(a):
    if i % 3 == 0:
        x0 = np.append(x0, item[0])
        y0 = np.append(y0, item[1])

mutation_rate = 0.2
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
        if parent.isdigit() or parent == 'x' or ('(-' in parent):
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

def makevaluate(heap):
    heap_n = []
    for i in heap:
        heap_n.append(i)
    #print(heap_n)
    for i in range(len(heap)-1, 0, -1):
        if heap_n[i] != '':
            #print(type(heap_n[i-1]))
            if heap_n[i//2] in ['sin', 'cos']:
                heap_n[i//2] = heap_n[i//2] + '(' + heap_n[i] + ')'

            else:
                heap_n[i//2] = heap_n[i-1] + heap_n[i//2] + heap_n[i]
                heap_n[i-1] = ''
    return heap_n[1]

def constant1():
    b = random.random()
    a = random.random()
    while a == 0:
        a = random.random()
    if b < 0.5:
        return str(10 * a)
    else:
        return '(' + '{}'.format(-10 * a) + ')'


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

def operator_generator(heap, i):
    op = operator[np.random.randint(0, len(operator)-2)]
    sc = operator[np.random.randint(len(operator)-2, len(operator))]
    changerate = np.random.rand()
    if heap[i] in '+-*/':
        if changerate < 0.5:
            heap[i] = op
        else:
            heap[i] = sc
            heap[i*2+1] = ''
        return heap
    if heap[i] == 'sin':
        if changerate < 0.5:
            heap[i] = 'cos'
        else:
            heap[i] = op
            heap[i*2+1] = num_genertator()
        return heap
    elif heap[i] == 'cos':
        if changerate < 0.3:
            heap[i] = 'sin'
        else:
            heap[i] = op
            heap[i*2+1] = num_genertator()
        return heap
    else:
        if i < 128:
            child1, child2 = i * 2, i * 2 + 1
            if changerate < 0.5:
                heap[i] = op
                heap[child1], heap[child2] = num_genertator(), num_genertator()
            else:
                heap[i] = sc
                heap[child1] = num_genertator()
        return heap


def num_genertator():
    rate = np.random.rand()
    if rate < 0.5:
        return constant1()
    else:
        return 'x'

def mutation(heap):
    indexmutation = []
    for i in range(1, len(heap)):
        if heap[i] != '':
            indexmutation = np.append(indexmutation, i)
    index = int(np.random.choice(indexmutation))
    #print(index)
    change = np.random.rand()
    if '.' in heap[index]:
        b = ''
        for j in heap[index]:
            b = b + j
        if change < 0.25:
            if '(' in heap[index]:
                heap[index] = '(' + str(float(b[1:-1]) * 1.1) + ')'
            else: heap[index] = str(float(heap[index])* 1.1)
        elif change < 0.5:
            if '(' in heap[index]:
                heap[index] = str(float(b[1:-1]) / 1.1)
            else: heap[index] = str(float(heap[index])/ 1.1)
        elif change < 0.75:
            heap[index] = 'x'
        else:
            if index < 128:
                heap = operator_generator(heap, index)
    elif heap[index] == 'x':
        if change < 0.2:
            heap[index] = str(constant1())
        elif change < 0.5:
            heap = operator_generator(heap, index)
    elif heap[index] in operator:
        if change < 0.5:
            if index < 128:
                heap = operator_generator(heap, index)
        if change < 0.75:
            heap[index] = constant1()
        else:
            heap[index] = 'x'
    for i in range(2, len(heap)):
        if heap[i//2] == '' or heap[i//2] not in operator:
            if i != 1:
                heap[i] = ''
    return heap

def run_hillclimber(times, a):
    heap = [''] * 256
    heap_1 = generateHeap(heap)
    goodfit = 100000
    funcgood = ''
    shortest_fitness = []
    for i in range(times):
        heap_2 = mutation(heap_1)
        func = makevaluate(heap_2)
        fitn = fitness(x0, func)
        if fitn < goodfit:
            goodfit = fitn
            funcgood= func
            heap_1 = heap_2
        shortest_fitness.append(goodfit)
        print('good:', goodfit)
        print('funcgood=', funcgood)
        print('times', i)
    np.savetxt('hc_best_fitness',shortest_fitness)
    #np.savetxt('hc_overall_fitness', over_fitness)
    np.savetxt('hc_bestfunction', funcgood, fmt='%s')
    return funcgood, fitn

func_good = run_hillclimber(time1, a)[0]
# func1111 = '(2.4983520005*((sin(sin(cos((-1.6981261635-4.1125370036))))+sin((x+(sin(2.9689779351)-2.8835261086))))*sin(x)))'
# func2222 = '2.5*(0.7+sin(x-2.711767314084740))*sin(x)'
#print(fitness(a, '2.5*(0.7+sin(x-2.711767314084740))*sin(x)'))
# print(fitness(a, func1111))
# print(fitness(a, func2222))
# func3 = 'sin(x)+cos((-9.15052698578442))'
# func4 = 'cos(9.984082701859851)+sin(x)'
# func5 = 'sin(x)-x/x'
# list_new_y = np.array([])
# for i in x0:
#     y = eval(func5.replace('x', str(i)))
#     list_new_y = np.append(list_new_y, y)
#     # print('i',i)
# # print('ly', list_new_y)
# # print(type(list_new_y))
# m = np.array([])
# for item in list_new_y:
#     m = np.append(m, item)


TOOLS = "pan,wheel_zoom,box_zoom,reset,save,box_select"
p1 = figure(title="Random Search", tools=TOOLS)
p1.circle(x0, y0)
p1.circle(x0, list_new_y, color = 'green', size = 0.1)
show(p1)


