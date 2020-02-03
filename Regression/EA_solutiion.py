import numpy as np
import numexpr as ne
from random import random
from numba import jit
from bokeh.plotting import figure, show
import math
import random
import time
import matplotlib.pyplot as plt
# from Assignment2.random_search.RandomSearch_hw2 import generateHeap, makevaluate, fitness, constant1

def cos(q):
    return math.cos(q)
def sin(q):
    return math.sin(q)
operator = ['+', '-', '*', '/', 'cos', 'sin']
generationtimes = 100
population_size = 20
group_size = 20
crossrate = 0.5
mutationrate = 0.1
selectrange = population_size/2

a = np.loadtxt('data.txt')
x0 = np.array([])
y0 = np.array([])
for item in a:
    x0 = np.append(x0, item[0])
    y0 = np.append(y0, item[1])

mutation_rate = 0.2
def generateHeap():
    heap = ['']*256
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

    #print('dddd',heap_n[1])
    return heap_n[1]

def constant1():
    a = 0
    b = np.random.rand()
    while a == 0:
        a = np.random.rand()
    if b < 0.5:
        return str(10 * a)
    else:
        return '(' + '{}'.format(-10 * a) + ')'


def fitness(x0, func):
    sums = 0
    for i in range(len(x0)):
        func1 = func.replace('x', str(x0[i]))
        fit_single = abs(y0[i] - eval(func1))
        sums += fit_single
    return sums/len(x0)

def combinefit(heap):
    # heap_2 = mutation(heap)
    func = makevaluate(heap)
    fitn = fitness(a, func)
    return fitn


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

def find_index(heap):
    indexmutation = []
    for i in range(1, len(heap)):
        if heap[i] != '':
            indexmutation.append(i)
    return int(indexmutation[random.randrange(len(indexmutation))])


def mutationforGA(heap1):
    heap = []
    for i in range(len(heap1)):
        heap.append(heap1[i])
    index = find_index(heap)
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
    # print('mutation')
    return heap
def modify_none(heap):
    for i in range(2, len(heap)):
        if heap[i//2] == '' or heap[i//2] not in operator:
            heap[i] = ''
    return heap

def find_spring_index(heap, index):
    index_list = [index]
    for i in range(len(heap)):
        if heap[i] != '':
            if i // 2 in index_list:
                index_list.append(i)
    return index_list

def switch(p1_switch, p_clist):
    for s in range(1, len(p_clist)):
        if p_clist[s] / p_clist[s-1] == 2:
            p1_switch.append(p1_switch[-1]*2)
        elif p_clist[s]-1 == p_clist[s-1]:
            p1_switch.append(p1_switch[-1]+1)
    return p1_switch

def init_pop(population_size):
    poplist = []
    for i in range(population_size):
        poplist.append(generateHeap())
    return poplist

def crossover(p1, p2):
    #print('make1',makevaluate(p1))
    #print('make2',makevaluate((p2)))
    n = 256

    p1_switch_1 = [0]
    p2_switch_2 = [0]
    p1_clist = [0]
    p2_clist = [0]
    while p1_switch_1[-1]+n > 256 or p2_switch_2[-1]-n > 256:
        index1 = find_index(p1)
        index2 = find_index(p2)
        #print(index1)
        #print(index2)
        place1, place2 = index1, index2
        p1_clist = find_spring_index(p1, index1)
        p2_clist = find_spring_index(p2, index2)
        # if index1 < index2:
        n = index2 - index1
        #print('n', n)
    #print('list1', p1_clist)
    #print('list2', p2_clist)
        mid = []
        for i in p1:
            mid.append(i)
        # print('mid', mid)
        p1_switch = [p1_clist[0]]
        # p2_switch = [p2_clist[0]]
        p1_switch_1 = switch(p1_switch, p2_clist)
        # p2_switch_2 = switch(p2_switch, p1_clist)
    #print('swithc1',p1_switch_1)
    #print('switch2',p2_switch_2)
    for j, char in enumerate(p1_switch_1):
        p1[char] = p2[p2_clist[j]]
        #print('p1j',p1[j])
    # for j2, char in enumerate(p2_switch_2):
    #     p2[char] = mid[p1_clist[j2]]
        #print('p2j',p2[j2])
    p1 = modify_none(p1)
    # p2 = modify_none(p2)
    # print('crossover')
    return p1


def select1(population):
    flag = np.random.rand()
    individual = population[int(flag * population_size)]
    return individual
def select2(population):
    new_people = sorted(population, key=combinefit)
    return new_people[:selectrange]
def select3(population100):
    poplist = ['']*5
    poplist[0] = sorted(population100, key=combinefit)[:20]
    poplist[1] = sorted(population100, key=combinefit)[20:40]
    poplist[2] = sorted(population100, key=combinefit)[40:60]
    poplist[3] = sorted(population100, key=combinefit)[60:80]
    poplist[4] = sorted(population100, key=combinefit)[80:100]
    return poplist



selectsize = int(population_size/2)

def fitnessall(population):
    best_individual=population[0]
    # print(type(best_individual))
    best_fitness=combinefit(best_individual)
    fitness_list=[0]*population_size
    sumfitness=0
    best_indlist = [0]*(selectsize)
    for i, individual in enumerate(population):
        fitness_list[i] = combinefit(individual)
        if fitness_list[i]<best_fitness:
            best_individual=individual
            best_fitness=fitness_list[i]
        sumfitness+=fitness_list[i]
    best_indlist1 = sorted(population, key=combinefit)
    print('33333',sorted(fitness_list))
    for i in range(selectsize):
        best_indlist[i] = best_indlist1[i]
    print('bf', best_fitness)
    return sumfitness,best_individual,fitness_list,best_fitness, best_indlist

# poplation = init_pop(population_size)
# for i in range(population_size-1):
#     print(fitness(a, makevaluate(crossover(poplation[i], poplation[i+1]))))
def generation(population, best_inlist, best_individual):
    # sumfitness,best_individual,fitness_list,best_fitness, best_inlist=fitnessall(population)
    #print('1')
    new_population=[0]*len(best_inlist)
    new_population[0] = best_individual
    for i in range(0, selectsize):
        new_population[i] = best_inlist[i]
    while len(new_population) < population_size:
        #print('3')
        new_population.append(born(population))
    #print(fitness_list)
    return new_population


def born(population):
    dad=select1(population)
    mom=select1(population)
    rate = random.uniform(0, 1)
    if rate < crossrate:
        child = crossover(dad, mom)
    if rate < mutationrate:
        child = mutationforGA(dad)
    else:
        child = generateHeap()
    return child

def run_GA1(i):
    population = init_pop(population_size)
    shortest_fitness = []
    over_fitness = []
    while i>0:
        sumfitness, best_individual, fitness_list,best_fitness, b = fitnessall(population)
        print("dis:%f,  generation:%d"%(best_fitness,generationtimes-i))
        shortest_fitness.append(best_fitness)
        over_fitness.append(fitness_list)
        population=generation(population, b, sumfitness)
        print('y=', makevaluate(best_individual))
        i-=1
    np.savetxt('GA1_best_fitness',shortest_fitness)
    np.savetxt('GA1_overall_fitness', over_fitness)
    np.savetxt('GA1_bestfunction', makevaluate(best_individual))


run_GA1(generationtimes)

# h1 = generateHeap()
# h2 = generateHeap()
# func1 = makevaluate(h1)
# func2 = makevaluate(h2)
# print('func1', func1)
# print('func2', func2)
# h11, h22 = crossover(h1, h2)
# func11 = makevaluate(h11)
# func22 = makevaluate(h22)
# print('11', func11)
# print('22', func22)
# poplation = init_pop(population_size)
# print(generation(poplation))

