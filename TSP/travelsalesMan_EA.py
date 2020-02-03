import numpy as np
from numba import njit
import random
import matplotlib.pyplot as plt
import time
CORSSRATE=0.6
MUTATIONRATE=0.1
GENERATION=50000
GENE_SIZE=1000
POPULATION_SIZE=20
a = np.loadtxt('/Users/chiqu/class2019fall/EA/Assignment1/tsp.txt')
CITY_POSITION = np.reshape(a, (1000, 2))


def make_population():
    population=np.vstack([np.random.permutation(GENE_SIZE) for _ in range(POPULATION_SIZE)])
    return population

#@njit
def fitness(population):
    best_individual=population[0]
    best_fitness=getfitness(best_individual)
    fitness_list=[0]*POPULATION_SIZE
    sumfitness=0
    for i,individual in enumerate(population):
        fitness_list[i]=getfitness(individual)
        if fitness_list[i]>best_fitness:
            best_individual=individual
            best_fitness=fitness_list[i]
        sumfitness+=fitness_list[i]
    return sumfitness,best_individual,fitness_list,best_fitness

@njit
def getfitness(individual):
    distance=0.000000
    for i in range(1, len(individual)):
        # distance += ((CITY_POSITION[individual[i]][0] - CITY_POSITION[individual[i-1]][0]) ** 2 +
        #              (CITY_POSITION[individual[i]][1] - CITY_POSITION[individual[i-1]][1]) ** 2)**0.5
        distance += np.linalg.norm(CITY_POSITION[individual[i]] - CITY_POSITION[individual[i-1]])
    return 1000/distance
#@jit
def cross(mom,dad):
    index1 = random.randint(0, GENE_SIZE - 1)  # 随机生成突变起始位置 #
    index2 = random.randint(index1, GENE_SIZE - 1)
    crossgene = []
    for i in range(index1, index2):
        crossgene.append(mom[i])
    child=[]
    for flag,item in enumerate(dad):
        if flag==index1:
            child.extend(crossgene)
        if item not in crossgene:
            child.append(item)
    return child

def mutaion(individual):
    index = random.sample(range(0, GENE_SIZE), 2)
    new_individual=individual[:]
    new_individual[index[0]],new_individual[index[1]]=new_individual[index[1]],new_individual[index[0]]
    return new_individual
'''def mutaion(individual):
    index0 = random.randint(0, GENE_SIZE - 2)
    index1=find_best_mutation(individual,index0)
    new_individual=individual[:]
    new_individual[index0+1],new_individual[index1]=new_individual[index1],new_individual[index0+1]
    return new_individual
@njit
def find_best_mutation(individual,index1):
    best_distance=1000
    best_index=index1
    for index,item in enumerate(individual):
        if index!=index1:
            distance=((CITY_POSITION[individual[index]][0] - CITY_POSITION[individual[index1]][0]) ** 2 +
                     (CITY_POSITION[individual[index]][1] - CITY_POSITION[individual[index1]][1]) ** 2)**0.5
            if distance<best_distance:
                best_distance=distance
                best_index=index
    return best_index'''
#@njit
def select(population,sumfitness):
    flag=np.random.rand()
    flag = flag * sumfitness
    for individual in population:
        flag = flag - getfitness(individual)
        if flag <= 25:
            return individual

def select2(population,sumfitness):
    flag = np.random.rand()
    individual = population[int(flag * POPULATION_SIZE)]
    return individual

def born(population,sumfitness):
    dad=select(population,sumfitness)
    mom=select(population,sumfitness)
    rate = random.uniform(0, 1)
    if rate < CORSSRATE:
        child = cross(dad, mom)
    else:
        child = dad
    rate = random.uniform(0, 1)
    if rate < MUTATIONRATE:
        child = mutaion(child)
    return child
#@njit
def generation(population):
    sumfitness,best_individual,fitness_list,best_fitness=fitness(population)
    new_population=[]
    new_population.append(best_individual)
    while len(new_population) < POPULATION_SIZE:
        new_population.append(born(population,sumfitness))
    return new_population
    #generation_times += 1

def draw(individual):
    x=[0]*GENE_SIZE
    y=[0]*GENE_SIZE
    z = CITY_POSITION
    for i in range(0,GENE_SIZE):
        x[i]=CITY_POSITION[individual[i]][0]
        y[i]=CITY_POSITION[individual[i]][1]
    fig_short = plt.figure()
    ax1 = fig_short.add_subplot(1, 1, 1)
    ax1.plot(x, y, lw=0.5, marker='o', markersize=1, mfc='w')
    ax1.set_title('EA_distance~~')
    np.savetxt('GA_short_x_point_s1.txt',x)
    np.savetxt('GA_short_y_point_s1.txt',y)
    plt.show()

start = time.time()
if __name__=='__main__':
    population=make_population()
    i=GENERATION
    shortest_distance=[]
    while i>0:
        sumfitness, best_individual, fitness_list,best_fitness=fitness(population)
        print("distance:%f,  generation:%d"%(1000/best_fitness,GENERATION-i))
        shortest_distance.append(1000/best_fitness)
        population=generation(population)
        i-=1
    draw(best_individual)
    np.savetxt('best_individual_shorts1.txt',best_individual)
    np.savetxt('short_distance_shorts1.txt', shortest_distance)
    fig_shortest = plt.figure()
    ax1 = fig_shortest.add_subplot(1, 1, 1)
    times=[i for i in range(GENERATION)]
    ax1.plot(times, shortest_distance)
    ax1.set_title('EA(shortest path)')
    ax1.set_xlabel('Times')
    ax1.set_ylabel('Distance')
    plt.show()
stop = time.time()
print('time=',stop-start)
