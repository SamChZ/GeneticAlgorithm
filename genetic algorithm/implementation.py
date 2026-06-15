'''This file contains a skeleton implementation for the practical assignment 
for the NACO 21/22 course. 

Your Genetic Algorithm should be callable like so:
    >>> problem = ioh.get_problem(...)
    >>> ga = GeneticAlgorithm(...)
    >>> ga(problem)

In order to ensure this, please inherit from the provided Algorithm interface,
in a similar fashion as the RandomSearch Example:
    >>> class GeneticAlgorithm(Algorithm):
    >>>     ...

Please be sure to use don't change this name (GeneticAlgorithm) for your implementation.

If you override the constructor in your algoritm code, please be sure to 
call super().__init__, in order to correctly setup the interface. 

Only use keyword arguments for your __init__ method. This is a requirement for the
test script to properly evaluate your algoritm.
'''

import ioh
import random
import numpy as np
from numpy.core.numeric import cross
from algorithm import Algorithm
import pandas as pd
import typing
import ast
import Levenshtein as lsh
import gc

currentRow = 0
rule = 0
t = 2
ct = 0
k = 2
dim = 0
data = pd.read_csv('ca_input.csv')


class RandomSearch(Algorithm):
    def __init__(self, max_iterations: int = 10000):
        super().__init__(max_iterations=max_iterations)
    '''An example of Random Search.'''
    def __call__(self, problem: ioh.problem.Integer, k) -> None:
        self.y_best: float = float("inf")
        for iteration in range(self.max_iterations):
            # Generate a random bit string
            x: list[int] = [random.randint(0, k-1) for _ in range(problem.meta_data.n_variables)]
            # Call the problem in order to get the y value    
            y: float = problem(x)
            # update the current state
            self.y_best = max(self.y_best, y)


class CellularAutomata:
    '''Skeleton CA, you should implement this.'''
    
    def __init__(self, rule_number: int):
        self.rule = rule_number #save the given rule
 

    def __call__(self, c0: typing.List[int], t: int) -> typing.List[int]:
        '''Evaluate for T timesteps. Return Ct for a given C0.'''
        rules=[] #initialaze an empty transition rules list
        global k
        if k == 2: #if we deal with a binairy representation of the problem
            binary = '{0:08b}'.format(self.rule) #transform the transition rule into an 8-bit binary number
            for j in range(8):
                rules.append(int(binary[j]))#append to our transition rules list bit by bit the binary number
            c0c = c0.copy()

            for i in range(t):
                
                c0c.insert(0, 0)#add zero-cell to the left border
                c0c.append(0)#add zero-cell to the right border

                ccopy = c0c.copy()#copy of the current grid to be able to perform updates on the grid without 
                                  #changing the values of the current grid
                for j in range(1,len(c0c)-1):#check for every cell, except for the zero-cells what their neighorhood
                                             #looks like
                    if c0c[j-1]==1 and c0c[j]==1 and c0c[j+1]==1:#update the cell according to the transition rules 
                        ccopy[j] = rules[0]
                    elif c0c[j-1]==1 and c0c[j]==1 and c0c[j+1]==0:
                        ccopy[j] = rules[1]
                    elif c0c[j-1]==1 and c0c[j]==0 and c0c[j+1]==1:
                        ccopy[j] =rules[2]
                    elif c0c[j-1]==1 and c0c[j]==0 and c0c[j+1]==0:
                        ccopy[j] =rules[3]
                    elif c0c[j-1]==0 and c0c[j]==1 and c0c[j+1]==1:
                        ccopy[j] =rules[4]
                    elif c0c[j-1]==0 and c0c[j]==1 and c0c[j+1]==0:
                        ccopy[j] =rules[5]
                    elif c0c[j-1]==0 and c0c[j]==0 and c0c[j+1]==1:
                        ccopy[j] =rules[6]
                    elif c0c[j-1]==0 and c0c[j]==0 and c0c[j+1]==0:
                        ccopy[j] =rules[7]
                        
                c0end=ccopy.copy()
                c0end.pop(0)#remove the left zero cell
                c0end.pop(-1)#remove the right zero cell
                c0c=c0end
            return c0end
        if k == 3:#if we deal with a ternary representation of the problem
            ternary=np.base_repr(self.rule,3)#transform the transition rule into a ternary number
            rules = [0]*27#make 27 entries with value 0
            for i in range(1,len(ternary)+1):
                rules[-i]=int(ternary[-i])#loop through our transition rules list and ternary number backwards
                                          #and update the entries of the rules list to the value of the corresponding bit
            c0c = c0.copy()
            for i in range(t):  
                c0c.insert(0, 0)#add zero-cell to the left border
                c0c.append(0)#add zero-cell to the right border
                ccopy = c0c.copy()#copy of the current grid to be able to perform updates on the grid without 
                                  #changing the values of the current grid
                for j in range(1,len(c0c)-1):#check for every cell, except for the zero-cells what their neighorhood
                                             #looks like
                    if c0c[j-1]==2 and c0c[j]==2 and c0c[j+1]==2:#update the cell according to the transition rules
                        ccopy[j] = rules[0]
                    elif c0c[j-1]==2 and c0c[j]==2 and c0c[j+1]==1:
                        ccopy[j] = rules[1]
                    elif c0c[j-1]==2 and c0c[j]==2 and c0c[j+1]==0:
                        ccopy[j] = rules[2]
                    elif c0c[j-1]==2 and c0c[j]==1 and c0c[j+1]==2:
                        ccopy[j] = rules[3]
                    elif c0c[j-1]==2 and c0c[j]==1 and c0c[j+1]==1:
                        ccopy[j] = rules[4]
                    elif c0c[j-1]==2 and c0c[j]==1 and c0c[j+1]==0:
                        ccopy[j] = rules[5]
                    elif c0c[j-1]==2 and c0c[j]==0 and c0c[j+1]==2:
                        ccopy[j] = rules[6]
                    elif c0c[j-1]==2 and c0c[j]==0 and c0c[j+1]==1:
                        ccopy[j] = rules[7]
                    elif c0c[j-1]==2 and c0c[j]==0 and c0c[j+1]==0:
                        ccopy[j] = rules[8]
                    elif c0c[j-1]==1 and c0c[j]==2 and c0c[j+1]==2:
                        ccopy[j] = rules[9]
                    elif c0c[j-1]==1 and c0c[j]==2 and c0c[j+1]==1:
                        ccopy[j] = rules[10]
                    elif c0c[j-1]==1 and c0c[j]==2 and c0c[j+1]==0:
                        ccopy[j] = rules[11]
                    elif c0c[j-1]==1 and c0c[j]==1 and c0c[j+1]==2:
                        ccopy[j] = rules[12]
                    elif c0c[j-1]==1 and c0c[j]==1 and c0c[j+1]==1:
                        ccopy[j] = rules[13]
                    elif c0c[j-1]==1 and c0c[j]==1 and c0c[j+1]==0:
                        ccopy[j] = rules[14]
                    elif c0c[j-1]==1 and c0c[j]==0 and c0c[j+1]==2:
                        ccopy[j] = rules[15]
                    elif c0c[j-1]==1 and c0c[j]==0 and c0c[j+1]==1:
                        ccopy[j] = rules[16]
                    elif c0c[j-1]==1 and c0c[j]==0 and c0c[j+1]==0:
                        ccopy[j] = rules[17]
                    elif c0c[j-1]==0 and c0c[j]==2 and c0c[j+1]==2:
                        ccopy[j] = rules[18]
                    elif c0c[j-1]==0 and c0c[j]==2 and c0c[j+1]==1:
                        ccopy[j] = rules[19]
                    elif c0c[j-1]==0 and c0c[j]==2 and c0c[j+1]==0:
                        ccopy[j] = rules[20]
                    elif c0c[j-1]==0 and c0c[j]==1 and c0c[j+1]==2:
                        ccopy[j] = rules[21]
                    elif c0c[j-1]==0 and c0c[j]==1 and c0c[j+1]==1:
                        ccopy[j] = rules[22]
                    elif c0c[j-1]==0 and c0c[j]==1 and c0c[j+1]==0:
                        ccopy[j] = rules[23]
                    elif c0c[j-1]==0 and c0c[j]==0 and c0c[j+1]==2:
                        ccopy[j] = rules[24]
                    elif c0c[j-1]==0 and c0c[j]==0 and c0c[j+1]==1:
                        ccopy[j] = rules[25]
                    elif c0c[j-1]==0 and c0c[j]==0 and c0c[j+1]==0:
                        ccopy[j] = rules[26]
                
                c0end=ccopy.copy()
                c0end.pop(0)#remove the left zero cell
                c0end.pop(-1)#remove the right zero cell
                c0c=c0end
            return c0end


def objective_function(c0_prime: typing.List[int]) -> float:
    ca = CellularAutomata(rule)
    ct_prime = ca(c0_prime, t)
    count = 0
    #compare ct with ct_prime
    for i in range (len(ct)):
        if ct_prime[i] == ct[i]:
            count +=1
    similarity = count/len(ct)
    return similarity


def LevenshteinDistance(c0_prime: typing.List[int]) -> float:
    ca = CellularAutomata(rule)
    ct_prime = ca(c0_prime, t)
    ct_prime_string = ''.join(str(e) for e in ct_prime)
    ct_string = ''.join(str(e) for e in ct)
    similarity = lsh.distance(ct_prime_string, ct_string)
    return  (dim - similarity)/60




            
#our implementation of a genetic algorithm
class GeneticAlgorithm(Algorithm):
    '''A skeleton (minimal) implementation of your Genetic Algorithm.'''
    
    #initialization of the algorithm's functions
    def __init__(self, init_method = 'Random', generation_method = 'Random', selection_method = 'Tournament', crossover_method = 'Uniform', mutation_method = 'Swap'):
        super().__init__(max_iterations = 10000)
        if init_method == 'Random':
            self.initialize_population = self._random_initialization
        else:
            self.initialize_population = self._half_initialization
        if generation_method == 'Random':
            self.generate_candidate = self._random_generation
        else:
            self.generate_candidate = self._bitflip_generation
        if selection_method == 'Rank':
            self.perform_selection = self._rank_selection
        elif selection_method == 'Tournament':
            self.perform_selection = self._tournament_selection
        if crossover_method == 'Uniform':
            self.perform_crossover = self._uniform_crossover
        elif crossover_method == 'OnePoint':
            self.perform_crossover = self._one_point_crossover
        elif crossover_method == 'KPoint':
            self.perform_crossover = self._k_point_crossover
        if mutation_method == 'Bitstring':
            self.perform_mutation = self._bitstring_mutation
        elif mutation_method == 'Bitflip':
            self.perform_mutation = self._bitflip_mutation
        elif mutation_method == 'Swap':
            self.perform_mutation = self._swap_mutation
        data = pd.read_csv('ca_input.csv', delimiter=',')
        row = 0
        self.ct = data.iloc[row,3]
        self.rule = data.iloc[row,1]
        self.t = data.iloc[row,2]

    

    #creates a single random bitstring
    def _random_initialization(self, dim):
        if k == 2:
            self.x = np.array(np.random.uniform(size=dim) > 0.5,dtype=int)
        elif k == 3:
            self.x = np.array(np.random.uniform(size=60)*3,dtype=int)

    #creates a single bitstring that's half 1s and half 0s
    def _half_initialization(self, dim):
        x = np.zeros(dim)
        idxs = np.random.choice(range(dim), int(dim/2), replace = False)
        x[idxs] = 1
        self.x = x

    def _random_generation(self, x):
        return np.array(np.random.uniform(size=x.shape) > 0.5, dtype=int)
    
    def _bitflip_generation(self, x):
        idx = np.random.randint(0, len(x))
        new_x = x.copy()
        new_x[idx] = 1-x[idx]
        return new_x

    #initializes multiple bitstrings
    #populationAmount = how many bits we want to generate
    def _initialize_multiple(self, dim, populationAmount):
        population = []
        for i in range(populationAmount):
            self.initialize_population(dim)
            population.append(self.x)
        return population


    #returns the indices of the arrays with the highest fitness
    #winnerAmount is the amount of indeces returned
    def _rank_selection(self, fitnessArray, winnerAmount):
        return np.argpartition(fitnessArray, -winnerAmount)[-winnerAmount:]

    def _tournament_selection(self, population, winnerAmount):
        tournamentGroups = []
        populationIndices = range(len(population))
        stepSize = int(len(population)/winnerAmount)
        for i in range(0, len(population), stepSize):
            #tournamentGroups.append(population[i:i + winnerAmount])
            tournamentGroups.append(populationIndices[i:i + winnerAmount])
        winners = np.zeros(shape=(len(tournamentGroups)))
        groupNumber = 0
        for group in tournamentGroups:
            groupFitness = [0]*len(group)
            for i in range(len(group)):
                groupFitness[i] = population[i + 10 * groupNumber]
            helpInt = (np.argpartition(groupFitness, -1)[-1:])
            winners[groupNumber] = helpInt[0] + 10 * groupNumber
            groupNumber = groupNumber + 1
        winners = [int(x) for x in winners]
        return winners

    
    #for every child, it loops over the genetic data of all parents bit by bit
    #every bit it chooses one of the parents to inherit from, with each of the 
    #parents having a 1/[amount of parents] chance to pass on their bit
    #amountOfChildren is the amount of children created from the parents
    def _uniform_crossover(self, parents, amountOfChildren):
        children = np.empty((amountOfChildren, (len(parents[0]))), int)
        for child in range(amountOfChildren):
            for i in range(len(parents[0])):
                children[child,i] = np.random.choice(parents[:,i], 1)
        return children


    def _one_point_crossover(self, parents, amountOfChildren):
        children = np.empty((amountOfChildren, (len(parents[0]))), int)
        crossoverPoint = np.random.randint(0, len(parents[0]))
        for child in range(amountOfChildren):
            for i in range(len(parents[0])):
                if i < crossoverPoint:
                    children[child,i] = parents[0,i]
                else:
                    children[child,i] = parents[1,i]
        return children


    def _k_point_crossover(self, parents, amountOfChildren):
        kpoints = 3
        children = np.empty((amountOfChildren, (len(parents[0]))), int)
        crossoverPoints = np.zeros(shape=kpoints)
        for i in range(len(crossoverPoints)):
            crossoverPoints[i] = np.random.randint(0, len(parents[0]))
        crossoverPoints = np.sort(crossoverPoints)
        currentPoint = 0
        for child in range(amountOfChildren):
            for i in range(len(parents[0])):
                if currentPoint < kpoints and i >= crossoverPoints[currentPoint]:
                    currentPoint = currentPoint + 1
                children[child,i] = parents[currentPoint%2,i]
        return children


    #loops over all the bits of all children and toggles between 0 and 1 with
    #a random chance. if no chance value is given, it generates one of
    #1/[length of genetic code] so approximately 1 mutation happens per child
    def _bitstring_mutation(self, children, chance = "no value"):
        if chance == "no value":
            chance = 2/len(children[0])
        for child in children:
            chance = np.random.randint(1,10)
            again = 1
            if chance <= 1:
                again = 2
            for i in range(again): #sometimes mutates twice
                for i in child:
                    if np.random.uniform(0,1) < chance:
                        if k == 2:
                            if child[i] == 1:
                                child[i] = 0
                            else:
                                child[i] = 1
                        else: #if k == 3
                            if child[i] == 0:
                                child[i] = 1
                            elif child[i] == 1:
                                child[i] = 2
                            else:
                                child[i] = 0
        return children

    #inverts all the bits of all children
    #incompatible with k == 3
    def _bitflip_mutation(self, children):
        for child in children:
            for i in child:
                if child[i] == 1:
                    child[i] = 0
                else:
                    child[i] = 1
        return children  

    #picks two random bits in the genome and swaps their values
    def _swap_mutation(self, children):
        rng = np.random.default_rng()
        for child in children:
            chance = np.random.randint(1,10)
            again = 1
            if chance <= 1:
                again = 2
            for i in range(again): #sometimes mutates twice
                swapLocations = rng.choice(len(children[0]), 2, replace = False)
                helpInt = child[swapLocations[0]]
                child[swapLocations[0]] = child[swapLocations[1]]
                child[swapLocations[1]] = helpInt
        return children
        
    #calculates the amount of 1s in an array
    def _calculate_fitness(self, array):
        return sum(array)

    #runs the problem until completion or until too many loops have been run
    def __call__(self, problem: ioh.problem.Integer, populationAmount = 100, amountOfParents = 'No Value') -> None:
        global dim
        dim = problem.meta_data.n_variables
        populationAmount = int(populationAmount)
        #must have at least two strings
        if populationAmount < 2:
            populationAmount = 2
        #if no parent value is given, default to the amount of strings in your
        #population divided by 10 with a minimum of 2
        if amountOfParents == 'No Value':
            amountOfParents = max(2,int(populationAmount/10))
        if self.perform_crossover == self._one_point_crossover: #only two parents allowed
            amountOfParents = 2
        amountOfParents = int(amountOfParents)
        #can't have more parents than strings
        if amountOfParents > populationAmount:
            amountOfParents = populationAmount
        population = self._initialize_multiple(dim, populationAmount)
        children = []
        while problem.state.evaluations < self.max_iterations and problem.state.current_best.y != 1:
            #empty the fitness array
            fitnessArray = []
            #fill the array with the fitnesses of all strings
            for i in range(populationAmount):
                #fitnessArray.append(problem(population[i]))
                fitnessArray.append(problem(population[i]))
            print(problem.state)
            #find the indices of the strings you want to keep
            mostFitIndices = self.perform_selection(fitnessArray, amountOfParents)
            parents = np.empty((amountOfParents, dim), int)
            #fill the parent array
            for i in range(len(mostFitIndices)):
                for j in range(dim):
                    parents[i,j] = population[mostFitIndices[i]][j]
            #create the children from crossing the parents
            children = self.perform_crossover(parents, populationAmount)
            #mutate the children
            children = self.perform_mutation(children)
            #enter these children into the problem to see if they solve it
            self.f = problem(children[self._rank_selection(fitnessArray,1)[0]])
            #if the problem hasn't been found yet, the new population 
            #will be the children of this generation
            population = children


    
            
def main():
    # Set a random seed in order to get reproducible results
    random.seed(42)
    global rule, t, ct, currentRow, k

    for simCalc in [objective_function, LevenshteinDistance]:
        if simCalc == objective_function:
            name = "objective_function"
        else:
            name = "LevenshteinDistance"
        problem = ioh.problem.wrap_integer_problem(
        simCalc,
        "objective_function_ca_1-{}".format(name),
        60, 
        ioh.OptimizationType.Maximization,
        ioh.IntegerConstraint([0]*60, [2]*60)
        )
        for currentRow in range(10):
            for selectionMethod in ['Rank', 'Tournament']:
                for crossoverMethod in ['Uniform', 'OnePoint']:
                    for mutationMethod in ['Swap', 'Bitstring']:
                        logger = ioh.logger.Analyzer(store_positions=True, algorithm_name="objective_function_ca_1-{}-{}-{}-{}".format(selectionMethod, crossoverMethod, mutationMethod, currentRow))
                        for trialNumber in range(5):
                            print("sm = {}, cm = {}, mm = {}, cr = {}, trialNumber = {}".format(selectionMethod, crossoverMethod, mutationMethod, currentRow, trialNumber))
                            if simCalc == LevenshteinDistance:
                                print("simcalc = Levenshtein")
                            else:
                                print("simcalc = Objective Function") 
                            problem.attach_logger(logger)
                            algorithm = GeneticAlgorithm('Random', 'Random', selectionMethod, crossoverMethod, mutationMethod)
                            # run your algoritm on the problem
                            k = data.iloc[currentRow,0]
                            rule = data.iloc[currentRow,1]
                            t = data.iloc[currentRow,2]
                            ct = data.iloc[currentRow,3]
                            ct = ast.literal_eval(ct)
                            ct = [int(i) for i in ct]
                            algorithm(problem, 60)
                            # Inspect the results
                            print("Best solution found:")
                            print("".join(map(str, problem.state.current_best.x)))
                            print("With an objective value of:", problem.state.current_best.y)
                            print()
                            problem.reset()

def randomMain():
    global rule, t, ct, currentRow, k
    for simCalc in [objective_function, LevenshteinDistance]:
        if simCalc == objective_function:
            name = "objective_function"
        else:
            name = "LevenshteinDistance"
        problem = ioh.problem.wrap_integer_problem(
        simCalc,
        "objective_function_ca_1-{}".format(name),
        60, 
        ioh.OptimizationType.Maximization,
        ioh.IntegerConstraint([0]*60, [2]*60)
        )
        for currentRow in range(10):
            logger = ioh.logger.Analyzer(store_positions=True, algorithm_name="randomSearch-{}-{}".format(name, currentRow))
            for trialNumber in range(5):
                print("problem = {}, trialNumber = {}".format(currentRow, trialNumber))
                if simCalc == LevenshteinDistance:
                    print("simcalc = Levenshtein")
                else:
                    print("simcalc = Objective Function") 
                problem.attach_logger(logger)
                algorithm = RandomSearch(10000)
                # run your algoritm on the problem
                k = data.iloc[currentRow,0]
                rule = data.iloc[currentRow,1]
                t = data.iloc[currentRow,2]
                ct = data.iloc[currentRow,3]
                ct = ast.literal_eval(ct)
                ct = [int(i) for i in ct]
                algorithm(problem, k)
                # Inspect the results
                print("Best solution found:")
                print("".join(map(str, problem.state.current_best.x)))
                print("With an objective value of:", problem.state.current_best.y)
                print()
                problem.reset()
if __name__ == '__main__':
    #randomMain() #main for the random search algorithm
    main()
