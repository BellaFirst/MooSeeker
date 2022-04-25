from pdb import pm
import string
import numpy as np
import copy
import random

from pymoo.core.problem import ElementwiseProblem
from pymoo.core.sampling import Sampling
from pymoo.core.crossover import Crossover
from pymoo.core.mutation import Mutation
from pymoo.core.duplicate import ElementwiseDuplicateElimination
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.optimize import minimize
from pymoo.visualization.scatter import Scatter

import os,inspect
from sympy import MutableMatrix, field_isomorphism, multinomial_coefficients
current_dir=os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
os.chdir(current_dir)
import sys
sys.path.append('../')

from lib.AChromesome import AChromesome
from lib.SingleNode import SingleLinkList, SingleReaction

class MyProblem(ElementwiseProblem):

    def __init__(self, S, P, abundant, pool_file, translator_file,):
        super().__init__(n_var=1, n_obj=3, n_constr=0)
        self.S = S
        self.P = P
        self.abundant = abundant
        self.pool_file = pool_file
        self.translator_file = translator_file
        
    def _evaluate(self, x, out, *args, **kwargs):
        F1 = self.func_1(x)
        F2 = self.func_2(x)
        F3 = self.func_3(x)

        out["F"] = np.array([F1, F2, F3], dtype=float)

    def func_1(self, chrom):
        # length of the pathway 
        # less is better, so get the opposite value of pathway, so more is better now
        return chrom.chrom.length()

    def func_2(self, chrom):
        # abs of Gib
        # More is better

        dGs = 0

        cur = chrom.chrom._head

        while(cur!=None):

            dGs = dGs + abs(cur.reaction['dG_Mean'])

        return - dGs
             
    def func_3(self, chrom):
        # yelid of pathway
        # more is better

        return - np.random.random()


class MySampling(Sampling):

    def _do(self, problem, n_samples, **kwargs):

        X = []

        count = 0

        while(count < n_samples):

            chrom = AChromesome(problem.S, problem.P, problem.abundant, 
                                problem.pool_file, problem.translator_file)

            if chrom.get_a_chrom():

                X.append(chrom)

                count = count + 1

        return X


class MyCrossover(Crossover):

    def __init__(self, Pc=0.8):

        # define the crossover: number of parents and number of offsprings
        super().__init__(2, 2)
        self.Pc = Pc

    def _do(self, problem, X, **kwargs):

        # The input of has the following shape (n_parents, n_matings, n_var)
        # _, n_matings, n_var = X.shape

        c = np.random.random()

        if c < self.Pc:

            try:

                crossA, crossB = self.MyCro(X)

                return[crossA, crossB]

            except:

                crossA = AChromesome(problem.S, problem.P, problem.abundant, 
                                problem.pool_file, problem.translator_file)

                crossB = AChromesome(problem.S, problem.P, problem.abundant, 
                                problem.pool_file, problem.translator_file)

                return[crossA, crossB]

        else:

            return X

    def MyCro(self, chromesomes):

        count = 0

        while(count < 100):

            print(count)

            i, j = random.sample(range(1,len(chromesomes)), 2)

            if self.has_same_node(chromesomes[i], chromesomes[j]):
                print(i, j)
                return self.cross(chromesomes[i], chromesomes[j])

            count = count + 1


    def cross(self, chromA, chromB):
        """cross chromA and chromB

        Args:
            chromA (class AChromsomes): 
            chromB (class AChromsomes): 

        Returns:
            newA (class AChromsomes): 
            newB (class AChromsomes): 
        """

        newA, newB = copy.deepcopy(chromA), copy.deepcopy(chromB)
        # newA, newB = chromA, chromB

        curA = newA.chrom._head

        while (curA!= None and curA.next!=None):

            curB = newB.chrom._head

            while (curB!= None and curB.next!=None): 

                if(curA.P==curB.P):
                    # find the same point and crossover
                    temp = curA.next
                    curA.next = curB.next
                    curB.next = temp

                    newA.update_sink_after_crossover_mutation()
                    newB.update_sink_after_crossover_mutation()

                    return  newA, newB

                else:
                    curB = curB.next

            curA = curA.next

    def has_same_node(self, chromA, chromB):
        """Judge chromA and chromB has the same product to crossover

        Args:
            chromA (class AChromsomes): _description_
            chromB (class AChromsomes): _description_

        Returns:
            bool: true or false
        """

        curA = chromA.chrom._head

        # 遍历到链表的前一个节点 除掉最后一个节点 
        while (curA != None and curA.next!=None):
            curB = chromB.chrom._head
            while (curB != None and curB.next!=None): 
                if(curA.P==curB.P):
                    # 如果都是头节点 则不满足有相同节点的条件
                    if(curA==chromA.chrom._head and curB==chromB.chrom._head):
                        curB = curB.next
                        continue

                    # 如果都是倒数第二节点，也不满足有相同节点的条件
                    elif(curA.next.next==None and curB.next.next==None): 
                        curB = curB.next
                        continue

                    else:
                        return True
                else:
                    curB = curB.next
            curA = curA.next
        
        return False



class MyMutation(Mutation):

    def __init__(self, Pm=0.1):
        super.__init__()
        self.Pm = Pm

    def _do(self, problem, X, **kwargs):

        m = np.random.random()

        if m < self.Pm:

            try:

                Mut_X = self.Muta(X, problem)

                return Mut_X

            except:

                mut_x = AChromesome(problem.s, problem.p, problem.abundant, 
                                        problem.pool_file, problem.translator_file)

                idx = np.random.randint(low=0, high=len(X))

                X[idx] = mut_x
                
                return X


    def Muta(self, chromesomes, problem):

        count = 0

        while(count < 100):

            idx = np.random.randint(len(chromesomes))

            cond, mut_cross = self.mut(chromesomes[idx], problem)

            if cond:

                chromesomes[idx] = mut_cross

                return chromesomes

    def mut(self, Chrome, problem):

        chrom = copy.deepcopy(Chrome)

        # 确定了底物和产物 
        m, n = random.sample(range(0, chrom.chrom.length()), 2)

        cur = chrom.chrom._head

        for i in range(0, chrom.chrom.length()):
            
            if i==min(m, n): 
                s = cur.P
                cur = cur.next
            elif i==max(m, n): 
                p = cur.P
                break
            else:
                cur = cur.next

        # generate mut_chrom
        while(1):

            mut_chrom = AChromesome(s, p, problem.abundant, problem.pool_file, problem.translator_file) 

            if mut_chrom.get_a_chrom(): 
                # mut_chrom.chrom.travel()
                break

        # connect chrom and mut_chrom
        cur = chrom.chrom._head 
        mut_cur = mut_chrom.chrom._head

        for i in range(1, max(m,n)+1):
            
            if i==min(m, n): 
                
                temp = cur.next
                cur.next = mut_cur

                while(mut_cur != None):
                    if mut_cur.next==None:
                        mut_cur.next = temp
                        # trim new chrom
                        chrom.trim_chrom()
                        # update chrom's sink
                        chrom.update_sink_after_crossover_mutation()

                        if chrom.chrom != Chrome.chrom:
                            return [True, chrom]
                        
                        else:
                            return [False, chrom]

                    else:
                        mut_cur = mut_cur.next
            else:
                cur = cur.next       


class MyDuplicateElimination(ElementwiseDuplicateElimination):

    def is_equal(self, a, b):

        return a.chrom==b.chrom




if __name__=='__main__':

    algorithm = NSGA2(pop_size=6,
                      sampling=MySampling(),
                      crossover=MyCrossover(),
                      mutation=MyMutation(),
                      eliminate_duplicates=MyDuplicateElimination())

    res = minimize(MyProblem(),
               algorithm,
               ('n_gen', 10),
               seed=1,
               verbose=False)

#%%

#%%
from pymoo.algorithms.moo.nsga2 import NSGA2

algorithm = NSGA2(pop_size=6,
                      sampling=MySampling(),
                      crossover=MyCrossover(),
                      mutation=MyMutation(),
                      eliminate_duplicates=MyDuplicateElimination())

# %%
