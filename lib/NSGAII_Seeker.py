import numpy as np
import copy
import json
import random
# random.seed(0)

# fast_non_dominatied_sort
from pymoo.util.nds.non_dominated_sorting import NonDominatedSorting


import os,inspect
from sympy import field_isomorphism, multinomial_coefficients
current_dir=os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
os.chdir(current_dir)
import sys
sys.path.append('../')

from lib.AChromesome import AChromesome
from lib.SingleNode import SingleLinkList, SingleReaction


class NSGAII_Seeker(object):
    def __init__(self, S, P, abundant, POOLFILE, translator_file,
                 G=10, NP=5, Pc=0.8, Pm=0.1) -> None:

        self.S = S
        self.P = P
        self.abundant = abundant
        self.POOLFILE = POOLFILE
        self.translator_file = translator_file

        self.G = G # 种群迭代次数
        self.NP = NP # 种群规模
        self.Pc = Pc # 交叉概率
        self.Pm = Pm # 变异概率
        

        self.chromesomes = []
        self.mut_chrom = None

        self.BestChrome = []
        self.BestFitness = []
        
    def init_chromesomes(self, NP=10):

        chromesomes = []

        count = 0

        while(count < NP):

            chrom = AChromesome(self.S, self.P, self.abundant, self.POOLFILE, 
                                translator_file=self.translator_file)    

            if chrom.get_a_chrom():
                
                chromesomes.append(chrom) # 直接把一个类添加到种群中

                print("-----Get the %dth chrom!----"%(count))
                chrom.chrom.travel()

                count += 1
 
        print("Init Population Successfully!")

        return chromesomes

    def save_chromesomes(self):
        FileDir = './'
        FileName = FileDir + 'ChromeSomes_class.npy'

        np.save(FileName, self.chromesomes)

        print("--------Save ChromeSomes Successfully!---------")

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

    def Crossover(self, chromesomes):

        count = 0

        while(count < 100):

            print(count)

            i, j = random.sample(range(1,len(chromesomes)), 2)

            if self.has_same_node(chromesomes[i], chromesomes[j]):
                print(i, j)
                return self.cross(chromesomes[i], chromesomes[j])

            count = count + 1

    def Mutation(self, chromesomes):

        count = 0

        while(count < 100):

            idx = np.random.randint(len(chromesomes))

            if self.mut(chromesomes[idx]):

                return self.mut_chrom

    def mut(self, Chrome):

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

            mut_chrom = AChromesome(s, p, self.abundant, self.POOLFILE, translator_file=self.translator_file) 

            if mut_chrom.get_a_chrom(): 
                mut_chrom.chrom.travel()
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
                            self.mut_chrom = chrom
                            return True
                        
                        else:
                            return False

                    else:
                        mut_cur = mut_cur.next
            else:
                cur = cur.next 

    def func_1(self, chrom):
        # length of the pathway 
        # less is better, so get the opposite value of pathway, so more is better now
        return -chrom.chrom.length()

    def func_2(self, chrom):
        # abs of Gib
        # More is better

        dGs = 0

        cur = chrom.chrom._head

        while(cur!=None):

            dGs = dGs + abs(cur.reaction['dG_Mean'])

        return dGs
             
    def func_3(self, chrom):
        # yelid of pathway
        # more is better

        return np.random.random()

    def get_fitness(self, Chrom):
        
        F1 = self.func_1(Chrom)
        F2 = self.func_2(Chrom)
        F3 = self.func_3(Chrom)
        F = [F1, F2, F3]
        
        return F

    def Fit(self, chromesomes):

        fitness = []

        if type(chromesomes).__name__=='list':

            for i in range(len(chromesomes)):

                fitness.append(self.get_fitness(chromesomes[i]))
        
        else:

            fitness.append(self.get_fitness(chromesomes))

        return fitness

    def Selection(self, chromesomes, fitness):

        # Tournament

        new_chromesomes = []
        new_fitness = []

        for i in range(self.NP):
            m, n, k = random.sample(range(0,len(chromesomes)), 3)

            max_idx = np.argmax([fitness[m], fitness[n], fitness[k]])

            if max_idx==0: 
                new_chromesomes.append(chromesomes[m])
                new_fitness.append(fitness[m])

            elif max_idx==1: 
                new_chromesomes.append(chromesomes[n])
                new_fitness.append(fitness[n])

            else: 
                new_chromesomes.append(chromesomes[k])
                new_fitness.append(fitness[k])


        return new_chromesomes, new_fitness

    def non_dominated_sorting(self, F):

        Front = NonDominatedSorting().do(F)

        return Front

    def cal_crowing_distance(self, F):

        from pymoo.algorithms.moo.nsga2 import calc_crowding_distance

        crowding = calc_crowding_distance(F=F)
        
        return crowding

    def Evolution(self):

        # Init Population
        pop = self.init_chromesomes(NP=self.NP)

        # Get Init_Pop Fitness
        fitness = self.Fit(pop)

        for i in range(self.G):

            print('-----This is the %d generation-----'%(i))

            # Crossover
            c = np.random.random()

            if c < self.Pc:

                print('-----Crossover-----')
                
                try:

                    crossA, crossB = self.Crossover(pop)
                    cross_fit = self.Fit([crossA, crossB])

                    pop.append(crossA)
                    pop.append(crossB)
                    fitness = fitness + cross_fit

                except:

                    cross_pop = self.init_chromesomes(NP=2)
                    cross_fit = self.Fit(cross_pop)

                    pop = pop + cross_pop
                    fitness = fitness + cross_fit


            # Mutation
            m = np.random.random()

            if m < self.Pm:
                print('-----Mutation-----')

                try:

                    mutA = self.Mutation(pop)
                    mut_fit = self.Fit(mutA)

                    pop.append(mutA)
                    fitness = fitness + mut_fit

                except:
                    mut_pop = self.init_chromesomes(NP=1)
                    mut_fit = self.Fit(mut_pop)

                    pop.append(mut_pop)
                    fitness = fitness + mut_fit

            # Selection
            print('-----Selection-----')
            new_pop, new_fitness = self.Selection(pop, fitness)

            idxes = np.argsort(new_fitness) # 根据适应度值从小到大排序
            sorted_pop = [new_pop[idx] for idx in idxes]
            sorted_fitness = [new_fitness[idx] for idx in idxes]

            # Save
            print('-----Save to best chrom-----')
            self.BestChrome.append(sorted_pop[-1])
            self.BestFitness.append(sorted_fitness[-1])

            print('The Best chrom is:')
            print(self.BestChrome[i].chrom.travel())
            print('The Best Fitness of the best chrom is %f' %(self.BestFitness[i]))

            pop = sorted_pop
            fitness = sorted_fitness
        
        print("-----Finally-----")

            









if __name__ == '__main__':

    file = '/home/caoyh/project/myseeker/db/KEGG_caoyh/MYPOOL/MYPOOL.npy'
    translator_file = '/home/caoyh/project/myseeker/db/KEGG_caoyh/CompDict_rn.json'
    mypool = np.load(file, allow_pickle=True).item() # 提示warning 特别慢
    abundant= ['C00001', 'C00002', 'C00003', 'C00004', 'C00005', 'C00006', 'C00007', 'C00008', 'C00009', 'C00010', 'C00080']

    ob_sustrate = 'C00103' #"C00022"
    ob_product = 'C00631' #"C07281"

    NR = 0

    mychromes = NSGAII_Seeker(S=ob_sustrate, P=ob_product, abundant=abundant, POOLFILE=mypool, translator_file=translator_file)
    mychromes.Evolution()












