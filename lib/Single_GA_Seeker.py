import numpy as np
import copy
import json
import random
# random.seed(0)

import os,inspect
from sympy import field_isomorphism, multinomial_coefficients
current_dir=os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
os.chdir(current_dir)
import sys
sys.path.append('../')

from lib.AChromesome import AChromesome
from lib.SingleNode import SingleLinkList, SingleReaction


class Single_GA_Seeker(object):
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

                # print(chromesomes[idx].chrom.travel())
                # print('Mutation Successfully!')
                # print(self.mut_chrom.travel())

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

    def get_fitness(self, Chrom):

        # return Chrom.chrom.length()
        return np.random.random()

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

    ob_sustrate = "C09871"#'C00103' 
    ob_product =  "C03069" #'C00631' 

    NR = 0

    mychromes = Single_GA_Seeker(S=ob_sustrate, P=ob_product, abundant=abundant, POOLFILE=mypool, translator_file=translator_file)
    mychromes.init_chromesomes(NP=5)

    print('--------------') 
    # mychromes.Mutation()
    print('--------------')


# -----Get the 0th chrom!----
# C09871  -->  C00005 , ['R08083'] : 1 C09871 + 1 C00006 <=> 1 C09847 + 1 C00005 + 1 C00080
# C00005  -->  C00006 , ['R08304'] : 1 C07501 + 1 C00005 + 1 C00080 + 1 C00007 <=> 1 C16582 + 1 C00006 + 1 C00001
# C00006  -->  C00332 , ['R05576'] : 1 C05116 + 1 C00006 <=> 1 C00332 + 1 C00005 + 1 C00080
# C00332  -->  C00356 , ['R01978'] : 1 C00356 + 1 C00010 <=> 1 C00024 + 1 C00001 + 1 C00332
# C00356  -->  C03231 , ['R02085'] : 1 C00356 <=> 1 C03231 + 1 C00001
# C03231  -->  C03069 , ['R04138'] : 1 C00002 + 1 C03069 + 1 C00288 <=> 1 C00008 + 1 C00009 + 1 C03231
# -----Get the 1th chrom!----
# C09871  -->  C02569 , ['R12672'] : 1 C02569 + 1 C00001 <=> 1 C09871 + 1 C00013
# C02569  -->  C00013 , ['R12638'] : 1 C02569 + 1 C00129 <=> 1 C00013 + 1 C19760
# C00013  -->  C00010 , ['R02194'] : 1 C00002 + 1 C01494 + 1 C00010 <=> 1 C00020 + 1 C00013 + 1 C00406
# C00010  -->  C03069 , ['R08095'] : 1 C16471 + 1 C00010 <=> 1 C03069 + 1 C00024
# -----Get the 2th chrom!----
# C09871  -->  C09847 , ['R08083'] : 1 C09871 + 1 C00006 <=> 1 C09847 + 1 C00005 + 1 C00080
# C09847  -->  C00341 , ['R08394'] : 1 C00341 <=> 1 C09847
# C00341  -->  C22152 , ['R12881'] : 1 C00758 + 1 C00341 <=> 1 C22152 + 1 C00013
# C22152  -->  C00758 , ['R12881'] : 1 C00758 + 1 C00341 <=> 1 C22152 + 1 C00013
# C00758  -->  C00021 , ['R02882'] : 1 C00019 + 1 C00758 <=> 1 C00021 + 1 C01557
# C00021  -->  C00019 , ['R07232'] : 1 C00019 + 1 C20648 <=> 1 C00021 + 1 C04154
# C00019  -->  C21591 , ['R11692'] : 1 C21590 + 1 C00019 <=> 1 C21591 + 1 C00021
# C21591  -->  C20297 , ['R12347'] : 1 C21591 + 1 C00001 <=> 1 C20297 + 1 C00033
# C20297  -->  C00004 , ['R10020'] : 1 C20297 + 1 C00003 <=> 1 C09592 + 1 C00004 + 1 C00080
# C00004  -->  C16469 , ['R08094'] : 1 C16469 + 1 C00003 <=> 1 C16471 + 1 C00004 + 1 C00080
# C16469  -->  C16471 , ['R08094'] : 1 C16469 + 1 C00003 <=> 1 C16471 + 1 C00004 + 1 C00080
# C16471  -->  C03069 , ['R08095'] : 1 C16471 + 1 C00010 <=> 1 C03069 + 1 C00024
# -----Get the 3th chrom!----
# C09871  -->  C09847 , ['R08083'] : 1 C09871 + 1 C00006 <=> 1 C09847 + 1 C00005 + 1 C00080
# C09847  -->  C09871 , ['R08083'] : 1 C09871 + 1 C00006 <=> 1 C09847 + 1 C00005 + 1 C00080
# C09871  -->  C00341 , ['R08393'] : 1 C00341 <=> 1 C09871
# C00341  -->  C20235 , ['R09940'] : 1 C00019 + 1 C00341 <=> 1 C00021 + 1 C20235
# C20235  -->  C00341 , ['R09940'] : 1 C00019 + 1 C00341 <=> 1 C00021 + 1 C20235
# C00341  -->  C00013 , ['R11087'] : 1 C00341 + 1 C21001 <=> 1 C00013 + 1 C21003
# C00013  -->  C00448 , ['R09627'] : 1 C00448 <=> 1 C19754 + 1 C00013
# C00448  -->  C00013 , ['R08373'] : 1 C00448 <=> 1 C09684 + 1 C00013
# C00013  -->  C00235 , ['R10365'] : 1 C20554 + 1 C00235 <=> 1 C20555 + 1 C00013
# C00235  -->  C16424 , ['R08051'] : 1 C00235 + 1 C00002 <=> 1 C16424 + 1 C00013
# C16424  -->  C16426 , ['R08061'] : 1 C16424 + 1 C00001 <=> 1 C16426 + 1 C00009
# C16426  -->  C00008 , ['R08052'] : 1 C00235 + 1 C00008 <=> 1 C16426 + 1 C00013
# C00008  -->  C00002 , ['R10953'] : 1 C00002 + 1 C01272 <=> 1 C00008 + 1 C01284
# C00002  -->  C00008 , ['R12853'] : 1 C22442 + 1 C00002 <=> 1 C22443 + 1 C00008
# C00008  -->  C00002 , ['R12293'] : 1 G13091 + 1 C00002 <=> 1 G13092 + 1 C00008
# C00002  -->  C00016 , ['R00161'] : 1 C00002 + 1 C00061 <=> 1 C00013 + 1 C00016
# C00016  -->  C16468 , ['R08092'] : 1 C16470 + 1 C00016 <=> 1 C16468 + 1 C01352
# C16468  -->  C16469 , ['R08093'] : 1 C16468 + 1 C00001 <=> 1 C16469
# C16469  -->  C16471 , ['R08094'] : 1 C16469 + 1 C00003 <=> 1 C16471 + 1 C00004 + 1 C00080
# C16471  -->  C03069 , ['R08095'] : 1 C16471 + 1 C00010 <=> 1 C03069 + 1 C00024
# -----Get the 4th chrom!----
# C09871  -->  C00341 , ['R08393'] : 1 C00341 <=> 1 C09871
# C00341  -->  C03190 , ['R02007'] : 1 C03190 <=> 1 C00341
# C03190  -->  C01765 , ['R08529'] : 1 C03190 + 1 C00001 <=> 1 C01765 + 1 C00013
# C01765  -->  C11338 , ['R08532'] : 1 C00024 + 1 C01765 <=> 1 C00010 + 1 C11338
# C11338  -->  C00024 , ['R08532'] : 1 C00024 + 1 C01765 <=> 1 C00010 + 1 C11338
# C00024  -->  C02232 , ['R00829'] : 1 C00091 + 1 C00024 <=> 1 C00010 + 1 C02232
# C02232  -->  C00091 , ['R02990'] : 1 C00091 + 1 C00846 <=> 1 C00042 + 1 C02232
# C00091  -->  C00356 , ['R02084'] : 1 C00091 + 1 C03761 <=> 1 C00042 + 1 C00356
# C00356  -->  C03231 , ['R02085'] : 1 C00356 <=> 1 C03231 + 1 C00001
# C03231  -->  C03069 , ['R04138'] : 1 C00002 + 1 C03069 + 1 C00288 <=> 1 C00008 + 1 C00009 + 1 C03231









