import numpy as np
import copy
import json
import sys

import os,inspect
current_dir=os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
os.chdir(current_dir)
import sys
sys.path.append('../')

from lib.SingleNode import SingleLinkList, SingleReaction

class AChromesome(object):
    def __init__(self, ob_sustrate, ob_product, abundant, PoolFile, translator_file,
                NrMax=30, 
                StrictInitialization=False):

        self.ob_sustrate = ob_sustrate
        self.ob_product = ob_product
        self.abundant = abundant
        self.PoolFile = PoolFile
        self.NrMax = NrMax

        self.StrictInitialization = StrictInitialization

        self.sink = dict()
        self.chrom = SingleLinkList()
        self.compounds = set()

        self._translator = self.init_translator(translator_file)

     
    def get_a_chrom(self):

        count = 0

        while(1):

            if count==10000: 
                
                return False

            elif self.get_a_chrom_roulette(): 

                return True

            else: count += 1


    def init_translator(self, translator_file):

        with open(translator_file, 'r', encoding='utf8') as fp:

            json_data = json.load(fp)

        return json_data

    def get_init_sink(self):
        """
        重新初始化 sink 
        :return:
        """
        # print("-----------Init SINK------------")

        compounds = self.abundant + [self.ob_sustrate]

        for c in compounds:

            temp = list(self.PoolFile[c])

            if self.StrictInitialization:

                for i in range(len(temp)):

                    if self._all_in(temp[i].reaction['substrates'], self.abundant):

                        self.sink.setdefault(c, set()).add(temp[i])

            else:

                self.sink[c] = temp

        # print("----------- Done! ------------")

    def update_init_sink(self, Node):
        """
        根据给定的化合物 更新 pool
        :param compound:
        :return:
        """
        if Node.reaction['reversibility'] == 0:

            compounds = Node.reaction['substrates'] + Node.reaction['products']

        if Node.reaction['reversibility'] == 1:

            compounds = Node.reaction['products']

        if Node.reaction['reversibility'] == -1:

            compounds = Node.reaction['substrates']


        for c in compounds:

            c = self._translator[c]

            temp = list(self.PoolFile[c])

            if self.StrictInitialization:

                for i in range(len(temp)):

                    if self._all_in(temp[i].reaction['substrates'], self.abundant):

                        self.sink.setdefault(c, set()).add(temp[i])

            else:

                if c in self.sink: continue

                else:

                    self.sink[c] = temp


    def update_sink_after_crossover_mutation(self):
        """ 
        After crossover and mutation
        Update sink of chrom
        return self.sink, self.compounds
        """

        self.get_init_sink()

        cur = self.chrom._head

        while(cur != None):

            self.update_init_sink(cur)

            cur = cur.next

        self.get_all_compounds()


    def roulette(self, NodeList):
        """轮盘赌算法
        :param NodeList: NodeList 是一个链表
        :return: 返回索引值
        """
        sum_val = 0
        for i in range(len(NodeList)):
            if NodeList[i].Tanimoto != 0 and NodeList[i].Tanimoto != None:
                sum_val = sum_val + NodeList[i].Tanimoto

        random_val = np.random.random()
        probability = 0 #累计概率

        for i in range(len(NodeList)):
            if NodeList[i].Tanimoto != 0 and NodeList[i].Tanimoto != None:
                probability = probability + NodeList[i].Tanimoto / sum_val
            if probability >= random_val:
                    return i
            else: continue

    def get_a_chrom_roulette(self):
        # 按照底物相似性 根据轮盘赌算法 生成个体

        self.get_init_sink()

        index = self.roulette(list(self.sink[self.ob_sustrate]))
        Node = copy.deepcopy(list(self.sink[self.ob_sustrate])[index])
        self.update_init_sink(Node)
        self.chrom.append(Node)

        while(Node.P!=self.ob_product and self.chrom.length() < self.NrMax):
            index = self.roulette(list(self.sink[Node.P]))
            Node = copy.deepcopy(list(self.sink[Node.P])[index])
            self.update_init_sink(Node)
            self.chrom.append(Node)

        if Node.P == self.ob_product:
            # Get a chromesome sucessfully and trim the chrom
            self.trim_chrom()

            return True
        else:
            # "Don't get a chromesome!
            self.chrom.clear()
            self.sink.clear()
            return False

    def trim_chrom(self):
        head = self.chrom._head
        while(head.next != None):
            if(head.S==head.next.P and head.P==head.next.S):
                head = head.next.next
            else:
                break
        self.chrom._head = head

        try:
            # 修剪情况1: a-->b, b-->a
            if(self.chrom!=None and self.chrom.length() > 2): 

                left = head.next
                right = left.next
            
                while(right!=None):
                    pre = head
                    if(left.S==right.P and left.P==right.S):
                        while(pre.next != left):
                            pre = pre.next
                        pre.next = right.next
                        left = pre
                        right = left.next
                    else:
                        left = left.next
                        right = right.next
            
            # 修剪情况2: a-->a
            cur = self.chrom._head.next
            while(cur!=None):
                pre = head
                if(cur.S == cur.P):
                    while(pre.next != cur):
                            pre = pre.next
                    pre.next = cur.next
                cur = cur.next                

        except:
            # print(self.chrom.travel())
            pass
            
        self.trim_sink

    def trim_sink(self):
        """Trim sink after trim self.chrom
        """

        self.get_all_compounds()

        for c in list(set(self.sink.keys())^set(self.compounds)):

                self.sink.pop(c)


    def get_all_compounds(self):

        for c in self.abundant:

            self.compounds.add(c)

        cur = self.chrom._head

        while(cur != None):
            
            for c in cur.reaction['substrates']:

                self.compounds.add(self._translator[c])
            
            for c in cur.reaction['products']:

                self.compounds.add(self._translator[c])

            cur = cur.next


    # def mut_update_a_chrom(self, Sustrate, Product):

    #     mut_chrom = SingleLinkList()
    #     mut_sink = dict()

    #     index = self.roulette(list(self.sink[Sustrate]))
    #     Node = copy.deepcopy(list(self.sink[Sustrate])[index])
    #     self.update_init_sink(Node)
    #     self.chrom.append(Node)

    #     while(Node.P!=self.ob_product and self.chrom.length() < self.NrMax):
    #         index = self.roulette(list(self.sink[Node.P]))
    #         Node = copy.deepcopy(list(self.sink[Node.P])[index])
    #         self.update_init_sink(Node)
    #         self.chrom.append(Node)

    #     if Node.P == self.ob_product:
    #         # Get a chromesome sucessfully and trim the chrom
    #         self.trim_chrom()

    #         return True
    #     else:
    #         # "Don't get a chromesome!
    #         self.chrom.clear()
    #         self.sink.clear()
    #         return False

    # def roulette_get_a_chrom(self, Sustrate, Product):
    #     # 按照底物相似性 根据轮盘赌算法 生成个体
    #     chrom = SingleLinkList()
    #     index = self.roulette(list(self.sink[Sustrate]))
    #     Node = copy.deepcopy(list(self.sink[Sustrate])[index])
    #     self.update_init_sink(Node)
    #     chrom.append(Node)

    #     while(Node.P!=Product and self.chrom.length() < self.NrMax):
    #         index = self.roulette(list(self.sink[Node.P]))
    #         Node = copy.deepcopy(list(self.sink[Node.P])[index])
    #         self.update_init_sink(Node)
    #         chrom.append(Node)

    #     if Node.P == Product:
    #         print("Get a chromesome sucessfully!")
    #         chrom.travel()
    #         return True
    #     else:
    #         print("Don't get a chromesome!")
    #         chrom.clear()
    #         return False

    # def random_get_a_chrom(self, Sustrate, Product):
    #     # 随机生成个体
    #     chrom = SingleLinkList()
    #     index = np.random.randint(len(self.sink[Sustrate]))
    #     Node = copy.deepcopy(list(self.sink[Sustrate])[index])
    #     self.update_init_sink(Node)
    #     chrom.append(Node)

    #     while(Node.P!=Product and self.chrom.length() < self.NrMax):
    #         index = np.random.randint(len(self.sink[Sustrate]))
    #         Node = copy.deepcopy(list(self.sink[Node.P])[index])
    #         self.update_init_sink(Node)
    #         chrom.append(Node)

    #     if Node.P == Product:
    #         print("Get a chromesome sucessfully!")
    #         chrom.travel()
    #         return True
    #     else:
    #         print("Don't get a chromesome!")
    #         chrom.clear()
    #         return False






if __name__ == '__main__':

    file = '/home/caoyh/project/myseeker/db/KEGG_caoyh/MYPOOL/MYPOOL.npy'
    translator_file = '/home/caoyh/project/myseeker/db/KEGG_caoyh/CompDict_rn.json'
    mypool = np.load(file, allow_pickle=True).item()
    abundant= ['C00001', 'C00002', 'C00003', 'C00004', 'C00005', 'C00006', 'C00007', 'C00008', 'C00009', 'C00010', 'C00080']

    ob_sustrate = "C00022"
    ob_product = "C07281"

    NR = 0
    temp_product = ob_sustrate
    cla = AChromesome(ob_sustrate, ob_product,abundant, mypool, translator_file=translator_file)
    cla.a_chrom()
    cla.chrom.travel()



