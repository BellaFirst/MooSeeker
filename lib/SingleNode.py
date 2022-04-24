# reaction = {
# 'name': 'alpha,alpha-trehalose glucohydrolase', 
# 'identifiers': {'KEGG': ['R00010']}, 
# 'enzymes': ['3.2.1.28'], 
# 'equation': 'alpha,alpha-Trehalose + H2O <=> 2 D-Glucose', 
# 'trans_equation': '1 C01083 + 1 C00001 <=> 2 C00031', '
# substrates': ['alpha,alpha-Trehalose', 'H2O'], 
# 'products': ['D-Glucose'], 
# 'Scoef': [1.0, 1.0], 
# 'Pcoef': [2.0], 
# 'reversibility': 0, 
# 'use': True}

# class SingleNode(object):
#     def __init__(self, Rname=None, Substrate=None, Product=None,
#                  Scoefficient=None, Pcoefficient=None,Tanimoto=0.0,
#                  Enzyme=None, next=None):
#         self.Rname = Rname
#         self.Substrate = Substrate
#         self.Product = Product
#         self.Scoefficients = Scoefficient
#         self.Pcoefficients = Pcoefficient
#         self.Tanimoto = Tanimoto
#         self.Enzyme = Enzyme
#         self.next = next


import copy


class SingleReaction(object):
    def __init__(self, S, P, reaction, Tanimoto=0.0, next=None):
        self.S = S
        self.P = P
        self.reaction = reaction
        self.Tanimoto = Tanimoto
        self.next = next
    # def __eq__(self, other) -> bool:
    #     return self.S == other.S and self.P == other.P


class SingleLinkList(object):
    def __init__(self):
        self._head = None

    def is_empty(self):
        return self._head == None

    def clear(self):
        self._head = None


    def travel(self):
        """遍历整个链表"""
        cur = self._head
        while cur != None:
            print(cur.S, " --> ", cur.P, ",",
                  cur.reaction['identifiers']['KEGG'], ":", cur.reaction['trans_equation'])
            # print(cur.item, end=" ")
            cur = cur.next

    def length(self):
        count = 0
        cur = self._head
        while cur != None:
            count += 1
            cur = cur.next
        return count

    def add(self, node):
        """在头部添加节点，节点是已知，更新头节点"""
        # node = SingleNode()
        # node = Node(item)
        node.next = self._head
        self._head = node

    def get_first_node(self):
        return self._head

    def append(self, node):
        """在尾部添加节点"""
        if self.is_empty():
            self._head = node
        else:
            cur = self._head
            # node = Node(item)
            while cur.next != None:
                cur = cur.next
            cur.next = node

    def insert(self, pos, node):
        """在选定的位置添加节点"""
        cur = self._head
        # node = Node(item)
        count = 0
        if pos <= 0:
            self.add(node)
        elif pos > (self.length() - 1):
            self.append(node)
        else:
            while count < (pos - 1):
                count += 1
                cur = cur.next
            node.next = cur.next
            cur.next = node

    # def exchange(self, chromosomeA, chromosomeB):
    #     a, b = chromosomeA, chromosomeB
    #     # Node_temp = SingleNode()
    #     Node_temp = a.next
    #     a.next = b.next
    #     b.next = Node_temp
    #     return a, b



    # def getInterSection(self, chromosomeA, chromosomeB):
    #     """
    #     交叉：染色体A 和 染色体B 在某个相同的位置进行交叉，如果有相同的位置则交叉，否则返回None
    #     :param chromosomeA:
    #     :param chromosomeB:
    #     :return: 返回交叉后的 染色体A 和 染色体B
    #     """
    #     a, b = chromosomeA, chromosomeB
    #     for i in range(self.length(chromosomeA)):
    #         for j in range(self.length(chromosomeB)):
    #             if a.Product == b.Product:
    #                 new_a, new_b = self.exchange(a, b)
    #                 return new_a, new_b
    #     return None

    # def getMutation(self, chromosome):
    #     """
    #     变异：
    #     染色体chromosome在某一个随机位置进行变异,
    #     利用 get_a_chromosome(ob_sustrate, ob_product)在变异位置进行交叉
    #     交叉成功后返回新的染色体
    #     :param chromosome:
    #     :return: new_chromosome
    #     """
    #     pass




