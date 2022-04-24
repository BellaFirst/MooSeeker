#!usr/bin/env python
# -*- coding:utf-8 -*-
'''
Created on 202101011
@author: think
'''


'''
Created on 20210111
@author: kakamesyli
'''
import random
import operator
import copy
import math
import matplotlib.pyplot as plt

class Pop(object):

    def __init__(self, pop_index, x_min, x_max, x_num, f_num, fun_name):
        self.index = pop_index;  # #code
        self.var = [];  #
        for i in range(x_num):
            self.cal_var(x_min, x_max)
        self.value = []  #
        #for j in range(f_num):
        self.cal_output(self.var, x_num, fun_name)
        self.pareto_rank = 0
        self.crowding = 0
        self.dom_num = 0
        self.dom_asm = []
        self.x_max = x_max
        self.x_min = x_min

    def cal_var(self, x_min, x_max):
        self.var.append(x_min + (x_max - x_min) * random.random())
    
    def cal_output(self, x, x_num, fun_name):
        self.value[:] = obj_fun(x, x_num, fun_name)
        #self.value.append(obj_fun(x, x_num, fun_name));
    
    def create_rank_asm(self, Pop):
        self.rank_asm.append(Pop);
    
    def refresh(self):
        self.pareto_rank = 0
        self.crowding = 0
        self.dom_num = 0
        self.dom_asm = []

        
class Chromo(object):

    def __init__(self, pop_num, x_min, x_max, x_num, f_num, fun_name):
        self.pop_num = pop_num;
        self.pop_asm = [];
        self.rank_asm = []; 
        for i in range(self.pop_num):
            self.create_pop_asm(i, x_min, x_max, x_num, f_num, fun_name);
    
    def create_pop_asm(self, pop_index, x_min, x_max, x_num, f_num, fun_name):
        self.pop_asm.append(Pop(pop_index, x_min, x_max, x_num, f_num, fun_name));
        
    def non_dominate_sort(self, pop_asm, f_num):
        pareto_rank = 1;
        for pop_i in pop_asm:
            for pop_j in pop_asm:
                less = 0;
                equal = 0;
                greater = 0;
                for k in range(f_num):
                    if (pop_i.value[k] < pop_j.value[k]):
                        less = less + 1;
                    elif (pop_i.value[k] == pop_j.value[k]):
                        equal = equal + 1;
                    else:
                        greater = greater + 1;
                if (less == 0 and equal != f_num):
                    pop_i.dom_num = pop_i.dom_num + 1;
                elif (greater == 0 and equal != f_num):
                    pop_i.dom_asm.append(pop_j);
                if pop_i.dom_num == 0:
                    pop_i.pareto_rank = pareto_rank;
                    self.rank_asm[pareto_rank].append(pop_i);
        while (len(self.rank_asm[pareto_rank]) != 0):
            for rank_ele in self.rank_asm[pareto_rank]:
                if (len(rank_ele.dom_asm) != 0):
                    for dom_asm_ele in rank_ele.dom_asm:
                        dom_asm_ele.dom_num = dom_asm_ele.dom_num - 1;
                        if dom_asm_ele.dom_num == 0:
                            dom_asm_ele.pareto_rank = pareto_rank + 1;
                            self.rank_asm[pareto_rank + 1].append(dom_asm_ele);
            pareto_rank = pareto_rank + 1;
            
    def crowding_distance_sort(self, pop_asm, f_num):
        # temp = pop_asm.sorted(key=lambda x:x.pareto_rank, reverse=False);
        for pareto_rank in range(len(self.rank_asm)):
            for func in range(f_num):
                func_temp = []
                func_temp[pareto_rank] = sorted(self.rank_asm[pareto_rank], key=lambda x:x.value[func], reverse=False)
                current_index = 0;
                fmin = func_temp[pareto_rank][0].value[func]
                fmax = func_temp[pareto_rank][-1].value[func]
                func_temp[pareto_rank][0] = float('Inf')
                func_temp[pareto_rank][-1] = float('Inf')
                for i in range(1, len(func_temp[pareto_rank]) - 2):  # for pop in func_temp[pareto_rank][1, -2]:
                    next_fun_value = func_temp[pareto_rank][i - 1].value[func]  # next_fun_value = func_temp[pareto_rank][pop.index()-1].value[func]
                    pre_fun_value = func_temp[pareto_rank][i + 1].value[func]  # pre_fun_value = func_temp[pareto_rank][pop.index()+1].value[func]
                    if fmax == fmin:
                        func_temp[pareto_rank][i].crowding = float('Inf')
                    else:
                        func_temp[pareto_rank][i].crowding = func_temp[pareto_rank][i].crowding + (next_fun_value - pre_fun_value) / (fmax - fmin)
            func_temp[pareto_rank] = sorted(func_temp[pareto_rank], key=lambda x:x.index, reverse=False)
            for i in range(len(func_temp[pareto_rank])):
                self.rank_asm[pareto_rank][i].crowding = func_temp[pareto_rank][i].crowding
                
                '''
            for pop in temp[current_index, len(self.rank_asm[pareto_rank]) - 1]:
                for func in range(f_num):
                    fun_value_temp[func] = 
                    fun_value_temp[func].append(pop.value[func])
            
                y_asm.append(pop)
            for f in range(f_num):
                f_sort[f] !! = y_asm.sorted(key=lambda x:x.value[f], reverse=False)
                fmin = f_sort[1];
                fmax = f_sort[-1]
                f_sort[1].crowding = float('inf')
                f_sort[-1].crowding = float('inf')
                for f_s in range(2, len(f_sort) - 1):
                    pre_f = f_sort[f_s - 1]
                    next_f = f_sort[f_s + 1]
                    if fmin == fmax:
                        f_sort[f_s] = float('inf')
                    else:
                        f_sort[f_s] = (next_f - pre_f) / (fmax - fmin)
            for f in range(f_num):
                nd = 
'''


def non_dominate_sort(pop_asm, f_num):
    rank_asm = []
    rank_asm_temp = []
    pareto_rank = 1
    for pop_i in pop_asm:
        pop_i.dom_asm = []
        pop_i.dom_num = 0
        for pop_j in pop_asm:
            #dom_flag = 0
            less = 0;
            equal = 0;
            greater = 0;
            for k in range(f_num):
                if (pop_i.value[k] < pop_j.value[k]):
                    less = less + 1
                elif (pop_i.value[k] == pop_j.value[k]):
                    equal = equal + 1
                else:
                    greater = greater + 1
                    
            if (less == 0 and equal != f_num):
                pop_i.dom_num = pop_i.dom_num + 1
            elif (greater == 0 and equal != f_num):
                pop_i.dom_asm.append(pop_j)

        if pop_i.dom_num == 0:
            pop_i.pareto_rank = pareto_rank
            rank_asm_temp.append(pop_i)
    #rank_asm.append(copy.deepcopy(rank_asm_temp, None, []))
    rank_asm.append(rank_asm_temp)
            #rank_asm.append(copy.deepcopy(pop_i, None, []))
        
    while (len(rank_asm[pareto_rank-1]) != 0):
        rank_asm_temp = []
        for rank_ele in rank_asm[pareto_rank-1]:
            if (len(rank_ele.dom_asm) != 0):
                for dom_asm_ele in rank_ele.dom_asm:
                    dom_asm_ele.dom_num = dom_asm_ele.dom_num - 1
                    if dom_asm_ele.dom_num == 0:
                        dom_asm_ele.pareto_rank = pareto_rank+1
                        rank_asm_temp.append(dom_asm_ele)
        pareto_rank += 1
        if len(rank_asm_temp) == 0:
            break
        #rank_asm.append(copy.deepcopy(rank_asm_temp,None,[]))
        rank_asm.append(rank_asm_temp)
    return rank_asm
def crowding_distance_sort(rank_asm, f_num):
    # temp = pop_asm.sorted(key=lambda x:x.pareto_rank, reverse=False);
    for pareto_rank in range(len(rank_asm)):
        for func in range(f_num):
            #func_temp = copy.deepcopy(rank_asm) #dont change origin
            func_temp = rank_asm
            #func_temp.append(sorted(rank_asm[pareto_rank], key=lambda x:x.value[func], reverse=False))
            func_temp[pareto_rank] = sorted(func_temp[pareto_rank], key=lambda x:x.value[func], reverse=False)
            current_index = 0
            fmin = func_temp[pareto_rank][0].value[func]
            fmax = func_temp[pareto_rank][-1].value[func]
            func_temp[pareto_rank][0].crowding = float('Inf')
            func_temp[pareto_rank][-1].crowding = float('Inf')
            for i in range(1, len(func_temp[pareto_rank]) - 1):  # for pop in func_temp[pareto_rank][1, -2]:
                next_fun_value = func_temp[pareto_rank][i - 1].value[func]  # next_fun_value = func_temp[pareto_rank][pop.index()-1].value[func]
                pre_fun_value = func_temp[pareto_rank][i + 1].value[func]  # pre_fun_value = func_temp[pareto_rank][pop.index()+1].value[func]
                if fmax == fmin:
                    func_temp[pareto_rank][i].crowding = float('Inf')
                else:
                    func_temp[pareto_rank][i].crowding = func_temp[pareto_rank][i].crowding + (next_fun_value - pre_fun_value) / (fmax - fmin)
        func_temp[pareto_rank] = sorted(func_temp[pareto_rank], key=lambda x:x.index, reverse=False)
        for i in range(len(func_temp[pareto_rank])):
            #rank_asm_crowding[pareto_rank][i].crowding = func_temp[pareto_rank][i].crowding
            rank_asm[pareto_rank][i].crowding = func_temp[pareto_rank][i].crowding
    return rank_asm
        
def elitsm(pop_num, chromo2):
    chromo_temp = copy.deepcopy(chromo2, None, [])
    pre_ind = -1#match the length and index
    current_rank = 0
    current_ind = 0
    ind = []
    current_rank_temp = []
    chromo_elit = []
    
    current_ind = len(chromo_temp[current_rank]) - 1
    
    while current_ind < pop_num:
        #chromo_elit.append(chromo_temp[current_rank][:]) #create list in list
        chromo_elit += chromo_temp[current_rank]
        current_rank += 1
        pre_ind = current_ind
        current_ind = current_ind + len(chromo_temp[current_rank])
    
    elit_temp = sorted(chromo_temp[current_rank],key = lambda x:x.crowding,reverse = False)
    #chromo_elit.append(elit_temp[0:(pop_num-1) - pre_ind]) #create list in list
    chromo_elit += elit_temp[0:(pop_num-1)-pre_ind]
    '''
    for pareto_rank in range(len(chromo_temp)):
        current_ind = current_ind + len(chromo_temp[pareto_rank])
    for pop in chromo_temp:       
        if pop.pareto_rank == current_rank:
            current_ind = current_ind + 1
        elif pop.pareto_rank > current_rank:
            ind[current_rank] = current_ind
            current_rank = current_rank + 1
            current_ind = current_ind + 1
        ind[current_rank] = current_ind
            
    for i in range(len(ind)):
        if ind[i] < pop_num or ind[i] == pop_num:
            chromo_elit.append(chromo_temp[i - 1])
        elif ind[i] > pop_num:
            for j in range(ind[i - 1], pop_num - ind[i - 1]):
                chromo_elit.append(chromo_temp[j - 1])
           
           
           
           
            current_rank_temp[current_rank].append(pop)
            if current_ind < pop_num and current_ind == pop_num:
                chromo_elit.append(current_rank_temp[current_rank])
            elif current_ind > pop_num:
                crowding_distance_sort(current_rank_temp[current_rank])
        elif pop.pareto_rank > current_rank:
            current_rank = current_rank + 1
            current_ind = current_ind + 1
            current_rank_temp[current_rank].append(pop)
        else:
            print('error')
        if pop.pareto_rank == current_rank:
            if current_ind < pop_num:
                chromo_elit.append(pop)
                current_ind = current_ind + 1
            elif current_ind > pop_num and current_ind == pop_num:
                
        elif pop.pareto_rank > current_rank:
            
        
            
            chromo_temp = chromo_temp + 1
        elif current_ind > pop_num and current_ind == pop_num:
            
    for i in range(1,max(chromo_temp.pareto_rank)):
        current_index_pareto = 
        '''
    return chromo_elit


def tournament_selection(chromo):
    tournament = 2
    #k = round(len(chromo)/2)
    chromo_len = len(chromo)
    tournament_index_temp = [None for x in range(tournament)]
    tournament_chromo = []
    chromo_temp = []
    chromo_rank_temp = []
    for i in range(chromo_len):
        for j in range(tournament):
            tournament_index_temp[j] = int(round((chromo_len-1) * random.random()))
            if j>0:
                while tournament_index_temp[j] == tournament_index_temp[j-1]:
                    tournament_index_temp[j] = int(round((chromo_len-1) * random.random()))
            chromo_temp.append(copy.deepcopy(chromo[tournament_index_temp[j]]))
        
        min_rank_chromo = min(chromo_temp,key = lambda x:x.pareto_rank)
        min_rank = min_rank_chromo.pareto_rank
        for k in range(tournament):
            if chromo_temp[k].pareto_rank == min_rank:
                chromo_rank_temp.append(chromo_temp[k])
        del chromo_temp[:]
        if len(chromo_rank_temp) == 1:
            tournament_chromo.append(chromo_rank_temp[0])
        else:
            max_crowding_index = chromo_rank_temp.index(max(chromo_rank_temp,key = lambda x:x.crowding))
            tournament_chromo.append(chromo_rank_temp[max_crowding_index])
        del chromo_rank_temp[:]
    return tournament_chromo

def cross_mutation(chromo,x_num,x_max,x_min):
    pc = 0.9
    pm = 1/30
    n = 1
    fun_name = 'ZDT1'
    chromo_len = len(chromo)
    x1_c_chromo = []
    x2_c_chromo = []
    v1_c_chromo = []
    v2_c_chromo = []
    chromo_cross_mutation = copy.deepcopy(chromo, None, [])
    for i in range(chromo_len//2):
        x1_index = int(round((chromo_len-1) * random.random()))
        x2_index = int(round((chromo_len-1) * random.random()))
        while x1_index == x2_index:
            x2_index = int(round((chromo_len-1) * random.random()))
        #x1_f_chromo = copy.deepcopy(chromo[x1_index], None, [])
        #x2_f_chromo = copy.deepcopy(chromo[x2_index], None, [])
        
        x1_f_chromo = chromo[x1_index].var
        x2_f_chromo = chromo[x2_index].var
        v1_c_chromo = chromo[x1_index].value
        v2_c_chromo = chromo[x2_index].value
        if random.random()<pc:
            for j in range(x_num):
                u1 = random.random()#u1=[0,1)
                if u1 <= 0.5:
                    beta = (2*u1) ** (1/n+1)
                else:
                    beta = (1/(2-2*u1)) ** (1/n+1)
                x1_c_chromo.append(0.5*(x2_f_chromo[j]+x1_f_chromo[j]) - 0.5*beta*(x2_f_chromo[j]-x1_f_chromo[j]))
                if x1_c_chromo[-1] > x_max:
                    x1_c_chromo[-1] = x_max
                elif x1_c_chromo[-1] < x_min:
                    x1_c_chromo[-1] = x_min
                x2_c_chromo.append(0.5*(x2_f_chromo[j]+x1_f_chromo[j]) + 0.5*beta*(x2_f_chromo[j]-x1_f_chromo[j]))
                if x2_c_chromo[-1] > x_max:
                    x2_c_chromo[-1] = x_max
                elif x2_c_chromo[-1] < x_min:
                    x2_c_chromo[-1] = x_min
            v1_c_chromo = obj_fun(x1_c_chromo, x_num, fun_name)
            v2_c_chromo = obj_fun(x2_c_chromo, x_num, fun_name)
        if random.random()<pm:
            for j in range(x_num):
                u2 = random.random()
                if u2 <= 0.5:
                    deta = (2*u2) ** (1/(n+1))
                else:
                    deta = (1-(2-2*u2)) ** (1/(n+1))
                x1_c_chromo[j] = x1_c_chromo[j] + deta
                if x1_c_chromo[j] > x_max:
                    x1_c_chromo[j] = x_max
                elif x1_c_chromo[j] < x_min:
                    x1_c_chromo[j] = x_min
            v1_c_chromo = obj_fun(x1_c_chromo, x_num, fun_name)
        if random.random()<pm:
            for j in range(x_num):
                u2 = random.random()
                if u2 <= 0.5:
                    deta = (2*u2) ** (1/(n+1))
                else:
                    deta = (1-(2-2*u2)) ** (1/(n+1))
                x2_c_chromo[j] = x2_c_chromo[j] + deta
                if x2_c_chromo[j] > x_max:
                    x2_c_chromo[j] = x_max
                elif x2_c_chromo[j] < x_min:
                    x2_c_chromo[j] = x_min
            v2_c_chromo = obj_fun(x2_c_chromo, x_num, fun_name)
        chromo_cross_mutation[x1_index].value[:] = v1_c_chromo
        chromo_cross_mutation[x2_index].value[:] = v2_c_chromo
    return chromo_cross_mutation

def obj_fun(x, x_num, fun_name):
    if operator.eq(fun_name, 'ZDT1'):
        f = []
        f.append(x[0])
        s = 0
        for i in range(1,len(x)):
            s = s + x[i]
        g = 1 + 9 * (s / (x_num - 1))
        f.append(g * (1 - (f[0] / g) ** 0.5))
        return f

    
class Person():

    def __init__(self, age, name):
        self.age = age
        self.name = name


def personSort():
    persons = [Person(age, name) for (age, name) in [(12, "lili"), (18, "lulu"), (16, "kaka"), (12, "xixi")]]
    persons.sort(key=lambda x:x.name, reverse=False)
    for element in persons:
        print (element.age, ":", element.name)


if __name__ == "__main__":
    x_min = 0
    x_max = 1
    x_num = 200
    f_num = 2
    fun_name = 'ZDT1'
    pop_num = 100
    gen = 10
    value_output = []
    
    chromo = Chromo(pop_num, x_min, x_max, x_num, f_num, fun_name)
    chromo_domi = non_dominate_sort(chromo.pop_asm, f_num)
    chromo_domi_crowding = crowding_distance_sort(chromo_domi, f_num)
    chromo_select = tournament_selection(chromo.pop_asm)
    chromo_result = cross_mutation(chromo_select, x_num, x_max, x_min)
    for i in range(gen):
        chromo_parent = tournament_selection(chromo_result)
        chromo_offspring = cross_mutation(chromo_parent, x_num, x_max, x_min)
        chromo_combine = chromo_parent + chromo_offspring
        chromo1_combine = non_dominate_sort(chromo_combine, f_num)
        chromo2_combine = crowding_distance_sort(chromo1_combine, f_num)
        chromo_result = elitsm(pop_num, chromo2_combine)
        if i % 10 == 0:
            print('--%d generation has completed!' % i)
    for f in range(f_num):
        value_output.append([chromo_result[i].value[f] for i in range(len(chromo_result))])
    plt.scatter(value_output[0], value_output[1])
    plt.show()