from NSGAII import MyCrossover,MyMutation, MyProblem, MySampling, MyDuplicateElimination
from lib.SingleNode import SingleReaction, SingleLinkList

import string
import numpy as np

from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.optimize import minimize


algorithm = NSGA2(pop_size=6,
                  sampling=MySampling(),
                  crossover=MyCrossover(),
                  mutation=MyMutation(),
                  eliminate_duplicates=MyDuplicateElimination())



file = '/home/caoyh/project/myseeker/db/KEGG_caoyh/MYPOOL/MYPOOL.npy'
translator_file = '/home/caoyh/project/myseeker/db/KEGG_caoyh/CompDict_rn.json'
mypool = np.load(file, allow_pickle=True).item() # 提示warning 特别慢
abundant= ['C00001', 'C00002', 'C00003', 'C00004', 'C00005', 'C00006', 'C00007', 'C00008', 'C00009', 'C00010', 'C00080']

ob_sustrate = 'C00103' #"C00022"
ob_product = 'C00631' #"C07281"


res = minimize(MyProblem(S=ob_sustrate, P=ob_product, abundant=abundant, pool_file=mypool, translator_file=translator_file),
               algorithm,
               ('n_gen', 10),
               seed=1,
               verbose=False)
print(res.F)



