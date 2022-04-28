# from NSGAII import MyCrossover,MyMutation, MyProblem, MySampling, MyDuplicateElimination
# from lib.SingleNode import SingleReaction, SingleLinkList

# import string
# import numpy as np

# from pymoo.algorithms.moo.nsga2 import NSGA2
# from pymoo.optimize import minimize


# algorithm = NSGA2(pop_size=6,
#                   sampling=MySampling(),
#                   crossover=MyCrossover(),
#                   mutation=MyMutation(),
#                   eliminate_duplicates=MyDuplicateElimination())



# file = '/home/caoyh/project/myseeker/db/KEGG_caoyh/MYPOOL/MYPOOL.npy'
# translator_file = '/home/caoyh/project/myseeker/db/KEGG_caoyh/CompDict_rn.json'
# mypool = np.load(file, allow_pickle=True).item() # 提示warning 特别慢
# abundant= ['C00001', 'C00002', 'C00003', 'C00004', 'C00005', 'C00006', 'C00007', 'C00008', 'C00009', 'C00010', 'C00080']

# ob_sustrate = 'C00103' #"C00022"
# ob_product = 'C00631' #"C07281"


# res = minimize(MyProblem(S=ob_sustrate, P=ob_product, abundant=abundant, pool_file=mypool, translator_file=translator_file),
#                algorithm,
#                ('n_gen', 10),
#                seed=1,
#                verbose=False)
# print(res.F)


#%%
from IPython.display import display, display_html
from cobramod import __version__
from cobra import __version__ as cobra_version
print(f'CobraMod version: {__version__}')
print(f'COBRApy version: {cobra_version}')
# From Escher:
# This option turns off the warning message if you leave or refresh this page
import escher
escher.rc['never_ask_before_quit'] = True
#%%
from pathlib import Path

from cobramod import add_pathway, add_reactions
from cobra.io import write_sbml_model, read_sbml_model

# %%
# Defining system path for data
dir_data = Path.cwd().resolve().joinpath("Theoretical Yields/data")
# %%
# Loading model of iJO1366
original = read_sbml_model("Theoretical Yields/data/iJO1366.xml")
model = original.copy()
model
# %%
wt_growth = model.optimize() 
wt_growth.objective_value  
# %%
min_growth=0.8*wt_growth.objective_value
min_growth
# %%
print(f'Number of reaction prior addition: {len(model.reactions)}')
# %%
reactions = ['R12672', 
             'R12638', 
             'R02194', 
             'R08095']

from cobramod.core.creation import _reaction_from_string

# for reaction in reacrions:

#     try:

#         _reaction_from_string()

# %%
reactions = ['R12672, c', 
             'R12638_c, test | 1 C02569_c + 1 C00129_c <=> 1 C00013_c + 1 C19760_c', 
             'R02194, c', 
             'R08095, c']

add_pathway(
    model=model,
    pathway=reactions,
    directory=dir_data,
    database="KEGG",
    compartment="c",
    genome="ecc",
    group="Propanoyl-CoA-PWY",
    filename="Propanoyl-CoA-PWY.txt",
    model_id="iJO1366"
)
# %%
print(f'Number of reactions after addition: {len(model.reactions)}')
