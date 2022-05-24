import json
import re

from pathlib import Path

from cobramod import add_pathway
from cobra.io import write_sbml_model, read_sbml_model
from cobra.flux_analysis import production_envelope

from cobramod import model_convert

dir_data = Path.cwd().resolve().joinpath("./data")

with open("data/RXN_DICT_20220414.json") as file_:
    RXN_DICT = json.load(file_)

compartment = 'c'


def add_c(match) -> str:
    cpd = match.group()
    return str(cpd + '_c')

def cpd_c(rxn):
    pat = r'C\d{5}'
    result = re.sub(pat, add_c, rxn)
    return result

def input(reactions):
    list = []
    for reaction in reactions:
        if reaction in RXN_DICT:
            RXN = RXN_DICT[reaction]
            if RXN["name"]==[""]:
                name = RXN["substrates"][0]
                equation = cpd_c(RXN["trans_equation"])
                reaction = reaction + '_' + compartment + ', ' + name + ' | ' + equation
                list.append(reaction)
            else:
                reaction = reaction + ', ' + compartment
                list.append(reaction)
        else:
            print("This pathway is not feasible!")
    return list


def get_yields(biomass_reactions, objective_id):
    original = read_sbml_model("./data/iJO1366.xml")
    model = original.copy()
    wt_growth = model.optimize()
    max_growth = wt_growth.objective_value
    min_growth = 0.8*max_growth
    reactions_list = input(biomass_reactions)
    objective = objective_id + '_' + compartment
    add_pathway(
        model=model,
        pathway= reactions_list,
        directory=dir_data,
        database="KEGG",
        compartment="c",
        genome="ecc",
        group="PWY",
        model_id="iJO1366"
    )
    model_convert(model)
    with model:
        medium = model.medium
        medium["EX_o2_e"] = 20.0
        medium["EX_glc__D_e"] = 20.0
        model.medium = medium
        growth_rate = model.optimize().objective_value
        if growth_rate < min_growth:
            print("This pathway is not feasible!")
        else:
            demand = model.add_boundary(model.metabolites.get_by_id(objective), type = "demand")
            model.objective = demand
            solution = model.optimize()
            production = solution.objective_value
            reactions = model.groups.get_by_id("PWY").members
            for reaction in reactions:
                reaction0 = str(reaction)
                fluxes = print(reaction0.split(':')[0], ':', reaction.flux)
            prod_env = production_envelope(model, ["EX_glc__D_e", "EX_o2_e"], objective=model.objective,
                                           carbon_sources = "EX_glc__D_e")
    return fluxes, print('Maximum productivity =', production, 'mmol/gDW*h'), \
           print('Maximum theoretical yield =', prod_env.loc[0, 'mass_yield_maximum'], 'g/g')

if __name__ == '__main__':
    get_yields(
        ["R12188",
        "R12514",
        "R04951",
        "R04951",
        "R10553",
        "R10617",
        "R10617",
        "R10617",
        "R10617",
        "R09591",
        "R09622",
        "R10301",
        "R12539",
        "R11156",
        "R11157",
        "R11156",
        "R11157",
        "R04789",
        "R01628",
        "R02193",
        "R00931"], 'C00100')
    get_yields(
        ["R03105",
        "R00896",
        "R04051",
        "R04051",
        "R04051",
        "R00355",
        "R00489",
        "R00904",
        "R03139",
        "R05329"],'C00986'
    )

