# %%[markdown]
# Author: Nelson Liu
# 
# Email: [nliu@uncharted.software](mailto:nliu@uncharted.software)

# %%[markdown]
# # H5N1 Model
# 
# Try to build a compartmental model for bird-cow-human infection model across US states.

# %%
import os
import json
import numpy
import pandas
import sympy
from typing import Optional, Iterable
import itertools

from mira.modeling.viz import GraphicalModel
from mira.metamodel import *
from mira.modeling import Model
from mira.modeling.amr.petrinet import template_model_to_petrinet_json
from mira.sources.amr.petrinet import template_model_from_amr_json

# %%
MIRA_REST_URL = 'http://34.230.33.149:8771/api'
PATH = "./data"

# %%
def generate_summary_table(model):

    data = {"name": [t.name for t in model.templates]}
    for k in ("subject", "outcome", "controller"):
        data[k] = [getattr(t, k).name if hasattr(t, k) else None for t in model.templates]

    data["controllers"] = [[c.name for c in getattr(t, k)] if hasattr(t, "controllers") else None for t in model.templates]
    data["controller(s)"] = [i if j == None else j for i, j in zip(data["controller"], data["controllers"])]
    __ = data.pop("controller")
    __ = data.pop("controllers")

    data["rate_law"] = [t.rate_law for t in model.templates]
    data["interactor_rate_law"] = [t.get_interactor_rate_law() for t in model.templates]

    df = pandas.DataFrame(data)

    return df

# %%[markdown]
# ## SIRHD Model

# %%
day_units = lambda: Unit(expression = sympy.Symbol("day"))
per_day_units = lambda: Unit(expression = 1 / sympy.Symbol("day"))

# Basic components
concepts = {
    c: Concept(name = c, display_name = c, description = f"Population of {d} agents")
    for c, d in zip(["S", "I", "R", "D"], ["susceptible", "infected", "recovered", "deceased"])
}

initials = {
    c: Initial(concept = concepts[c], expression = sympy.Float(0.0)) for c in concepts.keys()
}

# %%
model = TemplateModel(
    templates = [
        ControlledConversion(
            subject = concepts["S"],
            outcome = concepts["I"],
            controller = concepts["I"],
            name = "StoIviaI",
        ).with_mass_action_rate_law("a"),
        NaturalConversion(
            subject = concepts["I"],
            outcome = concepts["R"],
            name = "ItoR"
        ).with_mass_action_rate_law("b"),
        NaturalConversion(
            subject = concepts["I"],
            outcome = concepts["D"],
            name = "ItoD"
        ).with_mass_action_rate_law("c")
    ],
    initials = initials,
    observables = {},
    time = Time(name = "t", units = day_units()),
    annotations = Annotations(name = f"SIRD")
)

# %%
# Adjacency of US states
adj = []
with open("./data/us_states_adj.txt", "r") as f:
    while l := f.readline():
        adj.append(l.rstrip().split(","))

adj_edges = []
for a in adj:

    if len(a) < 2:
        continue

    adj_edges.extend(list(itertools.combinations(a, 2)))

list_us_states = []
for a in adj:
    list_us_states.extend(a)
list_states = list(set(list_us_states))

# %%
# Stratify by vaccination status
model_sp = stratify(
    model, 
    key = "species", 
    strata = ["B", "D", "H"], 
    cartesian_control = True, 
    structure = [["B", "D"], ["B", "H"], ["D", "H"]],
    concepts_to_stratify = None,
    params_to_stratify = None,
    param_renaming_uses_strata_names = True
)
model_sp.annotations.name = "SIRDwSp"

# %%
# Stratify by locations

model_sp_loc = stratify(
    model_sp,
    key = "locations",
    strata = list_us_states,
    cartesian_control = False,
    structure = None,
    concepts_to_stratify = None,
    params_to_stratify = None,
    param_renaming_uses_strata_names = True
)
model_sp_loc.annotations = "SIRDwSpwLoc"

# %%









