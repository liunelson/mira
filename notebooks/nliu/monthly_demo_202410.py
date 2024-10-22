# %%[markdown]
# Author: Nelson Liu
# 
# Email: [nliu@uncharted.software](mailto:nliu@uncharted.software)

# %%[markdown]
# # Monthly Demo (2024-10)

# %%
import os
import json
import numpy
import sympy
import pandas
import matplotlib.pyplot as plt
import matplotlib as mpl
import requests
import copy

from mira.sources import biomodels
from mira.modeling.viz import GraphicalModel
from mira.modeling.amr.petrinet import template_model_to_petrinet_json
from mira.sources.amr.petrinet import template_model_from_amr_json
from mira.metamodel import *
# from mira.modeling import Model

from mira.modeling.amr.ops import *
from mira.metamodel.io import expression_to_mathml

MIRA_REST_URL = "http://34.230.33.149:8771/api"

# %%
# Generate summary table of a template model
def generate_summary_table(model) -> pandas.DataFrame:

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

# Generate initial condition and parameter tables
def generate_init_param_tables(model) -> tuple[pandas.DataFrame, pandas.DataFrame]:

    data = {}
    data["name"] = [name for name, __ in model.initials.items()]
    data["expression"] = [init.expression for __, init in model.initials.items()]
    df_initials = pandas.DataFrame(data)

    data = {}
    data["name"] = [name for name, __ in model.parameters.items()]
    data["value"] = [param.value for __, param in model.parameters.items()]
    df_params = pandas.DataFrame(data)

    return (df_initials, df_params)

# %%
with open("./data/monthly_demo_202410/Prob 5 Model A.json", "r") as f:
    model = template_model_from_amr_json(json.load(f))

GraphicalModel.for_jupyter(model)

# %%[markdown]
# Problem 5
# 
# Test stratification to ~2,000 counties

# %%
model_stratified = stratify(
    model,
    key = 'location',
    strata = ['05001', '05003'],
    structure = [],
    directed = True,
    cartesian_control = False, 
    concepts_to_stratify = {'s', 'i', 'r'},
    params_to_stratify = {'y', 'R'},
    param_renaming_uses_strata_names = True
)

GraphicalModel.for_jupyter(model_stratified)

# %%
generate_summary_table(model_seirhd_effect_vacc)

# %%
generate_init_param_tables(model_seirhd_effect_vacc)[0]

# %%

