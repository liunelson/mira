# %%[markdown]
# Author: Nelson Liu
# 
# Email: [nliu@uncharted.software](mailto:nliu@uncharted.software)

# %%[markdown]
# # Epi Evaluation (2024-04-05)
# 
# Try out different double stratification schemes

# %%
import os
import json
import numpy
import pandas
import sympy
from typing import Optional, Iterable

from mira.modeling.viz import GraphicalModel
from mira.metamodel import *
from mira.modeling import Model
from mira.modeling.ode import OdeModel, simulate_ode_model
from mira.modeling.amr.petrinet import template_model_to_petrinet_json
from mira.sources.amr.regnet import template_model_from_amr_json

# %%
p = "scenario5_regnet.json"
with open(p, "r") as f:
    j = json.load(f)
    if "configuration" in j.keys():
        tm = template_model_from_amr_json(j["configuration"])
    else:
        tm = template_model_from_amr_json(j)

# %%
data = {"name": [], "value": [], "min": [], "max": []}
for __, p in tm.parameters.items():
    data["name"].append(p.name)
    data["value"].append(p.value)
    data["min"].append(p.distribution.parameters["minimum"])
    data["max"].append(p.distribution.parameters["maximum"])

df = pandas.DataFrame(data)

df

# %%
# df["min"] < df["max"]
# df["min"] <= df["value"]
df["value"] <= df["max"]

# %%
def gen_template_table(templates):

    props = {}
    for k in ["subject", "outcome", "controller", "controllers", "rate_law"]:
        if k == "controllers":
            props[k] = [[c.name for c in getattr(t, k)] if hasattr(t, k) else "" for t in templates]
        elif k == "rate_law":
            props[k] = [getattr(t, k) if hasattr(t, k) else "" for t in templates]
        else:
            props[k] = [getattr(t, k).name if hasattr(t, k) else "" for t in templates]

    return pandas.DataFrame(props)

gen_template_table(tm.templates)

# %%
om = OdeModel(Model(tm), initialized=True)

om.kinetics

# %%
