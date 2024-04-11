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
from mira.sources.amr.petrinet import template_model_from_amr_json

# %%
p = "../evaluation_2024.03/epi_scenario4/scenario4_petrinet.json"
with open(p, "r") as f:
    j = json.load(f)
    if "configuration" in j.keys():
        tm = template_model_from_amr_json(j["configuration"])
    else:
        tm = template_model_from_amr_json(j)

# %%

