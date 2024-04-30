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
from copy import deepcopy
from tqdm import tqdm

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
all_used_parameters = tm.get_all_used_parameters() - set(list(tm.parameters.keys()))- set(["time"])

for parameter in all_used_parameters:
    tm.add_parameter(parameter, name = parameter, value = 0.0)

# %%
model = Model(tm)

real_params = {k: v for k, v in model.parameters.items() if not v.placeholder}
pmap = {parameter.key: idx for idx, (pkey, parameter) in enumerate(real_params.items())}
vmap = {variable.key: idx for idx, variable in enumerate(model.variables.values())}

concept_map = {variable.concept.name: variable.key for variable in model.variables.values()}
parameter_map = {parameter.concept.name: parameter.key for parameter in real_params.values()}

y = sympy.MatrixSymbol('y', len(model.variables), 1)
p = sympy.MatrixSymbol('p', len(real_params), 1)

# %%
for transition in tqdm(model.transitions.values()):
    # Use rate if available which is a symbolic expression
    if transition.template.rate_law:
        rate = deepcopy(transition.template.rate_law.args[0])
        for symbol in rate.free_symbols:
            sym_str = str(symbol)
            if sym_str in concept_map:
                rate = rate.subs(symbol, y[vmap[concept_map[sym_str]]])
            elif sym_str in pmap:
                rate = rate.subs(symbol, p[pmap[parameter_map[sym_str]]])
            elif model.template_model.time and sym_str == model.template_model.time.name:
                rate = rate.subs(symbol, 't')
            else:
                assert False


# %%
om = OdeModel(Model(tm), initialized=True)
om.kinetics

# %%
