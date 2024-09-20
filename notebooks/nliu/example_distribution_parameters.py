# %%[markdown]
# Author: Nelson Liu
# 
# Email: [nliu@uncharted.software](mailto:nliu@uncharted.software)

# %%[markdown]
# # Test new MIRA features
# 
# 1. [Support Expressions in Parameters of Distributions](https://github.com/gyorilab/mira/pull/358)
# 2. [Implement adding observables by pattern](https://github.com/gyorilab/mira/pull/364)


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
from mira.modeling.amr.petrinet import template_model_to_petrinet_json
from mira.sources.amr.petrinet import template_model_from_amr_json

# %%
# Import example model
with open("./data/monthly_demo_202408/model_seirhd.json", "r") as f:
    model = template_model_from_amr_json(json.load(f))

GraphicalModel.for_jupyter(model)

# %%[markdown]
# ## 1. Expression for distribution parameters 

# %%
new_parameters = {
    'u': Parameter(
        name = 'u', 
        display_name = 'u', 
        description = 'Population count uncertainty', 
        value = 1.0
    ),
    'b_shape': Parameter(
        name = 'b_shape', display_name = 'b_shape', description = 'Beta distribution shape parameter for the parameter "b"', 
        distribution = Distribution(
            type = 'Beta1', 
            parameters = {
                'alpha': sympy.Integer(1),
                'beta': sympy.Integer(10)
            }
        )
    ),
    'p_u': Parameter(
        name = 'p_u', 
        display_name = 'p_u', 
        description = 'Percent uncertainty on probability transitions', 
        value = 0.10
    )
}

model.parameters = model.parameters | new_parameters

# %%
# Edit parameters to use new parameters in their distributions

model.parameters['b'] = Parameter(
    name = 'b',
    display_name = 'b',
    description = 'Infection rate',
    distribution = Distribution(
        type = 'InverseGamma1',
        parameters = {
            'shape': sympy.Symbol('b_shape') ** sympy.Float(1.5),
            'scale': 0.01
        }
    )
)

# for p in model.parameters.keys():
#     if p[0] == "r":
#         v = model.parameters[p].value
#         model.parameters[p].distribution.parameters = {
#             'minimum': sympy.Float(v) - sympy.Symbol('u'),
#             'maximum': sympy.Float(v) + sympy.Symbol('u')
#         }

# %%
with open('./data/example_models/seirhd_dist_params.json', 'w') as f:
    json.dump(template_model_to_petrinet_json(model), f, indent = 4)

# %%[markdown]
# ## 2. Assisted observable addition by patterns

# %%
# By vaccination status
model_vax = stratify(
    model, 
    key = "vaccination_status",
    strata = ["vaccinated", "unvaccinated"],
    directed = True,
    structure = [["unvaccinated", "vaccinated"]],
    cartesian_control = True, 
    concepts_to_stratify = ["S", "E", "I"],
    params_to_stratify = [],
    param_renaming_uses_strata_names = True
)

GraphicalModel.for_jupyter(model_vax)

# By vaccine
model_vax_vaccine = stratify(
    model_vax, 
    key = "vaccine",
    strata = ["moderna", "pfizer", "jj"],
    directed = False,
    structure = [],
    cartesian_control = False,
    concepts_to_stratify = ["S_vaccinated", "E_vaccinated", "I_vaccinated"],
    params_to_stratify = [],
    param_renaming_uses_strata_names = True
)

GraphicalModel.for_jupyter(model_vax_vaccine)

# By age
model_vax_vaccine_age = stratify(
    model_vax_vaccine,
    key = "age",
    strata = ["children", "adults"],
    directed = False,
    structure = [],
    cartesian_control = True,
    concepts_to_stratify = None,
    params_to_stratify = None,
    param_renaming_uses_strata_names = True
)

GraphicalModel.for_jupyter(model_vax_vaccine_age)

# %%[markdown]
# ### Example 1:
# Add an observable named TotalSusceptibleVaccinatedModernaChildren to the model that is the total number of (1) susceptible (2) vaccinated (3) with Moderna (4) children
# 
# ```python
# add_observable_pattern(
#     model_vax_vaccine_age,
#     'TotalSusceptibleVaccinatedModernaChildren',
#     identifiers = {'ido': '0000514'},
#     context = {
#         'vaccination_status': 'vaccinated', # stratification_key: strata_name
#         'vaccine': 'moderna',
#         'age': 'children'
#     }
# )
# ```
# 
# There'd be a new observable added to the model
# `TotalSusceptibleVaccinatedModernaChildren = S_vaccinated_moderna_child`
# 
# ### Example 2:
# Add an observable named TotalUnvaccinatedChildren that is the total number of (1) unvaccinated (2) children
# 
# ```python
# add_observable_pattern(
#     model_vax_vaccine_age,
#     'TotalunvaccinatedChildren',
#     context = {
#         'vaccination_status': 'unvaccinated',
#         'age': 'children'
#     }
# )
# ```
# 
# There'd be a new observable added to the model
# `TotalUnvaccinatedhildren = S_unvaccinated_children + I_unvaccinated_children + R_unvaccinated_children`

# %%
