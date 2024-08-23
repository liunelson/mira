# %%[markdown]
# Author: Nelson Liu
# 
# Email: [nliu@uncharted.software](mailto:nliu@uncharted.software)

# %%[markdown]
# # Monthly Demo (2024-08)

# %%
import os
import json
import mira.metamodel
import numpy
import sympy
import pandas
import matplotlib.pyplot as plt
import matplotlib as mpl
import requests
import copy

import mira
from mira.sources import biomodels
from mira.modeling.viz import GraphicalModel
from mira.modeling.amr.petrinet import template_model_to_petrinet_json
from mira.sources.amr.petrinet import template_model_from_amr_json
from mira.dkg.web_client import is_ontological_child_web
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
with open("./data/monthly_demo_202408/model_seirhd.json", "r") as f:
    model_seirhd = template_model_from_amr_json(json.load(f))

GraphicalModel.for_jupyter(model_seirhd)

# %%[markdown]
# ## Problem 3
# 
# 12-month evaluation, scenario 2
# 
# ### Question 1
# 
# 1. Start from the base SEIRHD model from Scenario 1
# 2. Edit the rate law of the exposure/infection process to include the effects of vaccination, masking, contact
# 3. Stratify only the state `S` and parameter `beta`` by vaccination status (vaccinated, unvaccinated)
# 4. Stratify only the vaccinated states and  `beta`` by vaccine used (pfizer, moderna, jj)
# 5. Stratify the mRNA-vaccinated states and beta by the number of doses (1, 2nd same, 2nd difference)
# 6. Remove all 0-2 dose transitions

# %%
# Add a masking-effect factor `m` to the exposure process (S -> E controlled by I)
model_seirhd_effect = copy.deepcopy(model_seirhd)
model_seirhd_effect.add_parameter(
    "v", "v", "Vaccine efficacy effect on the disease exposure process (0 = complete efficacy, 1 = no efficacy)", 1.0
)
model_seirhd_effect.add_parameter(
    "m", "m", "Masking effect on the disease exposure process (0 = unmasked, 1 = fully masked)", 0.0
)
model_seirhd_effect.add_parameter(
    "c", "c", "Contact effect on the disease exposure process (0 = no contact, 1 = full contact)", 1.0
)

# Find the exposure transition
id = [template.name for template in model_seirhd_effect.templates if (template.subject.name == "S") & (template.outcome.name == "E")][0]

# Change its rate law
amr = template_model_to_petrinet_json(model_seirhd_effect)
amr = replace_rate_law_sympy(amr, id, sympy.Symbol("c") * (1 - sympy.Symbol("m")) * sympy.Symbol("v") * sympy.Symbol("b") * sympy.Symbol("invN") * sympy.Symbol("S") * sympy.Symbol("I"))
model_seirhd_effect = template_model_from_amr_json(amr)

# Check
generate_summary_table(model_seirhd_effect)

# %%
model_seirhd_effect.annotations.name = "SEIRHD model with effects"
model_seirhd_effect.annotations.description = "Edit of the SEIRHD model with effect parameters on the exposure process"

# Save
with open("./data/monthly_demo_202408/model_seirhd_effect.json", "w") as f:
    json.dump(template_model_to_petrinet_json(model_seirhd_effect), f, indent = 4)

# %%
# Stratify by vaccination status
# * only unvaccinated -> vaccinated is permitted
# * only the infection rate `beta` and masking effect `m`
# * only for S, E, I

model_seirhd_effect_vacc = mira.metamodel.stratify(
    model_seirhd_effect,
    key = "vaccination",
    strata = ["unvaccinated", "vaccinated"],
    structure = [["unvaccinated", "vaccinated"]],
    directed = True,
    cartesian_control = False, 
    concepts_to_stratify = {"S"},
    params_to_stratify = {"v", "m"},
    param_renaming_uses_strata_names = True
)

GraphicalModel.for_jupyter(model_seirhd_effect_vacc)

# %%
generate_summary_table(model_seirhd_effect_vacc)

# %%
generate_init_param_tables(model_seirhd_effect_vacc)[0]

# %%
generate_init_param_tables(model_seirhd_effect_vacc)[1]

# %%
model_seirhd_effect_vacc.annotations.name = "SEIRHD model with effects (vaccination)"
model_seirhd_effect_vacc.annotations.description = "Edit of the SEIRHD model with effect parameters on the exposure process"

# Save
with open("./data/monthly_demo_202408/model_seirhd_effect_vacc.json", "w") as f:
    json.dump(template_model_to_petrinet_json(model_seirhd_effect_vacc), f, indent = 4)

# %%
# Stratify by vaccine used in vaccination
# * only vaccinated pops
# * only vaccine efficacy (v_*) and vaccination rate (p)

concepts_to_stratify = {
    c.name 
    for c in model_seirhd_effect_vacc.get_concepts_name_map().values() 
    if c.name.endswith("_vaccinated")
}
params_to_stratify = {
    p 
    for p in model_seirhd_effect_vacc.parameters.keys()
    if p.startswith("v_") or p.startswith("p")
}

model_seirhd_effect_vacc2 = mira.metamodel.stratify(
    model_seirhd_effect_vacc,
    key = "vaccine",
    strata = ["pfizer", "moderna", "jj"],
    structure = [],
    directed = True,
    cartesian_control = False,
    concepts_to_stratify = concepts_to_stratify,
    params_to_stratify = params_to_stratify,
    param_renaming_uses_strata_names = True
)

GraphicalModel.for_jupyter(model_seirhd_effect_vacc2)

# %%
generate_summary_table(model_seirhd_effect_vacc2)

# %%
generate_init_param_tables(model_seirhd_effect_vacc2)[0]

# %%
generate_init_param_tables(model_seirhd_effect_vacc2)[1]

# %%
model_seirhd_effect_vacc2.annotations.name = "SEIRHD model with effects (vaccination, vaccine brands)"
model_seirhd_effect_vacc2.annotations.description = "Edit of the SEIRHD model with effect parameters on the exposure process"

# Save
with open("./data/monthly_demo_202408/model_seirhd_effect_vacc2.json", "w") as f:
    json.dump(template_model_to_petrinet_json(model_seirhd_effect_vacc2), f, indent = 4)

# %%
# Stratify by the number of vaccine doses
# * only mRNA-vaccinated pops (S_vaccinated_moderna, S_vaccinated_pfizer)
# * only vaccine efficacy "v_vaccinated*" and vaccination rate "p_*"
# * one dose, two dose, different second dose

concepts_to_stratify = {
    c.name 
    for c in model_seirhd_effect_vacc2.get_concepts_name_map().values() 
    if c.name.startswith(("S_vaccinated_pfizer")) or c.name.startswith("S_vaccinated_moderna")
}

params_to_stratify = {
    p 
    for p in model_seirhd_effect_vacc2.parameters.keys()
    if p.startswith("v_vaccinated")
}

model_seirhd_effect_vacc2_dose = mira.metamodel.stratify(
    model_seirhd_effect_vacc2,
    key = "dose",
    strata = ["one", "two", "diff"],
    structure = [["one", "two"], ["one", "diff"]],
    directed = True,
    cartesian_control = False,
    concepts_to_stratify = concepts_to_stratify,
    params_to_stratify = params_to_stratify,
    param_renaming_uses_strata_names = True
)

# %%
# # Second dose of vaccination
# amr = template_model_to_petrinet_json(model_seirhd_effect_vacc2)
# amr = add_transition(
#     amr, 
#     new_transition_id = "VaccinationProcess2", 
#     src_id = "S_vaccinated", 
#     tgt_id = "S_vaccinated", 
#     rate_law_mathml = expression_to_mathml(sympy.Symbol("v") * sympy.Symbol("S_unvaccinated")),
#     params_dict = {
#         "v": {"value": 0.01, "description": "Vaccination rate of susceptible persons"}
#     }
# )
# model_seirhd_effect_vacc2 = template_model_from_amr_json(amr)

# %%
# Remove forbidden transitions:
# * unvaccinated -> 2 doses of any vaccine (*_two, *_diff)

forbidden_transition_ids = [t.name 
    for t in model_seirhd_effect_vacc2_dose.templates 
    if isinstance(t, NaturalConversion)
        and t.subject.name.endswith("unvaccinated") 
        and (t.outcome.name.endswith("two") or t.outcome.name.endswith("diff"))
]

amr = template_model_to_petrinet_json(model_seirhd_effect_vacc2_dose)
for t in model_seirhd_effect_vacc2_dose.templates:
    if t.name in forbidden_transition_ids:
        amr = remove_transition(amr, transition_id = t.name)

model_seirhd_effect_vacc2_dose = template_model_from_amr_json(amr)

GraphicalModel.for_jupyter(model_seirhd_effect_vacc2_dose)

# %%
generate_summary_table(model_seirhd_effect_vacc2_dose)

# %%
generate_init_param_tables(model_seirhd_effect_vacc2_dose)[0]

# %%
generate_init_param_tables(model_seirhd_effect_vacc2_dose)[1]

# %%
model_seirhd_effect_vacc2_dose.annotations.name = "SEIRHD model with effects (vaccination, vaccine brands, doses)"
model_seirhd_effect_vacc2_dose.annotations.description = "Edit of the SEIRHD model with effect parameters on the exposure process"

# Save
with open("./data/monthly_demo_202408/model_seirhd_effect_vacc2_dose.json", "w") as f:
    json.dump(template_model_to_petrinet_json(model_seirhd_effect_vacc2_dose), f, indent = 4)

# %%[markdown]
# ### Question 2
# 
# 1. Start from the final SEIRHD model from above
# 2. Stratify all states (except `D`) by three age groups ("0to9", "10to19", "20above")
# 3. Preserve the `v`, `invN` parameters
# 4. Remove all vaccination processes and vaccinated states associated with the age group "0to9"

# %%
params_to_preserve = {
    p 
    for p in model_seirhd_effect_vacc2_dose.parameters.keys()
    if (p in ("b", "invN")) or p.startswith("v_")
}

model_seirhd_effect_vacc2_dose_age = mira.metamodel.stratify(
    model_seirhd_effect_vacc2_dose,
    key = "age",
    strata = ["0to9", "10to19", "20above"],
    structure = [],
    directed = True,
    cartesian_control = True,
    concepts_to_preserve = {"D"},
    params_to_preserve = params_to_preserve,
    param_renaming_uses_strata_names = True
)

# %%
# Remove all vaccination processes associated with age group "0to9"
forbidden_transition_ids = [t.name 
    for t in model_seirhd_effect_vacc2_dose_age.templates 
    if isinstance(t, NaturalConversion)
        and t.subject.name.endswith("0to9") 
        and t.outcome.name.startswith("S_vaccinated")
]

amr = template_model_to_petrinet_json(model_seirhd_effect_vacc2_dose_age)
for t in model_seirhd_effect_vacc2_dose_age.templates:
    if t.name in forbidden_transition_ids:
        amr = remove_transition(amr, transition_id = t.name)

model_seirhd_effect_vacc2_dose_age = template_model_from_amr_json(amr)

GraphicalModel.for_jupyter(model_seirhd_effect_vacc2_dose_age)

# %%
generate_summary_table(model_seirhd_effect_vacc2_dose_age)

# %%
generate_init_param_tables(model_seirhd_effect_vacc2_dose_age)[0]

# %%
generate_init_param_tables(model_seirhd_effect_vacc2_dose_age)[1]

# %%
model_seirhd_effect_vacc2_dose_age.annotations.name = "SEIRHD model with effects (vaccination, vaccine brands, doses, age contacts)"
model_seirhd_effect_vacc2_dose_age.annotations.description = "Edit of the SEIRHD model with effect parameters on the exposure process"

# Save
with open("./data/monthly_demo_202408/model_seirhd_effect_vacc2_dose_age.json", "w") as f:
    json.dump(template_model_to_petrinet_json(model_seirhd_effect_vacc2_dose_age), f, indent = 4)

# %%
