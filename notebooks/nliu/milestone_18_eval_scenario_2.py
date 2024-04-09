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
from mira.modeling.amr.petrinet import template_model_to_petrinet_json

# %%
MIRA_REST_URL = 'http://34.230.33.149:8771/api'
PATH = "./data/milestone_18_hackathon/scenario_1"

# %%[markdown]
# ## SIRHD Model

# %%
day_units = lambda: Unit(expression = sympy.Symbol("day"))
per_day_units = lambda: Unit(expression = 1 / sympy.Symbol("day"))

# Basic components
concepts = {
    c: Concept(name = c, display_name = c, description = f"Population of {d} agents")
    for c, d in zip(["S", "I", "R", "D", "H"], ["susceptible", "infected", "recovered", "deceased", "hospitalized"])
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
            name = "S_to_I_via_I",
        ),
        NaturalConversion(
            subject = concepts["I"],
            outcome = concepts["H"],
            name = "I_to_H",
        ),
        NaturalConversion(
            subject = concepts["I"],
            outcome = concepts["R"],
            name = "I_to_R"
        ),
        NaturalConversion(
            subject = concepts["H"],
            outcome = concepts["R"],
            name = "H_to_R"
        ),
        NaturalConversion(
            subject = concepts["H"],
            outcome = concepts["D"],
            name = "H_to_D"
        )
    ],
    initials = initials,
    observables = {},
    time = Time(name = "t", units = day_units()),
    annotations = Annotations(name = f"SIRHD")
)

# %%
# model.draw_jupyter()
GraphicalModel.for_jupyter(model)

# %%
# Alphabet
alpha = list(map(chr, range(ord('a'), ord('z')+1)))

# Set mass action kinetics
for i, template in enumerate(model.templates):
    p = alpha[i]
    template.set_mass_action_rate_law(p)
    model.add_parameter(
        parameter_id = p,
        name = p,
        value = 1.0,
        description = f"Rate constant for {template.name}"
    )

# %%
# Customize the I-to-H rate law

vaxItoH, ageItoH = sympy.symbols('vaxItoH, ageItoH')
for t in model.templates:
    if t.name == "I_to_H":
        t.rate_law = t.rate_law.args[0] * vaxItoH * ageItoH

# %%
# Stratify by vaccination status
model_vax = stratify(
    model, 
    key = "vaxstatus", 
    strata = ["unvax", "vax"], 
    cartesian_control = True, 
    structure = [],
    concepts_to_stratify = None,
    params_to_stratify = ["vaxI_to_H"]
)

GraphicalModel.for_jupyter(model_vax)

# %%
# Add vaccination dynamics

t = NaturalConversion(
    subject = model_vax.get_concepts_name_map()["S_unvax"],
    outcome = model_vax.get_concepts_name_map()["S_vax"],
    name = "vaccination"
)
t.set_mass_action_rate_law(alpha[len(model_vax.parameters)])

model_vax = model_vax.add_template(t)

# %%
GraphicalModel.for_jupyter(model_vax)

# %%

########################



# %%
# Stratified by 16 age groups
num_ages = 16
model_age = stratify(
    model, 
    key = "age",
    strata = [f"A{i}" for i in range(num_ages)],
    cartesian_control = True, 
    structure = [], 
    directed = True,
    params_to_stratify = [p for t in model.templates if t.subject.name in ("S") for p in t.get_parameter_names()],
    concepts_to_stratify = ["S", "E", "M", "C"]
)

# %%
GraphicalModel.for_jupyter(model_age)

# %%
# Multi-stratify the base model
# 1. number of vaccine doses: zero, once, twice
# 2. age group: y, o
# 3. sex: male or female
# 4. ethnicity: white, black, asian, latino

model_vac = stratify(
    model, 
    key = "vac",
    strata = ["0", "1", "2"],
    cartesian_control = True, 
    structure = [["0", "1"], ["1", "2"]], 
    directed = True,
    params_to_stratify = None,
    concepts_to_stratify = None
)

# %%
model_vac_age = stratify(
    model_vac,
    key = "age",
    strata = ["y", "o"],
    cartesian_control = True,
    structure = [],
    directed = True,
    params_to_stratify = None,
    concepts_to_stratify = None
)

# %%
model_vac_age_sex = stratify(
    model_vac_age,
    key = "sex",
    strata = ["m", "f", "q"],
    cartesian_control = True,
    structure = [],
    directed = True,
    params_to_stratify = None,
    concepts_to_stratify = None
)

# %%
model_vac_age_sex_eth = stratify(
    model_vac_age_sex,
    key = "ethnicity",
    strata = ["w", "b", "l", "a"],
    cartesian_control = True,
    structure = [],
    directed = True,
    params_to_stratify = None,
    concepts_to_stratify = None
)

# %%
# GraphicalModel.for_jupyter(model_vac_age_sex_eth)

# %%
model_age.annotations.name = "SEIIRD_Age16"
model_vac.annotations.name = "SEIIRD_Vac3"
model_vac_age.annotations.name = "SEIIRD_Vac3_Age2"
model_vac_age_sex.annotations.name = "SEIIRD_Vac3_Age2_Sex3"
model_vac_age_sex_eth.annotations.name = "SEIIRD_Vac3_Age2_Sex3_Eth4"

# %%
for name, m in zip(("model", "model_age", "model_vac", "model_vac_age", "model_vac_age_sex", "model_vac_age_sex_eth"), (model, model_age, model_vac, model_vac_age, model_vac_age_sex, model_vac_age_sex_eth)):
    with open(os.path.join(PATH, f"{name}.json"), "w") as f:
        j = template_model_to_petrinet_json(m)
        json.dump(j, f, indent = 4)

# %%
