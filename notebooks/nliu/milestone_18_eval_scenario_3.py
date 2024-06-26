# %%[markdown]
# Author: Nelson Liu
# 
# Email: [nliu@uncharted.software](mailto:nliu@uncharted.software)

# %%[markdown]
# # 18-Month Epi Evaluation Scenario 3

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

import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcParams['text.usetex'] = True

# %%
MIRA_REST_URL = 'http://34.230.33.149:8771/api'
PATH = "./data/milestone_18_hackathon/scenario_1"

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

# %%
def generate_odesys(model, latex_align: bool = False) -> list:

    odeterms = {var: 0 for var in model.get_concepts_name_map().keys()}

    for t in model.templates:
        if hasattr(t, "subject"):
            var = t.subject.name
            odeterms[var] -= t.rate_law.args[0]
        
        if hasattr(t, "outcome"):
            var = t.outcome.name
            odeterms[var] += t.rate_law.args[0]


    # Convert to LaTeX
    if latex_align == True:
        odesys = [
            "\\frac{d " + f"{var}" + "}{d t}" + " &= " + sympy.latex(terms)
            for var, terms in odeterms.items()
        ]
        odesys = "\\begin{align*}\n" + "\n".join(["    " + expr + "\\\\" for expr in odesys]) + "\n\\end{align}"

    else:
        odesys = [
            "\\frac{d " + f"{var}" + "}{d t}" + " = " + sympy.latex(terms)
            for var, terms in odeterms.items()
        ]

    return odesys

# %%[markdown]
# ## Build Model

# %%
# Define units of measurement
person_units = lambda: Unit(expression = sympy.Symbol("person"))
day_units = lambda: Unit(expression = sympy.Symbol("day"))
per_day_units = lambda: Unit(expression = 1 / sympy.Symbol("day"))
dimensionless_units = lambda: Unit(expression=sympy.Integer('1'))
per_day_per_person_units = lambda: Unit(expression = 1 / (sympy.Symbol('day') * sympy.Symbol('person')))

# Basic components
# Not safe to use "I", might be interpreted as imaginary "i" by Sympy
concepts = {
    c: Concept(name = c, display_name = c, description = f"Population of {d} agents")
    for c, d in zip(["S", "E", "F", "R"], ["susceptible", "exposed", "infected", "recovered"])
}

# Add groundings
concepts["S"].identifiers = {"ido": "0000514"}
concepts["E"].identifiers = {"ido": "0000154"}
concepts["F"].identifiers = {"ido": "0000511"}
concepts["R"].identifiers = {"ido": "0000592"}

# Initial conditions
initials = {
    c: Initial(concept = concepts[c], expression = sympy.Float(0.0)) for c in concepts.keys()
}

# Parameters
parameters = {
    "b": Parameter(name = "beta", value = sympy.Float(0.2), units = per_day_units(), description = "Exposure factor"),
    "f": Parameter(name = "mew", value = sympy.Float(0.0), units = per_day_units(), description = "Mask efficacy coefficient (time-dependent)"),
    "c": Parameter(name = "mcw", value = sympy.Float(0.0), units = per_day_units(), description = "Mask compliance coefficient (time-dependent)"),
    "m": Parameter(name = "M", value = sympy.Float(10.0), units = per_day_units(), description = "Contact rate"),
    "n": Parameter(name = "N", value = sympy.Float(10305660 + 15281905 + 12154442 + 50*3 + 50*3), units = person_units(), description = "Total population"),
    "g": Parameter(name = "rEI", value = sympy.Float(0.08), units = per_day_units(), description = "Infection rate (from Exposed to Infected)"),
    "h": Parameter(name = "rIR", value = sympy.Float(0.06), units = per_day_units(), description = "Recovery rate (from Infected to Recovered)")
}

S, E, F, R, b, f, c, m, n, g, h = sympy.symbols("S E F R b f c m n g h")

# Add prior on the `b` parameter
parameters["b"].distribution = Distribution(type = "StandardUniform1", parameters = {"minimum": 0.01, "maximum": 0.99})

# %%
model = TemplateModel(
    templates = [
        ControlledConversion(
            subject = concepts["S"],
            outcome = concepts["E"],
            controller = concepts["F"],
            name = "ExposureProcess",
            rate_law = b * S * (1 - f * c ) * m * F / n,
        ),
        NaturalConversion(
            subject = concepts["E"],
            outcome = concepts["F"],
            name = "InfectionProcess",
        ).with_mass_action_rate_law("g"),
        NaturalConversion(
            subject = concepts["F"],
            outcome = concepts["R"],
            name = "RecoveryProcess"
        ).with_mass_action_rate_law("h"),
    ],
    parameters = parameters,
    initials = initials,
    observables = {},
    time = Time(name = "t", units = day_units()),
    annotations = Annotations(
        name = "SEFR model",
        description = "COVID-19 compartmental model with Susceptible, Exposed, inFected (avoid confusion with imaginary I), and Recovered states",
        authors = [{"name": "Nelson Liu"}, {"name": "Ben Gyori"}, {"name": "MITRE Team"}],
        license = "CC0",
        references = ["pubmed:32616574"],
        locations = ["geonames:5128581"],
        hosts = ["ncbitaxon:9606"],
        diseases = ["doid:0080600"],
        pathogens = ["ncbitaxon:2697049"],
        model_types = ["mamo:0000028", "mamo:0000046"],
        # time_start = "2020-03-01T00:00:00",
        # time_end = "2020-08-01T00:00:00"
    )
)

# %%
# model.draw_jupyter()
GraphicalModel.for_jupyter(model)

# %%
odesys = generate_odesys(model, latex_align = True)
print(odesys)

# %%
generate_summary_table(model)

# %%
with open("./data/milestone_18_evaluation/scenario3/SEFR.json", "w") as f:
    json.dump(template_model_to_petrinet_json(model), f, indent = 4)

# %%
# Stratify by three age groups
model_age = stratify(
    model, 
    key = "age", 
    strata = ["1", "2", "3"], 
    cartesian_control = True, 
    structure = [],
    concepts_to_stratify = ["S", "E", "F", "R"],
    params_to_stratify = ["m"],
    param_renaming_uses_strata_names = True
)
model_age.annotations.name = "SEFR model stratified by 3 age groups"
model_age.annotations.description = "COVID-19 compartmental model with Susceptible, Exposed, inFected (avoid confusion with imaginary I), and Recovered states; stratified by three age groups",

# %%
# Set initial conditions

initial_conditions = {
    "S_1": 10305660.0,
    "S_2": 15281905.0,
    "S_3": 12154442.0,
    "E_1": 50.0,
    "E_2": 50.0,
    "E_3": 50.0,
    "F_1": 50.0,
    "F_2": 50.0,
    "F_3": 50.0,
    "R_1": 0.0,
    "R_2": 0.0,
    "R_3": 0.0
}

for k, v in initial_conditions.items():
    model_age.initials[k].expression = sympy.Float(v)

# %% 
# Set contact matrix values
df_ContactMatrix = pandas.read_csv("./data/milestone_18_evaluation/scenario3/ContactMatrix.csv")
df_ContactMatrix = df_ContactMatrix.drop(df_ContactMatrix.columns[0], axis = 1)

for i in range(3):
    for j in range(3):
        name = f"m_{i + 1}_{j + 1}"
        model_age.parameters[name].value = df_ContactMatrix.iloc[i, j]
        model_age.parameters[name].description = f"Contact rate (between Age Group {i + 1} and {j + 1})"


# %%
GraphicalModel.for_jupyter(model_age)

# %%
generate_summary_table(model_age)

# %%
odesys = generate_odesys(model_age, latex_align = True)
print(odesys)

# %%
# Export to AMR
with open("./data/milestone_18_evaluation/scenario3/SEFR_age.json", "w") as f:
    json.dump(template_model_to_petrinet_json(model_age), f, indent = 4)

# %%