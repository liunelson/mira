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
from mira.sources.amr.petrinet import template_model_from_amr_json

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
            "\\frac{\\mathrm{d}}{\\mathrm{d} t}" + f"{var} & = " + sympy.latex(terms)
            for var, terms in odeterms.items()
        ]
        odesys = "\\begin{align*}\n" + "\n".join(["    " + expr + "\\\\" for expr in odesys]) + "\n\\end{align}"

    else:
        odesys = [
            "\\frac{\\mathrm{d}}{\\mathrm{d} t}" + f"{var} = " + sympy.latex(terms)
            for var, terms in odeterms.items()
        ]

    return odesys

# %%[markdown]
# ## SIRHD Model

# %%
day_units = lambda: Unit(expression = sympy.Symbol("day"))
per_day_units = lambda: Unit(expression = 1 / sympy.Symbol("day"))

# Basic components
concepts = {
    c: Concept(name = c, display_name = c, description = f"Population of {d} agents")
    for c, d in zip(["S", "C", "R", "D", "H"], ["susceptible", "infected", "recovered", "deceased", "hospitalized"])
}

initials = {
    c: Initial(concept = concepts[c], expression = sympy.Float(0.0)) for c in concepts.keys()
}

# %%
model = TemplateModel(
    templates = [
        ControlledConversion(
            subject = concepts["S"],
            outcome = concepts["C"],
            controller = concepts["C"],
            name = "S_to_C_via_C",
        ).with_mass_action_rate_law("a"),
        NaturalConversion(
            subject = concepts["C"],
            outcome = concepts["H"],
            name = "C_to_H",
        ).with_mass_action_rate_law("b"),
        NaturalConversion(
            subject = concepts["C"],
            outcome = concepts["R"],
            name = "C_to_R"
        ).with_mass_action_rate_law("c"),
        NaturalConversion(
            subject = concepts["H"],
            outcome = concepts["R"],
            name = "H_to_R"
        ).with_mass_action_rate_law("d"),
        NaturalConversion(
            subject = concepts["H"],
            outcome = concepts["D"],
            name = "H_to_D"
        ).with_mass_action_rate_law("e")
    ],
    initials = initials,
    observables = {},
    time = Time(name = "t", units = day_units()),
    annotations = Annotations(name = f"SCRHD")
)

# %%
# model.draw_jupyter()
GraphicalModel.for_jupyter(model)

# %%
# Customize the I-to-H rate law
# vaxItoH, ageItoH = sympy.symbols('vaxItoH, ageItoH')
for t in model.templates:
    if t.name == "C_to_H":
        # t.rate_law = t.rate_law.args[0] * vaxItoH * ageItoH
        t.set_rate_law("C * b * vaxCtoH * ageCtoH")
        # Not safe to "I", might be interpreted as imaginary "i"

# %%
# Initialize the parameters
model.parameters = {p: Parameter(name = p, value = 1.0) for p in model.get_all_used_parameters()}

# %%
generate_summary_table(model)

# %%
with open("./scenario2.json", "w") as f:
    json.dump(template_model_to_petrinet_json(model), f, indent = 4)

# %%
# Stratify by vaccination status
model_vax = stratify(
    model, 
    key = "vaxstatus", 
    strata = ["u", "v"], 
    cartesian_control = True, 
    structure = [],
    concepts_to_stratify = ["S", "C"],
    # concepts_to_preserve = ["H", "R", "D"],
    params_to_stratify = ["a", "vaxCtoH"]
)
model_vax.annotations.name = "SCRHD Vax"

# %%
# Add vaccination dynamics
model_vax = model_vax.add_template(
    NaturalConversion(
        subject = model_vax.get_concepts_name_map()["S_u"],
        outcome = model_vax.get_concepts_name_map()["S_v"],
        name = "vaccination"
    ).with_mass_action_rate_law("f")
)

model_vax.parameters["f"] = Parameter(name = "f", value = 1.0)

# %%
GraphicalModel.for_jupyter(model_vax)

# %%
generate_summary_table(model_vax)

# %%
# # Need to remove two redundant templates
# # H -> R, H -> D
# model_vax.templates = model_vax.templates[:9] + model_vax.templates[11:]
# generate_summary_table(model_vax)

# Fixed

# %%
# Export to AMR
with open("./scenario2_vax.json", "w") as f:
    json.dump(template_model_to_petrinet_json(model_vax), f, indent = 4)

# %%
# Stratified by 2 age groups
model_vax_age = stratify(
    model_vax, 
    key = "age",
    strata = ["y", "o"],
    cartesian_control = True, 
    structure = [], 
    params_to_preserve = ["b", "vaxCtoH_0", "vaxCtoH_1"],
    concepts_to_stratify = None
)
model_vax_age.annotations.name = "SCRHD Vax Age"

# %%
generate_summary_table(model_vax_age)

# %%
GraphicalModel.for_jupyter(model_vax_age)

# %%
# Export to AMR
with open("./scenario2_vax_age.json", "w") as f:
    json.dump(template_model_to_petrinet_json(model_vax_age), f, indent = 4)

# %%
# Review the `vaxCtoH_*` and `ageCtoH_*` parameters
# 
# vaxCtoH_* <- ['vaxCtoH_0', 'vaxCtoH_1']
# ageCtoH_*` <- ['ageCtoH_0', 'ageCtoH_1', 'ageCtoH_2', 'ageCtoH_3']


# %%
# ## Investigate Collapsible Equation Generation

# %%
odesys = generate_odesys(model, latex_align = True)
print(odesys)

# %%
m, a, i, n = sympy.symbols('m a i n')
x, y, z = sympy.symbols("x y z", cls = sympy.Function)
z = sympy.summation((m * x(i) + a - y(i)) ** 2, (i, 0, 3))
print(sympy.latex(z))

x = sympy.IndexedBase("x")
y = sympy.IndexedBase("y")
z = sympy.Sum((m * x[i] + a - y[i]) ** 2, (i, 0, 3))
print(sympy.latex(z.doit()))


# %%
# def generate_matrices(model, parameter_name):

#     parameters = list(model.get_all_used_parameters())
#     root_parameters = {p.split("_")[0]: [] for p in parameters}
#     for p in parameters:
#         root_parameters[p.split("_")[0]].extend(p.split("_"))

#     root_matrices = {p: {k: [] for k in ("subject_outcome", "subject_controller(s)", "outcome_controller(s)")} for p in root_parameters.keys()}

#     for t in model.templates:
#         if len(set(root_parameters[parameter_name]).intersection(t.get_parameter_names())) > 0:

#             x = []
#             for k in ("subject", "outcome"):
#                 if hasattr(t, k):
#                     x.append(getattr(t, k).name)
#                 else:
#                     None

#             if hasattr(t, "subject") & hasattr(t, "outcome")



#     root_matrices[parameter_name]["subject_outcome"] = [
#         [t.subject.name, t.outcome.name, set(root_parameters[parameter_name]).intersection(t.get_parameter_names())] if hasattr(t, "subject") & hasattr(t, "outcome") 

#         for t in model.templates if len(set(root_parameters[parameter_name]).intersection(t.get_parameter_names())) > 0
#     ]

#     return root_matrices

# matrices = generate_matrices(model_vax, "a")

# %%