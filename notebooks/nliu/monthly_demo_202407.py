# %%[markdown]
# Author: Nelson Liu
# 
# Email: [nliu@uncharted.software](mailto:nliu@uncharted.software)

# %%[markdown]
# # Monthly Demo (2024-07)

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

# %%[markdown]
# ## Q1
# 
# Define a model with vaccination by assuming vaccinated pop do not get infected

# %%
day_units = lambda: Unit(expression = sympy.Symbol("day"))
per_day_units = lambda: Unit(expression = 1 / sympy.Symbol("day"))

# Basic components
concepts = {
    c: Concept(name = c, display_name = c, description = f"Population of {d} agents")
    for c, d in zip(["S", "V", "C", "R", "D", "H"], ["susceptible", "vaccinated", "infected", "recovered", "deceased", "hospitalized"])
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
            name = "StoCviaC",
        ).with_mass_action_rate_law("a"),
        NaturalConversion(
            subject = concepts["C"],
            outcome = concepts["H"],
            name = "CtoH",
        ).with_mass_action_rate_law("b"),
        NaturalConversion(
            subject = concepts["C"],
            outcome = concepts["R"],
            name = "CtoR"
        ).with_mass_action_rate_law("c"),
        NaturalConversion(
            subject = concepts["H"],
            outcome = concepts["R"],
            name = "HtoR"
        ).with_mass_action_rate_law("d"),
        NaturalConversion(
            subject = concepts["H"],
            outcome = concepts["D"],
            name = "HtoD"
        ).with_mass_action_rate_law("e"),
        NaturalConversion(
            subject = concepts["S"],
            outcome = concepts["V"],
            name = "StoV",
        ).with_mass_action_rate_law("f"),
    ],
    initials = initials,
    observables = {},
    time = Time(name = "t", units = day_units()),
    annotations = Annotations(name = f"SVCRHD")
)

# %%
# model.draw_jupyter()
GraphicalModel.for_jupyter(model)

# %%
with open("./monthlydemo202407Q1.json", "w") as f:
    json.dump(template_model_to_petrinet_json(model), f, indent = 4)
# %%
