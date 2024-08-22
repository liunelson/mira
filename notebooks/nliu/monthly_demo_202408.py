# %%[markdown]
# Author: Nelson Liu
# 
# Email: [nliu@uncharted.software](mailto:nliu@uncharted.software)

# %%[markdown]
# # Monthly Demo (2024-08)

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
# ## Question 1
# 
# 1. Start with published SIDARTHE model
# 2. Do Unit Test 1: configure it as per the paper and compare the results
# 3. Implement 

# %%
person_units = lambda: Unit(expression = sympy.Symbol('person'))
day_units = lambda: Unit(expression = sympy.Symbol('day'))
per_day_units = lambda: Unit(expression = 1/sympy.Symbol('day'))
dimensionless_units = lambda: Unit(expression = sympy.Integer('1'))
per_day_per_person_units = lambda: Unit(expression = 1/(sympy.Symbol('day')*sympy.Symbol('person')))

# Symbols
s = lambda x: sympy.Symbol(x)

# Basic components
concepts = {
    c: Concept(name = c, display_name = c, description = f"Population of {d} agents", units = person_units())
    for c, d in zip(["S", "V", "C", "R", "D", "H"], ["susceptible", "vaccinated", "infected", "recovered", "deceased", "hospitalized"])
}

# Groundings
concepts["S"].identifiers = {'ido': '0000514'}
concepts["V"].identifiers = {'vo': '0001376'}
concepts["C"].identifiers = {'ido': '0000511'}
concepts["R"].identifiers = {'ido': '0000592'}
concepts["D"].identifiers = {'ncit': 'C28554'}
concepts["H"].identifiers = {'ncit': 'C25179'}

# Initial conditions
initials = {
    c: Initial(concept = concepts[c], expression = sympy.Float(0.0)) for c in concepts.keys()
}

# Parameters
parameters = {
    p: Parameter(name = p, display_name = p, value = 0.0, units = per_day_units())
    for p in ("a", "b", "c", "d", "f", "g")
}
parameters["a"].description = "Infection rate"
parameters["b"].description = "Hospitalization rate"
parameters["c"].description = "Recovery rate"
parameters["d"].description = "Hospitalized recovery rate"
parameters["f"].description = "Hospitalized death rate"
parameters["g"].description = "Vaccination rate"


# Observables
observables = {
    "TotalPop": Observable(name = "TotalPop", expression = s("S") + s("V") + s("C") + s("R") + s("D") + s("H")),
    "NonInfectable": Observable(name = "NonInfectable", expression = s("V") + s("R") + s("D"))

}

# %%
model = TemplateModel(
    templates = [
        ControlledConversion(
            subject = concepts["S"],
            outcome = concepts["C"],
            controller = concepts["C"],
            name = "Infection",
        ).with_mass_action_rate_law("a"),
        NaturalConversion(
            subject = concepts["C"],
            outcome = concepts["H"],
            name = "Hospitalization",
        ).with_mass_action_rate_law("b"),
        NaturalConversion(
            subject = concepts["C"],
            outcome = concepts["R"],
            name = "Recovery"
        ).with_mass_action_rate_law("c"),
        NaturalConversion(
            subject = concepts["H"],
            outcome = concepts["R"],
            name = "HospitalRecovery"
        ).with_mass_action_rate_law("d"),
        NaturalConversion(
            subject = concepts["H"],
            outcome = concepts["D"],
            name = "HospitalDeath"
        ).with_mass_action_rate_law("f"),
        NaturalConversion(
            subject = concepts["S"],
            outcome = concepts["V"],
            name = "Vaccination",
        ).with_mass_action_rate_law("g"),
    ],
    initials = initials,
    parameters = parameters,
    observables = observables,
    time = Time(name = "t", units = day_units()),
    annotations = Annotations(
        name = "SVCRDH", 
        description = "Compartmental model with vaccination and hospitalization for 6-month Milestone Epi Eval Scenario 2.",
        locations = ["geonames:5368381"],
        diseases = ["doid:0080600"],
        hosts = ["ncbitaxon:9606"],
        pathogens = ["ncbitaxon:2697049"],
        model_types = ["mamo:0000028", "mamo:0000046"]
    )
)

# %%
# model.draw_jupyter()
GraphicalModel.for_jupyter(model)

# %%
with open("./monthlydemo202407Q1.json", "w") as f:
    json.dump(template_model_to_petrinet_json(model), f, indent = 4)
    
# %%
