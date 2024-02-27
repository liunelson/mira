# %%[markdown]
# Author: Nelson Liu
# 
# Email: [nliu@uncharted.software](mailto:nliu@uncharted.software)

# %%[markdown]
# # Define Models for Epi Hackathon (2024-02)

# %%
import os
import json
import numpy
import pandas
import sympy
from typing import Optional, Iterable

from mira.metamodel import *
from mira.modeling import Model
from mira.modeling.amr.petrinet import template_model_to_petrinet_json
from mira.modeling.amr.regnet import template_model_to_regnet_json
from mira.sources.amr.regnet import template_model_from_amr_json
from mira.modeling.ode import OdeModel, simulate_ode_model

# %%
MIRA_REST_URL = 'http://34.230.33.149:8771/api'

# %%[markdown]
# N-species Lotke-Volterra model
# 
# For every species 
# X_i: d X_i / d t = r_i X_i + \sum_j a_ij X_i X_j, where j = 1...n
#
# Templates:
# * one production controlled by X_i
# * n production controlled by X_i and X_j

# %%
def GenerateLotkeVolterraModel(num_species: int = 2, config: dict = {}) -> TemplateModel:

    # Units
    cfu_units = lambda: Unit(expression = sympy.Symbol("CFU"))
    day_units = lambda: Unit(expression = sympy.Symbol("day"))
    cfu_per_day_units = lambda: Unit(expression = sympy.Symbol("CFU") / sympy.Symbol("day"))

    # Basic components
    concepts = {
        f"X_{i}": Concept(name = f"X_{i}", display_name = f"X_{i}", units = cfu_units(), description = f"Population of species {i} (CFU)") 
        for i in range(1, num_species + 1)
    }

    initials = {
        c: Initial(concept = concepts[c], expression = sympy.Float(0.5)) for c in concepts.keys()
    }

    parameters = {
        f"a_{i}_{j}": Parameter(name = f"a_{i}_{j}", display_name = f"a_{i}_{j}", value = 1.0, units = cfu_per_day_units(), description = f"Interaction matrix coefficient between species {i} and {j} (CFU/day)")
        for i in range(1, num_species + 1) for j in range(1, num_species + 1)
    }
    parameters = parameters | {
        f"r_{i}": Parameter(name = f"r_{i}", display_name = f"r_{i}", value = 1.0, units = cfu_per_day_units(), description = f"Growth rate of species {i} (CFU/day)")
        for i in range(1, num_species + 1)
    }

    # Initialize with config values if given
    if isinstance(config, dict):

        if "initials" in config.keys():
            if len(config["initials"]) >= num_species:
                for i, c in enumerate(concepts):
                    initials[c].expression = sympy.Float(config["initials"][i])

        if "parameters" in config.keys():
            if "r" in config["parameters"]:

                if len(config["parameters"]["r"]) >= num_species:
                    for i, c in enumerate(concepts):
                        r_i = f"r_{i + 1}"
                        parameters[r_i].value = config["parameters"]["r"][i]

            if "a" in config["parameters"]:
                a = numpy.array(config["parameters"]["a"])
                if numpy.shape(a) == (num_species, num_species):
                    for i in range(1, num_species + 1):
                        for j in range(1, num_species + 1):
                            a_ij = f"a_{i}_{j}"
                            parameters[a_ij].value = a[i - 1, j - 1]

    # Templates
    templates = []
    for i in range(1, num_species + 1):

        X_i = f"X_{i}"
        r_i = f"r_{i}"

        templates.append(ControlledProduction(
            outcome = concepts[X_i],
            controller = concepts[X_i],
            rate_law = sympy.parsing.sympy_parser.parse_expr(f"{r_i} * {X_i}"),
            name = f"ProductionOf{X_i}ControlledBy{X_i}"
        ))

        for j in range(1, num_species + 1):

            X_j = f"X_{j}"
            a_ij = f"a_{i}_{j}"

            templates.append(GroupedControlledProduction(
                outcome = concepts[X_i],
                controllers = [concepts[X_i], concepts[X_j]],
                rate_law = sympy.parsing.sympy_parser.parse_expr(f"{a_ij} * {X_i} * {X_j}"),
                name = f"ProductionOf{X_i}ControlledBy{X_i}And{X_j}"
            ))

    model = TemplateModel(
        templates = templates,
        parameters = parameters,
        initials = initials,
        observables = {},
        time = Time(name = "t", units = day_units()),
        annotations = Annotations(name = f"LotkeVolterraModelWith{num_species}Species")
    )

    return model

# %%
models = {}
for n in [2, 4, 6]:

    # Specify configuration from docs
    config = {}
    if n == 4:
        config = {
            "initials": [0.51, 0.39, 0.88, 0.4],
            "parameters": {
                "r": [0.53, 0.42, 0.49, 0.33],
                "a": [
                    [-0.5, -0.01, 0.002, -0.009], 
                    [0, -0.5, 0, -0.169], 
                    [-0.002, -0.003, -0.5, 0.02], 
                    [0, -0.226, -0.04, -0.5]
                ]
            }
        }

    if n == 6:
        config = {
            "initials": [0.51, 0.39, 0.88, 0.4, 0.2, 0.8],
            "parameters": {
                "r": [0.53, 0.42, 0.49, 0.33, 0.7, 0.3],
                "a": [
                    [-0.5, -0.01, 0.02, -0.009, -0.002, 0.01], 
                    [0, -0.5, 0, -0.169, 0, 0], 
                    [-0.002, -0.003, -0.5, 0.02, 0.03, -0.04],
                    [0, -0.226, -0.04, -0.5, 0, 0.01], 
                    [0, -0.1, -0.02, 0, -0.5, 0], 
                    [0, -0.04, -0.05, 0, 0, -0.5]
                ]
            }
        }

    if len(config) > 0:    
        __ = pandas.DataFrame(config["initials"], index = [f"X_{i}" for i in range(1, n + 1)], columns = ["Initials"]).to_csv(f"./data/milestone_18_hackathon/scenario_4_LV{n}_initials.csv")
        __ = pandas.DataFrame(config["parameters"]["r"], index = [f"X_{i}" for i in range(1, n + 1)], columns = ["GrowthRate"]).to_csv(f"./data/milestone_18_hackathon/scenario_4_LV{n}_r.csv")
        __ = pandas.DataFrame(config["parameters"]["a"], index = [f"X_{i}" for i in range(1, n + 1)], columns = [f"X_{i}" for i in range(1, n + 1)]).to_csv(f"./data/milestone_18_hackathon/scenario_4_LV{n}_a.csv")


    # Generate model
    models[n] = GenerateLotkeVolterraModel(n, config = config)

    # Plot model graph
    models[n].draw_graph(f"./data/milestone_18_hackathon/scenario_4_LV{n}.png")

    # Export to AMR JSON
    with open(f"./data/milestone_18_hackathon/scenario_4_LV{n}_regnet.json", "w") as f:
        j = template_model_to_regnet_json(models[n])
        json.dump(j, f, indent = 3)

    with open(f"./data/milestone_18_hackathon/scenario_4_LV{n}_petrinet.json", "w") as f:
        j = template_model_to_petrinet_json(models[n])
        json.dump(j, f, indent = 3) 

# %%[markdown]
# # Ben Gyori's Scenario 4 Models

# %%
MIRA_HACKATHON_PATH = "../hackathon_2024.02" 

for n in ("scenario4_4spec_regnet.json", "scenario4_6spec_regnet.json"):
    p = os.path.join(MIRA_HACKATHON_PATH, "scenario4", n)
    with open(p, "r") as f:
        models[n] = template_model_from_amr_json(json.load(f))

# %%
models_ode = {n: OdeModel(Model(models[n]), initialized=True) for n in models.keys()}

# %%
models_ode["scenario4_4spec_regnet.json"]

# %%
res = simulate_ode_model(models_ode["scenario4_4spec_regnet.json"], times = numpy.linspace(0, 30, 100))
plt.plot(res)
plt.ylim([0, 1.2])

# %%
