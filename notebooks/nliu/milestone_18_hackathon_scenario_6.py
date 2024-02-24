# %%[markdown]
# Author: Nelson Liu
# 
# Email: [nliu@uncharted.software](mailto:nliu@uncharted.software)

# %%[markdown]
# # Define Models for Milestone 18 Hackathon
#
# Scenario 6:
# RegNets for gene expression and regulation

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
from mira.modeling.amr.regnet import AMRRegNetModel

# %%
MIRA_REST_URL = 'http://34.230.33.149:8771/api'

# %%[markdown]
# Chen Model
# 
# d r_i / d t = C_i p_1 - V_i r_i
# d p_i / d t = L_i r_i - U_i p_i
# for i = 1, 2, 3
# 
# r_i -> ControlledProduction + ControlledDegradation
# p_i -> ControlledProduction + ControlledDegradation

# %%
def GenerateChenModel(config: dict = {}) -> TemplateModel:

  num_types = 3

  # Units
  nM_units = lambda: Unit(expression = sympy.Symbol("nM"))
  min_units = lambda: Unit(expression = sympy.Symbol("min"))
  per_min_units = lambda: Unit(expression = 1 / sympy.Symbol("min"))
  nM_per_min_units = lambda: Unit(expression = sympy.Symbol("nM") / sympy.Symbol("min"))

  # Basic components
  concepts = {}
  initials = {}
  for s, m in zip(("r", "p"), ("mRNA", "protein")):

    concepts = concepts | {
      f"{s}_{i}": Concept(name = f"{s}_{i}", display_name = f"{s}_{i}", units = nM_units(), description = f"Concentration of Type-{i} {m} in units of nM)")
      for i in range(1, num_types + 1)
    }

    initials = initials | {
      c: Initial(concept = concepts[c], expression = sympy.Float(0.5)) for c in concepts.keys()
    }

  parameters = {}
  for s, m, n, u in zip(("C", "L", "V", "U"), ("mRNA", "protein", "mRNA", "protein"), ("Transcription", "Translation", "Degradation", "Degradation"), (per_min_units(), per_min_units(), per_min_units(), per_min_units())):
    parameters = parameters | {
      f"{s}_{i}": Parameter(name = f"{s}_{i}", display_name = f"{s}_{i}", value = 1.0, units = u, description = f"{n} rate for type-{i} proteins")
      for i in range(1, num_types + 1)
    }

  # Initialize with config values if given
  if isinstance(config, dict):

    if "initials" in config.keys():
      for c in concepts:
        v, i = c.split("_")
        initials[c].expression = sympy.Float(config["initials"][v][int(i) - 1])

    if "parameters" in config.keys():
      for p in parameters:
        v, i = p.split("_")
        parameters[p].value = config["parameters"][v][int(i) - 1]

  # Templates
  templates = []
  for i in range(1, num_types + 1):
    
    r_i = f"r_{i}"
    p_1 = "p_1"
    p_i = f"p_{i}"
    C_i = f"C_{i}"
    L_i = f"L_{i}"
    V_i = f"V_{i}"
    U_i = f"U_{i}"

    # transcription of type-i mRNA
    templates.append(ControlledProduction(
      outcome = concepts[r_i],
      controller = concepts[p_1],
      rate_law = sympy.parsing.sympy_parser.parse_expr(f"{C_i} * {p_1}"),
      name = f"TranscriptionOf{r_i}ControlledBy{p_1}"
    ))

    # degradation of type-i mRNA
    templates.append(ControlledDegradation(
      subject = concepts[r_i],
      controller = concepts[r_i],
      rate_law = sympy.parsing.sympy_parser.parse_expr(f"{V_i} * {r_i}"),
      name = f"DegradationOf{r_i}"
    ))

    # translation of type-i protein
    templates.append(ControlledProduction(
      outcome = concepts[p_i],
      controller = concepts[r_i],
      rate_law = sympy.parsing.sympy_parser.parse_expr(f"{L_i} * {r_i}"),
      name = f"TranslationOf{p_i}ControlledBy{r_i}"
    ))

    # degradation of type-i protein
    templates.append(ControlledDegradation(
      subject = concepts[p_i],
      controller = concepts[p_i],
      rate_law = sympy.parsing.sympy_parser.parse_expr(f"{U_i} * {p_i}"),
      name = f"DegradationOf{p_i}"
    ))

  # Generate model
  model = TemplateModel(
    templates = templates,
    parameters = parameters,
    initials = initials,
    observables = {},
    time = Time(name = "t", units = min_units()),
    annotations = Annotations(name = f"ChenModelWith{num_types}Types")
  )

  return model

# %%
# Model configuration (Chen and Hunt)
config = {
  "initials": {
    "r": [3, 6, 5],
    "p": [100, 500, 1]
  },
  "parameters": {
    "C": [0.03, 0.03, 0.024],
    "L": [2, 2, 2],
    "V": [0.03, 0.03, 0.03],
    "U": [0.15, 0.15, 0.015],
    "a": [60, 140, 170],
    "b": [120, 140, 180],
    "d": [120, 150, 260]
  }
}
__ = pandas.DataFrame(list(config["initials"].values()), columns = [1, 2, 3], index = ["r_i (concentration of type-i mRNA in units of nM)", "p_i (concentration of type-i protein in units of nM)"]).transpose().to_csv(f"./data/milestone_18_hackathon/scenario_6/initials.csv")
__ = pandas.DataFrame(list(config["parameters"].values()), columns = [1, 2, 3], index = list(config["parameters"].keys())).transpose().to_csv(f"./data/milestone_18_hackathon/scenario_6/params.csv")

# %%
ChenModel = GenerateChenModel(config = config)

ChenModel.draw_jupyter("./data/milestone_18_hackathon/scenario_6/ChenModel.png")

# %%
with open("./data/milestone_18_hackathon/scenario_6/ChenModel_regnet_amr.json", "w") as f:
  model_amr = template_model_to_regnet_json(ChenModel)
  json.dump(model_amr, f, indent = 3)

with open("./data/milestone_18_hackathon/scenario_6/ChenModel_petrinet_amr.json", "w") as f:
  model_amr = template_model_to_petrinet_json(ChenModel)
  json.dump(model_amr, f, indent = 3)

# %%[markdown]
# Hunt Model
# 
# d r_1 / d t = (1 / (1 + p_1 ** 2 / a_1 ** 2)) C_1 p_1 - (1 / (1 + p_2 / b_1)) V_1 r_1
# 
# d r_i / d t = (1 / (1 + p_i / a_i)) C_i p_1 - (1 / (1 + p_2 / b_i)) V_i r_i
# 
# d p_i / d t = (1 / (1 + p_i / d_i)) L_i r_i - U_i p_i
# 
# for i = 1, 2, 3
# 
# r_i -> GroupedControlledProduction + GroupedControlledDegradation
# p_i -> GroupedControlledProduction + GroupedControlledDegradation
  
# %%
def GenerateHuntModel(config: dict = {}) -> TemplateModel:

  num_types = 3

  # Units
  nM_units = lambda: Unit(expression = sympy.Symbol("nM"))
  min_units = lambda: Unit(expression = sympy.Symbol("min"))
  per_min_units = lambda: Unit(expression = 1 / sympy.Symbol("min"))
  nM_per_min_units = lambda: Unit(expression = sympy.Symbol("nM") / sympy.Symbol("min"))

  # Basic components
  concepts = {}
  initials = {}
  for s, m in zip(("r", "p"), ("mRNA", "protein")):

    concepts = concepts | {
      f"{s}_{i}": Concept(name = f"{s}_{i}", display_name = f"{s}_{i}", units = nM_units(), description = f"Concentration of Type-{i} {m} in units of nM)")
      for i in range(1, num_types + 1)
    }

    initials = initials | {
      c: Initial(concept = concepts[c], expression = sympy.Float(0.5)) for c in concepts.keys()
    }

  parameters = {}
  for s, m, n, u in zip(("C", "L", "V", "U", "a", "b", "d"), ("mRNA", "protein", "mRNA", "protein", "mRNA", "mRNA", "protein"), ("Transcription", "Translation", "Degradation", "Degradation", "Effectiveness factor for feedback on transcription", "Effectiveness factor for feedback on degradation", "Effectiveness factor for feedback on translation"), (per_min_units(), per_min_units(), per_min_units(), per_min_units(), nM_units(), nM_units(), nM_units())):
    parameters = parameters | {
      f"{s}_{i}": Parameter(name = f"{s}_{i}", display_name = f"{s}_{i}", value = 1.0, units = u, description = f"{n} rate for type-{i} proteins")
      for i in range(1, num_types + 1)
    }

  # Initialize with config values if given
  if isinstance(config, dict):

    if "initials" in config.keys():
      for c in concepts:
        v, i = c.split("_")
        initials[c].expression = sympy.Float(config["initials"][v][int(i) - 1])

    if "parameters" in config.keys():
      for p in parameters:
        v, i = p.split("_")
        parameters[p].value = config["parameters"][v][int(i) - 1]

  # Templates
  templates = []
  for i in range(1, num_types + 1):
    
    r_i = f"r_{i}"
    p_1 = "p_1"
    p_2 = "p_2"
    p_i = f"p_{i}"
    C_i = f"C_{i}"
    L_i = f"L_{i}"
    V_i = f"V_{i}"
    U_i = f"U_{i}"
    a_i = f"a_{i}"
    b_i = f"b_{i}"
    d_i = f"d_{i}"

    if i == 1:

      # transcription of type-1 mRNA
      templates.append(GroupedControlledProduction(
        outcome = concepts[r_i],
        controllers = [concepts[p_1], concepts[p_i]],
        rate_law = sympy.parsing.sympy_parser.parse_expr(f"(1 / (1 + ({p_i})**2 / ({a_i})**2)) * {C_i} * {p_1}"),
        name = f"TranscriptionOf{r_i}ControlledBy{p_1}"
      ))

    else:

      # transcription of type-i mRNA, i = 2, 3
      templates.append(GroupedControlledProduction(
        outcome = concepts[r_i],
        controllers = [concepts[p_1], concepts[p_i]],
        rate_law = sympy.parsing.sympy_parser.parse_expr(f"(1 / (1 + {p_i} / {a_i})) * {C_i} * {p_1}"),
        name = f"TranscriptionOf{r_i}ControlledBy{p_1}"
      ))

    # degradation of type-i mRNA
    templates.append(GroupedControlledDegradation(
      subject = concepts[r_i],
      controllers = [concepts[r_i], concepts[p_2]],
      rate_law = sympy.parsing.sympy_parser.parse_expr(f"(1 / (1 + {p_2} / {b_i})) * {V_i} * {r_i}"),
      name = f"DegradationOf{r_i}"
    ))


    # translation of type-i protein
    templates.append(GroupedControlledProduction(
      outcome = concepts[p_i],
      controllers = [concepts[r_i], concepts[p_i]],
      rate_law = sympy.parsing.sympy_parser.parse_expr(f"(1 / (1 + {p_i} / {d_i})) * {L_i} * {r_i}"),
      name = f"TranslationOf{p_i}ControlledBy{r_i}And{p_i}"
    ))

    # degradation of type-i protein
    templates.append(ControlledDegradation(
      subject = concepts[p_i],
      controller = concepts[p_i],
      rate_law = sympy.parsing.sympy_parser.parse_expr(f"{U_i} * {p_i}"),
      name = f"DegradationOf{p_i}"
    ))

  # Generate model
  model = TemplateModel(
    templates = templates,
    parameters = parameters,
    initials = initials,
    observables = {},
    time = Time(name = "t", units = min_units()),
    annotations = Annotations(name = f"HuntModelWith{num_types}Types")
  )

  return model

# %%
HuntModel = GenerateHuntModel()

HuntModel.draw_jupyter("./data/milestone_18_hackathon/scenario_6/HuntModel.png")

# %%
with open("./data/milestone_18_hackathon/scenario_6/HuntModel_regnet_amr.json", "w") as f:
  model_amr = template_model_to_regnet_json(HuntModel)
  json.dump(model_amr, f, indent = 3)

# %%
with open("./data/milestone_18_hackathon/scenario_6/HuntModel_petrinet_amr.json", "w") as f:
  model_amr = template_model_to_petrinet_json(HuntModel)
  json.dump(model_amr, f, indent = 3)

# %%
