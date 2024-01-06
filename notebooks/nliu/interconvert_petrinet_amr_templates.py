# %%[markdown]
# Author: Nelson Liu
# 
# Email: [nliu@uncharted.software](mailto:nliu@uncharted.software)

# %%[markdown]
# # Interconvert between Flat Model and Decomposed Templates in PetriNet Framework
# 
# We want to convert a PetriNet model represented by a single AMR JSON to 
# multiple AMR JSONs, each containing a sub-model equivalent to a single MIRA template, 
# and vice versa.
# 
# [MIRA templates](https://github.com/gyorilab/mira/blob/main/mira/metamodel/templates.py):
# 1. natural conversion
# 2. natural production
# 3. natural degradation
# 4. controlled conversion
# 5. controlled production
# 6. controlled degradation
# 7. observable (not originally in MIRA)
#
# Example: 
# `SIR` = `ControlledConversion(S -> I, p = b, C = I)` + `NaturalConversion(I -> R, p = g)`


# %%
import os
import glob
import json
import sympy as sp
import tqdm
from typing import Optional

from mira.metamodel import *
from mira.metamodel import expression_to_mathml
from mira.modeling import Model
from mira.modeling.amr.petrinet import AMRPetriNetModel
from mira.sources.amr.petrinet import template_model_from_amr_json

# %%
PATH = "data/sir_templates_example"

# %%
# Define common objects

S, I, R, b, g, N = sp.symbols('S I R b g N')

day_units = lambda: Unit(expression = sp.Symbol("day"))
per_day_units = lambda: Unit(expression = 1 / sp.Symbol("day"))

concepts = {
    "S": Concept(name = "S", units = None, description = "SusceptiblePopulationVariable"),
    "I": Concept(name = "I", units = None, description = "InfectedPopulationVariable"),
    "R": Concept(name = "R", units = None, description = "RecoveredPopulationVariable")
}

initials = {
    c: Initial(concept = concepts[c], expression = sp.Float(1)) for c in concepts.keys()
}

parameters = {
    "b": Parameter(name = "b", value = 1.0, units = per_day_units(), description = "beta"),
    "g": Parameter(name = "g", value = 1.0, units = per_day_units(), description = "gamma"),
}

observables = {
    "N": Observable(name = "N", units = None, expression = S + I + R, description = "TotalPopulation")
}
# observables = {
#     "N": Observable(name = "N", units = None, expression = "S + I + R", description = "TotalPopulation")
# }

time = Time(name = "t", units = day_units())

# %%[markdown]
# ## Define SIR model

# %%
sir_tm = TemplateModel(
    templates = [
        ControlledConversion(
            subject = concepts["S"],
            outcome = concepts["I"],
            controller = concepts["I"],
            rate_law = b * S * I,
            name = "InfectionProcess"
        ),
        NaturalConversion(
            subject = concepts["I"],
            outcome = concepts["R"],
            rate_law = g * I,
            name = "RecoveryProcess"
        )
    ],
    parameters = parameters,
    initials = initials,
    observables = {"N": observables["N"]},
    time = time,
    annotations = Annotations(name = "SirModel")
)

sir_json = AMRPetriNetModel(Model(sir_tm)).to_json()


# Save as AMR JSON
with open(os.path.join(PATH, "sir.json"), "w") as f:
    json.dump(sir_json, f, indent = 4)

# %%[markdown]
# ## Convert a 'flat' model into templates

# %%
# Convert a model into a set of templates
def convert_model_into_templates(model_amr_path: str, save: Optional[bool] = True) -> tuple[dict, list[dict]]:

    """
    Convert a model in AMR repr into a list of models in AMR 
    that have only 1 MIRA TemplateModel inside.
    """

    # Load model AMR as MIRA TemplateModel
    with open(model_amr_path, "r") as f:
        model_json = json.load(f)
        model_tm = template_model_from_amr_json(model_json)

    # Define a new MIRA TemplateModel for each template in the model
    templates_tm = []
    templates_json = []
    for t in model_tm.templates:

        initials = {c.name: model_tm.initials[c.name] for c in t.get_concepts()}
        parameters = {p: model_tm.parameters[p] for p in t.get_parameter_names()}

        tm = TemplateModel(
            templates = [t],
            parameters = parameters,
            initials = {c.name: model_tm.initials[c.name] for c in t.get_concepts()},
            observables = {},
            time = model_tm.time,
            annotations = Annotations(name = t.name)
        )
        templates_tm.append(tm)
        templates_json.append(AMRPetriNetModel(Model(tm)).to_json())

    # Define a new MIRA TemplateModel for each observable in the model
    # concepts = model_tm.get_concepts_name_map()
    # staticconcepts = [str(symb): concepts[str(symb)] for symb in model_tm.observables["N"].expression.free_symbols]
    for obs_name, obs in model_tm.observables.items():
        tm = TemplateModel(
            templates = [
                StaticConcept(subject = model_tm.get_concepts_name_map()[str(symb)])
                for symb in obs.expression.free_symbols
            ],
            observables = {obs_name: obs},
            time = model_tm.time,
            annotations = Annotations(name = obs_name)
        )
        templates_tm.append(tm)
        templates_json.append(AMRPetriNetModel(Model(tm)).to_json())

    # Save templates as AMR
    if save == True:

        # Check if directory exists
        p = os.path.join(model_amr_path.rsplit("/", 1)[0], "templates")
        if not os.path.exists(p):
            os.makedirs(p)

        for j in templates_json:
            name = j["header"]["name"]
            with open(os.path.join(p, f"{name}.json"), "w") as f:
                json.dump(j, f, indent = 4)

    return model_json, templates_json

# %%
p = os.path.join(PATH, "sir.json")
sir_json, sir_templates_json = convert_model_into_templates(p, save = True)

# %%[markdown]
# ## Convert a set of templates into a flat model

# %%
def convert_templates_into_model(templates_amr_path: str, save: Optional[bool]) -> tuple[list[dict], dict]:

    """
    Convert a list of models in AMR repr into a single model in AMR repr 
    that is the concatenation of their MIRA TemplateModel templates.
    """
    
    # Load the list of templates in AMR repr as MIRA TemplateModel
    templates_json = []
    templates_tm = []
    fnames = glob.glob(os.path.join(templates_amr_path, "*.json"))
    for fname in fnames:
        with open(fname, "r") as f:
            j = json.load(f)
            templates_json.append(j)
            templates_tm.append(template_model_from_amr_json(j))

    # Define a new model from the union of all templates and observables in TemplateModel list
    model_tm = TemplateModel(
        templates = [t for tm in templates_tm for t in tm.templates],
        parameters = {p: v for tm in templates_tm for p, v in tm.parameters.items()},
        initials = {i: v for tm in templates_tm for i, v in tm.initials.items()},
        observables = {obs_name: obs for tm in templates_tm for (obs_name, obs) in tm.observables.items()},
        time = templates_tm[0].time,
        annotations = Annotations(name = " + ".join([tm.annotations.name for tm in templates_tm]))
    )

    if save == True:
        p = templates_amr_path.rsplit("/", 1)[0]
        with open(os.path.join(p, "model_templates.json"), "w") as f:
            j = AMRPetriNetModel(Model(model_tm)).to_json()
            json.dump(j, f, indent = 4)

    return (templates_tm, model_tm)

# %%
p = os.path.join(PATH, "templates")
templates_tm_, model_tm_ = convert_templates_into_model(p, save = True)

# %%
