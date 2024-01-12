# %%[markdown]
# Author: Nelson Liu
# 
# Email: [nliu@uncharted.software](mailto:nliu@uncharted.software)

# %%[markdown]
# # Define AMR of Template Models for PetriNet Framework
# 
# We define template models as the "basis vectors" of the PetriNet framework: 
# any PetriNet model should be a linear combination of such template models. 
# 
# Per the [MIRA package](https://github.com/gyorilab/mira/blob/main/mira/metamodel/templates.py):
# 1. natural conversion
# 2. natural production
# 3. natural degradation
# 4. controlled conversion
# 5. controlled production
# 6. controlled degradation
# 7. observable (not originally in MIRA)

# %%
import os
import json
import sympy as sp
import tqdm

from mira.metamodel import *
from mira.modeling import Model
from mira.modeling.amr.petrinet import AMRPetriNetModel

# %%
PATH = "data/petrinet_templates"

# %%
# Define common objects

day_units = lambda: Unit(expression = sp.Symbol("day"))
per_day_units = lambda: Unit(expression = 1 / sp.Symbol("day"))

concepts = {
    "A": Concept(name = "A", units = None, description = "Variable named 'A'"),
    "B": Concept(name = "B", units = None, description = "Variable named 'B'"),
    "C": Concept(name = "C", units = None, description = "Variable named 'C'")
}

initials = {c: Initial(concept = concepts[c], expression = sp.Float(1)) for c in concepts.keys()}

parameters = {
    "p": Parameter(name = "p", value = 1.0, units = per_day_units(), description = "Parameter named 'p'")
}

observables = {
    "O": Observable(name = "O", units = None, expression = sp.parsing.sympy_parser.parse_expr("0"), description = "Observable named 'O'")
}

time = Time(name = "t", units = day_units())

models = {}

# %%
models["NaturalConversion"] = TemplateModel(
    templates = [NaturalConversion(
        subject = concepts["A"],
        outcome = concepts["B"],
        rate_law = sp.parsing.sympy_parser.parse_expr("p * A * B"),
        name = "NaturalConversion"
    )],
    parameters = {"p": parameters["p"]},
    initials = {"A": initials["A"], "B": initials["B"]},
    observables = {},
    time = time,
    annotations = Annotations(name = "NaturalConversionTemplateModel")
)

models["NaturalProduction"] = TemplateModel(
    templates = [NaturalProduction(
        outcome = concepts["A"],
        rate_law = sp.parsing.sympy_parser.parse_expr("p"),
        name = "NaturalProduction"
    )],
    parameters = {"p": parameters["p"]},
    initials = {"A": initials["A"]},
    observables = {},
    time = time,
    annotations = Annotations(name = "NaturalProductionTemplateModel")
)


models["NaturalDegradation"] = TemplateModel(
    templates = [NaturalDegradation(
        subject = concepts["A"],
        rate_law = sp.parsing.sympy_parser.parse_expr("p * A"),
        name = "NaturalDegradation"
    )],
    parameters = {"p": parameters["p"]},
    initials = {"A": initials["A"]},
    observables = {},
    time = time,
    annotations = Annotations(name = "NaturalDegradationTemplateModel")
)


models["ControlledConversion"] = TemplateModel(
    templates = [ControlledConversion(
        subject = concepts["A"],
        outcome = concepts["B"],
        controller = concepts["C"],
        rate_law = sp.parsing.sympy_parser.parse_expr("p * A * B"),
        name = "ControlledConversion"
    )],
    parameters = {"p": parameters["p"]},
    initials = {"A": initials["A"], "B": initials["B"], "C": initials["C"]},
    observables = {},
    time = time,
    annotations = Annotations(name = "ControlledConversionTemplateModel")
)


models["ControlledProduction"] = TemplateModel(
    templates = [ControlledProduction(
        outcome = concepts["A"],
        controller = concepts["C"],
        rate_law = sp.parsing.sympy_parser.parse_expr("p * C"),
        name = "ControlledProduction"
    )],
    parameters = {"p": parameters["p"]},
    initials = {"A": initials["A"], "C": initials["C"]},
    observables = {},
    time = time,
    annotations = Annotations(name = "ControlledProductionTemplateModel")
)

models["ControlledDegradation"] = TemplateModel(
    templates = [ControlledDegradation(
        subject = concepts["A"],
        controller = concepts["C"],
        rate_law = sp.parsing.sympy_parser.parse_expr("p * A * C"),
        name = "ControlledDegradation"
    )],
    parameters = {"p": parameters["p"]},
    initials = {"A": initials["A"], "C": initials["C"]},
    observables = {},
    time = time,
    annotations = Annotations(name = "ControlledDegradationTemplateModel")
)


models["Observable"] = TemplateModel(
    templates = [],
    parameters = {},
    initials = {},
    observables = {"O": observables["O"]},
    annotations = Annotations(name = "ObservableTemplateModel")
)

# %%
# Save as AMR JSON
for k, tm in tqdm.tqdm(models.items()):
    with open(os.path.join(PATH, f"{k}.json"), "w") as f:
        j = AMRPetriNetModel(Model(tm)).to_json()
        json.dump(j, f, indent = 4)

# %%
