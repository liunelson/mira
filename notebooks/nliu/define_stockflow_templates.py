# %%[markdown]
# Author: Nelson Liu
# 
# Email: [nliu@uncharted.software](mailto:nliu@uncharted.software)

# %%[markdown]
# # Test MIRA Import and Export of StockFlow Models
# 
# MIRA supports ingestion of StockFlow models in the Vensim MDL and Stella XMILE file formats.
# * [https://github.com/liunelson/mira/tree/nliu/experiment/mira/sources/system_dynamics](https://github.com/liunelson/mira/tree/nliu/experiment/mira/sources/system_dynamics)
# * [https://github.com/liunelson/mira/blob/nliu/experiment/mira/modeling/amr/stockflow.py](https://github.com/liunelson/mira/blob/nliu/experiment/mira/modeling/amr/stockflow.py)
# 
# We get example models from the 
# * [SDXorg repository](https://exchange.iseesystems.com/directory/isee)
# * [ISEE Exchange repository](https://exchange.iseesystems.com/directory/isee)


# %%
import os
import glob
import json
import tqdm

from mira.metamodel import *
from mira.sources.system_dynamics.vensim import template_model_from_mdl_file
from mira.sources.system_dynamics.stella import template_model_from_stella_model_file
from mira.modeling.amr.stockflow import AMRStockFlowModel
from mira.modeling.amr.stockflow import template_model_to_stockflow_json

# %%
# MIRA_REST_URL = 'http://34.230.33.149:8771/api'
PATH = "data/stockflow_examples/SDXorg"

# %%[markdown]
# # Test Vensim MDL Models

# %%
models = []
for m in os.listdir(PATH):

    p = os.path.join(PATH, m, "*.mdl")
    fp = glob.glob(p)
    if len(fp) == 0:
        continue

    try:
        models[m]["tm"] = template_model_from_mdl_file(fp)
    except Exception as err:
        print(err)

    try:
        models[m]["amr"] = template_model_to_stockflow_json(models[m]["tm"])
    except Exception as er:
        print(err)


# %%
p = os.path.join(PATH, "SIR.mdl")
model = vensim.template_model_from_mdl_file(p)

model.draw_jupyter()

# %%
p = os.path.join(PATH, "SIR.xmile")
model = stella.template_model_from_stella_model_file(p)

model.draw_jupyter()

# %%
p = os.path.join(PATH, "roessler_chaos.mdl")
model = vensim.template_model_from_mdl_file(p)

model.draw_jupyter()

# %%
p = os.path.join(PATH, "workforce.mdl")
model = vensim.template_model_from_mdl_file(p)

model.draw_jupyter()

# %%
p = os.path.join(PATH, "COVID-19-Model.stmx")
model = stella.template_model_from_stella_model_file(p)

model.draw_jupyter()


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
