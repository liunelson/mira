# %%
import os
import sympy
from sympy.abc import _clash
from sympy.parsing.latex import parse_latex
# pip install antlr4-python3-runtime==4.11
import json
import pandas

from mira.metamodel import *
from mira.modeling import Model
from mira.sources.amr.petrinet import template_model_from_amr_json
from mira.modeling.amr.petrinet import template_model_to_petrinet_json
from mira.modeling.amr.regnet import template_model_to_regnet_json
from mira.modeling.amr.regnet import AMRRegNetModel
from mira.sources.sympy_ode import template_model_from_sympy_odes

from mira.modeling.viz import GraphicalModel

# %%
# %%
def generate_odesys(model, latex: bool = False, latex_align: bool = False) -> list:

    odeterms = {var: 0 for var in model.get_concepts_name_map().keys()}

    for t in model.templates:
        if hasattr(t, "subject"):
            var = t.subject.name
            odeterms[var] -= t.rate_law.args[0]
        
        if hasattr(t, "outcome"):
            var = t.outcome.name
            odeterms[var] += t.rate_law.args[0]

    # Time
    symb = lambda x: sympy.Symbol(x)
    try:
        time = model.time.name
    except:
        time = "t"
    finally:
        t = symb(time)

    # Construct Sympy equations
    odesys = [
        sympy.Eq(sympy.diff(sympy.Function(var)(t), t), terms) 
        if latex == False
        else sympy.latex(sympy.Eq(sympy.diff(sympy.Function(var)(t), t), terms))
        for var, terms in odeterms.items()
    ]
    
    if (latex == True) & (latex_align == True):
        odesys = "\\begin{align*} \n    " + " \\\\ \n    ".join([eq.replace(" = ", " &= ") for eq in odesys]) + "\n\\end{align*}"
        # odesys = "\\begin{align*}     " + " \\\\    ".join([eq.replace(" = ", " &= ") for eq in odesys]) + "\\end{align*}"

    return odesys

def generate_summary_table(model):

    data = {"name": [t.name for t in model.templates]}
    data['type'] = [t.type for t in model.templates]
    for k in ("subject", "outcome", "controller", "controllers"):
        if k != "controllers":
            data[k] = [getattr(t, k).name if hasattr(t, k) else None for t in model.templates]
        else:
            data[k] = [[c.name for c in t.controllers] if hasattr(t, k) else None for t in model.templates]

    data["rate_law"] = [t.rate_law for t in model.templates]
    data["interactor_rate_law"] = [t.get_interactor_rate_law() for t in model.templates]

    df = pandas.DataFrame(data)

    return df

# %%
odes_latex = [
    r"\frac{d S(t)}{d t} = -b * S(t) * I(t)",
    r"\frac{d I(t)}{d t} = b * S(t) * I(t) - g * I(t)", 
    r"\frac{d R(t)}{d t} = g * I(t)",
]

odes_sympy = [sympy.parsing.latex.parse_latex(ode) for ode in odes_latex]

__ = [print(ode) for ode in odes_sympy]

model = template_model_from_sympy_odes(odes_sympy)

generate_summary_table(model)

# %%
concepts_name_map = model.get_concepts_name_map()
initials = {}
controller_concept_list = []

# subject_concept = concepts_name_map.get("S")
# outcome_concept = concepts_name_map.get("I")
# controller_concept_list = [concepts_name_map.get(c) for c in ["S", "I", "R"]]

subject_concept = Concept(name = "Susceptible")
outcome_concept = Concept(name = "Infected")
controller_concept_list = [Concept(name = "Susceptible"), Concept(name = "Infected"), Concept(name = "Recovered")]


model1 = model.add_template(
    template = GroupedControlledConversion(
        subject = subject_concept,
        outcome = outcome_concept,
        controllers = controller_concept_list,
        # rate_law = safe_parse_expr("b * S * I / (S + I + R)", local_dict = _clash),
        rate_law = safe_parse_expr("b * Susceptible * Infected / (Susceptible + Infected + Recovered)", local_dict = _clash),
        name = "test"
    )
)

generate_summary_table(model1)

# %%
template_model_to_petrinet_json(model1)

# %%
concepts_name_map = model.get_concepts_name_map()
initials = {}
controller_concept_list = []

if "WastewaterViralLoad" in concepts_name_map:
    outcome_concept = concepts_name_map.get("WastewaterViralLoad")
else:
    outcome_concept = Concept(name = "WastewaterViralLoad")
    initials["WastewaterViralLoad"] = Initial(concept = outcome_concept, expression = sympy.Float(1))

for controller_name in ['I', 'H']:
    if controller_name in concepts_name_map:
        controller_concept_list.append(concepts_name_map.get(controller_name))
    else:
        controller_concept = Concept(name = controller_name)
        controller_concept_list.append(controller_concept)
        initials[controller_name] = Initial(concept = controller_concept, expression = sympy.Float(1))

# Define parameters
parameters = {}
if "shed_rate" in model.parameters: #note this is checks for paremeter's symbol
    parameters["shed_rate"] = model.parameters.get("shed_rate")
else: 
    parameters["shed_rate"] = Parameter(name = "shed_rate", value = 1, description = "")


# Add process as new template to the model
model2 = model.add_template(
    template = GroupedControlledProduction(
        outcome = outcome_concept,
        controllers = controller_concept_list,
        rate_law = safe_parse_expr("shed_rate * (I + H)", local_dict = _clash),
        name = "Production of WastewaterViralLoad"
    ),
    parameter_mapping = parameters,
    initial_mapping = initials
)

generate_summary_table(model2)

# %%
GraphicalModel.for_jupyter(model2)

# %%
template_model_to_petrinet_json(model2)

# %%