# %%[markdown]
# Author: Nelson Liu
# 
# Email: [nliu@uncharted.software](mailto:nliu@uncharted.software)

# %%[markdown]
# # Define Rename-State Function
# 

# %%
import os
import glob
import json
import sympy
import tqdm
from typing import Optional
import copy
import warnings

from mira.metamodel import *
from mira.metamodel import expression_to_mathml
from mira.modeling import Model
from mira.modeling.amr.petrinet import AMRPetriNetModel
from mira.sources.amr.petrinet import template_model_from_amr_json

# %%
PATH = "data/sir_templates_example"

# %%
# Define common objects

S, I, R, b, g, N = sympy.symbols('S I R b g N')

day_units = lambda: Unit(expression = sympy.Symbol("day"))
per_day_units = lambda: Unit(expression = 1 / sympy.Symbol("day"))

concepts = {
    "S": Concept(name = "S", units = None, description = "SusceptiblePopulationVariable"),
    "I": Concept(name = "I", units = None, description = "InfectedPopulationVariable"),
    "R": Concept(name = "R", units = None, description = "RecoveredPopulationVariable")
}

initials = {
    c: Initial(concept = concepts[c], expression = sympy.Float(1)) for c in concepts.keys()
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

# %%
# Add a new template "A -> B" with rate "r * A"
A, B, r = sympy.symbols("A B r")

sirab_tm = sir_tm.add_template(
    template = NaturalConversion(
        subject = Concept(name = "A", display_name = "A"), 
        outcome = Concept(name = "B", display_name = "B"), 
        rate_law = r * A, 
        name = "A->B"
    ),
    parameter_mapping = {
        "r": Parameter(name = "r", value = 1.0)
    }, 
    initial_mapping = {
        # "A": Initial(concept = Concept(name = "A"), expression = sympy.Float(1)), 
        "B": Initial(concept = Concept(name = "B"), expression = sympy.Float(1))
    }
)

# %%
def rename_state(tm: TemplateModel, old_name: str, new_name: str) -> TemplateModel:

    """
    Given a TemplateModel model and the old & new name of a state/concept therein, 
    1. if the old name doesn't exist, do nothing
    2. if the new name is unused, replace the name of that state/concept
    3. if the new name is used, replace the old-name concept by the new-name concept
    
    Apply same replacement to the rate laws, observable expressions, and initials.
    """

    # Check if concept with old name exists
    concepts_name_map = tm.get_concepts_name_map()
    if old_name not in concepts_name_map:
        raise ValueError(f"State with name {old_name} not found in model.")
    
    # Check if concept with new name exists
    if new_name in concepts_name_map:
        new_concept = concepts_name_map[new_name]
        warnings.warn(f"State with name {new_name} exists already in model and will replace state {old_name} everywhere.", UserWarning)

    # Check if old name is already used by a parameter or observable
    if old_name in tm.observables.keys():
        raise ValueError(f"Name {old_name} already used by a model observable.")
    if old_name in tm.parameters.keys():
        raise ValueError(f"Name {old_name} already used by a model parameter.")
        
    # Rename name of concept
    for template in tm.templates:
        concept_keys = template.concept_keys
        if old_name in template.get_concept_names():
            for concept in template.get_concepts():
                if concept.name == old_name: 

                    if new_name not in concepts_name_map:
                        concept.name = new_name
                    else:
                        for role in template.concept_keys:
                            if getattr(template, role).name == old_name:
                                setattr(template, role, new_concept)

            template.rate_law = SympyExprStr(
                template.rate_law.args[0].subs(
                    sympy.Symbol(old_name), sympy.Symbol(new_name))
                )
    
    # Update observable expressions with new state name
    for observable in tm.observables.values():
        if sympy.Symbol(old_name) in observable.expression.free_symbols:
            observable.expression = SympyExprStr(
                observable.expression.args[0].subs(
                    sympy.Symbol(old_name), sympy.Symbol(new_name))
                )
    
    # Ditto for initials
    if (old_name in tm.initials) & (new_name not in tm.initials):
        tm.initials[new_name] = tm.initials.pop(old_name)
        tm.initials[new_name].concept.name = new_name
    if (old_name in tm.initials) & (new_name in tm.initials):
        __ = tm.initials.pop(old_name)
    if (new_name not in tm.initials):
        concept = tm.get_concepts_name_map()[new_name]
        tm.initials[new_name] = Initial(concept = concept, expression = sympy.Float(1))

    return tm

# %%
tm = copy.deepcopy(sirab_tm)

tm_ = rename_state(tm = tm, old_name = "A", new_name = "C")
# tm_ = rename_state(tm = tm_, old_name = "C", new_name = "S")

# tm_ = rename_state(tm = tm, old_name = "S", new_name = "C")

# tm_ = rename_state(tm = tm, old_name = "S", new_name = "R")

print(f"concepts: {list(tm_.get_concepts_name_map().keys())}")
print(f"rate laws: {[t.rate_law for t in tm_.templates]}")
print(f"initials: {list(tm_.initials.keys())}")

# %%


