# %%[markdown]
# Author: Nelson Liu
# 
# Email: [nliu@uncharted.software](mailto:nliu@uncharted.software)

# %%
# Example stratification workflow:
#
# 1. Start with SIR model
# 2. Stratify by 2 age groups
# 3. Add a D compartment
# 4. stratify by 2 locations
# 5. Add a Z compartment

# %%
import sympy
import itertools
from mira.metamodel import *
# from mira.metamodel import ControlledConversion, NaturalConversion, Concept, Template, TemplateModel
from mira.modeling import Model
# from mira.modeling.viz import GraphicalModel
from mira.modeling.askenet.petrinet import AskeNetPetriNetModel

import requests
from PIL import Image
from typing import NoReturn

# %%
MIRA_REST_URL = 'http://34.230.33.149:8771/api'

# %%
def viz_mmt(mmt, filename) -> NoReturn:

    res = requests.post(url = f'{MIRA_REST_URL}/viz/to_image', data = mmt.json())

    with open(filename, 'wb') as f:
        f.write(res.content)

    with Image.open(filename).convert('RGB') as im:
        im.show()

# %%
# Units
person_units = lambda: Unit(expression=sympy.Symbol('person'))
day_units = lambda: Unit(expression=sympy.Symbol('day'))
per_day_units = lambda: Unit(expression=1/sympy.Symbol('day'))
dimensionless_units = lambda: Unit(expression=sympy.Integer('1'))
per_day_per_person_units = lambda: Unit(expression=1/(sympy.Symbol('day')*sympy.Symbol('person')))

# Concepts
concepts = {
    'S': Concept(name = 'S', units = person_units(), identifiers = {'ido': '0000514'}),
    'I': Concept(name = 'I', units = person_units(), identifiers = {'ido': '0000511'}),
    'R': Concept(name = 'R', units = person_units(), identifiers = {'ido': '0000592'})
}

# Parameters
parameters = {
    'b': Parameter(name = 'b', value = 0.4, units = per_day_units()),
    'g': Parameter(name = 'g', value = 1.0/11.0, units = per_day_units()),
}

# Initial conditions
initials = {
    'S': Initial(concept = Concept(name = 'S'), value = 1_000 - 1),
    'I': Initial(concept = Concept(name = 'I'), value = 1),
    'R': Initial(concept = Concept(name = 'R'), value = 0)
}

# Symbols
S, I, R, b, g = sympy.symbols('S I R b g')

# Observables
observables = {'total_population': Observable(name = 'Total Population', expression = S + I + R)}

# %%
# S -> I transition
t1 = ControlledConversion(
    subject = concepts['S'],
    outcome = concepts['I'],
    controller = concepts['I'],
    rate_law = S * I * b
)

# I -> R transition
t2 = NaturalConversion(
    subject = concepts['I'], 
    outcome = concepts['R'],
    rate_law = g * I
)

# S -> I transition where the rate law is correct but not the template type
t1_incorrect = NaturalConversion(
    subject = concepts['S'],
    outcome = concepts['I'],
    rate_law = S * I * b
) 

# %%
# Build SIR model from templates
sir_model = TemplateModel(
    templates = [t1, t2],
    parameters = parameters,
    initials = initials,
    time = Time(name = 't', units = day_units()),
    observables = observables,
    annotations = Annotations(name = 'SIR model')
)

# Visualize
viz_mmt(sir_model, 'sir_model.png')

# Generate AMR JSON
AskeNetPetriNetModel(Model(sir_model)).to_json_file('sir_model.json')

# %%
# Build SIR model from the "incorrect" templates
sir_model_incorrect = TemplateModel(
    templates = [t1_incorrect, t2],
    parameters = parameters,
    initials = initials,
    time = Time(name = 't', units = day_units()),
    observables = observables,
    annotations = Annotations(name = 'SIR model with incorrect templates')
)

# Visualize
viz_mmt(sir_model_incorrect, 'sir_model_incorrect.png')

# Generate AMR JSON
AskeNetPetriNetModel(Model(sir_model_incorrect)).to_json_file('sir_model_incorrect.json')

# %%
# Do stratification with 2 age groups

sir_2age_model = stratify(
    sir_model,
    key = 'age',
    strata = ['A1', 'A2'],
    structure = [],
    directed = False,
    cartesian_control = True,
    params_to_stratify = {'b', 'g'},
    concepts_to_stratify = {'S', 'I', 'R'}
)

sir_2age_model.annotations.name = 'SIR model + stratified by 2 age groups'

viz_mmt(sir_2age_model, 'sir_2age_model.png')

AskeNetPetriNetModel(Model(sir_2age_model)).to_json_file('sir_2age_model.json')

# %%
# Repeat with "incorrect" model

sir_2age_model_incorrect = stratify(
    sir_model_incorrect,
    key = 'age',
    strata = ['A1', 'A2'],
    structure = [],
    directed = False,
    cartesian_control = True,
    params_to_stratify = {'b', 'g'},
    concepts_to_stratify = {'S', 'I', 'R'}
)

sir_2age_model_incorrect.annotations.name = 'Incorrect SIR model + stratified by 2 age groups'

viz_mmt(sir_2age_model_incorrect, 'sir_2age_model_incorrect.png')

AskeNetPetriNetModel(Model(sir_2age_model_incorrect)).to_json_file('sir_2age_model_incorrect.json')

# %%
# Compare the resulting ODE systems

sir_2age_model_amr = AskeNetPetriNetModel(Model(sir_2age_model)).to_json()
__ = [
    print(f'{rate["target"]}: {rate["expression"].replace("*", " ")}') 
    for rate in sir_2age_model_amr["semantics"]["ode"]["rates"]
]

print("\n")

sir_2age_model_incorrect_amr = AskeNetPetriNetModel(Model(sir_2age_model_incorrect)).to_json()
__ = [
    print(f'{rate["target"]}: {rate["expression"].replace("*", " ")}') 
    for rate in sir_2age_model_incorrect_amr["semantics"]["ode"]["rates"]
]

# t1: I_A1 S_A1 b_0
# t2: I_A2 S_A1 b_1
# t3: I_A2 S_A2 b_2
# t4: I_A1 S_A2 b_3
# t5: I_A1 g_0
# t6: I_A2 g_1
# 
# 
# t1: I_A1 S_A1 b_0
# t2: I_A2 S_A2 b_1
# t3: I_A1 g_0
# t4: I_A2 g_1

# %%
# %%
# Repeat with `cartesian_control = False`

sir_2age_model = stratify(
    sir_model,
    key = 'age',
    strata = ['A1', 'A2'],
    structure = [],
    directed = False,
    cartesian_control = False,
    params_to_stratify = {'b', 'g'},
    concepts_to_stratify = {'S', 'I', 'R'}
)

sir_2age_model.annotations.name = 'SIR model + stratified by 2 age groups'

viz_mmt(sir_2age_model, 'sir_2age_model.png')

# AskeNetPetriNetModel(Model(sir_2age_model)).to_json_file('sir_2age_model.json')

# Repeat with "incorrect" model
sir_2age_model_incorrect = stratify(
    sir_model_incorrect,
    key = 'age',
    strata = ['A1', 'A2'],
    structure = [],
    directed = False,
    cartesian_control = False,
    params_to_stratify = {'b', 'g'},
    concepts_to_stratify = {'S', 'I', 'R'}
)

sir_2age_model_incorrect.annotations.name = 'Incorrect SIR model + stratified by 2 age groups'

viz_mmt(sir_2age_model_incorrect, 'sir_2age_model_incorrect.png')

# AskeNetPetriNetModel(Model(sir_2age_model_incorrect)).to_json_file('sir_2age_model_incorrect.json')


# Compare the resulting ODE systems
sir_2age_model_amr = AskeNetPetriNetModel(Model(sir_2age_model)).to_json()
__ = [
    print(f'{rate["target"]}: {rate["expression"].replace("*", " ")}') 
    for rate in sir_2age_model_amr["semantics"]["ode"]["rates"]
]

print("\n")

sir_2age_model_incorrect_amr = AskeNetPetriNetModel(Model(sir_2age_model_incorrect)).to_json()
__ = [
    print(f'{rate["target"]}: {rate["expression"].replace("*", " ")}') 
    for rate in sir_2age_model_incorrect_amr["semantics"]["ode"]["rates"]
]

# t1: I_A1 S_A1 b_0
# t2: I_A2 S_A2 b_1
# t3: I_A1 g_0
# t4: I_A2 g_1
#
#
# t1: I_A1 S_A1 b_0
# t2: I_A2 S_A2 b_1
# t3: I_A1 g_0
# t4: I_A2 g_1


# %%[markdown]
# Note:
# * When `cartesian_control = True`, not setting up the controller relationship of the template impacts both the structure and ODE system of the stratified model
# * However, when `cartesian_control = False`, the ODE systems are identical (missing cross-terms)
