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
def viz_mmt(mmt) -> NoReturn:

    res = requests.post(url = f'{MIRA_REST_URL}/viz/to_image', data = sir_2age_d_model.json())

    with open('sir_2age_d_model.png', 'wb') as f:
        f.write(res.content)

    im = Image.open('sir_2age_d_model.png')
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
# observables = {'total_population': Observable(name = 'Total Population', expression = S + I + R)}

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

# %%
# Build SIR model from templates
sir_model = TemplateModel(
    templates = [t1, t2],
    parameters = parameters,
    initials = initials,
    time = Time(name = 't', units = day_units()),
    # observables = observables,
    annotations = Annotations(name = 'SIR model')
)

# Generate AMR JSON
AskeNetPetriNetModel(Model(sir_model)).to_json_file('sir_model.json')

# %%
# Do stratification with 2 age groups

sir_2age_model = stratify(
    sir_model,
    key = 'age',
    strata = ['A1', 'A2'],
    structure = [['A1', 'A2'], ['A2', 'A1']],
    directed = False,
    cartesian_control = True,
    params_to_stratify = {'b', 'g'},
    concepts_to_stratify = {'S', 'I', 'R'}
)

sir_2age_model.annotations.name = 'SIR model + stratified by 2 age groups'

AskeNetPetriNetModel(Model(sir_2age_model)).to_json_file('sir_2age_model.json')

# %%
# Add D compartment

# New concept
concepts['D'] = Concept(name = 'D', units = person_units(), identifiers = {'ido': '0000509'})

# New `initials`
initials = sir_2age_model.initials
initials['D'] = Initial(concept = Concept(name = 'D'), value = 0)

# New parameters
parameters = {**parameters, **sir_2age_model.parameters, **{
    'd1': Parameter(name = 'd1', value = 1.0/20.0, units = per_day_units()),
    'd2': Parameter(name = 'd2', value = 1.0/20.0, units = per_day_units())
}}

# New symbols
I_A1, I_A2, D, d1, d2 = sympy.symbols('I_A1 I_A2 D d1 d2')

# New templates
# I_A1 -> D
t3 = NaturalConversion(
    subject = sir_2age_model.get_concept('I_A1'),
    outcome = concepts['D'],
    rate_law = d1 * I_A1
)

# I_A2 -> D
t4 = NaturalConversion(
    subject = sir_2age_model.get_concept('I_A2'),
    outcome = concepts['D'],
    rate_law = g * I_A2
)

sir_2age_d_model = TemplateModel(
    templates = sir_2age_model.templates + [t3, t4],
    parameters = parameters,
    initials = initials,
    time = Time(name = 't', units = day_units()),
    # observables = observables,
    annotations = Annotations(name = 'SIR model + stratified by 2 age groups + D')
)

AskeNetPetriNetModel(Model(sir_2age_d_model)).to_json_file('sir_2age_d_model.json')

# %%

# %%
# 'Z': Concept(name = 'D', units = person_units(), identifiers = {'ido': '0000511'})


