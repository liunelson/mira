# %%[markdown]
# Author: Nelson Liu
# 
# Email: [nliu@uncharted.software](mailto:nliu@uncharted.software)

# %%
# Test a stratification edge case for Dan:
#
# Start with SIR-type model where some parameters are shared across transitions and some states are not stratified
# 1. SIRD
# 3. S -(b*c*S*I)-> I
# 4. I -(g*I)-> R
# 5. I -((1-g)*I)-> D

# %%
import sympy
import itertools
from mira.metamodel import *
from mira.modeling import Model
# from mira.modeling.viz import GraphicalModel
from mira.modeling.askenet.petrinet import AMRPetriNetModel

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
    'R': Concept(name = 'R', units = person_units(), identifiers = {'ido': '0000592'}),
    'D': Concept(name = 'D', units = person_units(), identifiers = {'ido': '0000592'})
}

# Parameters
parameters = {
    'b': Parameter(name = 'b', value = 0.4, units = dimensionless_units()),
    'g': Parameter(name = 'g', value = 1.0/10.0, units = dimensionless_units()),
    'c': Parameter(name = 'c', value = 0.5, units = dimensionless_units()),
    'r': Parameter(name = 'r', value = 0.5, units = dimensionless_units())
}

# Initial conditions
initials = {
    'S': Initial(concept = Concept(name = 'S'), value = 1_000 - 1, expression = sympy.Float(1_000 - 1)),
    'I': Initial(concept = Concept(name = 'I'), value = 1, expression = sympy.Float(1)),
    'R': Initial(concept = Concept(name = 'R'), value = 0, expression = sympy.Float(0)),
    'D': Initial(concept = Concept(name = 'D'), value = 0, expression = sympy.Float(0)),
}

# Symbols
S, I, R, D, b, g, c, r = sympy.symbols('S I R D b g c r')

# Observables
observables = {'N': Observable(name = 'N', expression = S + I + R + D)}

# %%
# S -> I transition
t1 = ControlledConversion(
    subject = concepts['S'],
    outcome = concepts['I'],
    controller = concepts['I'],
    rate_law = b * c * S * I
)

# I -> R transition
t2 = NaturalConversion(
    subject = concepts['I'], 
    outcome = concepts['R'],
    rate_law = g * c * I
)

# I -> D transition
t3 = NaturalConversion(
    subject = concepts['I'],
    outcome = concepts['D'],
    rate_law = (1 - g) * r * I
) 

# %%
# Build SIRD model from templates
model = TemplateModel(
    templates = [t1, t2, t3],
    parameters = parameters,
    initials = initials,
    time = Time(name = 't', units = day_units()),
    observables = observables,
    annotations = Annotations(name = 'SIRD model')
)

# Visualize
viz_mmt(model, 'sird_model.png')

# Generate AMR JSON
AMRPetriNetModel(Model(model)).to_json_file('sird_model.json')

# %%
# Stratify with 2 strata

model_stratified = stratify(
    model,
    key = 'age',
    strata = ['A1', 'A2'],
    structure = [],
    directed = False,
    cartesian_control = True,
    params_to_stratify = {'b', 'g', 'r'},
    params_to_preserve = {'c'},
    concepts_to_stratify = {'S', 'I'},
    concepts_to_preserve = {'R', 'D'}
)

model_stratified.annotations.name = 'SIRD model stratified by 2 age groups'

viz_mmt(model_stratified, 'sird_stratified_model.png')

AMRPetriNetModel(Model(model_stratified)).to_json_file('sird_stratified_model.json')

# %%
# New rate laws

for i, t in enumerate(model_stratified.templates):
    if hasattr(t, "controller"):
        print(f"Template {i}: {t.subject.name} (controlled by {t.controller.name}) -> {t.outcome.name} via {str(t.rate_law)}")
    elif hasattr(t, "controllers"):
        print(f"Template {i}: {t.subject.name} (controlled by {[c.name for c in t.controllers]}) -> {t.outcome.name} via {str(t.rate_law)}")
    else:
        print(f"Template {i}: {t.subject.name} -> {t.outcome.name} via {str(t.rate_law)}")

# %%
# Template 0: S_A1 (controlled by I_A1) -> I_A1 via I_A1*S_A1*b_0*c
# Template 1: S_A1 (controlled by I_A2) -> I_A1 via I_A2*S_A1*b_1*c
# Template 2: S_A2 (controlled by I_A2) -> I_A2 via I_A2*S_A2*b_2*c
# Template 3: S_A2 (controlled by I_A1) -> I_A2 via I_A1*S_A2*b_3*c
# Template 4: I_A1 -> R via I_A1*c*g_0
# Template 5: I_A2 -> R via I_A2*c*g_1
# Template 6: I_A1 -> D via I_A1*r_0*(1 - g_2)
# Template 7: I_A2 -> D via I_A2*r_1*(1 - g_3)

# %%

