# %%[markdown]
# Author: Nelson Liu
# 
# Email: [nliu@uncharted.software](mailto:nliu@uncharted.software)

# %%
import sympy
import itertools
from mira.metamodel import *
# from mira.metamodel import ControlledConversion, NaturalConversion, Concept, Template, TemplateModel
from mira.modeling import Model
# from mira.modeling.viz import GraphicalModel
# from mira.modeling.askenet.petrinet import AskeNetPetriNetModel
from mira.modeling.amr.petrinet import AMRPetriNetModel

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
}

# Parameters
parameters = {
    'beta': Parameter(name = 'beta', value = 0.4, units = per_day_units()),
    'gamma': Parameter(name = 'gamma', value = 1.0/11.0, units = per_day_units()),
}

# Initial conditions
initials = {
    'S': Initial(concept = Concept(name = 'S'), value = 1_000 - 1, expression = sympy.Float(1_000 - 1)),
    'I': Initial(concept = Concept(name = 'I'), value = 1, expression = sympy.Float(1)),
    'R': Initial(concept = Concept(name = 'R'), value = 0, expression = sympy.Float(0))
}

# Symbols
S, I, R, beta, gamma = sympy.symbols('S I R beta gamma')

# Observables
observables = {
    'total_population': Observable(name = 'Total Population', expression = SympyExprStr(S) + SympyExprStr(I) + SympyExprStr(R))
}

# S -> I transition
t1 = ControlledConversion(
    subject = concepts['S'],
    outcome = concepts['I'],
    controller = concepts['I'],
    rate_law = S * I * beta
)

# I -> R transition
t2 = NaturalConversion(
    subject = concepts['I'], 
    outcome = concepts['R'],
    rate_law = gamma * I
)

# %%
# Build SIR model from templates
sir_model = TemplateModel(
    templates = [t1, t2],
    parameters = parameters,
    initials = initials,
    time = Time(name = 't', units = day_units()),
    observables = observables,
    annotations = Annotations(name = 'SIR model from MIRA')
)

# Generate AMR JSON
AMRPetriNetModel(Model(sir_model)).to_json_file('sir_model.json')

# %%
# Do stratification with 2 locations

sir_loc_model = stratify(
    sir_model,
    key = 'location',
    strata = ['TOR', 'MTL'],
    structure = [['TOR', 'MTL'], ['MTL', 'TOR']],
    directed = False,
    cartesian_control = False,
    params_to_stratify = {'beta', 'gamma'},
    concepts_to_stratify = {'S', 'I', 'R'}
)

sir_loc_model.annotations.name = 'SIR model from MIRA, 2-location stratified with MIRA'

AMRPetriNetModel(Model(sir_loc_model)).to_json_file('sir_loc_model.json')

# %%
