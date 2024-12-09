# %%
import os
import sympy
import json

from mira.metamodel import *
from mira.modeling import Model
from mira.modeling.amr.petrinet import template_model_to_petrinet_json
from mira.modeling.amr.regnet import template_model_to_regnet_json
from mira.modeling.amr.regnet import AMRRegNetModel

from mira.modeling.viz import GraphicalModel

# %%
MIRA_REST_URL = "http://mira-epi-dkg-lb-dc1e19b273dedaa2.elb.us-east-1.amazonaws.com"
os.environ["MIRA_REST_URL"] = MIRA_REST_URL

# %%
# Define concepts
# susceptible = Concept(name="S", identifiers={"ido": "0000514"})
# infected = Concept(name="I", identifiers={"ido": "0000511"})
# recovered = Concept(name="R", identifiers={"ido": "0000592"})

concepts = {k: Concept(name = k, display_name = k) for k in 'SEIVM'}
parameters = {p: Parameter(name = p, display_name = p, value = 1.0) for p in 'lkdabgmb'} | {f'{k}0': Parameter(name = f'{k}0', display_name = f'{k}0', value = 1.0) for k in concepts.keys()}
initials = {k: Initial(concept = c, expression = sympy.Symbol(f'{k}0')) for k, c in concepts.items()}

local_dict = {k: sympy.Symbol(k) for k in concepts.keys()} | {k: sympy.Symbol(k) for k in parameters.keys()}

# %%
# Define a non-mass-conserving model

model = TemplateModel(
    templates = [
        ControlledConversion(
            name = 'Infection',
            subject = concepts['S'],
            outcome = concepts['E'],
            controller = concepts['I'],
            rate_law = safe_parse_expr('l * S * I', local_dict = local_dict)
        ),
        NaturalConversion(
            name = 'Incubation',
            subject = concepts['E'],
            outcome = concepts['I'],
            rate_law = safe_parse_expr('k * E', local_dict = local_dict)
        ),
        NaturalDegradation(
            name = 'Recovery',
            subject = concepts['I'],
            rate_law = safe_parse_expr('d * I', local_dict = local_dict)
        ),
        ControlledProduction(
            name = 'Viral production',
            outcome = concepts['V'],
            controller = concepts['I'],
            rate_law = safe_parse_expr('a * b * (1 - g) * I', local_dict = local_dict)
        ),
        NaturalProduction(
            name = 'Natural birth',
            outcome = concepts['S'],
            rate_law = safe_parse_expr('b * S', local_dict = local_dict)
        ),
        GroupedControlledDegradation(
            name = 'Mask supply',
            subject = concepts['M'],
            controllers = [concepts['S'], concepts['E'], concepts['I']],
            rate_law = safe_parse_expr('m * (S + E + I)', local_dict = local_dict)
        )
    ],
    parameters = parameters,
    initials = initials,
    observables = {'N': Observable(name = 'N', expression = safe_parse_expr('S + E + I', local_dict = local_dict))},
    time = Time(name = 't', units = Unit(expression = safe_parse_expr('day')))
)

model.draw_jupyter()

# %%
with open('./data/example_models/non_mass_conserved.json', 'w') as fp:
    json.dump(template_model_to_petrinet_json(model), fp, indent = 4)

# %%
