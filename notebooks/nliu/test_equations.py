# %%
import os
import sympy
from sympy.parsing.latex import parse_latex
# pip install antlr4-python3-runtime==4.11
import json
import pandas

from mira.metamodel import *
from mira.modeling import Model
from mira.modeling.amr.petrinet import template_model_to_petrinet_json
from mira.modeling.amr.regnet import template_model_to_regnet_json
from mira.modeling.amr.regnet import AMRRegNetModel
from mira.sources.sympy_ode import template_model_from_sympy_odes

from mira.modeling.viz import GraphicalModel

# %%
MIRA_REST_URL = "http://mira-epi-dkg-lb-dc1e19b273dedaa2.elb.us-east-1.amazonaws.com"
os.environ["MIRA_REST_URL"] = MIRA_REST_URL

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
    for k in ("subject", "outcome", "controller"):
        data[k] = [getattr(t, k).name if hasattr(t, k) else None for t in model.templates]

    data["controllers"] = [[c.name for c in getattr(t, k)] if hasattr(t, "controllers") else None for t in model.templates]
    data["controller(s)"] = [i if j == None else j for i, j in zip(data["controller"], data["controllers"])]
    __ = data.pop("controller")
    __ = data.pop("controllers")

    data["rate_law"] = [t.rate_law for t in model.templates]
    data["interactor_rate_law"] = [t.get_interactor_rate_law() for t in model.templates]

    df = pandas.DataFrame(data)

    return df

# %%
# Define concepts
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
# Define model from equations
# 
# 0. Basic SIR
# 1. With birth, death, controlled production
# 2. Branching ratios 

# %%
# Model 0

t = sympy.Symbol('t')
S, E, I, R, V = sympy.symbols('S E I R V', cls = sympy.Function)
l, k, d, a, b, g, m, b = sympy.symbols('l k d a b g m b')

odes = [
    sympy.Eq(S(t).diff(t), - b * S(t) * I(t)),
    sympy.Eq(I(t).diff(t), b * S(t) * I(t) - g * I(t)),
    sympy.Eq(R(t).diff(t), g * I(t)),
]

model0 = template_model_from_sympy_odes(odes)

print(generate_odesys(model0, latex = True, latex_align = True))

with open('./data/model_equations/model0.json', 'w') as fp:
    json.dump(template_model_to_petrinet_json(model0), fp, indent = 4)

generate_summary_table(model0)

# %%
# Need note: 
# - missing initials
# - missing 't' unit
# - missing template names
# - missing parameter values

# %%
odes = [
    sympy.Eq(S(t).diff(t), - b * S(t) * I(t) + l - m * S(t)),
    sympy.Eq(I(t).diff(t), b * S(t) * I(t) - g * I(t) - m * I(t)),
    sympy.Eq(R(t).diff(t), g * I(t) - m * R(t)),
]

model1 = template_model_from_sympy_odes(odes)
model1.annotations = Annotations(name = 'SIR with births and deaths')

print(generate_odesys(model1, latex = True, latex_align = True))

with open('./data/model_equations/model1.json', 'w') as fp:
    json.dump(template_model_to_petrinet_json(model1), fp, indent = 4)

generate_summary_table(model1)
# %%
model1.draw_jupyter()

# %%
# Test with LaTeX-to-SymPy

odes_latex = [
    r"\frac{d S(t)}{d t} = -b * S(t) * I(t) + l - m * S(t)",
    r"\frac{d I(t)}{d t} = b * S(t) * I(t) - g * I(t) - m * I(t)", 
    r"\frac{d R(t)}{d t} = g * I(t) - m * R(t)",
]

odes_sympy = [sympy.parsing.latex.parse_latex(ode) for ode in odes_latex]

__ = [print(ode) for ode in odes_sympy]

model = template_model_from_sympy_odes(odes_sympy)
generate_summary_table(model)

with open('./data/model_equations/model.json', 'w') as fp:
    json.dump(template_model_to_petrinet_json(model), fp, indent = 4)

# %%
# Variation:
# All death processes depend on 'S'

odes = [
    sympy.Eq(S(t).diff(t), - b * S(t) * I(t) + l - m * S(t)),
    sympy.Eq(I(t).diff(t), b * S(t) * I(t) - g * I(t) - m * S(t)),
    sympy.Eq(R(t).diff(t), g * I(t) - m * S(t)),
]

model1a = template_model_from_sympy_odes(odes)

print(generate_odesys(model1a, latex = True, latex_align = True))

with open('./data/model_equations/model1a.json', 'w') as fp:
    json.dump(template_model_to_petrinet_json(model1a), fp, indent = 4)

generate_summary_table(model1a)

# Note:
# All 'm * S' go missing!

# %%
model1a.draw_jupyter()

# %%
# Variation: 
# Controlled production of viral particles

odes = [
    sympy.Eq(S(t).diff(t), - l * S(t) * I(t)),
    sympy.Eq(E(t).diff(t), l * S(t) * I(t) - k * E(t)),
    sympy.Eq(I(t).diff(t), k * E(t) - d * I(t)),
    sympy.Eq(V(t).diff(t), a * b * (1 - g) * I(t)),
]

model1b = template_model_from_sympy_odes(odes)

print(generate_odesys(model1b, latex = True, latex_align = True))

with open('./data/model_equations/model1b.json', 'w') as fp:
    json.dump(template_model_to_petrinet_json(model1b), fp, indent = 4)

generate_summary_table(model1b)

# %%
model1b.draw_jupyter()

# %%
# Test with LaTeX-to-SymPy

odes_latex = [
    r"\frac{d S(t)}{d t} = -l S(t) I(t)",
    r"\frac{d E(t)}{d t} = l S(t) I(t) - k E(t)", 
    r"\frac{d I(t)}{d t} = k E(t) - d * I(t)",
    r"\frac{d V(t)}{d t} = a * b * (1 - g) * I(t)"
]

odes_sympy = [sympy.parsing.latex.parse_latex(ode) for ode in odes_latex]

__ = [print(ode) for ode in odes_sympy]


# %%
# Model 2

odes = [
    sympy.Eq(S(t).diff(t), - b * S(t) * I(t)),
    sympy.Eq(I(t).diff(t), b * S(t) * I(t) - k * g * I(t) - (1 - k) * g * I(t)),
    sympy.Eq(R(t).diff(t), k * g * I(t)),
    sympy.Eq(V(t).diff(t), (1 - k) * g * I(t))
]

model2a = template_model_from_sympy_odes(odes)

print(generate_odesys(model2a, latex = True, latex_align = True))

with open('./data/model_equations/model2a.json', 'w') as fp:
    json.dump(template_model_to_petrinet_json(model2a), fp, indent = 4)

generate_summary_table(model2a)

# %%
# Test with LaTeX-to-SymPy

odes_latex = [
    r"\frac{d S(t)}{d t} = -b S(t) I(t)",
    r"\frac{d I(t)}{d t} = b * S(t) * I(t) - k * g * I(t) - (1 - k) * g * I(t)", 
    r"\frac{d R(t)}{d t} = k * g * I(t)",
    r"\frac{d V(t)}{d t} = (1 - k) * g * I(t)"
]

odes_sympy = [sympy.parsing.latex.parse_latex(ode) for ode in odes_latex]

__ = [print(ode) for ode in odes_sympy]

model = template_model_from_sympy_odes(odes_sympy)
generate_summary_table(model)

# %%
model2a.draw_jupyter()

# %%
# Model 2b (simplified/hidden branching)

odes = [
    sympy.Eq(S(t).diff(t), - b * S(t) * I(t)),
    sympy.Eq(I(t).diff(t), b * S(t) * I(t) - g * I(t)),
    sympy.Eq(R(t).diff(t), k * g * I(t)),
    sympy.Eq(V(t).diff(t), (1 - k) * g * I(t))
]

model2b = template_model_from_sympy_odes(odes)

print(generate_odesys(model2b, latex = True, latex_align = True))

with open('./data/model_equations/model2b.json', 'w') as fp:
    json.dump(template_model_to_petrinet_json(model2b), fp, indent = 4)

generate_summary_table(model2b)

# %%
model2b.draw_jupyter()

# %%
