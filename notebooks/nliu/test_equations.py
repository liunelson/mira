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
MIRA_REST_URL = "http://mira-epi-dkg-lb-dc1e19b273dedaa2.elb.us-east-1.amazonaws.com"
os.environ["MIRA_REST_URL"] = MIRA_REST_URL

# %%
def generate_odesys(model, latex: bool = False, latex_align: bool = False) -> list:

    # State variables
    odeterms = {var: [] for var in sorted(model.get_concepts_name_map().keys())}

    # ODE terms from template rate laws
    for template in model.templates:
        if hasattr(template, "subject"):
            var = template.subject.name
            if template.rate_law is not None:
                odeterms[var].append(-template.rate_law.args[0])
            else:
                odeterms[var].append(0)
        
        if hasattr(template, "outcome"):
            var = template.outcome.name
            if template.rate_law is not None:
                odeterms[var].append(template.rate_law.args[0])
            else:
                odeterms[var].append(0)

    # Sort the terms such that all negative ones come first
    # odeterms = {var: sorted(terms, key = lambda term: 0 if str(term)[0] == '-' else 1) for var, terms in odeterms.items()}
    odeterms = {var: sorted(terms, key = lambda term: str(term)) for var, terms in odeterms.items()}

    # Time variable
    symb = lambda x: sympy.Symbol(x)
    try:
        time = model.time.name
    except:
        time = "t"
    finally:
        t = symb(time)

    # Add "(t)" for all the state variables as time-dependent symbols
    for var, terms in odeterms.items():
        for i, term in enumerate(terms):
            if hasattr(term, 'atoms'):
                for atom in term.atoms(sympy.Symbol):
                    if str(atom) in odeterms.keys():
                        term = term.subs(atom, sympy.Function(str(atom))(t))
                terms[i] = term

    # Construct equations
    num_terms = 5 # Max number of terms per line in LaTeX align
    odesys = []
    exprs = ""
    for i, (var, terms) in enumerate(odeterms.items()):

        lhs = sympy.diff(sympy.Function(var)(t), t)
        rhs = sum(terms)

        if latex == False:
            odesys.append(sympy.Eq(lhs, rhs))
        
        elif (latex == True) & (latex_align == False):
            odesys.append(sympy.latex(sympy.Eq(lhs, rhs)))

        elif (latex == True) & (latex_align == True):
            
            exprs += sympy.latex(lhs) + " ={}& "

            # Few equation terms = no wrapping needed
            if len(terms) < num_terms:
                exprs += sympy.latex(rhs)

            else:
                rhs = [sympy.latex(sum(terms[j:(j + num_terms)])) for j in range(0, len(terms), num_terms)]
                rhs = [line if (j == 0) | (line[0] == '-') else "+ " + line for j, line in enumerate(rhs)] # Add '+ ' to all lines past the first if not start with '- '

                exprs += " \\\\ \n    &".join(rhs)

            if i < (len(odeterms) - 1):
                exprs += " \\\\ \n"

            odesys = [exprs]

        else:
            pass
        

    # Repeat for observables if present
    if len(model.observables) > 0:

        # Sort observables alphabetically
        observables = {obs: model.observables[obs].expression.args[0] for obs in sorted(model.observables.keys())}

        # Add "(t)" for all the state variables as time-dependent symbols
        for obs, expr in observables.items():
            if hasattr(expr, 'atoms'):
                for atom in expr.atoms(sympy.Symbol):
                    if str(atom) in odeterms.keys():
                        expr = expr.subs(atom, sympy.Function(str(atom))(t))
                observables[obs] = expr
        
        for i, (obs, expr) in enumerate(observables.items()):
            lhs = sympy.Function(obs)(t)
            rhs = expr
            if latex == False:
                odesys.append(sympy.Eq(lhs, rhs))
            
            elif (latex == True) & (latex_align == False):
                odesys.append(sympy.latex(sympy.Eq(lhs, rhs)))

            elif (latex == True) & (latex_align == True):
                exprs = "     " + sympy.latex(lhs) + " ={}& " + sympy.latex(rhs)
                if i == 0:
                    exprs = " \\\\ \n" + exprs
                if i < (len(observables) - 1):
                    exprs += " \\\\ \n"

                odesys[0] += exprs

    if (latex == True) & (latex_align == True):
        odesys[0] = "\\begin{align} \n    " + odesys[0] + "\n\\end{align}"

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
    data['display name'] = [t.display_name for t in model.templates]

    df = pandas.DataFrame(data)

    return df

def post_mira_cleanup(mmt: TemplateModel) -> TemplateModel:

    # Assume time unit = 'day'
    mmt.time.units = Unit(expression = sympy.Symbol('day'))

    # Assign default display names and units to all concepts in every template
    for t in mmt.templates:
        for i in t.get_concepts():
            i.display_name = f'{i.name}(t)'
            i.units = Unit(expression = sympy.Integer(1))

    # Ditto for all parameters
    for p, param in mmt.parameters.items():
        param.display_name = p
        param.units = Unit(expression = sympy.Integer(1))

    # Ditto for observables
    for obs in mmt.observables.values():
        obs.display_name = f'{obs.name}(t)'
        obs.units = Unit(expression = sympy.Integer(1))

    # Ensure model templates have unique names/ids
    if len(mmt.templates) > 0:
        # if len({t.name for t in mmt.templates}) < len(mmt.templates):
        for i, t in enumerate(mmt.templates):

            t.name = f't{i}:'

            t.display_name = f'{t.type}'
            if 'Production' in t.type:
                t.name += f'->{t.outcome.name}'
                t.display_name += f' of {t.outcome.display_name}'
            elif 'Degradation' in t.type:
                t.name += f'{t.subject.name}->'
                t.display_name += f' of {t.subject.display_name}'
            elif 'Conversion' in t.type:
                t.name += f'{t.subject.name}->{t.outcome.name}'
                t.display_name += f' from {t.subject.display_name} to {t.outcome.display_name}'

            if getattr(t, 'controller', False):
                t.display_name += f' controlled by {t.controller.display_name}'
            elif getattr(t, 'controllers', False):
                t.display_name += f' controlled by {" and ".join([c.display_name for c in t.controllers])}'
        
    # Set default values for the model parameters
    for param in mmt.parameters.values():
        param.value = 0.0

    # Ensure every state variable has an initial condition
    # mmt.initials = mmt.initials | {c: Initial(concept = concept, expression = sympy.Symbol(f'{c}0')) for c, concept in mmt.get_concepts_name_map().items()}
    mmt.initials = mmt.initials | {c: Initial(concept = concept, expression = sympy.Float(f'0')) for c, concept in mmt.get_concepts_name_map().items()}

    # Set default values for the initial condition parameters
    mmt.parameters = mmt.parameters | {f'{c}0': Parameter(name = f'{c}0', display_name = f'{c}0', description = f'Initial value of state variable "{c}"', value = 0.0) for c in mmt.get_concepts_name_map().keys()}
    
    return mmt

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

# %%
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

with open('./data/model_equations/model.json', 'w') as fp:
    json.dump(template_model_to_petrinet_json(model), fp, indent = 4)

generate_summary_table(model)

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

# %%
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
    r"\frac{d S(t)}{d t} = -b * S(t) * I(t)",
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
# Model 3 (grouped controlled conversion)

odes = [
    sympy.Eq(S(t).diff(t), - b * S(t) * I(t)),
    sympy.Eq(I(t).diff(t), b * S(t) * I(t) - g * (I(t) + S(t) + R(t))),
    sympy.Eq(R(t).diff(t), k * g * I(t)),
    sympy.Eq(V(t).diff(t), (1 - k) * g * I(t))
]

model3 = template_model_from_sympy_odes(odes)

model3.observables = {
    'Dummy': Observable(name = 'Dummy', expression = sympy.Symbol('S') + sympy.Symbol('I') + sympy.Symbol('R'))
}

model3 = post_mira_cleanup(model3)

print(generate_odesys(model3, latex = True, latex_align = True))

with open('./data/model_equations/model3.json', 'w') as fp:
    json.dump(template_model_to_petrinet_json(model3), fp, indent = 4)

generate_summary_table(model3)

# %%
model3.draw_jupyter()

# %%
def generate_model_latex(model):

    odeterms = {var: 0 for var in model.get_concepts_name_map().keys()}

    for t in model.templates:
        if hasattr(t, "subject"):
            var = t.subject.name
            odeterms[var] -= t.rate_law.args[0]

        if hasattr(t, "outcome"):
            var = t.outcome.name
            odeterms[var] += t.rate_law.args[0]

    # Time
    if model.time and model.time.name:
        time = model.time.name
    else:
        time = "t"

    t = sympy.Symbol(time)

    # Construct Sympy equations
    odesys = []
    for var, terms in odeterms.items():
        lhs = sympy.diff(sympy.Function(var)(t), t)
  
        # Write (time-dependent) symbols with "(t)"
        rhs = terms
        for atom in terms.atoms(sympy.Symbol):
            if str(atom) in odeterms.keys():
                rhs = rhs.subs(atom, sympy.Function(str(atom))(t))

        odesys.append(sympy.latex(sympy.Eq(lhs, rhs)))

    # Observables
    if len(model.observables) > 0:

        # Write (time-dependent) symbols with "(t)"
        obs_eqs = []
        for obs in model.observables.values():
            lhs = sympy.Function(obs.name)(t)
            terms = obs.expression.args[0]
            rhs = terms
            for atom in terms.atoms(sympy.Symbol):
                if str(atom) in odeterms.keys():
                    rhs = rhs.subs(atom, sympy.Function(str(atom))(t))
            obs_eqs.append(sympy.latex(sympy.Eq(lhs, rhs)))

        # Add observables
        odesys += obs_eqs

    # Reformat:
    odesys = "\\begin{align} \n    " + " \\\\ \n    ".join([eq for eq in odesys]) + "\n\\end{align}"

    return odesys

# %%
odesys = generate_model_latex(model3)
print(odesys)

# %%
# Model 4 (sympy functions)

t = sympy.Symbol('t')
S, I, R, V = sympy.symbols('S I R V', cls = sympy.Function)
a, b, c, d, g, k = sympy.symbols('a b c d g k')

odes = [
    sympy.Eq(S(t).diff(t), - S(t) * (a * I(t) + b * R(t) + c * R(t) + d * V(t)) - sympy.exp(-a * S(t))),
    sympy.Eq(I(t).diff(t), S(t) * (a * I(t) + b * R(t) + c * R(t) + d * V(t)) - g * I(t)),
    sympy.Eq(R(t).diff(t), k * g * I(t)),
    # sympy.Eq(V(t).diff(t), (1 - k) * g * I(t) - sympy.exp(a) * V(t))
    sympy.Eq(V(t).diff(t), (1 - k) * g * I(t))
]

model4 = template_model_from_sympy_odes(odes)

model4.observables = {
    'Dummy': Observable(name = 'Dummy', expression = safe_parse_expr('S + I + R + V + sympy.sin(t)', local_dict = _clash))
}

# model4 = post_mira_cleanup(model4)

GraphicalModel.for_jupyter(model4)

# %%
generate_odesys(model4)

# %%
generate_summary_table(model4)

# %%
check_simplify_rate_laws(model4)

# %%
with open('./data/model_equations/model4.json', 'w') as fp:
    j = template_model_to_petrinet_json(model4)
    json.dump(j, fp, indent = 4)

# %%
m = template_model_from_amr_json(j)

# %%
# Model 5 (Hepatitis A)

# odes_latex = [
#     r'\frac{d S(t)}{d t} = - \frac{\beta * S(t) * I(t)}{N} - \frac{\omega * \tau * \theta * S(t)}{N}',
#     r'\frac{d L(t)}{d t} = \frac{\beta * S(t) * I(t)}{N} - \alpha * L(t)',
#     r'\frac{d I(t)}{d t} = \alpha * L(t) + \sigma * R(t) - \gamma * I(t)',
#     r'\frac{d Z(t)}{d t} = \frac{\omega * \tau * \theta * S(t)}{N} + \eta * \gamma * I(t)',
#     r'\frac{d R(t)}{d t} = (1 - \eta) * \gamma * I(t) - \sigma * R(t)',
#     r'\beta = \beta_s + \frac{\beta_l - \beta_s}{1 + \exp{(-c * (t - t_s))}}'
# ]

odes_latex = [
    "\\frac{d S(t)}{d t} = - \\frac{\\beta * S(t) * I(t)}{N} - \\frac{\\omega * \\tau * \\theta * S(t)}{N}",
    "\\frac{d L(t)}{d t} = \\frac{\\beta * S(t) * I(t)}{N} - \\alpha * L(t)",
    "\\frac{d I(t)}{d t} = \\alpha * L(t) + \\sigma * R(t) - \\gamma * I(t)",
    "\\frac{d Z(t)}{d t} = \\frac{\\omega * \\tau * \\theta * S(t)}{N} + \\eta * \\gamma * I(t)",
    "\\frac{d R(t)}{d t} = (1 - \\eta) * \\gamma * I(t) - \\sigma * R(t)",
    "\\beta = \\beta_s + \\frac{\\beta_l - \\beta_s}{1 + \\exp{-c * (t - t_s)}}"
]


import sympy
from sympy.abc import _clash
t = sympy.symbols('t')
beta_s, beta_l, c, t_s, omega, tau, theta, alpha, sigma, gamma, eta, N = sympy.symbols('beta_s beta_l c t_s omega tau theta alpha sigma gamma eta N')
S, L, I, Z, R = sympy.symbols('S L I Z R', cls=sympy.Function)
beta = beta_s + (beta_l - beta_s) / (1 + 2.7183**(-c * (t - t_s)))
odes_sympy = [
    sympy.Eq(S(t).diff(t), -beta * S(t) * I(t) / N - omega * tau * theta * S(t) / N),
    sympy.Eq(L(t).diff(t), beta * S(t) * I(t) / N - alpha * L(t)),
    sympy.Eq(I(t).diff(t), alpha * L(t) + sigma * R(t) - gamma * I(t)),
    sympy.Eq(Z(t).diff(t), omega * tau * theta * S(t) / N + eta * gamma * I(t)),
    sympy.Eq(R(t).diff(t), (1 - eta) * gamma * I(t) - sigma * R(t))
]

model5 = template_model_from_sympy_odes(odes_sympy)
generate_summary_table(model5)

# %%
GraphicalModel.for_jupyter(model5)

# %%