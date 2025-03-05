# %%
import os
import sympy
from sympy.abc import _clash
from sympy.parsing.latex import parse_latex
# pip install antlr4-python3-runtime==4.11
import json
import pandas
from copy import deepcopy

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
        odesys[0] = "\\begin{align*} \n    " + odesys[0] + "\n\\end{align*}"

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
# Define base model

# Define time variable
t = sympy.symbols("t")

# Define variables with time derivative to be time-dependent functions
S, I, R = sympy.symbols("S I R", cls=sympy.Function)

# Define all parameters
a, b, m, g = sympy.symbols("a b m g")

# Define the equations without time-derivative on the left hand side
equation_output = [
    sympy.Eq(S(t).diff(t), -b * S(t) * (I(t) + a * R(t)) - m * S(t)),
    # sympy.Eq(I(t).diff(t), b * S(t) * I(t) - g * I(t) - m * I(t)),
    sympy.Eq(I(t).diff(t), b * S(t) * (I(t) + a * R(t)) - g * I(t) - m * I(t)),
    sympy.Eq(R(t).diff(t), g * I(t) - m * R(t))
]

model = template_model_from_sympy_odes(equation_output)

model = post_mira_cleanup(model)

generate_summary_table(model)

# %%
GraphicalModel.for_jupyter(model)

# %%
model_simp = simplify_rate_laws(model)
generate_summary_table(model_simp)

# %%
generate_odesys(model_simp)

# %%
print(f'Number of templates in complex model: {len(model.templates)}')
print(f'... in the simplified model: {len(model_simp.templates)}')

with open('./data/example_complex_model.json', 'w') as fp:
    json.dump(template_model_to_petrinet_json(model), fp, indent = 4)

with open('./data/example_simplified_model.json', 'w') as fp:
    json.dump(template_model_to_petrinet_json(model_simp), fp, indent = 4)

# %%
with open('./data/final_evaluation/Antagonistic Infection Model.json', 'r') as fp:
    large_model = template_model_from_amr_json(json.load(fp))

large_model_simp = simplify_rate_laws(large_model)

print(f'Number of templates in complex model: {len(large_model.templates)}')
print(f'... in the simplified model: {len(large_model_simp.templates)}')

with open('./data/final_evaluation/Antagonistic Infection Model - Simplified.json', 'w') as fp:
    json.dump(template_model_to_petrinet_json(large_model_simp), fp, indent = 4)

# %%


