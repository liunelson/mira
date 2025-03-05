# %%
import os
import sympy
from sympy.parsing.latex import parse_latex
# pip install antlr4-python3-runtime==4.11
import json
import pandas
import copy

from mira.metamodel import *
from mira.modeling import Model
from mira.modeling.amr.petrinet import template_model_to_petrinet_json
from mira.modeling.amr.regnet import template_model_to_regnet_json
from mira.modeling.amr.regnet import AMRRegNetModel
from mira.sources.sympy_ode import template_model_from_sympy_odes
from mira.sources.amr.petrinet import template_model_from_amr_json
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

# Generate initial condition and parameter tables
def generate_init_param_tables(model) -> tuple[pandas.DataFrame, pandas.DataFrame]:

    data = {}
    data["name"] = [name for name, __ in model.initials.items()]
    data["expression"] = [init.expression for __, init in model.initials.items()]
    df_initials = pandas.DataFrame(data)

    data = {}
    data["name"] = [name for name, __ in model.parameters.items()]
    data["value"] = [param.value for __, param in model.parameters.items()]
    df_params = pandas.DataFrame(data)

    return (df_initials, df_params)

# %%
# Load model
with open('./data/example_models/SIR (3).json', 'r') as fp:
    model = template_model_from_amr_json(json.load(fp))

GraphicalModel.for_jupyter(model)

# %%
generate_summary_table(model)

# %%
generate_odesys(model)

# %%
# Add contact factor "c_beta" to stratifying "beta" directly

add_param_factor = True

params_to_stratify = ['beta', 'gamma']

model_ = copy.deepcopy(model)
if add_param_factor:

    # Not f_ to avoid Terarium grouping issue
    new_params = {param: 'f' + param for param in params_to_stratify}

    for param, factor in new_params.items():

        # In case of parameter name clash
        while factor in model_.parameters.keys():
            factor += '_0'
        new_params[param] = factor


        # Replace 'param' with 'param * factor' in all rate laws
        for template in model_.templates:
            template.rate_law = SympyExprStr(
                template.rate_law.args[0].subs(
                    sympy.Symbol(param),
                    sympy.Symbol(param) * sympy.Symbol(factor)
                )
            )

        # Add the factor as a new model parameter
        model_.add_parameter(
            parameter_id = factor,
            name = factor,
            description = f'Stratification factor of the parameter {param}.',
            value = 1.0
        )

    params_to_stratify = list(new_params.values())

# %%
model_strat = stratify(
    template_model = model_,
    key = 'age',
    strata = ['y', 'o'],
    structure = [],
    directed = True,
    cartesian_control = True,
    concepts_to_stratify = ['S', 'I'],
    concepts_to_preserve = ['R'],
    params_to_stratify = params_to_stratify,
    params_to_preserve = None,
    modify_names = True,
    param_renaming_uses_strata_names = True
)

# %%
# Remove unused`f<param>` parameters (should be automatic?)
for __, factor in new_params.items():
    if factor in model_strat.parameters.keys():
        __ = model_strat.parameters.pop(factor)

# %%
generate_summary_table(model_strat)

# %%
generate_odesys(model_strat)

# %%
(df_initials, df_params) = generate_init_param_tables(model_strat)
df_initials

# %%
df_params


# %%
