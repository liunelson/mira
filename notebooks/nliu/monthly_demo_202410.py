# %%[markdown]
# Author: Nelson Liu
# 
# Email: [nliu@uncharted.software](mailto:nliu@uncharted.software)

# %%[markdown]
# # Monthly Demo (2024-10)

# %%
import json
import pandas
from sympy.parsing.sympy_parser import parse_expr
from tabulate import tabulate
from tqdm import tqdm
import sympy

from mira.modeling.viz import GraphicalModel
from mira.modeling.amr.petrinet import template_model_to_petrinet_json
from mira.sources.amr.petrinet import template_model_from_amr_json
from mira.metamodel import *
# from mira.modeling.amr.ops import *

# %%
# Generate summary table of a template model
def generate_summary_table(model) -> pandas.DataFrame:

    data = {"name": [t.name for t in model.templates]}
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

# Generate Sympy equations from a template model
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

# %%
with open("./data/monthly_demo_202410/Prob 5 Model A.json", "r") as f:
    model = template_model_from_amr_json(json.load(f))

GraphicalModel.for_jupyter(model)

# %%
# model.initials = {k: Initial(concept = Concept(name = k), expression = parse_expr(f'{k}0')) for k in ('s', 'i', 'r')}
# model.parameters = model.parameters | {f'{k}0': Parameter(name = f'{k}0', value = 0.0) for k in 'sir'}

# %%[markdown]
# ## Problem 5
# 
# Test stratification to ~2,000 counties

# %%
model_2 = stratify(
    model,
    key = 'location',
    strata = ['05001', '05003'],
    structure = [],
    directed = True,
    cartesian_control = False, 
    concepts_to_stratify = {'s', 'i', 'r'},
    params_to_stratify = {'y', 'R'},
    param_renaming_uses_strata_names = True
)

GraphicalModel.for_jupyter(model_2)

# %%
generate_summary_table(model_2)

# %%
# Load FIPS dataset from Josh Hewitt

df_hewitt = pandas.read_csv(
    './data/monthly_demo_202410/from_josh_hewitt/sir_parameters.csv',
    dtype = {
        'fips': str,
        't0': str,
        'r_init': float
    }
)

df_hewitt['t0'] = pandas.to_datetime(df_hewitt['t0'])

df_hewitt

# %%
# Load US DoT FIPS code table

df_dot = pandas.read_csv(
    './data/monthly_demo_202410/State__County_and_City_FIPS_Reference_Table.csv',
    dtype = str
)

# Note: in 2015, Shannon County -> Oglala Lakota County
# FIPS 46113 <-> FIPS 46102
df_dot_ = df_dot[df_dot['StCnty FIPS Code'] == '46113']
df_dot_.loc[:, 'StCnty FIPS Code'] = '46102'
df_dot = pandas.concat([df_dot, df_dot_]).reset_index()

# %%
x, y = generate_init_param_tables(model_2)

print(tabulate(x, headers = 'keys'))
print(tabulate(y, headers = 'keys'))

# %%
# Repeat for `n` counties and time it

def stratify_n(model, n):

    if isinstance(n, int) and (n < len(df_hewitt['fips'])):
        counties = list(df_hewitt['fips'].iloc[:n])
    elif n == 'all':
        counties = list(df_hewitt['fips'])

    model_n = stratify(
        model,
        key = 'location',
        strata = counties,
        structure = [],
        directed = True,
        cartesian_control = False, 
        concepts_to_stratify = {'s', 'i', 'r'},
        params_to_stratify = {'y', 'R'},
        param_renaming_uses_strata_names = True
    )

    return model_n

def configure_n(model, df_hewitt):

    for k, p in tqdm(model.parameters.items()):
        if '_' not in k:
            continue

        fips = k.split('_')[-1]
        if k[:2] == 'R_':
            p.value = float(df_hewitt[df_hewitt['fips'] == fips]['Rl'].values[0])
        elif k[:2] == 'y_':
            p.value = float(df_hewitt[df_hewitt['fips'] == fips]['gamma'].values[0])
        else:
            pass

    for k, i in model.initials.items():
        fips = k.split('_')[-1]
        if k[:2] == 's_':
            i.expression = parse_expr(str(float(df_hewitt[df_hewitt['fips'] == fips]['s_init'].values[0])))
        elif k[:2] == 'i_':
            i.expression = parse_expr(str(float(df_hewitt[df_hewitt['fips'] == fips]['i_init'].values[0])))
        elif k[:2] == 'r_':
            i.expression = parse_expr(str(float(df_hewitt[df_hewitt['fips'] == fips]['r_init'].values[0])))
        else:
            pass

    return model

def describe_n(model, df_dot):

    for k, p in tqdm(model.parameters.items()):
        if '_' not in k:
            continue

        fips = k.split('_')[-1]
        n = list(df_dot[df_dot['StCnty FIPS Code'] == fips].iloc[0, 1:3].values)
        if k[0] == 'R':
            s = 'Basic reproduction number'
        elif k[0] == 'y':
            s = 'Recovery rate'

        p.description = f'{s}, county of {n[1]}, state of {n[0]}'

    for k, c in tqdm(model.get_concepts_name_map().items()):
        if '_' not in k:
            continue
        
        fips = k.split('_')[-1]
        n = list(df_dot[df_dot['StCnty FIPS Code'] == fips].iloc[0, 1:3].values)

        if k[0] == 's':
            s = 'Susceptible population'
        elif k[0] == 'i':
            s = 'Infected populaton'
        elif k[0] == 'r':
            s = 'Recovered population'

        c.description = f'{s}, county of {n[1]}, state of {n[0]}'

    for t in tqdm(model.templates):

        if has_subject(t):
            k = t.subject.name
            if '_' not in k:
                continue
            
            fips = k.split('_')[-1]
            n = list(df_dot[df_dot['StCnty FIPS Code'] == fips].iloc[0, 1:3].values)

            if k[0] == 's':
                s = 'Susceptible population'
            elif k[0] == 'i':
                s = 'Infected populaton'
            elif k[0] == 'r':
                s = 'Recovered population'

            t.subject.description = f'{s}, county of {n[1]}, state of {n[0]}'

        if has_outcome(t):
            k = t.outcome.name
            if '_' not in k:
                continue
            
            fips = k.split('_')[-1]
            n = list(df_dot[df_dot['StCnty FIPS Code'] == fips].iloc[0, 1:3].values)

            if k[0] == 's':
                s = 'Susceptible population'
            elif k[0] == 'i':
                s = 'Infected populaton'
            elif k[0] == 'r':
                s = 'Recovered population'

            t.outcome.description = f'{s}, county of {n[1]}, state of {n[0]}'

        if has_controller(t):
            k = t.controller.name
            if '_' not in k:
                continue
            
            fips = k.split('_')[-1]
            n = list(df_dot[df_dot['StCnty FIPS Code'] == fips].iloc[0, 1:3].values)

            if k[0] == 's':
                s = 'Susceptible population'
            elif k[0] == 'i':
                s = 'Infected populaton'
            elif k[0] == 'r':
                s = 'Recovered population'

            t.controller.description = f'{s}, county of {n[1]}, state of {n[0]}'
    
    return model

# %%
%timeit -n 1 -r 1 model_n = stratify_n(model, 'all')
# 1.32 s

# %%
model_n = stratify_n(model, 'all')
model_n = configure_n(model_n, df_hewitt = df_hewitt)
model_n = describe_n(model_n, df_dot = df_dot)

# %%
# Bug: descriptions lost when using template_model_to_petrinet_json
model_to_json_file(model_n, './data/monthly_demo_202410/prob5_test.json')

j = template_model_to_petrinet_json(model_n)

for state in j['model']['states']:
    d = model_n.get_concepts_name_map()[state['id']].description
    state['description'] = d

with open('./data/monthly_demo_202410/Prob 5 Model A (1,992 counties).json', 'w') as f:
    json.dump(j, f, indent = 4)

# %%
