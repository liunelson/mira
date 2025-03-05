# %%
import os
import sympy
from sympy.abc import _clash
from sympy.parsing.latex import parse_latex
# pip install antlr4-python3-runtime==4.11
import json
import pandas
from copy import deepcopy
import itertools

from mira.metamodel import *
from mira.modeling import Model
from mira.sources.amr.petrinet import template_model_from_amr_json
from mira.modeling.amr.petrinet import template_model_to_petrinet_json
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

# %%
# Define base model

state_vars = {
    'GDP': Concept(name = 'GDP', display_name = 'GDP(t)', units = Unit(expression = sympy.Symbol('dollar')), description = 'Gross domestric product of a country.'),
    'R': Concept(name = 'R', display_name = 'R(t)', units = Unit(expression = SympyExprStr(safe_parse_expr("1", local_dict = _clash))), description = 'Natural resource reserve of a commodity in a country.'),
    'OS': Concept(name = 'OS', display_name = 'OS(t)', units = Unit(expression = SympyExprStr(safe_parse_expr("1", local_dict = _clash))), description = 'Stock of old scrap of a commodity in a country.'),
    'S': Concept(name = 'S', display_name = 'S(t)', units = Unit(expression = SympyExprStr(safe_parse_expr("1", local_dict = _clash))), description = 'Stock of a commodity in a country.')
}

initials = {
    v: Initial(concept = concept, expression = sympy.Float('0.0'))
    for v, concept in state_vars.items()
}

parameters = {
    'PPI': Parameter(name = 'PPI', display_name = 'PPI', description = 'Policy perception index of a country.', units = Unit(expression = SympyExprStr(safe_parse_expr("1", local_dict = _clash)))),
    'FWI': Parameter(name = 'FWI', display_name = 'FWI', description = 'Freedom House\'s Freedom in the World index of a country.', units = Unit(expression = SympyExprStr(safe_parse_expr("1", local_dict = _clash)))),
    'MC': Parameter(name = 'MC', display_name = 'MC', description = 'Military cooperation, refers to whether the country has a current collective defense arrangement with the US.'),
    'VA': Parameter(name = 'VA', display_name = 'VA', description = 'Value added (i.e. its contribution to a GDP) of an industry in the US.', units = Unit(expression = sympy.Symbol('dollar'))),
    'OP': Parameter(name = 'OP', display_name = 'OP', description = 'Operating profit of an industry in the US.', units = Unit(expression = sympy.Symbol('dollar'))),
    'EXP': Parameter(name = 'EXP', display_name = 'EXP', description = 'Expenditure of an industry on a commodity in the US.', units = Unit(expression = sympy.Symbol('dollar'))),
    'DeltaS': Parameter(name = 'DeltaS', display_name = 'DeltaS', description = 'Adjustments of industry and government stocks of a commodity in a country.', units = Unit(expression = SympyExprStr(safe_parse_expr("1", local_dict = _clash)))),
    'PP': Parameter(name = 'PP', display_name = 'PP', description = 'Primary production of a commodity in a country.', units = Unit(expression = SympyExprStr(safe_parse_expr("1", local_dict = _clash)))),
    'SP': Parameter(name = 'SP', display_name = 'SP', description = 'Secondary production of a commodity in a country.', units = Unit(expression = SympyExprStr(safe_parse_expr("1", local_dict = _clash)))),
    'g': Parameter(name = 'g', display_name = 'g', description = 'GDO growth rate of a country.', units = Unit(expression = SympyExprStr(safe_parse_expr("1", local_dict = _clash))))
}

# Define observables
# observables = {
#     'SR': Observable(name = 'SR', display_name = 'SR(t)',  expression = safe_parse_expr('(DP * TE * EV) ** (1/3)', local_dict = _clash), description = 'Supply risk of a commodity for a country relative to all others.'),
#     'DP': Observable(name = 'DP', display_name = 'DP(t)', expression = safe_parse_expr('(PS ** 2) * ASI * WSI', local_dict = _clash), description = 'Disruption potential of a commodity for a country relative to all others.'),
#     'TE': Observable(name = 'TE', display_name = 'TE(t)', expression = safe_parse_expr('(I - E + DeltaS) / (PP + SP + I - E + DeltaS)', local_dict = _clash), description = 'Trade exposure of a commodity for a country relative to all others.'),
#     'AC': Observable(name = 'AC', display_name = 'AC(t)', expression = safe_parse_expr('PP + SP + I - E + DeltaS', local_dict = _clash), description = 'Apparent consumption of a commodity by a country.'),
#     'EV': Observable(name = 'EV', display_name = 'EV(t)', expression = safe_parse_expr('(VA / GDP) * (EXP / OP)', local_dict = _clash), description = 'Economic vulnerability of a commodity for a country relative to all others.'),
#     'PS': Observable(name = 'PS', display_name = 'PS(t)', expression = safe_parse_expr('P / P', local_dict = _clash), description = 'Production share of a commodity by a country'),
#     'P': Observable(name = 'P', display_name = 'P(t)', expression = safe_parse_expr('PP + SP', local_dict = _clash), description = 'Production of a commodity by a country'),
#     'ASI': Observable(name = 'ASI', display_name = 'ASI(t)', expression = safe_parse_expr('(100 - PPI) / 100', local_dict = _clash), description = 'Ability of Supply Index (ASI), an index measuing factors affecting a country\'s ability to supply commodities.'),
#     'WSI': Observable(name = 'WSI', display_name = 'WSI(t)', expression = safe_parse_expr('(TT + SV) / 2', local_dict = _clash), description = 'Willingness to Supply Index (WSI), an index measuing factors affecting a country\'s willingness to supply commodities.'),
#     'TT': Observable(name = 'TT', display_name = 'TT(t)', expression = safe_parse_expr('(I + E) / GDP', local_dict = _clash), description = 'Trade ties, the amount of trade between country A and B and is measured as the monetary sum of A\'s imports and exports with A relative to the GDP of B'),
#     'SV': Observable(name = 'SV', display_name = 'SV(t)', expression = safe_parse_expr('FWI', local_dict = _clash), description = 'Shared values, the extent to which the ideological values of countries A, B align and is measured as the Euclidean distance across indicators of political rights and civil liberties (electoral process, political pluralism and participation, functioning of government, freedom of expression and belief, associational and organizational rights, rule of law, and personal autonomy and individual rights), as quantified by Freedom House\'s Freedom in the World (FIW) index')
# }

# Define model
model = TemplateModel(
    templates = [
        NaturalConversion(
            name = 'PrimaryProduction',
            display_name = 'R->S',
            rate_law = safe_parse_expr('PP', local_dict = {'PP': sympy.Symbol('PP')}),
            subject = state_vars['R'],
            outcome = state_vars['S']
        ),
        NaturalConversion(
            name = 'SecondaryProduction',
            display_name = 'OS->S',
            rate_law = safe_parse_expr('SP', local_dict = {'SP': sympy.Symbol('SP')}),
            subject = state_vars['OS'],
            outcome = state_vars['S']
        ),
        NaturalProduction(
            name = 'StockAdjustments',
            display_name = '->S',
            rate_law = safe_parse_expr('DeltaS', local_dict = {'DeltaS': sympy.Symbol('DeltaS')}),
            outcome = state_vars['S']
        ),
        NaturalProduction(
            name = 'GDPGrowth',
            display_name = '->GDP',
            rate_law = safe_parse_expr('g * GDP', local_dict = {'g': sympy.Symbol('g'), 'GDP': sympy.Symbol('GDP')}),
            outcome = state_vars['GDP']
        )
    ],
    initials = initials,
    parameters = parameters,
    # observables = observables,
    annotations = Annotations(name = 'CriticalMineralsModel', description = 'Based on Nassar et al. 2020.'),
    time = Time(name = 't', units = Unit(expression = sympy.Symbol('year')))
)

generate_summary_table(model)

# %%
GraphicalModel.for_jupyter(model)

# %%
with open('./data/example_models/criticalminerals.json', 'w') as fp:
    json.dump(template_model_to_petrinet_json(model), fp, indent = 4)

# %%
# Selectively stratify by country

# countries = ['us', 'ua', 'cn']
countries = ['US', 'UA', 'CN']

model_country = stratify(
    template_model = model,
    key = 'country',
    strata = countries,
    # structure = [['ua', 'us'], ['cn', 'us'], ['ua', 'cn']],
    structure = [],
    directed = False,
    cartesian_control = False,
    # params_to_stratify = ['PPI', 'FWI', 'MC', 'DeltaS', 'PP', 'SP', 'g'],
    params_to_stratify = ['DeltaS', 'PP', 'SP', 'g'],
    concepts_to_stratify = ['GDP', 'R', 'OS', 'S'],
    param_renaming_uses_strata_names = True
)

# Add missing stratified PPI, FWI, MC
for p in ('PPI', 'FWI', 'MC'):
    for c in countries:
        model_country.parameters[f'{p}_{c}'] = parameters[p]
        model_country.parameters[f'{p}_{c}'].name = f'{p}_{c}'
        model_country.parameters[f'{p}_{c}'].display_name = f'{p}_{c}(t)'
        model_country.parameters[f'{p}_{c}'].description = parameters[p].description[:-1] + f' {c}.'

    __ = model_country.parameters.pop(p)

# Import/export trade between country stocks `S_*`
# IE_i_j = trade from country i to country j
for c1, c2 in itertools.permutations(countries, 2):
    model_country.templates.append(
        NaturalConversion(
            subject = model_country.get_concepts_name_map()[f'S_{c1}'],
            outcome = model_country.get_concepts_name_map()[f'S_{c2}'],
            rate_law = safe_parse_expr(f'IE_{c1}_{c2}', local_dict = _clash),
            name = f'Trade_{c1}_{c2}',
            display_name = f'Trade_{c1}_{c2}'
        )
    )
    model_country.parameters[f'IE_{c1}_{c2}'] = Parameter(
        name = f'IE_{c1}_{c2}',
        display_name = f'IE_{c1}_{c2}',
        value = 0.0,
        description = f'Trade from {c1} to {c2}',
        units = Unit(expression = SympyExprStr(safe_parse_expr("1", local_dict = _clash)))
    )

GraphicalModel.for_jupyter(model_country)

# %%
generate_summary_table(model_country)

# %%
# Selectively stratify by industry/application (for trade exposure)

# industries = ['331313', '332431', '336111'] # household/industrial foil, containers/packaging, transportation (passengers and light trucks), 
industries = ['foil', 'containers', 'cars']

# model_country_industry = stratify(
#     template_model = model_country,
#     key = 'industry',
#     strata = industries, 
#     structure = [],
#     directed = False,
#     cartesian_control = False,
#     params_to_stratify = ['VA', 'OP', 'EXP'],
#     concepts_to_stratify = [],
#     param_renaming_uses_strata_names = True
# )

model_country_industry = deepcopy(model_country)
for p in ('VA', 'OP'):
    for ind in industries:
        model_country_industry.parameters[f'{p}_{ind}'] = parameters[p]
        model_country_industry.parameters[f'{p}_{ind}'].name = f'{p}_{ind}'
        model_country_industry.parameters[f'{p}_{ind}'].display_name = f'{p}_{ind}(t)'

    __ = model_country_industry.parameters.pop(p)

GraphicalModel.for_jupyter(model_country_industry)

# %% 
# Selectively stratify by commodity

commodities = ['bauxite', 'alumina', 'aluminum']

model_country_industry_commodity = stratify(
    template_model = model_country_industry,
    key = 'commodity',
    strata = commodities,
    structure = [],
    directed = False,
    cartesian_control = False,
    params_to_stratify = [p for p in model_country_industry.parameters.keys() if p.split('_')[0] in ('DeltaS', 'SP', 'PP', 'IE')],
    concepts_to_stratify = [v for v in model_country_industry.get_concepts_name_map().keys() if v.split('_')[0] in ('R', 'OS', 'S')],
    param_renaming_uses_strata_names = True
)

# Add doubly stratified parameter `EXP`
p = 'EXP'
for ind in industries:
    for com in commodities:
        name = f'{p}_{ind}_{com}'
        model_country_industry_commodity .parameters[name] = parameters[p]
        model_country_industry_commodity .parameters[name].name = name
        model_country_industry_commodity .parameters[name].display_name = f'{name}(t)'

__ = model_country_industry_commodity .parameters.pop(p)


GraphicalModel.for_jupyter(model_country_industry_commodity)

# %%
generate_summary_table(model_country_industry_commodity)

# %%
# Define observables

# Trade ties (TT)
for c1, c2 in itertools.permutations(countries, 2):
    name = f'TT_{c1}_{c2}'
    str = ' + '.join([f'IE_{c1}_{c2}_{com} + IE_{c2}_{c1}_{com}' for com in commodities])
    str = '(' + str + f') / GDP_{c1}'
    model_country_industry_commodity.observables[name] = Observable(
        name = name, 
        display_name = f'{name}(t)', 
        expression = safe_parse_expr(str, local_dict = _clash), 
        description = f'Trade ties, the amount of trade between countries {c1} and {c2} and is measured as the monetary sum of {c1}\'s imports and exports with {c2} relative to the GDP of {c1}.'
    )

# Shared values (SV)
for c1, c2 in itertools.permutations(countries, 2):
    name = f'SV_{c1}_{c2}'
    model_country_industry_commodity.observables[name] = Observable(
        name = name, 
        display_name = f'{name}(t)', 
        expression = safe_parse_expr(f'((FWI_{c1} - FWI_{c2}) ** 2) ** (1/2)', local_dict = _clash), 
        description = 'Shared values, the extent to which the ideological values of countries A, B align and is measured as the Euclidean distance across indicators of political rights and civil liberties (electoral process, political pluralism and participation, functioning of government, freedom of expression and belief, associational and organizational rights, rule of law, and personal autonomy and individual rights), as quantified by Freedom House\'s Freedom in the World (FIW) index'
    )

# Willingness to Supply index (WSI)
for c1, c2 in itertools.permutations(countries, 2):
    name = f'WSI_{c1}_{c2}'
    model_country_industry_commodity.observables[name] = Observable(
        name = name,
        display_name = f'{name}(t)',
        expression = safe_parse_expr(f'(TT_{c1}_{c2} + SV_{c1}_{c2}) / 2', local_dict = _clash), 
        description = 'Willingness to Supply Index (WSI), an index measuing factors affecting a country\'s willingness to supply commodities.'
    )

# Ability to Supply index (ASI)
for c in countries:
    name = f'ASI_{c}'
    model_country_industry_commodity.observables[name] = Observable(
        name = name,
        display_name = name,
        expression = safe_parse_expr(f'1 - (PPI_{c} / 100)', local_dict = _clash), 
        description = 'Ability of Supply Index (ASI), an index measuing factors affecting a country\'s ability to supply commodities.'
    )

# Production Share (PS)
for c in countries:
    name = f'PS_{c}_{com}'
    num = f'PP_{c}_{com} + SP_{c}_{com}'
    denom = ' + '.join([f'PP_{c}_{com} + SP_{c}_{com}' for c in countries])
    expr = safe_parse_expr(f'({num}) / ({denom})', local_dict = _clash)
    model_country_industry_commodity.observables[name] = Observable(
        name = name,
        display_name = name,
        expression = expr, 
        description = 'Production share of a commodity by a country'
    )

# Disruption Potential (DP)
c1 = 'us'
for com in commodities:
    name = f'DP_{com}'
    str = ' + '.join([f'PS_{c2}_{com} * PS_{c2}_{com} * ASI_{c} * WSI_{c1}_{c2}' for c2 in countries])
    expr = safe_parse_expr(str, local_dict = _clash)
    model_country_industry_commodity.observables[name] = Observable(
        name = name,
        display_name = name,
        expression = expr, 
        description = f'Disruption potential of a commodity {com} to a country {c1} relative to all others.'
    )

# Trade Exposure (TE)
for c1 in countries:
    for com in commodities:
        name = f'TE_{c1}_{com}'
        num = f'PP_{c1}_{com} + SP_{c1}_{com}'
        denom = f'DeltaS_{c1}_{com}' + ' + '.join([f'IE_{c1}_{c2}_{com} - IE_{c2}_{c1}_{com}' for c2 in countries if c2 != c1])
        expr = safe_parse_expr(f'1 / (1 + ({num}) / ({denom}))', local_dict = _clash)
        model_country_industry_commodity.observables[name] = Observable(
            name = name,
            display_name = name,
            expression = expr, 
            description = f'Trade exposure of a commodity {com} for a country {c1} relative to all others.'
        )

# Economic Vulnerability (EV)
c = 'us'
for com in commodities:
    name = f'EV_{c}_{com}'
    num = ' + '.join([f'VA_{ind} * EXP_{ind}_{com} / OP_{ind}' for ind in industries])
    denom = f'GDP_{c}'
    expr = safe_parse_expr(f'({num}) / ({denom})', local_dict = _clash)
    model_country_industry_commodity.observables[name] = Observable(
        name = name,
        display_name = name,
        expression = expr, 
        description = f'Economic vulnerability of a commodity {com} for a country {c} relative to all others.'
    )

# Supply Risk (SR)
c = 'us'
for com in commodities:
    name = f'SR_{c}_{com}'
    expr = safe_parse_expr(f'(DP_{com} * TE_{c}_{com} * EV_{c}_{com}) ** (1/3)', local_dict = _clash)
    model_country_industry_commodity.observables[name] = Observable(
        name = name,
        display_name = name,
        expression = expr,
        description = f'Supply risk of a commodity {com} for a country {c} relative to all others.'
    )

# %%
model_country_industry_commodity.observables

# %%
with open('./data/example_models/CriticalMineralsModel_gollmcard.json', 'r') as fp:
    j_gollmcard = json.load(fp)

j = template_model_to_petrinet_json(model_country_industry_commodity)
j['metadata'] = j_gollmcard['metadata']

with open('./data/example_models/criticalminerals_stratified.json', 'w') as fp:
    json.dump(j, fp, indent = 4)

# %%
