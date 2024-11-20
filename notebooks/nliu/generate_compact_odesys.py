# %%[markdown]
# Author: Nelson Liu
# 
# Email: [nliu@uncharted.software](mailto:nliu@uncharted.software)

# %%[markdown]
# # Generate ODE System of Equations from PetriNet Model (Compact)
# 
# Use MIRA to generate template model with Sympy rate laws 
# and assemble into ODE system of equations.
# 
# Use stratification metadata to use summation notation.

# %%
import os
import json
import sympy
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcParams['text.usetex'] = True

# from mira.modeling.viz import GraphicalModel
# from mira.metamodel import *
# from mira.modeling import Model
from mira.sources.amr.petrinet import template_model_from_amr_json
# from mira.modeling.amr.petrinet import template_model_to_petrinet_json

# %%
PATH = "./data/stratify_test_cases"

# %%
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
# Load example model
models_amr = []
models = []
for model_name in ('SIR.json', 'SIR (age).json', 'SIR (age, vax).json'):
    with open(os.path.join(PATH, model_name), "r") as fp:
        models_amr.append(json.load(fp))
        models.append(template_model_from_amr_json(models_amr[-1]))

# %%
# Generate LaTeX code (one per equation)
for model in models:
    odesys = generate_odesys(model, latex = True, latex_align = False)

    for eq in odesys:
        print(eq)
    print('\n')

# %%
model = models[1]














# %%
# Generate LaTeX code (one for entire set of equations)
print(generate_odesys(model, observables = True, latex = True, latex_align = True))

# %%
odesys = generate_odesys(model, observables = True, latex = True, latex_align = False)

# Remove `\\operatorname for plot
odesys = [re.sub(r"(\\operatorname)", "", eq) for eq in odesys]
odesys

# %%
fig, ax = plt.subplots(1, 1, figsize = (8, 8))
for y, eq in zip(np.flip(np.linspace(0, 1, len(odesys))), odesys):
    __ = ax.text(0.0, y, "$" + eq + "$", fontsize = 12)
__ = ax.axis("off")

# %%

