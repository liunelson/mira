# %%[markdown]
# Author: Nelson Liu
# 
# Email: [nliu@uncharted.software](mailto:nliu@uncharted.software)

# %%[markdown]
# Test Importing SBML Models with MIRA
# 
# 

# %%
import os
import glob
import json
import tqdm

from mira.metamodel.ops import simplify_rate_laws
from mira.modeling import Model
from mira.modeling.amr.petrinet import AMRPetriNetModel
from mira.sources.sbml import template_model_from_sbml_file
# from mira.sources.biomodels import query_biomodels, get_template_model

# %%
PATH = "data/biomodels/biomodels_sbml_popmodel" # 25/25
# PATH = "data/biomodels/biomodels_sbml_ode" # 61/100
# PATH = "data/biomodels/biomodels_sbml_constraint" # 0
# PATH = "data/biomodels/biomodels_sbml_pn" # 0
# PATH = "data/biomodels/biomodels_sbml_phys" # 6/8
# PATH = "data/biomodels/biomodels_sbml_meta" # 3/10
# PATH = "data/biomodels/biomodels_sbml_ppint" # 1/1
# PATH = "data/biomodels/biomodels_sbml_math" # 7/8

fnames = glob.glob(os.path.join(PATH, "*.*ml"))

fnames_succ = []
fnames_fail = []
for fname in tqdm.tqdm(fnames):
    try:
        model_tm = template_model_from_sbml_file(fname)
        model_tm_ = simplify_rate_laws(model_tm)
        model_pn = AMRPetriNetModel(Model(model_tm_))
        model_pn_json = model_pn.to_json()

        with open(".".join([fname.split(".")[0], "json"]), "w") as f:
            json.dump(model_pn_json, f, indent = 4)

        fnames_succ.append(fname)

    except:
        fnames_fail.append(fname)

print(f"{len(fnames_succ)} successes and {len(fnames_fail)} fails")

# %%
# https://www.ebi.ac.uk/biomodels/search?query=*%3A*%20AND%20modelformat%3A%22SBML%22&domain=biomodels&offset=0&numResults=10

