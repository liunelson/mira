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
import unittest
import tqdm

from mira.metamodel.ops import simplify_rate_laws
from mira.modeling import Model
from mira.modeling.acsets.petri import PetriNetModel
from mira.sources.sbml import template_model_from_sbml_file

# %%
from mira.sources.biomodels import query_biomodels, get_template_model

# %%







