# %%[markdown]
# Author: Nelson Liu
# 
# Email: [nliu@uncharted.software](mailto:nliu@uncharted.software)

# %%[markdown]
# # Test MIRA Utilities
# 
# Try to build a compartmental model for bird-cow-human infection model across US states.

# %%
import os
import json
import numpy
import pandas
import sympy
from typing import Optional, Iterable
import itertools

from mira.modeling.viz import GraphicalModel
from mira.metamodel import *
from mira.modeling import Model
from mira.modeling.amr.petrinet import template_model_to_petrinet_json
from mira.sources.amr.petrinet import template_model_from_amr_json

# %%


