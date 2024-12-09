# %%
import os

from mira.metamodel import *
from mira.modeling import Model
from mira.modeling.amr.petrinet import template_model_to_petrinet_json
from mira.modeling.amr.regnet import template_model_to_regnet_json
from mira.modeling.amr.regnet import AMRRegNetModel

from mira.modeling.viz import GraphicalModel

# %%
MIRA_REST_URL = "http://mira-epi-dkg-lb-dc1e19b273dedaa2.elb.us-east-1.amazonaws.com"
os.environ["MIRA_REST_URL"] = MIRA_REST_URL

# %%
# Define concepts
# susceptible = Concept(name="S", identifiers={"ido": "0000514"})
# infected = Concept(name="I", identifiers={"ido": "0000511"})
# recovered = Concept(name="R", identifiers={"ido": "0000592"})

susceptible = Concept(name="S")
infected = Concept(name="I")
recovered = Concept(name="R")


# %%
model_A = TemplateModel(
    templates = [
        ControlledConversion(
            name = 'Infection',
            subject = susceptible,
            outcome = infected,
            controller = infected
        )
    ],
    parameters = {
        'b': Parameter(name = 'b', value = 1.0)
    }
)

model_A.templates[0].set_mass_action_rate_law('b')

model_A.draw_jupyter()

# %%
model_B = TemplateModel(
    templates = [
        NaturalConversion(
            name = 'Recovery',
            subject = infected,
            outcome = recovered
        )
    ],
    parameters = {
        'g': Parameter(name = 'g', value = 1.0)
    }
)

model_B.templates[0].set_mass_action_rate_law('g')

model_B.draw_jupyter()

# %%
# Compose the two models together

model_AB = compose([model_A, model_B])

model_AB.draw_jupyter()
