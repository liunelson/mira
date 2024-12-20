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
S1 = Concept(name="Susceptible", identifiers={"ido": "0000514"})
I1 = Concept(name="Infected", identifiers={"ido": "0000511"})
R1 = Concept(name="Recovery", identifiers={"ido": "0000592"})

S2 = Concept(name="Susceptible", identifiers={"ido": "0000513"})
I2 = Concept(name="Infected", identifiers={"ido": "0000512"})
R2 = Concept(name="Recovery", identifiers={"ido": "0000593"})

S3 = Concept(name="S", identifiers={"ido": "0000514"})
I3 = Concept(name="I", identifiers={"ido": "0000511"})
R3 = Concept(name="R", identifiers={"ido": "0000592"})

S4 = Concept(name="S")
I4 = Concept(name="I")
R4 = Concept(name="R")

# %%
model_A1 = TemplateModel(
    templates = [
        ControlledConversion(
            name = 'Infection',
            subject = S1,
            outcome = I1,
            controller = I1
        )
    ],
    parameters = {
        'b': Parameter(name = 'b', value = 1.0)
    }
)

model_A1.templates[0].set_mass_action_rate_law('b')

model_A1.draw_jupyter()

# %%
model_B1 = TemplateModel(
    templates = [
        NaturalConversion(
            name = 'Recovery',
            subject = I1,
            outcome = R1
        )
    ],
    parameters = {
        'g': Parameter(name = 'g', value = 1.0)
    }
)

model_B1.templates[0].set_mass_action_rate_law('g')

model_B1.draw_jupyter()

# %%
model_B2 = TemplateModel(
    templates = [
        NaturalConversion(
            name = 'Recovery',
            subject = I2,
            outcome = R2
        )
    ],
    parameters = {
        'g': Parameter(name = 'g', value = 1.0)
    }
)

model_B2.templates[0].set_mass_action_rate_law('g')

model_B2.draw_jupyter()

# %%
model_B3 = TemplateModel(
    templates = [
        NaturalConversion(
            name = 'Recovery',
            subject = I3,
            outcome = R3
        )
    ],
    parameters = {
        'g': Parameter(name = 'g', value = 1.0)
    }
)

model_B3.templates[0].set_mass_action_rate_law('g')

model_B3.draw_jupyter()

# %%
model_B4 = TemplateModel(
    templates = [
        NaturalConversion(
            name = 'Recovery',
            subject = I4,
            outcome = R4
        )
    ],
    parameters = {
        'g': Parameter(name = 'g', value = 1.0)
    }
)

model_B4.templates[0].set_mass_action_rate_law('g')

model_B4.draw_jupyter()

# %%
model_A3 = TemplateModel(
    templates = [
        ControlledConversion(
            name = 'Infection',
            subject = S3,
            outcome = I3,
            controller = I3
        )
    ],
    parameters = {
        'b': Parameter(name = 'b', value = 1.0)
    }
)

model_A3.templates[0].set_mass_action_rate_law('b')
model_A3.draw_jupyter()

# %%
model_A4 = TemplateModel(
    templates = [
        ControlledConversion(
            name = 'Infection',
            subject = S4,
            outcome = I4,
            controller = I4
        )
    ],
    parameters = {
        'b': Parameter(name = 'b', value = 1.0)
    }
)

model_A4.templates[0].set_mass_action_rate_law('b')
model_A4.draw_jupyter()

# %%
# Compose the two models together

# (matching name and ids)
model_AB11 = compose([model_A1, model_B1])
model_AB11.draw_jupyter()

# composed

# %%
# matching names, mismatched ids
model_AB12 = compose([model_A1, model_B2])
model_AB12.draw_jupyter()

# no composition

# %%
# mismatched names, matching ids
model_AB13 = compose([model_A1, model_B3])
model_AB13.draw_jupyter()

# no composition

# %%
# matching names, no id + yes id
model_AB34 = compose([model_A3, model_B4])
model_AB34.draw_jupyter()

# no composition

# %%
# matching names, no id + no id
model_AB44 = compose([model_A4, model_B4])
model_AB44.draw_jupyter()

# composed

# %%
