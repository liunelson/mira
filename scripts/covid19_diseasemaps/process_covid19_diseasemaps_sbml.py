import json
import tqdm
import requests

from mira.sources.sbml import template_model_from_sbml_qual_string
from mira.modeling.amr.regnet import template_model_to_regnet_json
from mira.sources.biomodels import get_sbml_model



models = [
    "Apoptosis",
    "Coagulation-pathway",
    "ER_Stress",
    "ETC",
    "E_protein",
    "HMOX1_Pathway",
    "IFN-lambda",
    "Interferon1",
    "JNK_pathway",
    "Kynurenine_pathway",
    "NLRP3_Activation",
    "Nsp14",
    "Nsp4_Nsp6",
    "Nsp9_protein",
    "Orf10_Cul2_pathway",
    "Orf3a",
    "PAMP_signaling",
    "Pyrimidine_deprivation",
    "RTC-and-transcription",
    "Renin_angiotensin",
    "TGFB_pathway",
    "Virus_replication_cycle",
]


SBML_URL_BASE = (
    "https://git-r3lab.uni.lu/covid/models/-/raw/master/"
    "Executable%20Modules/SBML_qual_build/sbml"
)


if __name__ == "__main__":
    for model in tqdm.tqdm(models):
        url = f"{SBML_URL_BASE}/{model}_stable.sbml"
        model_text = requests.get(url).text
        tm = template_model_from_sbml_qual_string(model_text)
        regnet = template_model_to_regnet_json(tm)
        with open(f"regnet_amr_sbml/{model}.json", "w") as fh:
            json.dump(regnet, fh, indent=1)

    model_id = "BIOMD0000000562"
    model_text = get_sbml_model(model_id)
    tm = template_model_from_sbml_qual_string(model_text)
    regnet = template_model_to_regnet_json(tm)
    with open(f"regnet_amr_sbml/{model_id}.json", "w") as fh:
        json.dump(regnet, fh, indent=1)
