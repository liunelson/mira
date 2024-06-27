# %%[markdown]
# Author: Nelson Liu
# 
# Email: [nliu@uncharted.software](mailto:nliu@uncharted.software)

# %%[markdown]
# Convert All SBML Models in BioModels
# 
# * Scrapped all BioModels models found by filtering for SBML format
# * Convert from SBML to AMR JSON
# * Download paper PDFs
# * Upload to Terarium

# %%
import os
import glob
import json
from tqdm import tqdm
import requests
import time
import shutil

from mira.metamodel.ops import simplify_rate_laws
from mira.modeling import Model
from mira.modeling.amr.petrinet import AMRPetriNetModel
from mira.sources.sbml import template_model_from_sbml_file

from mira.modeling.amr.petrinet import template_model_to_petrinet_json
from mira.sources.amr.petrinet import model_from_json_file

# %%
PATH = "./_data/biomodels/biomodels_covid19"
URL_OA = "https://api.openaccessbutton.org"
URL_BIOMODELS = "https://www.ebi.ac.uk/biomodels"

# %%
# Search for all COVID-19 population models in SBML format
# https://www.ebi.ac.uk/biomodels/search?offset=0&numResults=50&sort=relevance-desc&query=*%3A*+AND+modelformat%3A%22SBML%22+AND+modellingapproach%3A%22population+model%22&domain=biomodels

url = f"{URL_BIOMODELS}/search"
params = {
    "query": '*:* AND modellingapproach:"population model" AND modelformat:"SBML"', 
    "domain": "biomodels",
    "offset": 0,
    "numResults": 100,
    "format": "json"
}
r = requests.get(url, params = params)
if r.ok:
    print(f"Number of matches: {r.json()['matches']}")
    models = {m["id"]: {} for m in r.json()["models"]}

# %%
# Get model metadata
for id in tqdm(models.keys()):
    url = f"{URL_BIOMODELS}/{id}"
    params = {"format": "json"}
    r = requests.get(url, params = params)
    if r.ok:
        models[id] = r.json()

# %%
with open(os.path.join(PATH, "model_data.json"), "w") as f:
    json.dump(models, f, indent = 4)

# %%
# Pull model files from previous cache
for id in tqdm(models.keys()):

    # Create directory
    p = os.path.join(PATH, id)
    if not os.path.exists(p):
        os.makedirs(p)

    for ext in ("xml", "json"):
        dst = os.path.join(PATH, id, f"{id}.{ext}")
        src = os.path.join("./_data/biomodels/biomodels_sbml_all", id, f"{id}.{ext}")
        try:
            shutil.copy(src, dst)
        except:
            print(f"{src} not found")

# %%
# Update model metadata
for id in tqdm(models):

    p = os.path.join(PATH, id, f"{id}.json")
    with open(p, "r") as fp:
        m = json.load(fp)

    m["header"]["description"] = models[id]["description"]
    m["metadata"]["annotations"]["authors"] = models[id]["publication"]["authors"]

    with open(p, "w") as fp:
        json.dump(m, fp, indent = 4)


# %%
# Get the open access URL if available
num_oa = 0
for id in tqdm(models.keys()):

    url = models[id]["publication"]["link"]

    if len(url) > 0:

        r = requests.get(url = URL_OA + "/find", params = {"id": url})

        if r.ok:
            if "url" in r.json().keys():
                num_oa += 1
                models[id]["publication"]["link_oa"] = r.json()["url"]

            else:
                print(f"Model {id} missing OA link")
                models[id]["publication"]["link_oa"] = None

    time.sleep(1.001)

print(f"Number of models with OA publication: {num_oa / len(models) * 100:.2f} % of models")

# %%
# Download publication PDF if open access URL is available

headers = {
    "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36 Edg/120.0.0.0",
    "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.7"
}

num_pdf = 0
pdf_fail = []
for id in tqdm(models):

    url = models[id]["publication"]["link_oa"]
    fp = os.path.join(PATH, id, f"{id}.pdf")

    try:
        r = requests.get(url, headers = headers)
        if r.ok:
            num_pdf += 1
            with open(fp, "wb") as f:
                f.write(r.content)
        else:
            pdf_fail.append(id)
    
    except:
        pdf_fail.append(id)

print(f"Number of publication PDF downloaded: {num_pdf} of {len(models)} of models ({num_pdf / len(models) * 100})")


# %%
# Manual link search
# pdf_fail = ['BIOMD0000001045', 'BIOMD0000000970', 'MODEL2212310001', 'BIOMD0000000956']


# %%
# Stats

num = {k: 0 for k in ["models", "xml", "json", "pub", "oa", "pdf"]}
for m in tqdm(models.keys()):
    for k in num.keys():

        if k == "models":
            num[k] += 1

        if (k == "pub") & (models[m]["publication"]["link"] != None):
            num[k] += 1

        if (k == "oa") & (models[m]["publication"]["link_oa"] != None):
            num[k] += 1

        fp = os.path.join(PATH, m, f"{m}.{k}")
        if os.path.isfile(fp):
            num[k] += 1


print(f"Number of models:\t\t{num['models']}")
print(f"\twith SBML:\t\t{num['xml'] / num['models'] * 100:.1f}%")
print(f"\tconvertable to AMR:\t{num['json'] / num['models'] * 100:.1f}%")
print(f"\twith PDF link:\t\t{num['pub'] / num['models'] * 100:.1f}%")
print(f"\twith OA PDF link:\t{num['oa'] / num['models'] * 100:.1f}%")
print(f"\twith downloaded PDF:\t{num['pdf'] / num['models'] * 100:.1f}%")

# %%
