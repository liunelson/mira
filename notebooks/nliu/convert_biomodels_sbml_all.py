# %%[markdown]
# Author: Nelson Liu
# 
# Email: [nliu@uncharted.software](mailto:nliu@uncharted.software)

# %%[markdown]
# Convert All SBML Models in BioModels
# 
# * Scrapped all BioModels models found by filtering for SBML format
# * Convert from SBML to AMR JSON
# * Upload to Terarium

# %%
import os
import glob
import json
import tqdm
import requests
import time

from mira.metamodel.ops import simplify_rate_laws
from mira.modeling import Model
from mira.modeling.amr.petrinet import AMRPetriNetModel
from mira.sources.sbml import template_model_from_sbml_file

# %%
PATH = "data/biomodels/biomodels_sbml_all"
URL_OA = "https://api.openaccessbutton.org"

# %%
# Load all model metadata
models = {}
# fp = "model_data.json"
fp = "model_data_oa.json"
with open(os.path.join(PATH, fp), "r") as f:
    models = json.load(f)

# %%
# Download the SBML file
num_sbml = 0
for m in tqdm.tqdm(models.keys()):

    fp = os.path.join(PATH, m, f"{m}.xml")
    if ~os.path.isfile(fp):
        
        if len(models[m]["model_files"]) > 0:

            url = models[m]["model_files"][0]
            r =requests.get(url)
            if r.ok:

                p = os.path.join(PATH, m)
                if not os.path.exists(p):
                    os.makedirs(p)

                with open(fp, "wb") as f:
                    f.write(r.content)
                
                num_sbml += 1

print(f"Number of downloaded model files: {num_sbml / len(models) * 100:.2f} % of models")
# X % of models

# %%
# Convert SBML to AMR

sbml_amr_succ = []
sbml_amr_fail = []
for m in tqdm.tqdm(models.keys()):
    
    fp = os.path.join(PATH, m, f"{m}.xml")
    fp_amr = os.path.join(PATH, m, f"{m}.json")

    if os.path.isfile(fp) & ~os.path.isfile(fp_amr):

        try:
            model_tm = template_model_from_sbml_file(fp)
            model_tm_ = simplify_rate_laws(model_tm)
            model_pn = AMRPetriNetModel(Model(model_tm_))
            model_pn_json = model_pn.to_json()

            with open(fp_amr, "w") as f:
                json.dump(model_pn_json, f, indent = 4)

            sbml_amr_succ.append(m)

            models[m]["amr"] = True

        except:
            sbml_amr_fail.append(m)
            models[m]["amr"] = False

print(f"Number of SBML-to-AMR conversion: {len(sbml_amr_succ) / len(models) * 100:.2f} % of models")
# 50.39 % of models

# %%
# Get the open access URL and metadata of publication 
# if available

num_oa = 0
for m in tqdm.tqdm(models.keys()):

    models[m]["publication_link_oa"] = ""
    models[m]["publication_metadata"] = {}

    url = models[m]["publication_link"]

    if len(url) > 0:

        r = requests.get(url = URL_OA + "/find", params = {"id": url})

        if r.ok:
            if "url" in r.json().keys():
                num_oa += 1
                models[m]["publication_link_oa"] = r.json()["url"]

            if "metadata" in r.json().keys():
                models[m]["publication_metadata"] = r.json()["metadata"]

    time.sleep(1.001)

print(f"Number of models with OA publication: {num_oa / len(models) * 100:.2f} % of models")

# %%
# Save open access data
with open(os.path.join(PATH, "model_data_oa.json"), "w") as f:
    json.dump(models, f, indent = 4)

# %%
# Download publication PDF if open access URL is available

headers = {
    "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36 Edg/120.0.0.0",
    "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.7"
}

# %%
num_pdf = 0
pdf_fail = []
for m in tqdm.tqdm(models.keys()):

    url = models[m]["publication_link_oa"]
    fp = os.path.join(PATH, m, f"{m}.pdf")

    if len(url) == 0:
        continue

    if os.path.isfile(fp):
        continue

    try:
        r = requests.get(url, headers = headers)
        if ~r.ok:
            continue

        num_pdf += 1
        with open(fp, "wb") as f:
            f.write(r.content)
    
    except:
        pdf_fail.append(m)

print(f"Number of publication PDF downloaded: {num_pdf / len(models) * 100} % of models")

# Failed: 'MODEL2003050002', 'BIOMD0000000904'

# %%
# Stats

num = {k: 0 for k in ["models", "xml", "json", "pub", "oa", "pdf"]}
for m in tqdm.tqdm(models.keys()):
    for k in num.keys():

        if k == "models":
            num[k] += 1

        if (k == "pub") & (len(models[m]["publication_link"]) > 0):
            num[k] += 1

        if (k == "oa") & (len(models[m]["publication_link_oa"]) > 0):
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
# models with OA PDF link but no downloaded PDF
x = [m for m in models if (len(models[m]["publication_link_oa"]) > 0) & (~os.path.isfile(os.path.join(PATH, m, f"{m}.pdf")))]

# # %%
# from mira.sources.amr.petrinet import template_model_from_amr_json

# PATH = "./data/biomodels/biomodels_sbml_all"

# models_amr = []

# for d in tqdm.tqdm(os.listdir(PATH)):
#     if os.path.isdir(os.path.join(PATH, d)):
#         l = os.listdir(os.path.join(PATH, d))
#         for f in l:
#             if f.split(".")[1] == "json":
#                 with open(os.path.join(PATH, d, f), "r") as f:
#                     models_amr.append(json.load(f))

# models_tm = []
# for m in tqdm.tqdm(models_amr):
#     try:
#         models_tm.append(template_model_from_amr_json(m))
#     except:
#         models_tm.append(None)


# %%



