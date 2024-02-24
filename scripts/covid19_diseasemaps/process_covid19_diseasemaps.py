import json
import tqdm
from mira.metamodel import Annotations
from mira.modeling.amr.regnet import template_model_to_regnet_json
from mira.sources.sif import template_model_from_sif_url


models = ['Apoptosis', 'Coagulation-pathway', 'ER_Stress', 'ETC', 'E_protein',
          'HMOX1_Pathway', 'IFN-lambda', 'Interferon1', 'JNK_pathway',
          'Kynurenine_pathway', 'NLRP3_Activation', 'Nsp14', 'Nsp4_Nsp6',
          'Nsp9_protein', 'Orf10_Cul2_pathway', 'Orf3a', 'PAMP_signaling',
          'Pyrimidine_deprivation', 'RTC-and-transcription',
          'Renin_angiotensin', 'TGFB_pathway', 'Virus_replication_cycle']


SIF_URL_BASE = ('https://git-r3lab.uni.lu/covid/models/-/raw/master/'
                'Executable%20Modules/SBML_qual_build/sif')


if __name__ == "__main__":
    for model in tqdm.tqdm(models):
        url = f'{SIF_URL_BASE}/{model}_stable.sif'
        tm = template_model_from_sif_url(url)
        tm.annotations = Annotations(
            name=model,
            description=("A submodel of the COVID-10 Disease Map "
                         f"focusing on {model}"),
            pathogens=["ncbitaxon:2697049"],
            diseases=["doid:0080600"],
            hosts=["ncbitaxon:9606"],
            references=["pubmed:34664389"]
        )
        regnet = template_model_to_regnet_json(tm)
        with open(f'regnet_amr/{model}.json', 'w') as fh:
            json.dump(regnet, fh, indent=1)
