# %%[markdown]
# Author: Nelson Liu
# 
# Email: [nliu@uncharted.software](mailto:nliu@uncharted.software)

# %%[markdown]
# # Demonstrate Geo-stratification of a Model
# 
# Reference: [Capturing the Effects of Transportation on the Spread of COVID-19 With a Multi-Networked SEIR Model](https://ieeexplore.ieee.org/document/9319714)

# %%
import os
import json
import numpy
import pandas
import sympy
from tqdm import tqdm
from datetime import date
import requests
import copy

from mira.modeling.viz import GraphicalModel
from mira.metamodel import *
from mira.modeling import Model
from mira.modeling.amr.petrinet import template_model_to_petrinet_json
from mira.sources.amr.petrinet import template_model_from_amr_json

# %%
MIRA_DKG_REST_URL = 'http://mira-epi-dkg-lb-c7b58edea41524e6.elb.us-east-1.amazonaws.com/'

# %%[markdown]
# ## Adjacency Matrix

# %%
# County adjacency
adj_county = pandas.read_csv(
    './data/COVID_19_aviation_data/Pre_processing/COVID19_datasets/county_adjacency2010.csv',
    dtype = {
        'countyname': str,
        'fipscounty': str,
        'neighborname': str,
        'fipsneighbor': str
    }
)

# %%
# State FIPS codes
x = {}
for __, row in adj_county.iterrows():
    x[row['countyname'].split(',')[-1].strip()] = row['fipscounty'][:2]
state_fips = pandas.DataFrame({'stateAbbr': x.keys(), 'stateFIPS': x.values()})

# %%
# State adjacency
adj_state = numpy.zeros((len(state_fips), len(state_fips)))
for __, row in adj_county.iterrows():
    i = pandas.Index(state_fips['stateFIPS']).get_loc(row['fipscounty'][:2])  
    j = pandas.Index(state_fips['stateFIPS']).get_loc(row['fipsneighbor'][:2])
    adj_state[i, j] += 1

# %%[markdown]
# ## State adjacency by flight network

# %%
p = './data/COVID_19_aviation_data/Pre_processing/Transportation_air_network_datasets'
flight_filenames = [os.path.join(dirpath, f) for [dirpath, __, filenames] in os.walk(p) for f in filenames if f.endswith('0.csv')]

dtype = {
    'ORIGIN_STATE_FIPS': str,
    'DEST_STATE_FIPS': str
}

df_flight = pandas.DataFrame()
for f in tqdm(flight_filenames):
    df_flight = pandas.concat([df_flight, pandas.read_csv(f, dtype = dtype)], ignore_index = True)

# %%
# Correction:
# GUM = Guam = = 66, SPN = Saipan = Northern Mariana Islands = 69, PPG = American Samoa = 60

for c, fips, name in [('GUM', '66', 'Guam'), ('SPN', '69', 'Northern Mariana Islands'), ('PPG', '60', 'American Samoa')]:
    i = df_flight['ORIGIN'] == c
    df_flight.loc[i, 'ORIGIN_STATE_FIPS'] = fips
    df_flight.loc[i, 'ORIGIN_STATE_NM'] = name

    i = df_flight['DEST'] == c
    df_flight.loc[i, 'DEST_STATE_FIPS'] = fips
    df_flight.loc[i, 'DEST_STATE_NM'] = name

# %%
adj_state_flight = numpy.zeros((len(state_fips), len(state_fips)))
for __, row in tqdm(df_flight.iterrows(), total = df_flight.shape[0]):
    i = pandas.Index(state_fips['stateFIPS']).get_loc(row['ORIGIN_STATE_FIPS'])  
    j = pandas.Index(state_fips['stateFIPS']).get_loc(row['DEST_STATE_FIPS'])
    adj_state_flight[i, j] += 1

# Normalize to number of flights per day
num_days = (date(2020, 7, 31) - date(2020, 1, 1)).days
adj_state_flight = adj_state_flight / num_days

# %%
import networkx as nx
G = nx.Graph()
G.add_nodes_from(state_fips['stateAbbr'].values)
for i in range(len(state_fips)):
    for j in range(len(state_fips)):
        if adj_state_flight[i, j] > 0:
            G.add_edge(state_fips.loc[i, 'stateAbbr'], state_fips.loc[j, 'stateAbbr'], weight = adj_state_flight[i, j])

nx.draw_networkx(G, pos = nx.kamada_kawai_layout(G, weight = 'weight'), with_labels = True)

# %%
# Get FIPS-Name map
x = {fips: name for __, (fips, name) in df_flight.loc[:, ['ORIGIN_STATE_FIPS', 'ORIGIN_STATE_NM']].iterrows()}
x = x | {fips: name for __, (fips, name) in df_flight.loc[:, ['DEST_STATE_FIPS', 'DEST_STATE_NM']].iterrows()}

# %%
# Missing
x['10'] = 'Delaware'
x['11'] = 'District of Columbia'

state_fips['stateName'] = [x[fips] for fips in state_fips['stateFIPS']]

# %%
y = []
for name in tqdm(state_fips['stateName']):
    r = requests.get(f'{MIRA_DKG_REST_URL}/api/ground/{name}')
    if r.status_code == 200:
        if len(r.json()['results']) > 0:
            y.append(r.json()['results'][0]['curie'])
        else:
            print(f'{name} - {r.json()["results"]}')
            y.append('')

    else:
        print(f'{r.status_code}\t{name}')

state_fips['stateCurie'] = y

# %%
# Missing states in Epi DKG
for name, curie in (('Delaware', 'geonames:4142224'), ('Vermont', 'geonames:5242283'), ('Maine', 'geonames:4971068'), ('Wyoming', 'geonames:5843591'), ('West Virginia', 'geonames:4826850')):
    i = pandas.Index(state_fips['stateName']).get_loc(name)
    state_fips.loc[i, 'stateCurie'] = curie

# %%
state_fips.to_csv('./data/example_models/state_fips.csv')

df = pandas.DataFrame(adj_state, index = state_fips['stateCurie'].values, columns = state_fips['stateCurie'].values)
df.to_csv('./data/example_models/adj_state.csv')

df = pandas.DataFrame(adj_state_flight, index = state_fips['stateCurie'].values, columns = state_fips['stateCurie'].values)
df.to_csv('./data/example_models/adj_state_flight.csv')

# %%
# Load a SEIRHD model
with open("./data/monthly_demo_202408/model_seirhd.json", "r") as f:
    model = template_model_from_amr_json(json.load(f))

GraphicalModel.for_jupyter(model)

# %%
model.annotations = Annotations(
    name = 'SEIRHD model',
    description = 'Edit of SIDARTHE model from Giordano 2020, stratified by U.S. states and territories to consider geographical adjacency and Jan-Jul 2020 flight network',
    authors = [{'name': 'Nelson Liu'}],
    references = ['10.1038/s41591-020-0883-7', '10.1109/LCSYS.2021.3050954'],
    time_scale = 'day',
    time_start = '2020-01-01T00:00:00',
    time_end = '2020-07-01T00:00:00',
    locations = [],
    pathogens = ["ncbitaxon:2697049"],
    diseases = ["doid:0080600"],
    hosts = ["ncbitaxon:9606"],
    model_types = ["mamo:0000028", "mamo:0000046"]
)

# %%
# Stratify the model by the entities in `state_fips`

model_travel_stateName = stratify(
    model,
    key = 'state',
    strata = [name.replace(' ', '').replace('.', '') for name in state_fips['stateName'].values],
    cartesian_control = True,
    directed = False,
    structure = None,
    concepts_to_stratify = ['S', 'E', 'I'],
    params_to_stratify =['b'],
    param_renaming_uses_strata_names = True
)

# %%
# Reset the initial expressions
for k, v in model_travel_stateName.initials.items():
    if len(k.split('_')) == 2:
        x, name = k.split('_')
        v.expression = sympy.Symbol(f'{x}0_{name}')

        # Add parameters
        model_travel_stateName.parameters = model_travel_stateName.parameters | {
            f'{x}0_{name}': Parameter(
                name = f'{x}0_{name}', 
                display_name = f'{x}0_{name}', 
                description = f'Initial {x} pop in {name}', 
                identifiers = v.concept.identifiers, 
                context = v.concept.context, 
                units = v.concept.units, 
                value = 1.0, 
            )
        }
    
# %%
# With Geonames curies
model_travel_stateCurie = copy.deepcopy(model_travel_stateName)

# %%
# Update annotations.locations
model_travel_stateCurie.annotations.locations = ['geonames:6252001'] + list(state_fips['stateCurie'].values),

# %%
m = {s['stateName'].replace(' ', '').replace('.', ''): s['stateCurie'] for __, s in state_fips.iterrows()}

# Change concept context from name to curie
for k, v in model_travel_stateCurie.get_concepts_name_map().items():
    if len(k) > 1:
        x, name = k.split('_')
        curie = m[name]
        v.context['location'] = curie

        # update initials concept
        model_travel_stateCurie.initials[k].concept.context['location'] = curie

# %%
# Update the initial parameters
for k, v in model_travel_stateCurie.initials.items():
    if len(k.split('_')) == 2:
        x, name = k.split('_')
        p = f'{x}0_{name}'
        curie = m[name]
        model_travel_stateCurie.parameters[p].context['location'] = curie
        
# %%
# Add curie to other parameteres
for k, p in model_travel_stateCurie.parameters.items():
    if k[:2] == 'b_':
        __, name_0, name_1 = k.split('_')
        p.context['location_0'] = m[name_0]
        p.context['location_1'] = m[name_1]

# %%
# GraphicalModel.for_jupyter(model_travel)

# %%
with open("./data/example_models/model_seirhd_travel_stateName.json", "w") as f:
    json.dump(template_model_to_petrinet_json(model_travel_stateName), f, indent = 4)

with open("./data/example_models/model_seirhd_travel_stateCurie.json", "w") as f:
    json.dump(template_model_to_petrinet_json(model_travel_stateCurie), f, indent = 4)

# %%
