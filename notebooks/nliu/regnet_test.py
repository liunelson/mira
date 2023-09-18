# %%[markdown]
# Author: Nelson Liu
# 
# Email: [nliu@uncharted.software](mailto:nliu@uncharted.software)

# %%
# Testing MIRA RegNet interfaces

# %%
import requests
import json
import sympy
import numpy as np
import scipy as sp
from PIL import Image
from typing import Optional, NoReturn
import networkx as nx
import matplotlib as mpl
import matplotlib.pyplot as plt

from mira.metamodel import *
import mira.sources.askenet as askenet

# %%
# Get different model AMRs from model representation repo:
# * SIR model in Petri-net framework 
# * Lotka-Volterra model in the gene-regulatory network framework

models = {}
models["url"] = [
    "https://raw.githubusercontent.com/DARPA-ASKEM/Model-Representations/main/petrinet/examples/sir.json",
    "https://raw.githubusercontent.com/DARPA-ASKEM/Model-Representations/main/regnet/examples/lotka_volterra.json",
]

models["amr"] = [requests.get(url).json() for url in models["url"]]

# %%
# Convert model AMR to a MIRA MMT
# Switch depending on `header.schema_name``

models["mmt"] = []
for amr in models["amr"]:
    if amr["header"]["schema_name"] == "petrinet":
        mmt = askenet.petrinet.template_model_from_askenet_json(amr)
    elif amr["header"]["schema_name"] == "regnet":
        mmt = askenet.regnet.template_model_from_askenet_json(amr)
    else:
        mmt = None

    models["mmt"].append(mmt)

# %%
def plot_model(model: dict, model_type: str, ax: Optional[mpl.axes.Axes] = None) -> nx.MultiDiGraph:

    if ax == None:
        fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (8, 8))

    G = nx.MultiDiGraph()
    pos = nx.kamada_kawai_layout(G)

    if model_type == "MMT":

        G.graph = {k: v for k, v in model.annotations}
        G.graph["schema_name"] = None

        for i, template in enumerate(model.templates):

            node_sub = None
            node_out = None
            nodes_con = []
            if hasattr(template, "subject"):
                node_sub = template.subject
            if hasattr(template, "outcome"):   
                node_out = template.outcome
            if hasattr(template, "controller"):
                nodes_con.append(template.controller)
            if hasattr(template, "controllers"):
                nodes_con.extend(template.controllers)

            # Add interactor nodes
            G.add_nodes_from([
                (node.name, {"type": "interactor", "id": node.name})
                for node in ([node_sub, node_out] + nodes_con) if node != None
            ])
            
            # Add template node
            G.add_nodes_from([(i, {"type": "template", "id": f"{str(template.rate_law)}"})])

            # Add subject-template-outcome edges
            if node_sub != None:
                G.add_edges_from([(node_sub.name, i, {"type": "subject", "id": i})])
            if node_out != None:
                G.add_edges_from([(i, node_out.name, {"type": "outcome", "id": i})])

            # Add controller edges
            G.add_edges_from([
                (node.name, i, {"type": "controller", "id": i})
                for node in nodes_con
            ])    

    
        # Graph layout
        pos = nx.kamada_kawai_layout(G)

        # Draw
        h = []
        for node_type, shape in zip(("interactor", "template"), ("o", "s")):
            nodelist = [node for node, data in G.nodes(data = True) if data["type"] == node_type]
            k = nx.draw_networkx_nodes(G, pos, ax = ax, nodelist = nodelist, node_shape = shape, label = f"{node_type} nodes")
            h.append(k)

        for edge_type, style in zip(("subject", "outcome", "controller"), ("-", "-", "--")):
            edgelist = [(src, tgt) for (src, tgt, data) in G.edges(data = True) if data["type"] == edge_type]
            k = nx.draw_networkx_edges(G, pos, ax = ax, edgelist = edgelist, style = style, connectionstyle = "arc3,rad=0.2", label = f"{edge_type}")
            h.append(k)

        __ = ax.legend(
            [h[0], h[1], mpl.lines.Line2D([0, 1], [0, 1], color = 'k', linestyle = '-'), mpl.lines.Line2D([0, 1], [0, 1], color = 'k', linestyle = '--')], 
            ["Concept", "Template", "Direct", "Control"], 
            loc = "lower right"
        )

    if model_type == "AMR":

        G.graph = {k: v for k, v in model["header"].items()}

        if model["header"]["schema_name"] == "petrinet":

            # Add state and transition nodes
            G.add_nodes_from([
                (node["id"], node | {"type": node_type})
                for node_type in ("states", "transitions") for node in model["model"][node_type]
            ])

            # Map transition array index to rate expression
            # m_i_id = {i: transition["id"] for i, transition in enumerate(model["model"]["transitions"])}
            m_id_expr = {rate["target"]: rate["expression"] for rate in model["semantics"]["ode"]["rates"]}

            # Replace transition node ID with rate expression
            for node, data in G.nodes(data = True):
                if data["type"] == "transitions":
                    i = data["id"]
                    G.nodes(data = True)[node]["id"] = m_id_expr[i]


            # Add edges
            G.add_edges_from([
                (src, transition["id"], {"type": "input", "id": i})
                for i, transition in enumerate(model["model"]["transitions"]) for src in transition["input"]
            ])
            G.add_edges_from([
                (transition["id"], tgt, {"type": "output", "id": i})
                for i, transition in enumerate(model["model"]["transitions"]) for tgt in transition["output"]
            ])
            
            # Graph layout
            pos = nx.kamada_kawai_layout(G)
            
            # Draw
            h = []
            for node_type, shape in zip(("states", "transitions"), ("o", "s")):
                nodelist = [node for node, data in G.nodes(data = True) if data["type"] == node_type]
                k = nx.draw_networkx_nodes(G, pos, ax = ax, nodelist = nodelist, node_shape = shape, label = f"{node_type} nodes")
                h.append(k)

            for edge_type, style in zip(("input", "output"), ("-", "-")):
                edgelist = [(src, tgt) for (src, tgt, data) in G.edges(data = True) if data["type"] == edge_type]
                k = nx.draw_networkx_edges(G, pos, ax = ax, edgelist = edgelist, style = style, connectionstyle = "arc3,rad=0.2", label = f"{edge_type}")
                h.append(k)

            __ = ax.legend(
                [h[0], h[1], mpl.lines.Line2D([0, 1], [0, 1], color = 'k', linestyle = '-'), mpl.lines.Line2D([0, 1], [0, 1], color = 'k', linestyle = '--')], 
                ["State", "Transition", "All"], 
                loc = "lower right"
            )

        if model["header"]["schema_name"] == "regnet":

            # Add vertex nodes
            G.add_nodes_from([
                (node["id"], node | {"type": node_type})
                for node_type in ("vertices", ) for node in model["model"][node_type]
            ])

            # Add edges
            G.add_edges_from([
                (edge["source"], edge["target"], edge | {"type": "edge"})
                for edge in model["model"]["edges"]
            ])

            # Add edges implied by vertex sign
            for vertex in model["model"]["vertices"]:

                # Natural production
                if vertex["sign"] == True:
                    G.add_edges_from([
                        (vertex["id"], vertex["id"], {"type": "sign edge", "id": f'+{vertex["rate_constant"]}'})
                    ])

                if vertex["sign"] == False:
                    G.add_edges_from([
                        (vertex["id"], vertex["id"], {"type": "sign edge", "id": f'-{vertex["rate_constant"]}'})
                    ])


            # Graph layout
            pos = nx.kamada_kawai_layout(G)
            
            # Draw
            h = []
            for node_type, shape in zip(("vertices", ), ("o", )):
                nodelist = [node for node, data in G.nodes(data = True) if data["type"] == node_type]
                k = nx.draw_networkx_nodes(G, pos, ax = ax, nodelist = nodelist, node_shape = shape, label = f"{node_type} nodes")
                h.append(k)

            for edge_type, style in zip(("edge", "sign edge"), ("-", "-")):
                edgelist = [(src, tgt) for (src, tgt, data) in G.edges(data = True) if (data["type"] == edge_type) and (src != tgt)]
                k = nx.draw_networkx_edges(G, pos, ax = ax, edgelist = edgelist, style = style, connectionstyle = "arc3,rad=0.2", label = f"{edge_type}")

                edgelist = [(src, tgt) for (src, tgt, data) in G.edges(data = True) if (data["type"] == edge_type) and (src == tgt)]
                k = nx.draw_networkx_edges(G, pos, ax = ax, edgelist = edgelist, style = style, label = f"{edge_type}")

            __ = ax.legend(
                [h[0], mpl.lines.Line2D([0, 1], [0, 1], color = 'k', linestyle = '-')], 
                ["Vertex", "Edge"], 
                loc = "lower right"
            )

    # Node labels
    labels = {node: data["id"] for node, data in G.nodes(data = True)}
    nx.draw_networkx_labels(G, pos, ax = ax, labels = labels)

    # Edge labels
    edge_labels = {(src, tgt): data["id"] for src, tgt, data in G.edges(data = True)}
    nx.draw_networkx_edge_labels(G, pos, ax = ax, edge_labels = edge_labels, font_size = 8, label_pos = 0.3, rotate = False, bbox = {"boxstyle": "round", "fc": (0, 0, 0), "alpha": 0.0})

    __ = plt.setp(ax, title = f'{G.graph["name"]} ({G.graph["schema_name"] if not None else ""}) ({model_type})')
    
    return G

# %%
fig, axes = plt.subplots(2, 2, figsize = (10, 10))

G = plot_model(models["amr"][0], model_type = "AMR", ax = fig.axes[0])
G = plot_model(models["amr"][1], model_type = "AMR", ax = fig.axes[1])
G = plot_model(models["mmt"][0], model_type = "MMT", ax = fig.axes[2])
G = plot_model(models["mmt"][1], model_type = "MMT", ax = fig.axes[3])


# %%
MIRA_REST_URL = 'http://34.230.33.149:8771/api'

def viz_mmt(mmt, filename) -> NoReturn:

    res = requests.post(url = f'{MIRA_REST_URL}/viz/to_image', data = mmt.json())

    with open(filename, 'wb') as f:
        f.write(res.content)

    with Image.open(filename).convert('RGB') as im:
        im.show()

# %%
viz_mmt(models["mmt"][1], 'lotka_volterra_model.png')

# %%






