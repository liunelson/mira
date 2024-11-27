from pyobo.api.utils import get_version
from pyobo.getters import _ensure_ontology_path
from pyobo.utils.path import prefix_directory_join
from obonet import read_obo
import networkx
import pickle


def download_convert_ncbitaxon_obo_to_graph():
    resource_prefix = "ncbitaxon"
    version = get_version(resource_prefix)

    # Checks to see if the pickled ncbitaxon obo graph exists in the container
    cached_relabeled_obo_graph_path = prefix_directory_join(resource_prefix,
                                                            name="relabeled_obo_graph.pkl",
                                                            version=version)
    if not cached_relabeled_obo_graph_path.exists():
        _, obo_path = _ensure_ontology_path(
            resource_prefix, force=False, version=version
        )
        obo_graph = read_obo(obo_path)

        # Normalize node indices
        relabeled_graph = networkx.relabel_nodes(
            obo_graph, lambda node_index: node_index.lower()
        )
        with open(
            cached_relabeled_obo_graph_path, "wb"
        ) as relabeled_graph_file:
            pickle.dump(relabeled_graph, relabeled_graph_file)


if __name__ == "__main__":
    download_convert_ncbitaxon_obo_to_graph()
