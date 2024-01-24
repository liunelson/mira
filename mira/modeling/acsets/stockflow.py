"""This module implements generation into stock and flow models which are
defined through a set of stocks, flows, links, and the input and output
connections between flows.
"""

__all__ = ["ACSetsStockFlowModel", "template_model_to_stockflow_ascet_json"]

from mira.modeling import Model
from mira.metamodel import *


class ACSetsStockFlowModel:
    """A class representing a stock and flow model."""

    def __init__(self, model: Model):
        """Instantiate a stock and flow model from a generic transition model.

        Parameters
        ----------
        model :
            The pre-compiled transition model.
        """
        self.properties = {}
        self.stocks = []
        self.flows = []
        self.links = []
        self.model_name = "Model"

        for idx, flow in enumerate(model.transitions.values()):
            fid = flow.template.name
            fname = flow.template.display_name

            input = flow.consumed[0].key
            output = flow.produced[0].key

            rate_law_str = (
                str(flow.template.rate_law) if flow.template.rate_law else None
            )

            flow_dict = {
                "_id": fid,
                "u": input,
                "d": output,
                "fname": fname,
                "ϕf": rate_law_str,
            }

            self.flows.append(flow_dict)

        vmap = {}
        for var_key, var in model.variables.items():
            vmap[var_key] = name = var.concept.name or str(var_key)
            display_name = var.concept.display_name or name

            stocks_dict = {
                "_id": name,
                "sname": display_name,
            }
            self.stocks.append(stocks_dict)

            # Declare 's' and 't' field of a link before assignment,
            # this is because if a stock is found to be
            # a target for a flow before it is a source, then the 't'
            # field of the link associated with a stock
            # will be displayed first
            links_dict = {"_id": name}
            links_dict["s"] = None
            links_dict["t"] = None
            for flow in model.transitions.values():
                if flow.consumed[0].concept.name == name:
                    links_dict["s"] = flow.template.name
                if flow.produced[0].concept.name == name:
                    links_dict["t"] = flow.template.name

            if not links_dict.get("s") and links_dict.get("t"):
                links_dict["s"] = links_dict.get("t")
            elif links_dict.get("s") and not links_dict.get("t"):
                links_dict["t"] = links_dict.get("s")

            self.links.append(links_dict)

    def to_json(self):
        """
        Return a JSON dict structure of the Petri net model.

        Returns
        -------
        : JSON
            The JSON dict structure of the stock and flow model.
        """
        return {"Flow": self.flows, "Stock": self.stocks, "Link": self.links}


def template_model_to_stockflow_ascet_json(tm: TemplateModel):
    """
    Convert a TemplateModel into stock flow JSON and return the converted
    model.

    Parameters
    ----------
    tm :
        The TemplateModel to be converted.

    Returns
    -------
    : JSON
        The JSON structure representing the stock flow model.
    """
    return ACSetsStockFlowModel(Model(tm)).to_json()
