"""This module implements an API interface for retrieving Vensim models by Ventana Systems
denoted by the .mdl extension through a locally downloaded file or URL. We then
convert the Vensim model into a generic pysd model object that will be parsed and converted to an
equivalent MIRA template model. We preprocess the vensim file to extract variable expressions.

Vensim model documentation:https://www.vensim.com/documentation/sample_models.html

Repository of sample Vensim models: https://github.com/SDXorg/test-models/tree/master/samples
"""

import tempfile
import re

import pysd
from pysd.translators.vensim.vensim_file import VensimFile
import requests

from mira.metamodel import TemplateModel
from mira.sources.system_dynamics.pysd import (
    template_model_from_pysd_model,
)

__all__ = ["template_model_from_mdl_file", "template_model_from_mdl_url"]

NEW_CONTROL_DELIMETER = (
    " ******************************************************** .Control "
    "********************************************************"
)
UTF_ENCODING = "{UTF-8} "


def template_model_from_mdl_file(fname) -> TemplateModel:
    """Return a template model from a local Vensim file

    Parameters
    ----------
    fname : str or pathlib.Path
        The path to the local Vensim file

    Returns
    -------
    :
        A MIRA template model
    """
    pysd_model = pysd.read_vensim(fname)
    vensim_file = VensimFile(fname)
    expression_map = extract_vensim_variable_expressions(vensim_file.model_text)

    return template_model_from_pysd_model(pysd_model, expression_map)


def template_model_from_mdl_url(url) -> TemplateModel:
    """Return a template model from a Vensim file provided by an url

    Parameters
    ----------
    url : str
        The url to the mdl file

    Returns
    -------
    :
        A MIRA Template Model
    """
    data = requests.get(url).content
    temp_file = tempfile.NamedTemporaryFile(
        mode="w+b", suffix=".mdl", delete=False
    )

    with temp_file as file:
        file.write(data)

    pysd_model = pysd.read_vensim(temp_file.name)
    vensim_file = VensimFile(temp_file.name)
    expression_map = extract_vensim_variable_expressions(vensim_file.model_text)

    return template_model_from_pysd_model(pysd_model, expression_map)


def extract_vensim_variable_expressions(model_text):
    """Method that extracts expressions for each variable in a Vensim file

    Parameters
    ----------
    model_text : str
        The plain-text information about the Vensim file

    Returns
    -------
    : dict[str,str]
        Mapping of variable name to string variable expression
    """
    expression_map = {}
    model_split_text = model_text.split("|")

    for text in model_split_text:
        if NEW_CONTROL_DELIMETER in text:
            break
        if "=" not in text:
            continue

        # first entry usually has encoding type
        if UTF_ENCODING in text:
            text = text.replace(UTF_ENCODING, "")

        var_declaration = text.split("~")[0].split("=")
        old_var_name = var_declaration[0].strip()
        text_expression = var_declaration[1]
        if "INTEG" in text_expression:
            text_expression = re.search(r"\(([^,]+),", text_expression).group(1)
        expression_map[old_var_name] = text_expression

    return expression_map
