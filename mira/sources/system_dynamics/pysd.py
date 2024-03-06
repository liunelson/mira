"""This module implements parsing of a generic pysd model irrespective of source and source type
and extracting its contents to create an equivalent MIRA template model.
"""

__all__ = ["template_model_from_pysd_model"]

import copy
import re
import typing as t
import logging
from more_itertools import chunked

import pandas as pd
import sympy

from mira.metamodel import *
from mira.metamodel.utils import safe_parse_expr
from mira.metamodel import Concept, TemplateModel
from mira.sources.util import (
    parameter_to_mira,
    transition_to_templates,
    get_sympy,
)

logger = logging.getLogger(__name__)

CONTROL_VARIABLE_NAMES = {
    "FINALTIME",
    "INITIALTIME",
    "SAVEPER",
    "TIMESTEP",
    "FINAL TIME",
    "INITIAL TIME",
    "TIME STEP",
}
UNITS_MAPPING = {
    sympy.Symbol("Person"): sympy.Symbol("person"),
    sympy.Symbol("Persons"): sympy.Symbol("person"),
    sympy.Symbol("Day"): sympy.Symbol("day"),
    sympy.Symbol("Days"): sympy.Symbol("day"),
}
SYMPY_FLOW_RATE_PLACEHOLDER = safe_parse_expr("xxplaceholderxx")


def template_model_from_pysd_model(
    pysd_model,
    expression_map,
    *,
    grounding_map=None,
    initials_map: t.Optional[t.Dict[str, float]] = None,
) -> TemplateModel:
    """Given a model and its accompanying expression_map, extract information from the arguments
    to create an equivalent MIRA template model.

    Parameters
    ----------
    pysd_model : Model
        The pysd model object
    expression_map : dict[str,str]
        Map of variable name to expression
    grounding_map: dict[str, Concept]
        A grounding map, a map from label to Concept
    initials_map: dict[str, float]
        A pre-populated mapping of initial values

    Returns
    -------
    :
        MIRA template model
    """
    # dataframe that contains information about each variable in the model
    model_doc_df = pysd_model.doc

    # mapping of expressions after they have been processed to be sympy compatible
    processed_expression_map = {}

    # Mapping of variable name in vensim model to variable python-equivalent name
    old_name_new_pyname_map = dict(
        zip(model_doc_df["Real Name"], model_doc_df["Py Name"])
    )

    # preprocess expression text to make it sympy parseable
    for var_name, var_expression in expression_map.items():
        new_var_name = old_name_new_pyname_map[var_name]
        processed_expression_map[new_var_name] = preprocess_expression_text(
            var_expression
        )

    # Mapping of variable's python name to symbol for expression parsing in sympy
    symbols = dict(
        zip(
            model_doc_df["Py Name"],
            list(map(lambda x: sympy.Symbol(x), model_doc_df["Py Name"])),
        )
    )

    # Retrieve the states
    model_states = model_doc_df.loc[
        (model_doc_df["Type"] == "Stateful")
        & (model_doc_df["Subtype"] == "Integ")
    ]

    concepts = {}
    all_states = set()

    # map between states and incoming/out-coming rate laws
    state_rate_map = {}

    # process states and build mapping of state to input rate laws and output rate laws
    for index, state in model_states.iterrows():
        concept_state = state_to_concept(state, grounding_map=grounding_map)
        concepts[concept_state.name] = concept_state
        all_states.add(concept_state.name)
        symbols[concept_state.name] = sympy.Symbol(concept_state.name)

        state_name = state["Py Name"]
        state_rate_map[state_name] = {"input_rates": [], "output_rates": []}
        state_expr_text = processed_expression_map[state_name]

        # retrieve the expression of inflows and outflows for the state
        state_arg_sympy = safe_parse_expr(state_expr_text, symbols)

        # Create a map of states and whether the flows involved with a state are going in
        # or out of the state
        if state_arg_sympy.args:
            # if it's just the negation of a single symbol
            if (
                sympy.core.numbers.NegativeOne() in state_arg_sympy.args
                and len(state_arg_sympy.args) == 2
            ):
                str_symbol = str(state_arg_sympy)
                state_rate_map[state_name]["output_rates"].append(
                    str_symbol[1:]
                )
            else:
                for rate_free_symbol in state_arg_sympy.args:
                    str_rate_free_symbol = str(rate_free_symbol)
                    if "-" in str_rate_free_symbol:
                        # If the symbol representing the flow has a negative sign, it is an
                        # outgoing flow
                        # Add the symbol to outputs symbol without the negative sign
                        state_rate_map[state_name]["output_rates"].append(
                            str_rate_free_symbol[1:]
                        )
                    else:
                        # else it is an incoming flow
                        state_rate_map[state_name]["input_rates"].append(
                            str_rate_free_symbol
                        )
        else:
            # if it's just a single symbol (i.e. no negation), args property will be empty
            state_rate_map[state_name]["input_rates"].append(
                str(state_arg_sympy)
            )

    # process initials, currently we use the value of the state at timestamp 0
    # state initial values are listed in the same order as states are for the pysd model

    # TODO: Current assumption is true for other models but not the Vensim hackathon model.
    # We have 19 states but 44 initial values so we do not assign the right initial values.
    # For the hackathon, the only initial that has a value other than 0 is "Susceptibles" with a
    # value of "1.3392e+09.
    # One solution for this through Vensim is to text parse the Integ function which is in the form
    # INTEG(<expression>,<value>) and extract the value as well and keep it in a separate mapping
    # for initials where we keep track of the value for each initial. Then pass the initial mapping
    # to "template_model_from_pysd_model" method. So this would be done in "vensim.py"
    mira_initials = {}
    for state_initial_value, (state_name, state_concept) in zip(
        pysd_model.state, concepts.items()
    ):
        if initials_map and (mapped_value := initials_map.get(state_name)):
            initial = Initial(
                concept=concepts[state_name].copy(deep=True),
                expression=SympyExprStr(sympy.Float(mapped_value)),
            )
        # if the state value is not a number
        elif not isinstance(state_initial_value, int) and not isinstance(
            state_initial_value, float
        ):
            logger.warning(f"got non-numeric state value for {state_name}: {state_initial_value}")
            initial = Initial(
                concept=concepts[state_name].copy(deep=True),
                expression=SympyExprStr(sympy.Float("0")),
            )
        else:
            initial = Initial(
                concept=concepts[state_name].copy(deep=True),
                expression=SympyExprStr(sympy.Float(state_initial_value)),
            )
        mira_initials[initial.concept.name] = initial

    # Stella issue
    # Currently cannot capture all parameters as some of them cannot have their expressions parsed
    # Additionally, some auxiliary variables are added as parameters due to inability to parse
    # the auxiliary expression and differentiate it from a parameter
    mira_parameters = {}
    # process parameters
    # No discernible way to identify only parameters in model_doc_df, so go through all variables
    # in the processed expression map

    for name, expression in processed_expression_map.items():
        # evaluate the expression
        # try catch exists for stella models, safe_parse_expr will error for incorrectly
        # constructed parameters in Stella
        # eval_expression returns a sympy object
        try:
            eval_expression = safe_parse_expr(expression).evalf()
        except TypeError:
            eval_expression = safe_parse_expr("0")

        # convert the sympy object to a string
        str_eval_expression = str(eval_expression)
        value = None
        is_initial = False

        # Stella issue
        # Sometimes parameter values reference a stock rather than a number
        if str_eval_expression in mira_initials:
            value = float(str(mira_initials[str_eval_expression].expression))
            is_initial = True

        # Stella issue
        # Replace negative signs for placeholder parameter value which is "-1" for Aux structures
        # that cannot be parsed from stella models. isdecimal() interprets the "-" as a dash.

        # If the sympy expression object isn't equal to the placeholder
        # If the string is evaluated to be a number after removing decimals, dashes and empty spaces
        # then create a parameter
        if str_eval_expression in mira_initials or (
            eval_expression != SYMPY_FLOW_RATE_PLACEHOLDER
            and str_eval_expression.replace(".", "")
            .replace(" ", "")
            .replace("-", "")
            .isdecimal()
        ):
            if not is_initial:
                value = float(str_eval_expression)
            model_parameter_info = model_doc_df[model_doc_df["Py Name"] == name]

            # if units exist
            if (
                model_parameter_info["Units"].values[0]
                and model_parameter_info["Units"].values[0] != "dimensionless"
            ):
                unit_text = (
                    model_parameter_info["Units"].values[0].replace(" ", "")
                )
                parameter = {
                    "id": name,
                    "value": value,
                    "description": model_parameter_info["Comment"].values[0],
                    "units": {"expression": unit_text},
                }
            else:
                parameter = {
                    "id": name,
                    "value": value,
                    "description": model_parameter_info["Comment"].values[0],
                }

            mira_parameters[name] = parameter_to_mira(parameter)

            # standardize parameter units if they exist
            if mira_parameters[name].units:
                param_unit = mira_parameters[name].units
                for old_unit_symbol, new_unit_symbol in UNITS_MAPPING.items():
                    param_unit.expression = param_unit.expression.subs(
                        old_unit_symbol, new_unit_symbol
                    )
    # construct transitions mapping that determine inputs and outputs states to a rate-law
    transition_map = {}

    # List of rates
    # Currently duplicate rates are added
    # Using a set breaks the current iteration of tests for system dynamics in
    # "tests/system_dynamics.py" as the tests there test for hard-coded values in the list of
    # templates

    # For example, we test that the first template in the list of templates associated with a
    # template model is of type ControlledConversion for the SIR model; however, using sets,
    # the first type of template is of type NaturalConversion. Would require a lot of rewriting
    # of tests.
    rates = []
    for state_rates in state_rate_map.values():
        for input_rates in state_rates["input_rates"]:
            rates.append(input_rates)
        for output_rates in state_rates["output_rates"]:
            rates.append(output_rates)

    # create map of transitions
    for rate_name in rates:
        # because we use a list of rates to not change the template order for testing purposes
        # skip adding previously parsed rates to the transition_map
        if rate_name in transition_map:
            continue
        inputs, outputs, controllers, rate_expr = [], [], [], None
        for state_name, in_out_rate_map in state_rate_map.items():
            # if a rate cannot be parsed, assign None
            # This is a Stella issue as some expressions aren't created properly for rates
            try:
                rate_expr = safe_parse_expr(
                    processed_expression_map[rate_name],
                    symbols,
                )
            except TypeError:
                rate_expr = SYMPY_FLOW_RATE_PLACEHOLDER

            # if a rate is leaving a state, then that state is an input to the rate
            if rate_name in in_out_rate_map["output_rates"]:
                inputs.append(state_name)

            # if a rate is going into a state, then that state is an output to the rate
            if rate_name in in_out_rate_map["input_rates"]:
                outputs.append(state_name)

                # if a state is present in a rate law, and the state isn't an input to the rate
                # law, then that state is a controller of the rate law
                if (
                    sympy.Symbol(state_name) in rate_expr.free_symbols
                    and rate_name
                    not in state_rate_map[state_name]["output_rates"]
                ):
                    controllers.append(state_name)

        transition_map[rate_name] = {
            "name": rate_name,
            "expression": rate_expr,
            "inputs": inputs,
            "outputs": outputs,
            "controllers": controllers,
        }

    used_states = set()

    # Create templates from transitions
    templates_ = []
    for template_id, (transition_name, transition) in enumerate(
        transition_map.items()
    ):
        input_concepts, input_names = [], []
        output_concepts, output_names = [], []
        controller_concepts = []

        for input_name in transition.get("inputs"):
            input_concepts.append(concepts[input_name].copy(deep=True))
            input_names.append(input_name)

        for output_name in transition.get("outputs"):
            output_concepts.append(concepts[output_name].copy(deep=True))
            output_names.append(output_name)

        for controller_name in transition.get("controllers"):
            controller_concepts.append(
                concepts[controller_name].copy(deep=True)
            )

        used_states.update(input_names, output_names)

        templates_.extend(
            transition_to_templates(
                input_concepts=input_concepts,
                output_concepts=output_concepts,
                controller_concepts=controller_concepts,
                transition_rate=transition["expression"]
                if transition["expression"] != SYMPY_FLOW_RATE_PLACEHOLDER
                else None,
                transition_id=str(template_id + 1),
                transition_name=transition_name,
            )
        )

    return TemplateModel(
        templates=templates_,
        parameters=mira_parameters,
        initials=mira_initials,
    )


def state_to_concept(state, grounding_map=None) -> Concept:
    """Create a MIRA Concept from a state

    Parameters
    ----------
    state : pd.Series
        The series that contains state data
    grounding_map: dict[str, Concept]
        A grounding map, a map from label to Concept

    Returns
    -------
    :
        The MIRA concept created from the state
    """
    name = state["Py Name"]
    description = state["Comment"]
    unit_dict = {
        "expression": state["Units"].replace(" ", "")
        if state["Units"]
        else None
    }
    unit_expr = get_sympy(unit_dict, UNIT_SYMBOLS)
    units_obj = Unit(expression=unit_expr) if unit_expr else None

    # If there's something hacked in, use that directly,
    # keeping the units and description from the model
    if grounding_map is not None and name in grounding_map:
        concept = copy.deepcopy(grounding_map[name])
        concept.name = name
        concept.description = description
        concept.units = units_obj
        return concept
    else:
        print("could not ground", name)

    return Concept(name=name, units=units_obj, description=description)


def preprocess_expression_text(expr_text):
    """Preprocess a string expression to convert the expression into sympy parseable string

    Parameters
    ----------
    expr_text : str
        The string expression

    Returns
    -------
    : str
        The processed string expression
    """

    if not expr_text:
        return expr_text

    # This regex removes spaces between operands and operators. Also removes spaces between
    # parenthesis and operators or operands. Works for symbols as well.
    expr_text = re.sub(
        r"(?<=[^\w\s])\s+(?=[^\w\s])|(?<=[^\w\s])\s+(?=\w)|(?<=\w)\s+(?=[^\w\s])",
        "",
        expr_text,
    )

    # TODO: Use regular expressions for all text preprocessing rather than using string replace
    # strip leading and trailing white spaces
    # replace space between two words that makeup a variable name with "_"'
    # replace single and doubel quotation marks
    # replace ampersand & with and
    expr_text = (
        expr_text.strip()
        .replace("^", "**")
        .replace(" ", "_")
        .replace("'", "")
        .replace('"', "")
        .replace("&", "_")
        .lower()
    )
    return expr_text


def ifthenelse_to_piecewise(expr_text):
    """Convert Vensim if then else expression to sympy Piecewise string

    Parameters
    ----------
    expr_text : str
        The string expression

    Returns
    -------
    : str
        The sympy Piecewise expression as a string
    """
    # We find the three components and rearrange them for piecewise
    if_then_else_regex = r"IF THEN ELSE\((.*),(.*)\,(.*)\)"
    match = re.search(if_then_else_regex, expr_text)
    if match:
        condition, then, else_ = match.groups()
        piecewise = f"Piecewise(({then}, {condition}), ({else_}, True))"
        # We also need to replace single = with == but make sure if there is a ==, we don't replace
        # it and also, we need to make sure we don't replace <= or >=
        piecewise = re.sub(r"(?<!<|>)=(?!=)", "==", piecewise)
        return piecewise


def with_lookup_to_piecewise(expr_text: str) -> str:
    """Convert a Vensim WITH LOOKUP expression to a piecewise function.

    The semantics of the ``WITH LOOKUP`` element are documented at
    https://www.vensim.com/documentation/fn_with_lookup.html.

    For example, this could come from:

    .. code-block::

        Infection Rate new arrivals= WITH LOOKUP (
            Time,
                (
                [(0,0)-(500,100)],(0,0),(1,2),(2,1),(3,0),(4,2),(5,1),(6,2),(7,3),(8,6),(9,2),(10,7\
                ),(11,10),(12,4),(13,10),(14,5),(15,
                11),(16,14),(17,14),(18,26),(19,34),(20,35),(21,45),(22,55),(23,38),(24,34),(25,24)\
                ,(26,40),(27,16),(28,20),(29,12),(30,
                23),(31,14),(32,8),(33,14),(34,12),(35,5),(36,9),(37,6),(38,0),(39,0),(40,0),(1000,\
                0)
                ))
            ~	Persons/Day
            ~		|

    Which ends up being the text (after normalization from vensim)

    .. code-block::

        with_lookup(time,([(0,0)-(500,100)],(0,0),(1,2),(2,1),(3,0),(4,2),(5,1),(6,2),\
        (7,3),(8,6),(9,2),(10,7),(11,10),(12,4),(13,10),(1000,0)))
    """
    # there's a variety of ways this is written WITH LOOKUP, with_lookup, withlookup
    # so just normalize it all out
    expr_text = expr_text.strip().replace(" ", "").replace("_", "").lower()
    if not expr_text.lower().startswith("withlookup"):
        raise ValueError(expr_text)
    expr_text = expr_text[len("withlookup") :].lstrip("(").rstrip(")")
    # The first input is either a value or an input variable name
    # The second input is a list of X,Y pairs
    variable, second = (x.strip() for x in expr_text.split(",", 1))
    second: str = second.lstrip("(").rstrip(")")

    # This is an undocumented part of this function, but it appears that the
    # list begins with some kind of definition of the form [(a,b)-(c,d)]
    if second.startswith("["):
        box, rest = second.lstrip("[").split("]", 1)
        x1, y1, x2, y2 = [
            float(x.strip())
            for x in box.replace("(", "")
            .replace(")", "")
            .replace("-", ",")
            .split(",")
        ]
    else:
        rest = second
        x1, y1, x2, y2 = [None] * 4

    # TODO how to use x1, y1, x2, and y2? It's not clear if/how these are used
    # to fill in gaps in the lookup table

    # This gets the list of x,y pairs in order from the lookup table
    pairs = list(
        chunked(
            (
                float(x.strip())
                for x in rest.strip()
                .replace("(", "")
                .replace(")", "")
                .split(",")
                if x.strip()
            ),
            2,
        )
    )
    # print(variable, x1, y1, x2, y2, pairs)

    # construct the sympy string as a simple lookup, where the y
    # value gets returned if the target variable is the given x value
    # for each x,y pair
    # FIXME what's the right way to write the conditional here
    conditions = ",".join(f"({y}, {variable} >= {x})" for x, y in pairs)
    sympy_str = f"Piecewise({conditions})"
    return sympy_str
