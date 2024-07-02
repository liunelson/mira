import requests
import json

import sympy

from mira.metamodel import *
from mira.sources.amr import model_from_url
from mira.sources.amr import petrinet
from mira.sources.amr import regnet
from mira.sources.amr import stockflow
from mira.modeling.amr.regnet import template_model_to_regnet_json
from mira.modeling.amr.petrinet import template_model_to_petrinet_json
from mira.modeling.amr.stockflow import template_model_to_stockflow_json

petrinet_example = 'https://raw.githubusercontent.com/DARPA-ASKEM/' \
                   'Model-Representations/main/petrinet/examples/sir.json'
regnet_example = 'https://raw.githubusercontent.com/DARPA-ASKEM/' \
                 'Model-Representations/main/regnet/examples/lotka_volterra.json'
stockflow_example = 'https://raw.githubusercontent.com/DARPA-ASKEM/' \
                    'Model-Representations/7f5e377225675259baa6486c64102f559edfd79f/stockflow/examples/sir.json'

template_model = TemplateModel(
    templates=[
        ControlledReplication(
            name='replication',
            controller=Concept(name='A'),
            subject=Concept(name='B'),
            rate_law=SympyExprStr(sympy.sympify('k * A * B / (1 + B)'))
        ),
        NaturalDegradation(
            name='degradation',
            subject=Concept(name='B'),
            rate_law=SympyExprStr(sympy.sympify('k * B'))
        )
    ],
    observables={
        'obs1': Observable(
            name='obs1',
            expression=SympyExprStr(sympy.sympify('A + B')),
            display_name='obs1'
        )
    },
    annotations=Annotations(
        name="test_name",
        description="test_description",
        diseases=["test_disease"],
        hosts=["test_host"],
        license="test_license",
        locations=["test_location"],
        model_types=["test_model_type"],
        pathogens=["test_pathogen"],
        references=["test_reference"],
        authors=[Author(name="test_author")],
        time_start="2020-03-01T00:00:00",
        time_end="2020-08-01T00:00:00",
        time_scale="days"
    ),
    time=Time(name='timexx'),
)


def stockflow_set_up_file():
    return requests.get(stockflow_example).json()


def test_petrinet_model_from_url():
    template_model = petrinet.model_from_url(petrinet_example)
    assert len(template_model.templates) == 2
    assert isinstance(template_model.templates[0], ControlledConversion)
    assert isinstance(template_model.templates[1], NaturalConversion)
    assert template_model.templates[0].controller.display_name == 'Infected'
    assert template_model.templates[0].controller.name == 'I'
    assert template_model.templates[0].subject.display_name == 'Susceptible'
    assert template_model.templates[0].outcome.display_name == 'Infected'
    assert template_model.templates[1].subject.display_name == 'Infected'
    assert template_model.templates[1].outcome.display_name == 'Recovered'


def test_regnet_model_from_url():
    template_model = regnet.model_from_url(regnet_example)
    assert len(template_model.templates) == 4
    assert isinstance(template_model.templates[0], NaturalReplication)
    assert isinstance(template_model.templates[1], NaturalDegradation)
    assert isinstance(template_model.templates[2], ControlledDegradation)
    assert isinstance(template_model.templates[3], ControlledReplication)


def test_model_from_url_generic():
    tm = model_from_url(regnet_example)
    assert len(tm.templates) == 4

    tm = model_from_url(petrinet_example)
    assert len(tm.templates) == 2


def test_stockflow_flow_to_template():
    tm = stockflow.model_from_url(stockflow_example)
    assert len(tm.templates) == 2
    assert isinstance(tm.templates[0], ControlledConversion)
    assert isinstance(tm.templates[1], NaturalConversion)
    assert tm.templates[0].name == 'flow1'
    assert tm.templates[1].name == 'flow2'
    assert tm.templates[0].subject.name == 'S'
    assert tm.templates[0].outcome.name == 'I'
    assert tm.templates[0].controller.name == 'I'
    assert tm.templates[1].subject.name == 'I'
    assert tm.templates[1].outcome.name == 'R'

    assert sympy.Symbol('S') in tm.templates[0].rate_law.free_symbols
    assert sympy.Symbol('p_cbeta') in tm.templates[0].rate_law.free_symbols
    assert sympy.Symbol('p_N') in tm.templates[0].rate_law.free_symbols
    assert sympy.Symbol('I') in tm.templates[0].rate_law.free_symbols
    assert sympy.Symbol('I') in tm.templates[1].rate_law.free_symbols
    assert sympy.Symbol('p_tr') in tm.templates[1].rate_law.free_symbols

    assert len(tm.templates[0].rate_law.free_symbols) == 4
    assert len(tm.templates[1].rate_law.free_symbols) == 2


def test_stockflow_parameter_to_mira():
    sfamr = stockflow_set_up_file()
    tm = stockflow.model_from_url(stockflow_example)
    for amr_param, tm_param in zip(sfamr['semantics']['ode']['parameters'],
                                   tm.parameters.values()):
        assert amr_param['id'] == tm_param.name
        assert amr_param['name'] == tm_param.display_name
        assert amr_param['description'] == tm_param.description
        assert amr_param['value'] == tm_param.value
        assert amr_param.get('units', {}) == (
            tm_param.units if tm_param.units else {})


def test_stockflow_initial_to_mira():
    sfamr = stockflow_set_up_file()
    tm = stockflow.model_from_url(stockflow_example)
    for amr_initial, tm_initial in zip(sfamr['semantics']['ode']['initials'],
                                       tm.initials.values()):
        assert amr_initial['target'] == tm_initial.concept.name
        assert amr_initial['expression'] == str(tm_initial.expression)
        assert amr_initial['expression_mathml'] == expression_to_mathml(
            tm_initial.expression)


def test_stockflow_stock_to_concept():
    stock = {
        "id": "S",
        "name": "Susceptible",
        "grounding": {
            "identifiers": {
                "ido": "0000514"
            },
            "modifiers": {
                "test_key": "test_value"
            }
        }
    }

    concept = stockflow.stock_to_concept(stock)
    assert concept.name == str(stock['id'])
    assert concept.display_name == stock['name']
    assert concept.context == stock['grounding']['modifiers']
    assert concept.identifiers == stock['grounding']['identifiers']


def test_regnet_rate_laws():
    # Make a simple template model with rate laws, then export
    # into AMR, ingest, then make sure we get proper rate laws back
    # out.
    amr_json = template_model_to_regnet_json(template_model)
    tm = regnet.template_model_from_amr_json(amr_json)
    assert isinstance(tm.templates[0].rate_law, SympyExprStr)
    assert isinstance(tm.templates[1].rate_law, SympyExprStr)
    assert tm.templates[1].rate_law.args[0].equals(
        sympy.sympify('k * A * B / (1 + B)'))
    assert tm.time.name == 'timexx'
    assert isinstance(tm.observables['obs1'].expression, SympyExprStr)


def test_annotation_serialization_ingestion():
    """Test to see if we can serialize template models that contain
    datetime in their annotations into different frameworks.

    Also test to see if we can extract all annotation related attributes
    when we ingest an amr
    """
    amrs = [template_model_to_regnet_json(template_model),
            template_model_to_petrinet_json(template_model),
            template_model_to_stockflow_json(template_model)
            ]
    for amr in amrs:
        json.dumps(amr)

    # test to see if we can extract all annotation attributes during ingestion
    # for each framework
    regnet_tm = regnet.template_model_from_amr_json(amrs[0])
    petrinet_tm = petrinet.template_model_from_amr_json(amrs[1])
    stockflow_tm = stockflow.template_model_from_amr_json(amrs[2])

    zipped_annotations = zip(regnet_tm.annotations.dict().values(),
                              petrinet_tm.annotations.dict().values(),
                              stockflow_tm.annotations.dict().values())

    for annotation_attribute_tuple in zipped_annotations:
        assert annotation_attribute_tuple[0]
        assert annotation_attribute_tuple[1]
        assert annotation_attribute_tuple[2]
