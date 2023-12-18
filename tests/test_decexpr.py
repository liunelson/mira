import sympy

from mira.examples.decapodes.decapodes_examples import (
    get_oscillator_decaexpr,
    get_friction_decaexpr,
)
from mira.metamodel.decapodes import RootVariable


def test_oscillator_decaexpr():
    # Check that variable types are correct (Form0,  Form1, Literal, infer,
    # Constant, Dualform_1, Dualform_2, etc.)
    # Check that we have the correct number of op1, op2, tangent variables,
    # and summations
    # Check that the op1, op2, tangent variables, and summations are the
    # expected ones (check source, target etc.)
    # Check that expressions are correct...

    oscillator_decapode = get_oscillator_decaexpr()
    # Original equation:
    # ∂ₜ(X) = V
    # ∂ₜ(∂ₜ(X)) = ∂ₜ(V) = -1*k*X

    variable_name_set = {
        "V",
        "X",
        "k",
        "-1",
        "∂ₜ(V)",
        "mult_1",
    }

    assert variable_name_set == {
        v.name for v in oscillator_decapode.variables.values()
    }
    name_to_variable = {
        v.name: v for v in oscillator_decapode.variables.values()
    }
    assert len(oscillator_decapode.variables) == 6
    assert name_to_variable["X"].type == "Form0"
    assert name_to_variable["V"].type == "Form0"
    assert name_to_variable["k"].type == "Constant"
    assert name_to_variable["-1"].type == "Literal"
    assert name_to_variable["∂ₜ(V)"].type == "infer"
    assert name_to_variable["mult_1"].type == "infer"

    # Check there is at least one RootVariable
    assert any(
        isinstance(v, RootVariable)
        for v in oscillator_decapode.variables.values()
    )

    assert len(oscillator_decapode.op1s) == 2  # Have dX/dt and dV/dt
    unary_targets = {op.tgt.name for op in oscillator_decapode.op1s.values()}
    assert unary_targets == {"V", "∂ₜ(V)"}
    unary_args = {op.src.name for op in oscillator_decapode.op1s.values()}
    assert unary_args == {"X", "V"}

    assert (
        len(oscillator_decapode.op2s) == 2
    )  # Have -1*k=mult_1, mult_1*X=mult_2
    for op2 in oscillator_decapode.op2s.values():
        assert {op2.proj1.name, op2.proj2.name} == {"-1", "k"} or {
            op2.proj1.name,
            op2.proj2.name,
        } == {"mult_1", "X"}
    assert {op2.res.name for op2 in oscillator_decapode.op2s.values()} == {
        "mult_1",
        "∂ₜ(V)",
    }

    assert len(oscillator_decapode.summations) == 0  # No summations

    assert (
        len(oscillator_decapode.tangent_variables) == 2
    )  # Have dX/dt and dV/dt
    tangent_variable_names = {
        v.incl_var.name for v in oscillator_decapode.tangent_variables.values()
    }
    assert tangent_variable_names == {"V", "∂ₜ(V)"}

    # Check that expressions are correct
    dt = sympy.Function("∂ₜ")
    X, V, k, minus_one = sympy.symbols("X V k -1")
    assert name_to_variable["k"].expression == k
    assert name_to_variable["-1"].expression == minus_one
    assert name_to_variable["X"].expression == X
    assert name_to_variable["V"].expression == dt(X)
    assert isinstance(name_to_variable["∂ₜ(V)"], RootVariable)
    assert {
        name_to_variable["∂ₜ(V)"].expression[0],
        name_to_variable["∂ₜ(V)"].expression[1],
    } == {dt(dt(X)), minus_one * k * X}
    assert name_to_variable["mult_1"].expression == minus_one * k


def test_friction_decaexpr():
    friction_decapode = get_friction_decaexpr()
    # Original equation:
    # ∂ₜ(Q) = κ*V + λ*(Q-Q₀)

    variable_set = {
        "V",
        "Q",
        "κ",
        "λ",
        "Q₀",
        "∂ₜ(Q)",
        "mult_1",
        "sub_1",
        "mult_2",
    }
    assert variable_set == {
        v.name for v in friction_decapode.variables.values()
    }
    name_to_variable = {v.name: v for v in friction_decapode.variables.values()}
    assert name_to_variable["V"].type == "Form0"
    assert name_to_variable["Q"].type == "Form0"
    assert name_to_variable["κ"].type == "Constant"
    assert name_to_variable["λ"].type == "Constant"
    assert name_to_variable["Q₀"].type == "Parameter"
    assert name_to_variable["∂ₜ(Q)"].type == "infer"
    assert name_to_variable["mult_1"].type == "infer"
    assert name_to_variable["sub_1"].type == "infer"
    assert name_to_variable["mult_2"].type == "infer"

    # Check there is at least one RootVariable
    assert any(
        isinstance(v, RootVariable)
        for v in friction_decapode.variables.values()
    )

    assert len(friction_decapode.op1s) == 1  # Only have dQ/dt
    assert friction_decapode.op1s[0].tgt.name == "∂ₜ(Q)"
    assert friction_decapode.op1s[0].src.name == "Q"
    assert friction_decapode.op1s[0].function_str == "∂ₜ"

    # κ*V=mult_1, Q-Q₀=sub_1, and λ*sub_1=mult_2
    assert len(friction_decapode.op2s) == 3
    for op2 in friction_decapode.op2s.values():
        assert (
            {op2.proj1.name, op2.proj2.name} == {"κ", "V"}
            or {op2.proj1.name, op2.proj2.name} == {"Q", "Q₀"}
            or {op2.proj1.name, op2.proj2.name} == {"λ", "sub_1"}
        )
    assert {op2.res.name for op2 in friction_decapode.op2s.values()} == {
        "mult_1",
        "sub_1",
        "mult_2",
    }

    assert (
        len(friction_decapode.summations) == 1
    )  # Only have mult_1+mult_2=sum_1
    assert friction_decapode.summations[0].sum.name == "∂ₜ(Q)"
    assert {v.name for v in friction_decapode.summations[0].summands} == {
        "mult_1",
        "mult_2",
    }

    assert len(friction_decapode.tangent_variables) == 1  # Only have dQ/dt
    assert friction_decapode.tangent_variables[0].incl_var.name == "∂ₜ(Q)"

    # Check that expressions are correct
    dt = sympy.Function("∂ₜ")
    Q, V, kappa, _lambda, Q_0 = sympy.symbols("Q V κ λ Q₀")
    assert name_to_variable["κ"].expression == kappa
    assert name_to_variable["λ"].expression == _lambda
    assert name_to_variable["Q"].expression == Q
    assert name_to_variable["V"].expression == V
    assert name_to_variable["Q₀"].expression == Q_0
    assert isinstance(name_to_variable["∂ₜ(Q)"], RootVariable)
    root_variable = name_to_variable["∂ₜ(Q)"]
    assert {root_variable.expression[0], root_variable.expression[1]} == {
        dt(Q), kappa * V + _lambda * (Q - Q_0)
    }
    assert name_to_variable["mult_1"].expression == kappa * V
    assert name_to_variable["sub_1"].expression == Q - Q_0
    assert (
        name_to_variable["mult_2"].expression
        == _lambda * name_to_variable["sub_1"].expression
    )
