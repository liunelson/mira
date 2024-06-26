"""
Python model 'IndiaNonSubscriptedPulsed.py'
Translated using PySD
"""

from pathlib import Path
import numpy as np

from pysd.py_backend.functions import pulse, xidz, if_then_else, ramp
from pysd.py_backend.statefuls import Integ, Delay
from pysd import Component

__pysd_version__ = "3.13.2"

__data = {"scope": None, "time": lambda: 0}

_root = Path(__file__).parent


component = Component()

#######################################################################
#                          CONTROL VARIABLES                          #
#######################################################################

_control_vars = {
    "initial_time": lambda: 0,
    "final_time": lambda: 500,
    "time_step": lambda: 0.0625,
    "saveper": lambda: 1,
}


def _init_outer_references(data):
    for key in data:
        __data[key] = data[key]


@component.add(name="Time")
def time():
    """
    Current time of the model.
    """
    return __data["time"]()


@component.add(
    name="FINAL TIME", units="Day", comp_type="Constant", comp_subtype="Normal"
)
def final_time():
    """
    The final time for the simulation.
    """
    return __data["time"].final_time()


@component.add(
    name="INITIAL TIME", units="Day", comp_type="Constant", comp_subtype="Normal"
)
def initial_time():
    """
    The initial time for the simulation.
    """
    return __data["time"].initial_time()


@component.add(
    name="SAVEPER",
    units="Day",
    limits=(0.0, np.nan),
    comp_type="Constant",
    comp_subtype="Normal",
)
def saveper():
    """
    The frequency with which output is stored.
    """
    return __data["time"].saveper()


@component.add(
    name="TIME STEP",
    units="Day",
    limits=(0.0, np.nan),
    comp_type="Constant",
    comp_subtype="Normal",
)
def time_step():
    """
    The time step for the simulation.
    """
    return __data["time"].time_step()


#######################################################################
#                           MODEL VARIABLES                           #
#######################################################################


@component.add(
    name="Open Duration",
    units="Day",
    limits=(0.0, 60.0, 1.0),
    comp_type="Constant",
    comp_subtype="Normal",
)
def open_duration():
    return 10


@component.add(
    name="End Lockdown Time",
    units="Day",
    limits=(0.0, 1000.0, 20.0),
    comp_type="Constant",
    comp_subtype="Normal",
)
def end_lockdown_time():
    return 400


@component.add(
    name="Lockdown Period",
    units="Dmnl",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={
        "lock_down_start_day": 1,
        "lock_down_duration": 2,
        "open_duration": 1,
        "end_lockdown_time": 1,
        "time": 1,
    },
)
def lockdown_period():
    return pulse(
        __data["time"],
        lock_down_start_day(),
        repeat_time=open_duration() + lock_down_duration(),
        width=lock_down_duration(),
        end=end_lockdown_time(),
    )


@component.add(
    name='"Fr symptomatic tested & isolated"',
    units="Dmnl",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={
        "low_fr_iso_sym": 2,
        "high_fr_iso_sym": 1,
        "high_contact_tracing_and_isolation": 1,
    },
)
def fr_symptomatic_tested_isolated():
    return low_fr_iso_sym() + high_contact_tracing_and_isolation() * (
        high_fr_iso_sym() - low_fr_iso_sym()
    )


@component.add(
    name="High Contact Tracing and Isolation",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={
        "time": 1,
        "april_14": 1,
        "policy_high_contact_tracing_and_isolation": 1,
    },
)
def high_contact_tracing_and_isolation():
    return if_then_else(
        time() < april_14(),
        lambda: 0,
        lambda: policy_high_contact_tracing_and_isolation(),
    )


@component.add(
    name="Contact Tracing Policy",
    units="Dmnl",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={
        "low_contact_tracing_policy": 2,
        "high_contact_tracing_policy": 1,
        "high_contact_tracing_and_isolation": 1,
    },
)
def contact_tracing_policy():
    """
    1 if contact tracing of Covid positive patients are done; 0 otherwise. Fractional values indicate the successfulness of the contact tracing.
    """
    return low_contact_tracing_policy() + high_contact_tracing_and_isolation() * (
        high_contact_tracing_policy() - low_contact_tracing_policy()
    )


@component.add(
    name="Net Fr requiring ICU",
    units="Dmnl",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"fr_multiplier_icu": 1, "fr_requiring_icu": 1},
)
def net_fr_requiring_icu():
    return np.minimum(0.4, fr_multiplier_icu() * fr_requiring_icu())


@component.add(
    name="Infection Rate new arrivals",
    units="Persons/Day",
    comp_type="Auxiliary",
    comp_subtype="with Lookup",
    depends_on={"time": 1},
)
def infection_rate_new_arrivals():
    return np.interp(
        time(),
        [
            0,
            1,
            2,
            3,
            4,
            5,
            6,
            7,
            8,
            9,
            10,
            11,
            12,
            13,
            14,
            15,
            16,
            17,
            18,
            19,
            20,
            21,
            22,
            23,
            24,
            25,
            26,
            27,
            28,
            29,
            30,
            31,
            32,
            33,
            34,
            35,
            36,
            37,
            38,
            39,
            40,
            1000,
        ],
        [
            0,
            2,
            1,
            0,
            2,
            1,
            2,
            3,
            6,
            2,
            7,
            10,
            4,
            10,
            5,
            11,
            14,
            14,
            26,
            34,
            35,
            45,
            55,
            38,
            34,
            24,
            40,
            16,
            20,
            12,
            23,
            14,
            8,
            14,
            12,
            5,
            9,
            6,
            0,
            0,
            0,
            0,
        ],
    )


@component.add(
    name="Contacts per day",
    units="1/Day",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={
        "lockdown_period": 1,
        "contact_under_lockdown_per_day": 1,
        "net_contacts_sum": 1,
    },
)
def contacts_per_day():
    return if_then_else(
        lockdown_period() == 1,
        lambda: contact_under_lockdown_per_day(),
        lambda: net_contacts_sum(),
    )


@component.add(
    name="Hygiene Mask Impact",
    units="Dmnl",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={
        "time": 1,
        "april_14": 1,
        "regular_hygiene_level": 1,
        "high_hygiene_level": 1,
    },
)
def hygiene_mask_impact():
    return if_then_else(
        time() < april_14(),
        lambda: regular_hygiene_level(),
        lambda: high_hygiene_level(),
    )


@component.add(
    name="Avg contact",
    units="1/Day",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"normal_contacts_per_day": 1},
)
def avg_contact():
    return normal_contacts_per_day()


@component.add(
    name="Net Fr fatality",
    units="Dmnl",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"fr_multiplier_fatality": 1, "fr_fatality": 1},
)
def net_fr_fatality():
    return np.minimum(0.4, fr_multiplier_fatality() * fr_fatality())


@component.add(
    name="Interaction Intensity",
    units="Dmnl",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={
        "hygiene_mask_impact": 1,
        "physical_distance_impact": 1,
        "lack_of_awareness_multiplier": 2,
        "april_14": 1,
        "time": 1,
    },
)
def interaction_intensity():
    return (
        hygiene_mask_impact()
        * physical_distance_impact()
        * if_then_else(
            time() < april_14(),
            lambda: lack_of_awareness_multiplier(),
            lambda: lack_of_awareness_multiplier(),
        )
    )


@component.add(
    name="Contact under lockdown per day",
    units="1/Day",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"normal_contacts_per_day": 1, "lockdown_effect_on_contacts": 1},
)
def contact_under_lockdown_per_day():
    return normal_contacts_per_day() * lockdown_effect_on_contacts()


@component.add(
    name="Contacts Iso Asym per day",
    units="1/Day",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"net_contacts_home": 1, "isolation_effect_on_asym_contacts": 1},
)
def contacts_iso_asym_per_day():
    return net_contacts_home() * isolation_effect_on_asym_contacts()


@component.add(
    name="Normal Contacts per day",
    units="1/Day",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={
        "base_contact_home": 1,
        "base_contact_school": 1,
        "base_contact_work": 1,
        "base_contacts_other": 1,
        "effect_of_density_on_contacts": 1,
    },
)
def normal_contacts_per_day():
    return (
        base_contact_home()
        + base_contact_school()
        + base_contact_work()
        + base_contacts_other()
    ) * effect_of_density_on_contacts()


@component.add(
    name="April 14", units="Day", comp_type="Constant", comp_subtype="Normal"
)
def april_14():
    return 50


@component.add(
    name="Policy High Contact Tracing and Isolation",
    units="Dmnl",
    limits=(0.0, 1.0, 1.0),
    comp_type="Constant",
    comp_subtype="Normal",
)
def policy_high_contact_tracing_and_isolation():
    return 0


@component.add(
    name="Contacts Q Asym per day",
    units="1/Day",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"quarantine_effect_on_asym_contacts": 1, "net_contacts_home": 1},
)
def contacts_q_asym_per_day():
    return quarantine_effect_on_asym_contacts() * net_contacts_home()


@component.add(
    name="Contacts Q Sym per day",
    units="1/Day",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"net_contacts_home": 1, "quarantine_effect_on_sym_contacts": 1},
)
def contacts_q_sym_per_day():
    return net_contacts_home() * quarantine_effect_on_sym_contacts()


@component.add(
    name="Contacts Iso Sym per day",
    units="1/Day",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"net_contacts_home": 1, "isolation_effect_on_sym_contacts": 1},
)
def contacts_iso_sym_per_day():
    return net_contacts_home() * isolation_effect_on_sym_contacts()


@component.add(
    name="Net Contacts Other",
    units="1/Day",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={
        "base_contacts_other": 1,
        "policy_multiplier_other": 1,
        "effect_of_density_on_contacts": 1,
    },
)
def net_contacts_other():
    return (
        base_contacts_other()
        * policy_multiplier_other()
        * effect_of_density_on_contacts()
    )


@component.add(
    name="Net Contacts School",
    units="1/Day",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={
        "base_contact_school": 1,
        "policy_multiplier_school": 1,
        "effect_of_density_on_contacts": 1,
    },
)
def net_contacts_school():
    return (
        base_contact_school()
        * policy_multiplier_school()
        * effect_of_density_on_contacts()
    )


@component.add(
    name="Base Contacts Other",
    units="1/Day",
    comp_type="Constant",
    comp_subtype="Normal",
)
def base_contacts_other():
    return 6


@component.add(
    name="Physical Distance Impact",
    units="Dmnl",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"density_ratio": 2},
)
def physical_distance_impact():
    return if_then_else(
        density_ratio() <= 1,
        lambda: 1,
        lambda: (np.log(density_ratio()) / np.log(50)) + 1,
    )


@component.add(
    name="Density of state",
    units="Persons/squareKM",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"total_population": 1, "area_of_states": 1},
)
def density_of_state():
    return total_population() / area_of_states()


@component.add(
    name="Density ratio",
    units="Dmnl",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"density_of_state": 1, "reference_population_density": 1},
)
def density_ratio():
    """
    Area and Density has to be properly defined
    """
    return (density_of_state() / reference_population_density()) * 0 + 1


@component.add(
    name="Effect of density on contacts",
    units="Dmnl",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"density_ratio": 2, "contacts_calibration_multiplier": 1},
)
def effect_of_density_on_contacts():
    return (
        if_then_else(
            density_ratio() <= 1,
            lambda: 1,
            lambda: (np.log(density_ratio()) / np.log(50)) + 1,
        )
        * contacts_calibration_multiplier()
    )


@component.add(
    name="Policy Multiplier Home",
    units="1/Day",
    limits=(0.0, 1.5, 0.1),
    comp_type="Constant",
    comp_subtype="Normal",
)
def policy_multiplier_home():
    return 1


@component.add(name="Base Contact Home", comp_type="Constant", comp_subtype="Normal")
def base_contact_home():
    return 8


@component.add(
    name="Base Contact School",
    units="1/Day",
    comp_type="Constant",
    comp_subtype="Normal",
)
def base_contact_school():
    return 5


@component.add(name="Base Contact Work", comp_type="Constant", comp_subtype="Normal")
def base_contact_work():
    return 7


@component.add(
    name="Net Contacts Sum",
    units="1/Day",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={
        "net_contacts_home": 1,
        "net_contacts_other": 1,
        "net_contacts_school": 1,
        "net_contacts_work": 1,
    },
)
def net_contacts_sum():
    return (
        net_contacts_home()
        + net_contacts_other()
        + net_contacts_school()
        + net_contacts_work()
    )


@component.add(
    name="Net Contacts Work",
    units="1/Day",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={
        "base_contact_work": 1,
        "policy_multiplier_work": 1,
        "effect_of_density_on_contacts": 1,
    },
)
def net_contacts_work():
    return (
        base_contact_work() * policy_multiplier_work() * effect_of_density_on_contacts()
    )


@component.add(
    name="Policy Multiplier School",
    units="Dmnl",
    limits=(0.0, 1.0, 0.1),
    comp_type="Constant",
    comp_subtype="Normal",
)
def policy_multiplier_school():
    return 1


@component.add(
    name="Policy Multiplier Work",
    units="Dmnl",
    limits=(0.0, 1.0, 0.1),
    comp_type="Constant",
    comp_subtype="Normal",
)
def policy_multiplier_work():
    return 1


@component.add(
    name="Reference Population Density",
    units="Persons/squareKM",
    comp_type="Constant",
    comp_subtype="Normal",
)
def reference_population_density():
    return 1


@component.add(
    name="Policy Multiplier Other",
    units="Dmnl",
    limits=(0.0, 1.0, 0.1),
    comp_type="Constant",
    comp_subtype="Normal",
)
def policy_multiplier_other():
    return 1


@component.add(
    name="Net Contacts Home",
    units="1/Day",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={
        "base_contact_home": 1,
        "policy_multiplier_home": 1,
        "effect_of_density_on_contacts": 1,
    },
)
def net_contacts_home():
    return (
        base_contact_home() * policy_multiplier_home() * effect_of_density_on_contacts()
    )


@component.add(
    name="Sum Infections from Symptomatics mild",
    units="1/Day",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"contacts_with_sypmtomatics_mild": 1},
)
def sum_infections_from_symptomatics_mild():
    return contacts_with_sypmtomatics_mild()


@component.add(
    name="Fraction Can't Work",
    units="Dmnl",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"sum_susc_asymp_recov_can_work": 1, "total_population": 1},
)
def fraction_cant_work():
    return 1 - sum_susc_asymp_recov_can_work() / total_population()


@component.add(
    name="Total Population in AgeGroup",
    units="Persons",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"total_population": 1},
)
def total_population_in_agegroup():
    return total_population()


@component.add(
    name="Contacts total per Susceptible",
    units="1/Day",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={
        "contacts_with_asymptomatics": 1,
        "contacts_with_asym_in_q": 1,
        "contacts_with_asym_in_isolation": 1,
        "contacts_with_sypmtomatics_mild": 1,
        "contacts_with_sym_in_q": 1,
        "contacts_with_sym_in_isolation": 1,
    },
)
def contacts_total_per_susceptible():
    return (
        contacts_with_asymptomatics()
        + contacts_with_asym_in_q()
        + contacts_with_asym_in_isolation()
        + contacts_with_sypmtomatics_mild()
        + contacts_with_sym_in_q()
        + contacts_with_sym_in_isolation()
    )


@component.add(
    name="Recovered Isolated Symptomatic Mild",
    comp_type="Stateful",
    comp_subtype="Integ",
    depends_on={"_integ_recovered_isolated_symptomatic_mild": 1},
    other_deps={
        "_integ_recovered_isolated_symptomatic_mild": {
            "initial": {},
            "step": {"recoveries_iso_sym": 1},
        }
    },
)
def recovered_isolated_symptomatic_mild():
    return _integ_recovered_isolated_symptomatic_mild()


_integ_recovered_isolated_symptomatic_mild = Integ(
    lambda: recoveries_iso_sym(),
    lambda: 0,
    "_integ_recovered_isolated_symptomatic_mild",
)


@component.add(
    name="lack of Awareness multiplier",
    units="Dmnl",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={
        "loastartlevel": 2,
        "loastarttime": 2,
        "loaendlevel": 1,
        "time": 1,
        "loaendtime": 2,
    },
)
def lack_of_awareness_multiplier():
    return loastartlevel() + ramp(
        __data["time"],
        (loaendlevel() - loastartlevel()) / (loaendtime() - loastarttime()),
        loastarttime(),
        loaendtime(),
    )


@component.add(name="LoAendtime", comp_type="Constant", comp_subtype="Normal")
def loaendtime():
    return 50


@component.add(name="LoAstarttime", comp_type="Constant", comp_subtype="Normal")
def loastarttime():
    return 15


@component.add(
    name="Regular Hygiene Level", comp_type="Constant", comp_subtype="Normal"
)
def regular_hygiene_level():
    return 1


@component.add(name="High Fr Iso Sym", comp_type="Constant", comp_subtype="Normal")
def high_fr_iso_sym():
    return 0.8


@component.add(
    name="High Contact Tracing Policy",
    limits=(0.0, 1.0, 0.1),
    comp_type="Constant",
    comp_subtype="Normal",
)
def high_contact_tracing_policy():
    return 0.8


@component.add(
    name="Low Contact Tracing Policy", comp_type="Constant", comp_subtype="Normal"
)
def low_contact_tracing_policy():
    return 0.1


@component.add(
    name="High Hygiene Level",
    limits=(0.0, 1.0, 0.05),
    comp_type="Constant",
    comp_subtype="Normal",
)
def high_hygiene_level():
    return 0.8


@component.add(name="Low Fr Iso Sym", comp_type="Constant", comp_subtype="Normal")
def low_fr_iso_sym():
    return 0.25


@component.add(name="LoAstartlevel", comp_type="Constant", comp_subtype="Normal")
def loastartlevel():
    return 1.2


@component.add(name="LoAendlevel", comp_type="Constant", comp_subtype="Normal")
def loaendlevel():
    return 1


@component.add(
    name="Default Ro",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={
        "avg_contact": 1,
        "default_infectivity": 1,
        "typical_duration_when_infection_is_spread": 1,
    },
)
def default_ro():
    return (
        avg_contact()
        * default_infectivity()
        * typical_duration_when_infection_is_spread()
    )


@component.add(
    name="Fr becoming serious",
    units="Dmnl",
    comp_type="Constant",
    comp_subtype="Normal",
)
def fr_becoming_serious():
    """
    0.02,0.02,0.02,0.02,0.17,0.17,0.17,0.17,0.17,0.28,0.28,0.27,0.27,0.4,0.4,0. 47
    """
    return 0.15


@component.add(
    name="Infectivity",
    units="Dmnl",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"default_infectivity": 1, "interaction_intensity": 1},
)
def infectivity():
    return default_infectivity() * interaction_intensity()


@component.add(
    name="Fr fatality", units="Dmnl", comp_type="Constant", comp_subtype="Normal"
)
def fr_fatality():
    return 0.23


@component.add(
    name="Default Infectivity",
    limits=(0.0, 0.03, 0.005),
    comp_type="Constant",
    comp_subtype="Normal",
)
def default_infectivity():
    return 0.015


@component.add(
    name="default delay disease diagnosis",
    units="Day",
    limits=(0.0, 20.0),
    comp_type="Constant",
    comp_subtype="Normal",
)
def default_delay_disease_diagnosis():
    return 5


@component.add(
    name="testing impact on delay",
    limits=(0.0, 1.0, 0.05),
    comp_type="Constant",
    comp_subtype="Normal",
)
def testing_impact_on_delay():
    return 1


@component.add(
    name="delay disease diagnosis",
    units="Day",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"default_delay_disease_diagnosis": 1, "testing_impact_on_delay": 1},
)
def delay_disease_diagnosis():
    return default_delay_disease_diagnosis() * testing_impact_on_delay()


@component.add(
    name="delay Asymp Noninfective to Infective",
    units="Day",
    comp_type="Constant",
    comp_subtype="Normal",
)
def delay_asymp_noninfective_to_infective():
    return 3


@component.add(
    name="Contacts Calibration Multiplier",
    units="Dmnl",
    limits=(0.0, 3.0, 0.1),
    comp_type="Constant",
    comp_subtype="Normal",
)
def contacts_calibration_multiplier():
    return 1.1


@component.add(
    name="Area of States", units="squareKM", comp_type="Constant", comp_subtype="Normal"
)
def area_of_states():
    return 1


@component.add(
    name="Base Infectivity",
    units="Dmnl",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={
        "default_ro_base_infectivity_estimation": 1,
        "typical_duration_when_infection_is_spread": 1,
        "avg_contact": 1,
    },
)
def base_infectivity():
    return default_ro_base_infectivity_estimation() / (
        avg_contact() * typical_duration_when_infection_is_spread()
    )


@component.add(
    name="PPE availability",
    units="Dmnl",
    limits=(0.0, 1.0, 0.1),
    comp_type="Constant",
    comp_subtype="Normal",
)
def ppe_availability():
    return 0.8


@component.add(
    name="Fraction Healthcare Can't Work",
    units="Dmnl",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"fraction_cant_work": 2, "ppe_availability": 1},
)
def fraction_healthcare_cant_work():
    return fraction_cant_work() + (1 - fraction_cant_work()) * (1 - ppe_availability())


@component.add(
    name="Sum Susc Asymp Recov Can Work",
    units="Persons",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"sum_asymtomatics": 1, "sum_recovered": 1, "sum_susceptibles": 1},
)
def sum_susc_asymp_recov_can_work():
    return sum_asymtomatics() + sum_recovered() + sum_susceptibles()


@component.add(
    name="Fr Q Asym via awareness",
    units="Dmnl",
    limits=(0.0, 1.0, 0.05),
    comp_type="Constant",
    comp_subtype="Normal",
)
def fr_q_asym_via_awareness():
    return 0


@component.add(
    name="Fr Q Asym via tracing",
    units="Dmnl",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"fr_symptomatic_tested_isolated": 1, "contact_tracing_policy": 1},
)
def fr_q_asym_via_tracing():
    return fr_symptomatic_tested_isolated() * contact_tracing_policy()


@component.add(
    name="fraction Q Asym",
    units="Dmnl",
    limits=(0.0, 1.0),
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"fr_q_asym_via_awareness": 1, "fr_q_asym_via_tracing": 1},
)
def fraction_q_asym():
    return fr_q_asym_via_awareness() + fr_q_asym_via_tracing()


@component.add(
    name="delay for death", units="Day", comp_type="Constant", comp_subtype="Normal"
)
def delay_for_death():
    return 5


@component.add(
    name="delay Asymp Infective to Symp",
    units="Day",
    comp_type="Constant",
    comp_subtype="Normal",
)
def delay_asymp_infective_to_symp():
    return 2


@component.add(
    name="delay Iso Asym", units="Day", comp_type="Constant", comp_subtype="Normal"
)
def delay_iso_asym():
    return 1


@component.add(
    name="delay Iso Sym",
    units="Day",
    limits=(0.0, 5.0, 0.5),
    comp_type="Constant",
    comp_subtype="Normal",
)
def delay_iso_sym():
    return 2.5


@component.add(
    name="delay Q Asym", units="Day", comp_type="Constant", comp_subtype="Normal"
)
def delay_q_asym():
    return 1


@component.add(
    name="delay Q Sym", units="Day", comp_type="Constant", comp_subtype="Normal"
)
def delay_q_sym():
    return 2.5


@component.add(
    name="delay Recovery Asym", units="Day", comp_type="Constant", comp_subtype="Normal"
)
def delay_recovery_asym():
    return 11


@component.add(
    name="delay Recovery Sym extreme",
    units="Day",
    comp_type="Constant",
    comp_subtype="Normal",
)
def delay_recovery_sym_extreme():
    return 14


@component.add(
    name="delay Recovery Sym mild",
    units="Day",
    comp_type="Constant",
    comp_subtype="Normal",
)
def delay_recovery_sym_mild():
    return 14


@component.add(
    name="delay Recovery Sym serious",
    units="Day",
    comp_type="Constant",
    comp_subtype="Normal",
)
def delay_recovery_sym_serious():
    return 14


@component.add(
    name="delay Worsening", units="Day", comp_type="Constant", comp_subtype="Normal"
)
def delay_worsening():
    return 5


@component.add(
    name="Fr developing symptoms",
    units="Dmnl",
    comp_type="Constant",
    comp_subtype="Normal",
)
def fr_developing_symptoms():
    return 0.7


@component.add(
    name="fraction Iso Asym",
    units="Dmnl",
    limits=(0.0, 1.0),
    comp_type="Constant",
    comp_subtype="Normal",
)
def fraction_iso_asym():
    return 0


@component.add(
    name='"Fr Iso+Q Asym"',
    units="Dmnl",
    limits=(0.0, 1.0),
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"fraction_iso_asym": 1, "fraction_q_asym": 1},
)
def fr_isoq_asym():
    return fraction_iso_asym() + fraction_q_asym()


@component.add(
    name='"Fr Iso+Q Sym"',
    units="Dmnl",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"fr_symptomatic_tested_isolated": 1, "fraction_q_sym": 1},
)
def fr_isoq_sym():
    return fr_symptomatic_tested_isolated() + fraction_q_sym()


@component.add(
    name="Fr multiplier fatality",
    units="Dmnl",
    comp_type="Constant",
    comp_subtype="Normal",
)
def fr_multiplier_fatality():
    return 2


@component.add(
    name="Fr multiplier ICU", units="Dmnl", comp_type="Constant", comp_subtype="Normal"
)
def fr_multiplier_icu():
    return 2


@component.add(
    name="fraction Q Sym",
    units="Dmnl",
    limits=(0.0, 1.0),
    comp_type="Constant",
    comp_subtype="Normal",
)
def fraction_q_sym():
    return 0


@component.add(
    name="Fr requiring ICU", units="Dmnl", comp_type="Constant", comp_subtype="Normal"
)
def fr_requiring_icu():
    return 0.22


@component.add(
    name="Healthcare Capacity Multiplier",
    units="Dmnl",
    limits=(0.0, 500.0, 10.0),
    comp_type="Constant",
    comp_subtype="Normal",
)
def healthcare_capacity_multiplier():
    return 200


@component.add(
    name="Infectivity reduction for asym",
    units="Dmnl",
    limits=(0.0, 1.0, 0.1),
    comp_type="Constant",
    comp_subtype="Normal",
)
def infectivity_reduction_for_asym():
    return 0.5


@component.add(
    name="Isolation effect on Asym contacts",
    units="Dmnl",
    comp_type="Constant",
    comp_subtype="Normal",
)
def isolation_effect_on_asym_contacts():
    return 0.5


@component.add(
    name="Isolation effect on Sym contacts",
    units="Dmnl",
    comp_type="Constant",
    comp_subtype="Normal",
)
def isolation_effect_on_sym_contacts():
    return 0.5


@component.add(
    name="Lock down duration",
    units="Day",
    limits=(0.0, 60.0, 1.0),
    comp_type="Constant",
    comp_subtype="Normal",
)
def lock_down_duration():
    """
    Indicates number of days of lockdown. Keep as 0 for no lockdown
    """
    return 14


@component.add(
    name="Lock down start day",
    units="Day",
    limits=(0.0, 30.0, 5.0),
    comp_type="Constant",
    comp_subtype="Normal",
)
def lock_down_start_day():
    return 30


@component.add(
    name="Lockdown effect on contacts",
    units="Dmnl",
    limits=(0.0, 8.0),
    comp_type="Constant",
    comp_subtype="Normal",
)
def lockdown_effect_on_contacts():
    return 0.3


@component.add(
    name="net delay Dprogress Iso Sym",
    units="Day",
    limits=(0.0, 100.0),
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"delay_disease_diagnosis": 1, "delay_iso_sym": 1},
)
def net_delay_dprogress_iso_sym():
    return np.maximum(delay_disease_diagnosis() - delay_iso_sym(), 1)


@component.add(
    name="net delay Dprogress Q Sym",
    units="Day",
    limits=(0.0, 100.0),
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"delay_disease_diagnosis": 1, "delay_q_sym": 1},
)
def net_delay_dprogress_q_sym():
    return np.maximum(delay_disease_diagnosis() - delay_q_sym(), 1)


@component.add(
    name="net delay incubation Iso ASym",
    units="Day",
    limits=(0.0, 100.0),
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"delay_asymp_infective_to_symp": 1, "delay_iso_asym": 1},
)
def net_delay_incubation_iso_asym():
    return np.maximum(delay_asymp_infective_to_symp() - delay_iso_asym(), 1)


@component.add(
    name="net delay incubation Q Asym",
    units="Day",
    limits=(0.0, 100.0),
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"delay_asymp_infective_to_symp": 1, "delay_q_asym": 1},
)
def net_delay_incubation_q_asym():
    return np.maximum(delay_asymp_infective_to_symp() - delay_q_asym(), 1)


@component.add(
    name="net delay recovery Iso Asym",
    units="Day",
    limits=(0.0, 100.0),
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"delay_recovery_asym": 1, "delay_iso_asym": 1},
)
def net_delay_recovery_iso_asym():
    return np.maximum(delay_recovery_asym() - delay_iso_asym(), 1)


@component.add(
    name="net delay recovery Iso Sym",
    units="Day",
    limits=(0.0, 100.0),
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"delay_recovery_sym_mild": 1, "delay_iso_sym": 1},
)
def net_delay_recovery_iso_sym():
    return np.maximum(delay_recovery_sym_mild() - delay_iso_sym(), 1)


@component.add(
    name="net delay recovery Q Asym",
    units="Day",
    limits=(0.0, 100.0),
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"delay_recovery_asym": 1, "delay_q_asym": 1},
)
def net_delay_recovery_q_asym():
    return np.maximum(delay_recovery_asym() - delay_q_asym(), 1)


@component.add(
    name="net delay recovery Q Sym",
    units="Day",
    limits=(0.0, 100.0),
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"delay_recovery_sym_mild": 1, "delay_q_sym": 1},
)
def net_delay_recovery_q_sym():
    return np.maximum(delay_recovery_sym_mild() - delay_q_sym(), 1)


@component.add(
    name="Normal Hospital Beds per person",
    units="Units/ Persons",
    comp_type="Constant",
    comp_subtype="Normal",
)
def normal_hospital_beds_per_person():
    return 0.9 / 1000


@component.add(
    name="Normal ICU capacity per person",
    units="Units/Persons",
    comp_type="Constant",
    comp_subtype="Normal",
)
def normal_icu_capacity_per_person():
    return 2.3 / 100000


@component.add(
    name="Patients per unit",
    units="Persons/Units",
    comp_type="Constant",
    comp_subtype="Normal",
)
def patients_per_unit():
    return 1


@component.add(
    name="Quarantine effect on Asym contacts",
    units="Dmnl",
    comp_type="Constant",
    comp_subtype="Normal",
)
def quarantine_effect_on_asym_contacts():
    return 1


@component.add(
    name="Quarantine effect on Sym contacts",
    units="Dmnl",
    comp_type="Constant",
    comp_subtype="Normal",
)
def quarantine_effect_on_sym_contacts():
    return 0.75


@component.add(
    name="Default Ro Base Infectivity Estimation",
    units="Dmnl",
    limits=(0.0, 8.0, 0.2),
    comp_type="Constant",
    comp_subtype="Normal",
)
def default_ro_base_infectivity_estimation():
    return 4.6


@component.add(
    name="Sum Infections from Asymtomatics",
    units="1/Day",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"contacts_with_asymptomatics": 1},
)
def sum_infections_from_asymtomatics():
    return contacts_with_asymptomatics()


@component.add(
    name="Total Population",
    units="Persons",
    comp_type="Constant",
    comp_subtype="Normal",
)
def total_population():
    """
    1.3392e+09
    """
    return 1339200000.0


@component.add(
    name="Typical Duration when infection is spread",
    units="Day",
    limits=(5.0, 30.0, 1.0),
    comp_type="Constant",
    comp_subtype="Normal",
)
def typical_duration_when_infection_is_spread():
    return 7


@component.add(
    name="Typical worsening rate",
    units="Persons/Day",
    comp_type="Stateful",
    comp_subtype="Delay",
    depends_on={"fr_requiring_icu": 1, "_delay_typical_worsening_rate": 1},
    other_deps={
        "_delay_typical_worsening_rate": {
            "initial": {"hospital_admissions": 1, "delay_worsening": 1},
            "step": {"hospital_admissions": 1, "delay_worsening": 1},
        }
    },
)
def typical_worsening_rate():
    return fr_requiring_icu() * _delay_typical_worsening_rate()


_delay_typical_worsening_rate = Delay(
    lambda: hospital_admissions(),
    lambda: delay_worsening(),
    lambda: 0 * hospital_admissions(),
    lambda: 3,
    time_step,
    "_delay_typical_worsening_rate",
)


@component.add(
    name="Recoveries seriousH",
    units="Persons/Day",
    comp_type="Stateful",
    comp_subtype="Delay",
    depends_on={"fr_requiring_icu": 1, "_delay_recoveries_serioush": 1},
    other_deps={
        "_delay_recoveries_serioush": {
            "initial": {"hospital_admissions": 1, "delay_recovery_sym_serious": 1},
            "step": {"hospital_admissions": 1, "delay_recovery_sym_serious": 1},
        }
    },
)
def recoveries_serioush():
    return (1 - fr_requiring_icu()) * _delay_recoveries_serioush()


_delay_recoveries_serioush = Delay(
    lambda: hospital_admissions(),
    lambda: delay_recovery_sym_serious(),
    lambda: 0 * hospital_admissions(),
    lambda: 3,
    time_step,
    "_delay_recoveries_serioush",
)


@component.add(
    name="Recoveries ICU overflow",
    units="Persons/Day",
    comp_type="Stateful",
    comp_subtype="Delay",
    depends_on={"_delay_recoveries_icu_overflow": 1, "net_fr_fatality": 1},
    other_deps={
        "_delay_recoveries_icu_overflow": {
            "initial": {
                "into_inf_sym_icu_overflow": 1,
                "delay_recovery_sym_extreme": 1,
            },
            "step": {"into_inf_sym_icu_overflow": 1, "delay_recovery_sym_extreme": 1},
        }
    },
)
def recoveries_icu_overflow():
    return _delay_recoveries_icu_overflow() * (1 - net_fr_fatality())


_delay_recoveries_icu_overflow = Delay(
    lambda: into_inf_sym_icu_overflow(),
    lambda: delay_recovery_sym_extreme(),
    lambda: 0 * into_inf_sym_icu_overflow(),
    lambda: 3,
    time_step,
    "_delay_recoveries_icu_overflow",
)


@component.add(
    name="Deaths overflow",
    units="Persons/Day",
    comp_type="Stateful",
    comp_subtype="Delay",
    depends_on={"_delay_deaths_overflow": 1, "net_fr_fatality": 1},
    other_deps={
        "_delay_deaths_overflow": {
            "initial": {"into_inf_sym_icu_overflow": 1, "delay_for_death": 1},
            "step": {"into_inf_sym_icu_overflow": 1, "delay_for_death": 1},
        }
    },
)
def deaths_overflow():
    return _delay_deaths_overflow() * net_fr_fatality()


_delay_deaths_overflow = Delay(
    lambda: into_inf_sym_icu_overflow(),
    lambda: delay_for_death(),
    lambda: 0 * into_inf_sym_icu_overflow(),
    lambda: 3,
    time_step,
    "_delay_deaths_overflow",
)


@component.add(
    name="Incubation Rate Quarantines",
    units="Persons/Day",
    comp_type="Stateful",
    comp_subtype="Delay",
    depends_on={"fr_developing_symptoms": 1, "_delay_incubation_rate_quarantines": 1},
    other_deps={
        "_delay_incubation_rate_quarantines": {
            "initial": {"quarantining_rate_asym": 1, "net_delay_incubation_q_asym": 1},
            "step": {"quarantining_rate_asym": 1, "net_delay_incubation_q_asym": 1},
        }
    },
)
def incubation_rate_quarantines():
    return fr_developing_symptoms() * _delay_incubation_rate_quarantines()


_delay_incubation_rate_quarantines = Delay(
    lambda: quarantining_rate_asym(),
    lambda: net_delay_incubation_q_asym(),
    lambda: 0 * quarantining_rate_asym(),
    lambda: 3,
    time_step,
    "_delay_incubation_rate_quarantines",
)


@component.add(
    name="Infected Sym Extreme ICU",
    units="Persons",
    comp_type="Stateful",
    comp_subtype="Integ",
    depends_on={"_integ_infected_sym_extreme_icu": 1},
    other_deps={
        "_integ_infected_sym_extreme_icu": {
            "initial": {},
            "step": {"icu_admissions": 1, "deaths": 1, "recoveries_icu": 1},
        }
    },
)
def infected_sym_extreme_icu():
    return _integ_infected_sym_extreme_icu()


_integ_infected_sym_extreme_icu = Integ(
    lambda: icu_admissions() - deaths() - recoveries_icu(),
    lambda: 0,
    "_integ_infected_sym_extreme_icu",
)


@component.add(
    name="Recoveries Q Asym",
    units="Persons/Day",
    comp_type="Stateful",
    comp_subtype="Delay",
    depends_on={"fr_developing_symptoms": 1, "_delay_recoveries_q_asym": 1},
    other_deps={
        "_delay_recoveries_q_asym": {
            "initial": {"quarantining_rate_asym": 1, "net_delay_recovery_q_asym": 1},
            "step": {"quarantining_rate_asym": 1, "net_delay_recovery_q_asym": 1},
        }
    },
)
def recoveries_q_asym():
    return (1 - fr_developing_symptoms()) * _delay_recoveries_q_asym()


_delay_recoveries_q_asym = Delay(
    lambda: quarantining_rate_asym(),
    lambda: net_delay_recovery_q_asym(),
    lambda: 0 * quarantining_rate_asym(),
    lambda: 3,
    time_step,
    "_delay_recoveries_q_asym",
)


@component.add(
    name="Recoveries QSym",
    units="Persons/Day",
    comp_type="Stateful",
    comp_subtype="Delay",
    depends_on={"fr_becoming_serious": 1, "_delay_recoveries_qsym": 1},
    other_deps={
        "_delay_recoveries_qsym": {
            "initial": {"into_q_sym": 1, "net_delay_recovery_q_sym": 1},
            "step": {"into_q_sym": 1, "net_delay_recovery_q_sym": 1},
        }
    },
)
def recoveries_qsym():
    return (1 - fr_becoming_serious()) * _delay_recoveries_qsym()


_delay_recoveries_qsym = Delay(
    lambda: into_q_sym(),
    lambda: net_delay_recovery_q_sym(),
    lambda: 0 * into_q_sym(),
    lambda: 3,
    time_step,
    "_delay_recoveries_qsym",
)


@component.add(
    name="Contacts with Sym in Q",
    units="1/Day",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={
        "contacts_q_sym_per_day": 1,
        "quarantined_symtomatics_mild": 1,
        "total_population_in_agegroup": 1,
    },
)
def contacts_with_sym_in_q():
    return (
        contacts_q_sym_per_day()
        * quarantined_symtomatics_mild()
        / total_population_in_agegroup()
    )


@component.add(
    name="Isolated Symtomatics mild",
    units="Persons",
    comp_type="Stateful",
    comp_subtype="Integ",
    depends_on={"_integ_isolated_symtomatics_mild": 1},
    other_deps={
        "_integ_isolated_symtomatics_mild": {
            "initial": {},
            "step": {
                "incubation_rate_iso": 1,
                "isolation_rate_sym": 1,
                "iso_disease_progress_rate": 1,
                "recoveries_iso_sym": 1,
            },
        }
    },
)
def isolated_symtomatics_mild():
    return _integ_isolated_symtomatics_mild()


_integ_isolated_symtomatics_mild = Integ(
    lambda: incubation_rate_iso()
    + isolation_rate_sym()
    - iso_disease_progress_rate()
    - recoveries_iso_sym(),
    lambda: 0,
    "_integ_isolated_symtomatics_mild",
)


@component.add(
    name="Disease Progression",
    units="Persons/Day",
    limits=(0.0, np.nan),
    comp_type="Stateful",
    comp_subtype="Delay",
    depends_on={
        "fr_isoq_sym": 1,
        "fr_becoming_serious": 1,
        "_delay_disease_progression": 1,
    },
    other_deps={
        "_delay_disease_progression": {
            "initial": {"symptom_setting": 1, "delay_disease_diagnosis": 1},
            "step": {"symptom_setting": 1, "delay_disease_diagnosis": 1},
        }
    },
)
def disease_progression():
    return (
        (1 - fr_isoq_sym())
        * np.minimum(0.9, fr_becoming_serious())
        * _delay_disease_progression()
    )


_delay_disease_progression = Delay(
    lambda: symptom_setting(),
    lambda: delay_disease_diagnosis(),
    lambda: 0 * symptom_setting(),
    lambda: 3,
    time_step,
    "_delay_disease_progression",
)


@component.add(
    name="Dead",
    units="Persons",
    comp_type="Stateful",
    comp_subtype="Integ",
    depends_on={"_integ_dead": 1},
    other_deps={
        "_integ_dead": {"initial": {}, "step": {"deaths": 1, "deaths_overflow": 1}}
    },
)
def dead():
    return _integ_dead()


_integ_dead = Integ(lambda: deaths() + deaths_overflow(), lambda: 0, "_integ_dead")


@component.add(
    name="Deaths",
    units="Persons/Day",
    comp_type="Stateful",
    comp_subtype="Delay",
    depends_on={"fr_fatality": 1, "_delay_deaths": 1},
    other_deps={
        "_delay_deaths": {
            "initial": {"icu_admissions": 1, "delay_for_death": 1},
            "step": {"icu_admissions": 1, "delay_for_death": 1},
        }
    },
)
def deaths():
    return fr_fatality() * _delay_deaths()


_delay_deaths = Delay(
    lambda: icu_admissions(),
    lambda: delay_for_death(),
    lambda: 0 * icu_admissions(),
    lambda: 3,
    time_step,
    "_delay_deaths",
)


@component.add(
    name="Sum Deaths reported daily",
    units="Persons/Day",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"deaths": 1, "deaths_overflow": 1},
)
def sum_deaths_reported_daily():
    return deaths() + deaths_overflow()


@component.add(
    name="Recoveries Iso Sym",
    comp_type="Stateful",
    comp_subtype="Delay",
    depends_on={"fr_becoming_serious": 1, "_delay_recoveries_iso_sym": 1},
    other_deps={
        "_delay_recoveries_iso_sym": {
            "initial": {"into_iso_sym": 1, "net_delay_recovery_iso_sym": 1},
            "step": {"into_iso_sym": 1, "net_delay_recovery_iso_sym": 1},
        }
    },
)
def recoveries_iso_sym():
    return (1 - np.minimum(0.9, fr_becoming_serious())) * _delay_recoveries_iso_sym()


_delay_recoveries_iso_sym = Delay(
    lambda: into_iso_sym(),
    lambda: net_delay_recovery_iso_sym(),
    lambda: 0 * into_iso_sym(),
    lambda: 3,
    time_step,
    "_delay_recoveries_iso_sym",
)


@component.add(
    name="Recoveries Iso Asym",
    units="Persons/Day",
    comp_type="Stateful",
    comp_subtype="Delay",
    depends_on={"fr_developing_symptoms": 1, "_delay_recoveries_iso_asym": 1},
    other_deps={
        "_delay_recoveries_iso_asym": {
            "initial": {"isolation_rate_asym": 1, "net_delay_recovery_iso_asym": 1},
            "step": {"isolation_rate_asym": 1, "net_delay_recovery_iso_asym": 1},
        }
    },
)
def recoveries_iso_asym():
    return (1 - fr_developing_symptoms()) * _delay_recoveries_iso_asym()


_delay_recoveries_iso_asym = Delay(
    lambda: isolation_rate_asym(),
    lambda: net_delay_recovery_iso_asym(),
    lambda: 0 * isolation_rate_asym(),
    lambda: 3,
    time_step,
    "_delay_recoveries_iso_asym",
)


@component.add(
    name="Symptom Setting",
    units="Persons/Day",
    comp_type="Stateful",
    comp_subtype="Delay",
    depends_on={
        "fr_isoq_asym": 1,
        "fr_developing_symptoms": 1,
        "_delay_symptom_setting": 1,
    },
    other_deps={
        "_delay_symptom_setting": {
            "initial": {"infection_rate_sum": 1, "delay_asymp_infective_to_symp": 1},
            "step": {"infection_rate_sum": 1, "delay_asymp_infective_to_symp": 1},
        }
    },
)
def symptom_setting():
    return (1 - fr_isoq_asym()) * fr_developing_symptoms() * _delay_symptom_setting()


_delay_symptom_setting = Delay(
    lambda: infection_rate_sum(),
    lambda: delay_asymp_infective_to_symp(),
    lambda: 0 * infection_rate_sum(),
    lambda: 3,
    time_step,
    "_delay_symptom_setting",
)


@component.add(
    name="Isolated Asymptomatics",
    units="Persons",
    comp_type="Stateful",
    comp_subtype="Integ",
    depends_on={"_integ_isolated_asymptomatics": 1},
    other_deps={
        "_integ_isolated_asymptomatics": {
            "initial": {},
            "step": {
                "isolation_rate_asym": 1,
                "recoveries_iso_asym": 1,
                "incubation_rate_iso": 1,
            },
        }
    },
)
def isolated_asymptomatics():
    return _integ_isolated_asymptomatics()


_integ_isolated_asymptomatics = Integ(
    lambda: isolation_rate_asym() - recoveries_iso_asym() - incubation_rate_iso(),
    lambda: 0,
    "_integ_isolated_asymptomatics",
)


@component.add(
    name="Iso Disease progress rate",
    units="Persons/Day",
    limits=(0.0, np.nan),
    comp_type="Stateful",
    comp_subtype="Delay",
    depends_on={"fr_becoming_serious": 1, "_delay_iso_disease_progress_rate": 1},
    other_deps={
        "_delay_iso_disease_progress_rate": {
            "initial": {"into_iso_sym": 1, "net_delay_dprogress_iso_sym": 1},
            "step": {"into_iso_sym": 1, "net_delay_dprogress_iso_sym": 1},
        }
    },
)
def iso_disease_progress_rate():
    return fr_becoming_serious() * _delay_iso_disease_progress_rate()


_delay_iso_disease_progress_rate = Delay(
    lambda: into_iso_sym(),
    lambda: net_delay_dprogress_iso_sym(),
    lambda: 0 * into_iso_sym(),
    lambda: 3,
    time_step,
    "_delay_iso_disease_progress_rate",
)


@component.add(
    name="Incubation Rate Iso",
    units="Persons/Day",
    comp_type="Stateful",
    comp_subtype="Delay",
    depends_on={"fr_developing_symptoms": 1, "_delay_incubation_rate_iso": 1},
    other_deps={
        "_delay_incubation_rate_iso": {
            "initial": {"isolation_rate_asym": 1, "net_delay_incubation_iso_asym": 1},
            "step": {"isolation_rate_asym": 1, "net_delay_incubation_iso_asym": 1},
        }
    },
)
def incubation_rate_iso():
    return fr_developing_symptoms() * _delay_incubation_rate_iso()


_delay_incubation_rate_iso = Delay(
    lambda: isolation_rate_asym(),
    lambda: net_delay_incubation_iso_asym(),
    lambda: 0 * isolation_rate_asym(),
    lambda: 3,
    time_step,
    "_delay_incubation_rate_iso",
)


@component.add(
    name="INTO Iso Sym",
    units="Persons/Day",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"incubation_rate_iso": 1, "isolation_rate_sym": 1},
)
def into_iso_sym():
    return incubation_rate_iso() + isolation_rate_sym()


@component.add(
    name="Q Disease Progress rate",
    units="Persons/Day",
    limits=(0.0, np.nan),
    comp_type="Stateful",
    comp_subtype="Delay",
    depends_on={"fr_becoming_serious": 1, "_delay_q_disease_progress_rate": 1},
    other_deps={
        "_delay_q_disease_progress_rate": {
            "initial": {"into_q_sym": 1, "net_delay_dprogress_q_sym": 1},
            "step": {"into_q_sym": 1, "net_delay_dprogress_q_sym": 1},
        }
    },
)
def q_disease_progress_rate():
    return np.minimum(0.9, fr_becoming_serious()) * _delay_q_disease_progress_rate()


_delay_q_disease_progress_rate = Delay(
    lambda: into_q_sym(),
    lambda: net_delay_dprogress_q_sym(),
    lambda: 0 * into_q_sym(),
    lambda: 3,
    time_step,
    "_delay_q_disease_progress_rate",
)


@component.add(
    name="Recoveries ICU",
    units="Persons/Day",
    comp_type="Stateful",
    comp_subtype="Delay",
    depends_on={"fr_fatality": 1, "_delay_recoveries_icu": 1},
    other_deps={
        "_delay_recoveries_icu": {
            "initial": {"icu_admissions": 1, "delay_recovery_sym_extreme": 1},
            "step": {"icu_admissions": 1, "delay_recovery_sym_extreme": 1},
        }
    },
)
def recoveries_icu():
    return (1 - np.minimum(0.9, fr_fatality())) * _delay_recoveries_icu()


_delay_recoveries_icu = Delay(
    lambda: icu_admissions(),
    lambda: delay_recovery_sym_extreme(),
    lambda: 0 * icu_admissions(),
    lambda: 3,
    time_step,
    "_delay_recoveries_icu",
)


@component.add(
    name="Recoveries Infected Sym mild",
    units="Persons/Day",
    comp_type="Stateful",
    comp_subtype="Delay",
    depends_on={
        "fr_isoq_sym": 1,
        "fr_becoming_serious": 1,
        "_delay_recoveries_infected_sym_mild": 1,
    },
    other_deps={
        "_delay_recoveries_infected_sym_mild": {
            "initial": {"symptom_setting": 1, "delay_recovery_sym_mild": 1},
            "step": {"symptom_setting": 1, "delay_recovery_sym_mild": 1},
        }
    },
)
def recoveries_infected_sym_mild():
    return (
        (1 - fr_isoq_sym())
        * (1 - np.minimum(0.9, fr_becoming_serious()))
        * _delay_recoveries_infected_sym_mild()
    )


_delay_recoveries_infected_sym_mild = Delay(
    lambda: symptom_setting(),
    lambda: delay_recovery_sym_mild(),
    lambda: 0 * symptom_setting(),
    lambda: 3,
    time_step,
    "_delay_recoveries_infected_sym_mild",
)


@component.add(
    name="Recoveries Asym",
    units="Persons/Day",
    comp_type="Stateful",
    comp_subtype="Delay",
    depends_on={
        "fr_isoq_asym": 1,
        "fr_developing_symptoms": 1,
        "_delay_recoveries_asym": 1,
    },
    other_deps={
        "_delay_recoveries_asym": {
            "initial": {"infection_rate_sum": 1, "delay_recovery_asym": 1},
            "step": {"infection_rate_sum": 1, "delay_recovery_asym": 1},
        }
    },
)
def recoveries_asym():
    return (
        (1 - fr_isoq_asym()) * (1 - fr_developing_symptoms()) * _delay_recoveries_asym()
    )


_delay_recoveries_asym = Delay(
    lambda: infection_rate_sum(),
    lambda: delay_recovery_asym(),
    lambda: 0 * infection_rate_sum(),
    lambda: 3,
    time_step,
    "_delay_recoveries_asym",
)


@component.add(
    name="Contacts with Sypmtomatics mild",
    units="1/Day",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={
        "contacts_per_day": 1,
        "infectious_symptomatics_mild": 1,
        "total_population_in_agegroup": 1,
    },
)
def contacts_with_sypmtomatics_mild():
    return (
        contacts_per_day()
        * infectious_symptomatics_mild()
        / total_population_in_agegroup()
    )


@component.add(
    name="Contacts with Asym in Q",
    units="1/Day",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={
        "contacts_q_asym_per_day": 1,
        "quarantined_asymptomatics": 1,
        "total_population_in_agegroup": 1,
        "infectivity_reduction_for_asym": 1,
    },
)
def contacts_with_asym_in_q():
    return (
        contacts_q_asym_per_day()
        * quarantined_asymptomatics()
        / total_population_in_agegroup()
    ) * (1 - infectivity_reduction_for_asym())


@component.add(
    name="Contacts with Asymptomatics",
    units="1/Day",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={
        "contacts_per_day": 1,
        "infectious_asymptomatics": 1,
        "total_population_in_agegroup": 1,
        "infectivity_reduction_for_asym": 1,
    },
)
def contacts_with_asymptomatics():
    return (
        contacts_per_day() * infectious_asymptomatics() / total_population_in_agegroup()
    ) * (1 - infectivity_reduction_for_asym())


@component.add(
    name="Contacts with Asym in Isolation",
    units="1/Day",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={
        "contacts_iso_asym_per_day": 1,
        "isolated_asymptomatics": 1,
        "total_population_in_agegroup": 1,
        "infectivity_reduction_for_asym": 1,
    },
)
def contacts_with_asym_in_isolation():
    return (
        contacts_iso_asym_per_day()
        * isolated_asymptomatics()
        / total_population_in_agegroup()
    ) * (1 - infectivity_reduction_for_asym())


@component.add(
    name="Contacts with Sym in Isolation",
    units="1/Day",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={
        "contacts_iso_sym_per_day": 1,
        "isolated_symtomatics_mild": 1,
        "total_population_in_agegroup": 1,
    },
)
def contacts_with_sym_in_isolation():
    return (
        contacts_iso_sym_per_day()
        * isolated_symtomatics_mild()
        / total_population_in_agegroup()
    )


@component.add(
    name="Cumulative Cases Reported",
    units="Persons",
    comp_type="Stateful",
    comp_subtype="Integ",
    depends_on={"_integ_cumulative_cases_reported": 1},
    other_deps={
        "_integ_cumulative_cases_reported": {
            "initial": {},
            "step": {"new_cases_reported": 1},
        }
    },
)
def cumulative_cases_reported():
    return _integ_cumulative_cases_reported()


_integ_cumulative_cases_reported = Integ(
    lambda: new_cases_reported(), lambda: 0, "_integ_cumulative_cases_reported"
)


@component.add(
    name="Exposed",
    units="Persons",
    comp_type="Stateful",
    comp_subtype="Integ",
    depends_on={"_integ_exposed": 1},
    other_deps={
        "_integ_exposed": {
            "initial": {},
            "step": {"exposure": 1, "infectivity_setting": 1},
        }
    },
)
def exposed():
    return _integ_exposed()


_integ_exposed = Integ(
    lambda: exposure() - infectivity_setting(), lambda: 0, "_integ_exposed"
)


@component.add(
    name="Exposure",
    units="Persons/Day",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={
        "contacts_total_per_susceptible": 1,
        "susceptibles": 1,
        "infectivity": 1,
    },
)
def exposure():
    return contacts_total_per_susceptible() * susceptibles() * infectivity()


@component.add(
    name="Hospital Admissions",
    units="Persons/Day",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={
        "incoming_demand_on_hospital": 1,
        "time_step": 1,
        "managable_hospital_inflow": 1,
    },
)
def hospital_admissions():
    return np.minimum(
        incoming_demand_on_hospital() / time_step(), managable_hospital_inflow()
    )


@component.add(
    name="Hospital Bed Capacity",
    units="Persons",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={
        "normal_hospital_beds_per_person": 1,
        "total_population": 1,
        "healthcare_capacity_multiplier": 1,
        "patients_per_unit": 1,
        "fraction_healthcare_cant_work": 1,
    },
)
def hospital_bed_capacity():
    return (
        normal_hospital_beds_per_person()
        * total_population()
        * healthcare_capacity_multiplier()
        * patients_per_unit()
        * (1 - fraction_healthcare_cant_work())
    )


@component.add(
    name="Hospital Overflow",
    units="Persons/Day",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={
        "incoming_demand_on_hospital": 1,
        "time_step": 1,
        "managable_hospital_inflow": 1,
    },
)
def hospital_overflow():
    return np.maximum(
        0, incoming_demand_on_hospital() / time_step() - managable_hospital_inflow()
    )


@component.add(
    name="ICU Admissions",
    units="Persons/Day",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"typical_worsening_rate": 1, "managable_icu_inflow": 1},
)
def icu_admissions():
    return np.minimum(typical_worsening_rate(), managable_icu_inflow())


@component.add(
    name="ICU Capacity",
    units="Persons",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={
        "normal_icu_capacity_per_person": 1,
        "total_population": 1,
        "healthcare_capacity_multiplier": 1,
        "patients_per_unit": 1,
        "fraction_healthcare_cant_work": 1,
    },
)
def icu_capacity():
    return (
        normal_icu_capacity_per_person()
        * total_population()
        * healthcare_capacity_multiplier()
        * patients_per_unit()
        * (1 - fraction_healthcare_cant_work())
    )


@component.add(
    name="ICU Overflow",
    units="Persons/Day",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"typical_worsening_rate": 1, "managable_icu_inflow": 1},
)
def icu_overflow():
    return np.maximum(0, typical_worsening_rate() - managable_icu_inflow())


@component.add(
    name="Incoming Demand on Hospital",
    units="Persons",
    comp_type="Stateful",
    comp_subtype="Integ",
    depends_on={"_integ_incoming_demand_on_hospital": 1},
    other_deps={
        "_integ_incoming_demand_on_hospital": {
            "initial": {},
            "step": {
                "disease_progression": 1,
                "iso_disease_progress_rate": 1,
                "q_disease_progress_rate": 1,
                "hospital_admissions": 1,
                "hospital_overflow": 1,
            },
        }
    },
)
def incoming_demand_on_hospital():
    return _integ_incoming_demand_on_hospital()


_integ_incoming_demand_on_hospital = Integ(
    lambda: disease_progression()
    + iso_disease_progress_rate()
    + q_disease_progress_rate()
    - hospital_admissions()
    - hospital_overflow(),
    lambda: 0,
    "_integ_incoming_demand_on_hospital",
)


@component.add(
    name="Infected Sym Hospital Overflow",
    units="Persons",
    comp_type="Stateful",
    comp_subtype="Integ",
    depends_on={"_integ_infected_sym_hospital_overflow": 1},
    other_deps={
        "_integ_infected_sym_hospital_overflow": {
            "initial": {},
            "step": {
                "hospital_overflow": 1,
                "recoveries_serioush_overflow": 1,
                "worsening_rate_overflow": 1,
            },
        }
    },
)
def infected_sym_hospital_overflow():
    return _integ_infected_sym_hospital_overflow()


_integ_infected_sym_hospital_overflow = Integ(
    lambda: hospital_overflow()
    - recoveries_serioush_overflow()
    - worsening_rate_overflow(),
    lambda: 0,
    "_integ_infected_sym_hospital_overflow",
)


@component.add(
    name="Infected Sym ICU Overflow",
    units="Persons",
    comp_type="Stateful",
    comp_subtype="Integ",
    depends_on={"_integ_infected_sym_icu_overflow": 1},
    other_deps={
        "_integ_infected_sym_icu_overflow": {
            "initial": {},
            "step": {
                "icu_overflow": 1,
                "worsening_rate_overflow": 1,
                "deaths_overflow": 1,
                "recoveries_icu_overflow": 1,
            },
        }
    },
)
def infected_sym_icu_overflow():
    return _integ_infected_sym_icu_overflow()


_integ_infected_sym_icu_overflow = Integ(
    lambda: icu_overflow()
    + worsening_rate_overflow()
    - deaths_overflow()
    - recoveries_icu_overflow(),
    lambda: 0,
    "_integ_infected_sym_icu_overflow",
)


@component.add(
    name="Infected Sym Serious Hospital",
    units="Persons",
    comp_type="Stateful",
    comp_subtype="Integ",
    depends_on={"_integ_infected_sym_serious_hospital": 1},
    other_deps={
        "_integ_infected_sym_serious_hospital": {
            "initial": {},
            "step": {
                "hospital_admissions": 1,
                "icu_admissions": 1,
                "icu_overflow": 1,
                "recoveries_serioush": 1,
            },
        }
    },
)
def infected_sym_serious_hospital():
    return _integ_infected_sym_serious_hospital()


_integ_infected_sym_serious_hospital = Integ(
    lambda: hospital_admissions()
    - icu_admissions()
    - icu_overflow()
    - recoveries_serioush(),
    lambda: 0,
    "_integ_infected_sym_serious_hospital",
)


@component.add(
    name="Infection Rate SUM",
    units="Persons/Day",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"infectivity_setting": 1, "infection_rate_new_arrivals": 1},
)
def infection_rate_sum():
    return infectivity_setting() + infection_rate_new_arrivals()


@component.add(
    name="Infectious Asymptomatics",
    units="Persons",
    comp_type="Stateful",
    comp_subtype="Integ",
    depends_on={"_integ_infectious_asymptomatics": 1},
    other_deps={
        "_integ_infectious_asymptomatics": {
            "initial": {},
            "step": {
                "infection_rate_new_arrivals": 1,
                "infectivity_setting": 1,
                "isolation_rate_asym": 1,
                "quarantining_rate_asym": 1,
                "recoveries_asym": 1,
                "symptom_setting": 1,
            },
        }
    },
)
def infectious_asymptomatics():
    return _integ_infectious_asymptomatics()


_integ_infectious_asymptomatics = Integ(
    lambda: infection_rate_new_arrivals()
    + infectivity_setting()
    - isolation_rate_asym()
    - quarantining_rate_asym()
    - recoveries_asym()
    - symptom_setting(),
    lambda: 0,
    "_integ_infectious_asymptomatics",
)


@component.add(
    name="Infectious Symptomatics Mild",
    units="Persons",
    comp_type="Stateful",
    comp_subtype="Integ",
    depends_on={"_integ_infectious_symptomatics_mild": 1},
    other_deps={
        "_integ_infectious_symptomatics_mild": {
            "initial": {},
            "step": {
                "symptom_setting": 1,
                "disease_progression": 1,
                "isolation_rate_sym": 1,
                "quarantine_rate_sym": 1,
                "recoveries_infected_sym_mild": 1,
            },
        }
    },
)
def infectious_symptomatics_mild():
    return _integ_infectious_symptomatics_mild()


_integ_infectious_symptomatics_mild = Integ(
    lambda: symptom_setting()
    - disease_progression()
    - isolation_rate_sym()
    - quarantine_rate_sym()
    - recoveries_infected_sym_mild(),
    lambda: 0,
    "_integ_infectious_symptomatics_mild",
)


@component.add(
    name="Infectivity Setting",
    units="Persons/Day",
    comp_type="Stateful",
    comp_subtype="Delay",
    depends_on={"_delay_infectivity_setting": 1},
    other_deps={
        "_delay_infectivity_setting": {
            "initial": {"exposure": 1, "delay_asymp_noninfective_to_infective": 1},
            "step": {"exposure": 1, "delay_asymp_noninfective_to_infective": 1},
        }
    },
)
def infectivity_setting():
    return _delay_infectivity_setting()


_delay_infectivity_setting = Delay(
    lambda: exposure(),
    lambda: delay_asymp_noninfective_to_infective(),
    lambda: 0 * exposure(),
    lambda: 3,
    time_step,
    "_delay_infectivity_setting",
)


@component.add(
    name="Sum Recovered",
    units="Persons",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={
        "recovered_from_hospital": 1,
        "recovered_asymtomatics": 1,
        "recovered_mild": 1,
        "recovered_isolated_symptomatic_mild": 1,
    },
)
def sum_recovered():
    return (
        recovered_from_hospital()
        + recovered_asymtomatics()
        + recovered_mild()
        + recovered_isolated_symptomatic_mild()
    )


@component.add(
    name="INTO Inf Sym ICU Overflow",
    units="Persons/Day",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"icu_overflow": 1, "worsening_rate_overflow": 1},
)
def into_inf_sym_icu_overflow():
    return icu_overflow() + worsening_rate_overflow()


@component.add(
    name="INTO Q Sym",
    units="Persons/Day",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"incubation_rate_quarantines": 1, "quarantine_rate_sym": 1},
)
def into_q_sym():
    return incubation_rate_quarantines() + quarantine_rate_sym()


@component.add(
    name="INTO Sym Serious",
    units="Persons/Day",
    limits=(0.0, np.nan),
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={
        "disease_progression": 1,
        "iso_disease_progress_rate": 1,
        "q_disease_progress_rate": 1,
    },
)
def into_sym_serious():
    return (
        disease_progression() + iso_disease_progress_rate() + q_disease_progress_rate()
    )


@component.add(
    name="Sum Recoveries reported daily",
    units="Persons/Day",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={
        "recoveries_icu": 1,
        "recoveries_icu_overflow": 1,
        "recoveries_iso_asym": 1,
        "recoveries_iso_sym": 1,
        "recoveries_serioush": 1,
        "recoveries_serioush_overflow": 1,
    },
)
def sum_recoveries_reported_daily():
    return (
        recoveries_icu()
        + recoveries_icu_overflow()
        + recoveries_iso_asym()
        + recoveries_iso_sym()
        + recoveries_serioush()
        + recoveries_serioush_overflow()
    )


@component.add(
    name="Isolation Rate Asym",
    units="Persons/Day",
    limits=(0.0, np.nan),
    comp_type="Stateful",
    comp_subtype="Delay",
    depends_on={"_delay_isolation_rate_asym": 1, "fraction_iso_asym": 1},
    other_deps={
        "_delay_isolation_rate_asym": {
            "initial": {"infection_rate_sum": 1, "delay_iso_asym": 1},
            "step": {"infection_rate_sum": 1, "delay_iso_asym": 1},
        }
    },
)
def isolation_rate_asym():
    return _delay_isolation_rate_asym() * fraction_iso_asym()


_delay_isolation_rate_asym = Delay(
    lambda: infection_rate_sum(),
    lambda: delay_iso_asym(),
    lambda: 0 * infection_rate_sum(),
    lambda: 3,
    time_step,
    "_delay_isolation_rate_asym",
)


@component.add(
    name="Isolation rate Sym",
    units="Persons/Day",
    limits=(0.0, np.nan),
    comp_type="Stateful",
    comp_subtype="Delay",
    depends_on={"fr_symptomatic_tested_isolated": 1, "_delay_isolation_rate_sym": 1},
    other_deps={
        "_delay_isolation_rate_sym": {
            "initial": {"symptom_setting": 1, "delay_iso_sym": 1},
            "step": {"symptom_setting": 1, "delay_iso_sym": 1},
        }
    },
)
def isolation_rate_sym():
    return fr_symptomatic_tested_isolated() * _delay_isolation_rate_sym()


_delay_isolation_rate_sym = Delay(
    lambda: symptom_setting(),
    lambda: delay_iso_sym(),
    lambda: 0 * symptom_setting(),
    lambda: 3,
    time_step,
    "_delay_isolation_rate_sym",
)


@component.add(
    name="Managable hospital inflow",
    units="Persons/Day",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={
        "incoming_demand_on_hospital": 2,
        "hospital_bed_capacity": 1,
        "infected_sym_serious_hospital": 1,
        "time_step": 1,
    },
)
def managable_hospital_inflow():
    return np.maximum(
        0,
        xidz(incoming_demand_on_hospital(), incoming_demand_on_hospital(), 1)
        * ((hospital_bed_capacity() - infected_sym_serious_hospital()) / time_step()),
    )


@component.add(
    name="managable ICU inflow",
    units="Persons/Day",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={
        "typical_worsening_rate": 2,
        "infected_sym_extreme_icu": 1,
        "icu_capacity": 1,
        "time_step": 1,
    },
)
def managable_icu_inflow():
    return np.maximum(
        0,
        xidz(typical_worsening_rate(), typical_worsening_rate(), 1)
        * ((icu_capacity() - infected_sym_extreme_icu()) / time_step()),
    )


@component.add(
    name="New Cases Reported",
    units="Persons/Day",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={
        "disease_progression": 1,
        "q_disease_progress_rate": 1,
        "isolation_rate_asym": 1,
        "isolation_rate_sym": 1,
    },
)
def new_cases_reported():
    return (
        disease_progression()
        + q_disease_progress_rate()
        + isolation_rate_asym()
        + isolation_rate_sym()
    )


@component.add(
    name="Quarantine Rate Sym",
    units="Persons/Day",
    comp_type="Stateful",
    comp_subtype="Delay",
    depends_on={"fraction_q_sym": 1, "_delay_quarantine_rate_sym": 1},
    other_deps={
        "_delay_quarantine_rate_sym": {
            "initial": {"symptom_setting": 1, "delay_q_sym": 1},
            "step": {"symptom_setting": 1, "delay_q_sym": 1},
        }
    },
)
def quarantine_rate_sym():
    return fraction_q_sym() * _delay_quarantine_rate_sym()


_delay_quarantine_rate_sym = Delay(
    lambda: symptom_setting(),
    lambda: delay_q_sym(),
    lambda: 0 * symptom_setting(),
    lambda: 3,
    time_step,
    "_delay_quarantine_rate_sym",
)


@component.add(
    name="Quarantined Asymptomatics",
    units="Persons",
    comp_type="Stateful",
    comp_subtype="Integ",
    depends_on={"_integ_quarantined_asymptomatics": 1},
    other_deps={
        "_integ_quarantined_asymptomatics": {
            "initial": {},
            "step": {
                "quarantining_rate_asym": 1,
                "recoveries_q_asym": 1,
                "incubation_rate_quarantines": 1,
            },
        }
    },
)
def quarantined_asymptomatics():
    return _integ_quarantined_asymptomatics()


_integ_quarantined_asymptomatics = Integ(
    lambda: quarantining_rate_asym()
    - recoveries_q_asym()
    - incubation_rate_quarantines(),
    lambda: 0,
    "_integ_quarantined_asymptomatics",
)


@component.add(
    name="Quarantined Symtomatics mild",
    units="Persons",
    comp_type="Stateful",
    comp_subtype="Integ",
    depends_on={"_integ_quarantined_symtomatics_mild": 1},
    other_deps={
        "_integ_quarantined_symtomatics_mild": {
            "initial": {},
            "step": {
                "incubation_rate_quarantines": 1,
                "quarantine_rate_sym": 1,
                "q_disease_progress_rate": 1,
                "recoveries_qsym": 1,
            },
        }
    },
)
def quarantined_symtomatics_mild():
    return _integ_quarantined_symtomatics_mild()


_integ_quarantined_symtomatics_mild = Integ(
    lambda: incubation_rate_quarantines()
    + quarantine_rate_sym()
    - q_disease_progress_rate()
    - recoveries_qsym(),
    lambda: 0,
    "_integ_quarantined_symtomatics_mild",
)


@component.add(
    name="Quarantining Rate Asym",
    units="Persons/Day",
    comp_type="Stateful",
    comp_subtype="Delay",
    depends_on={"_delay_quarantining_rate_asym": 1, "fraction_q_asym": 1},
    other_deps={
        "_delay_quarantining_rate_asym": {
            "initial": {"infection_rate_sum": 1, "delay_q_asym": 1},
            "step": {"infection_rate_sum": 1, "delay_q_asym": 1},
        }
    },
)
def quarantining_rate_asym():
    return _delay_quarantining_rate_asym() * fraction_q_asym()


_delay_quarantining_rate_asym = Delay(
    lambda: infection_rate_sum(),
    lambda: delay_q_asym(),
    lambda: 0 * infection_rate_sum(),
    lambda: 3,
    time_step,
    "_delay_quarantining_rate_asym",
)


@component.add(
    name="Sum Isolation Sym mild",
    units="Persons",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"isolated_symtomatics_mild": 1},
)
def sum_isolation_sym_mild():
    return isolated_symtomatics_mild()


@component.add(
    name="Recovered from Hospital",
    units="Persons",
    comp_type="Stateful",
    comp_subtype="Integ",
    depends_on={"_integ_recovered_from_hospital": 1},
    other_deps={
        "_integ_recovered_from_hospital": {
            "initial": {},
            "step": {
                "recoveries_icu_overflow": 1,
                "recoveries_icu": 1,
                "recoveries_serioush_overflow": 1,
                "recoveries_serioush": 1,
            },
        }
    },
)
def recovered_from_hospital():
    return _integ_recovered_from_hospital()


_integ_recovered_from_hospital = Integ(
    lambda: recoveries_icu_overflow()
    + recoveries_icu()
    + recoveries_serioush_overflow()
    + recoveries_serioush(),
    lambda: 0,
    "_integ_recovered_from_hospital",
)


@component.add(
    name="Recovered mild",
    units="Persons",
    comp_type="Stateful",
    comp_subtype="Integ",
    depends_on={"_integ_recovered_mild": 1},
    other_deps={
        "_integ_recovered_mild": {
            "initial": {},
            "step": {"recoveries_infected_sym_mild": 1, "recoveries_qsym": 1},
        }
    },
)
def recovered_mild():
    return _integ_recovered_mild()


_integ_recovered_mild = Integ(
    lambda: recoveries_infected_sym_mild() + recoveries_qsym(),
    lambda: 0,
    "_integ_recovered_mild",
)


@component.add(
    name="Sum Quarantine Sym mild",
    units="Persons",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"quarantined_symtomatics_mild": 1},
)
def sum_quarantine_sym_mild():
    return quarantined_symtomatics_mild()


@component.add(
    name="Sum Recovered via Hosp",
    units="Persons",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"recovered_from_hospital": 1},
)
def sum_recovered_via_hosp():
    return recovered_from_hospital()


@component.add(
    name="Sum Susceptibles",
    units="Persons",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"susceptibles": 1},
)
def sum_susceptibles():
    return susceptibles()


@component.add(
    name="Sum Symtomatics",
    units="Persons",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"infectious_symptomatics_mild": 1},
)
def sum_symtomatics():
    return infectious_symptomatics_mild()


@component.add(
    name="Sum Asymtomatics",
    units="Persons",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"infectious_asymptomatics": 1},
)
def sum_asymtomatics():
    return infectious_asymptomatics()


@component.add(
    name="Sum Cumulative Cases Reported",
    units="Persons",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"cumulative_cases_reported": 1},
)
def sum_cumulative_cases_reported():
    return cumulative_cases_reported()


@component.add(
    name="Sum Deaths",
    units="Persons",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"dead": 1},
)
def sum_deaths():
    return dead()


@component.add(
    name="Sum Isolation Asym",
    units="Persons",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"isolated_asymptomatics": 1},
)
def sum_isolation_asym():
    return isolated_asymptomatics()


@component.add(
    name="Sum In Hospital",
    units="Persons",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"infected_sym_serious_hospital": 1},
)
def sum_in_hospital():
    return infected_sym_serious_hospital()


@component.add(
    name="Sum in Hospital Overflow",
    units="Persons",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"infected_sym_hospital_overflow": 1},
)
def sum_in_hospital_overflow():
    return infected_sym_hospital_overflow()


@component.add(
    name="Sum in ICU",
    units="Persons",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"infected_sym_extreme_icu": 1},
)
def sum_in_icu():
    return infected_sym_extreme_icu()


@component.add(
    name="Sum in ICU Overflow",
    units="Persons",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"infected_sym_icu_overflow": 1},
)
def sum_in_icu_overflow():
    return infected_sym_icu_overflow()


@component.add(
    name="Sum of Infection Rate",
    units="Persons/Day",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"exposure": 1},
)
def sum_of_infection_rate():
    return exposure()


@component.add(
    name="Sum New Cases Reported",
    units="Persons/Day",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"new_cases_reported": 1},
)
def sum_new_cases_reported():
    return new_cases_reported()


@component.add(
    name="Susceptibles",
    units="Persons",
    comp_type="Stateful",
    comp_subtype="Integ",
    depends_on={"_integ_susceptibles": 1},
    other_deps={
        "_integ_susceptibles": {
            "initial": {"total_population_in_agegroup": 1},
            "step": {"exposure": 1},
        }
    },
)
def susceptibles():
    return _integ_susceptibles()


_integ_susceptibles = Integ(
    lambda: -exposure(), lambda: total_population_in_agegroup(), "_integ_susceptibles"
)


@component.add(
    name="Sum Quarantine Asym",
    units="Persons",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"quarantined_asymptomatics": 1},
)
def sum_quarantine_asym():
    return quarantined_asymptomatics()


@component.add(
    name="Sum Recovered after symptoms",
    units="Persons",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"recovered_from_hospital": 1, "recovered_isolated_symptomatic_mild": 1},
)
def sum_recovered_after_symptoms():
    return recovered_from_hospital() + recovered_isolated_symptomatic_mild()


@component.add(
    name="Sum Recovered Asym",
    units="Persons",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"recovered_asymtomatics": 1},
)
def sum_recovered_asym():
    return recovered_asymtomatics()


@component.add(
    name="Sum Recovered Mild",
    units="Persons",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"recovered_mild": 1},
)
def sum_recovered_mild():
    return recovered_mild()


@component.add(
    name="SumNoninfectives",
    comp_type="Auxiliary",
    comp_subtype="Normal",
    depends_on={"exposed": 1},
)
def sumnoninfectives():
    return exposed()


@component.add(
    name="Recovered Asymtomatics",
    units="Persons",
    comp_type="Stateful",
    comp_subtype="Integ",
    depends_on={"_integ_recovered_asymtomatics": 1},
    other_deps={
        "_integ_recovered_asymtomatics": {
            "initial": {},
            "step": {
                "recoveries_iso_asym": 1,
                "recoveries_q_asym": 1,
                "recoveries_asym": 1,
            },
        }
    },
)
def recovered_asymtomatics():
    return _integ_recovered_asymtomatics()


_integ_recovered_asymtomatics = Integ(
    lambda: recoveries_iso_asym() + recoveries_q_asym() + recoveries_asym(),
    lambda: 0,
    "_integ_recovered_asymtomatics",
)


@component.add(
    name="Recoveries seriousH overflow",
    units="Persons/Day",
    comp_type="Stateful",
    comp_subtype="Delay",
    depends_on={"_delay_recoveries_serioush_overflow": 1, "net_fr_requiring_icu": 1},
    other_deps={
        "_delay_recoveries_serioush_overflow": {
            "initial": {"hospital_overflow": 1, "delay_recovery_sym_serious": 1},
            "step": {"hospital_overflow": 1, "delay_recovery_sym_serious": 1},
        }
    },
)
def recoveries_serioush_overflow():
    return _delay_recoveries_serioush_overflow() * (1 - net_fr_requiring_icu())


_delay_recoveries_serioush_overflow = Delay(
    lambda: hospital_overflow(),
    lambda: delay_recovery_sym_serious(),
    lambda: 0 * hospital_overflow(),
    lambda: 3,
    time_step,
    "_delay_recoveries_serioush_overflow",
)


@component.add(
    name="Worsening Rate Overflow",
    units="Persons/Day",
    comp_type="Stateful",
    comp_subtype="Delay",
    depends_on={"_delay_worsening_rate_overflow": 1, "net_fr_requiring_icu": 1},
    other_deps={
        "_delay_worsening_rate_overflow": {
            "initial": {"hospital_overflow": 1, "delay_worsening": 1},
            "step": {"hospital_overflow": 1, "delay_worsening": 1},
        }
    },
)
def worsening_rate_overflow():
    return _delay_worsening_rate_overflow() * net_fr_requiring_icu()


_delay_worsening_rate_overflow = Delay(
    lambda: hospital_overflow(),
    lambda: delay_worsening(),
    lambda: 0 * hospital_overflow(),
    lambda: 3,
    time_step,
    "_delay_worsening_rate_overflow",
)
