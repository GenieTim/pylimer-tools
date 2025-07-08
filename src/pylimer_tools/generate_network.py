#!/usr/bin/env python
"""
Generate MC Network
===================

Use and MC procedure to generate a crosslinked polymer network.
"""

import math
import time
from typing import Union

import numpy as np
from scipy.stats import anderson, normaltest
from termcolor import colored

from pylimer_tools.calc.structure_analysis import \
    compute_crosslinker_conversion
from pylimer_tools.io.bead_spring_parameter_provider import (
    Parameters, get_parameters_for_polymer, get_supported_polymer_names)
from pylimer_tools_cpp import (DataFileWriter, MCUniverseGenerator,
                               MEHPForceBalance2, SLESolver, Universe)

###
# Configuration
###
normal_atom_type = 1
crosslink_type = 2
solvent_chain_type = 3
monofunctional_chains_type = 4

run_force_relaxation = False
###


def print_with_time(message: str) -> None:
    print(
        message
        + " at "
        + time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(time.time())),
        flush=True,
    )


def prepare_structure_generation(
    params: Parameters,
    target_p: Union[float, None] = None,
    target_wsol: Union[float, None] = None,
    n_chains_crosslinks: int = 0,
    target_f: int = 4,
    n_solvent_chains: int = 0,
    n_beads_per_solvent_chain: int = 0,
    n_chains_1: int = 1000,
    n_beads_per_chain_1: int = 0,
    n_chains_2: int = 0,
    n_beads_per_chain_2: int = 0,
    n_mono_chains: int = 0,
    n_mono_beads_per_chain: int = 0,
    n_beads_per_xlink: int = 1,
    remove_wsol: bool = False,
    functionalize_discrete: bool = False,
    disable_primary_loops: bool = False,
    disable_secondary_loops: bool = False,
    z_score_std_mult: float = 3.0,
) -> MCUniverseGenerator:
    if n_beads_per_xlink < 1:
        raise ValueError("n_beads_per_xlink must be at least 1")
    if target_p is None and target_wsol is None:
        raise ValueError("target_p or target_wsol must be provided")
    if target_p is not None and target_wsol is not None:
        raise ValueError("Only one of target_p and target_wsol should be provided")

    n_crosslinks = n_chains_crosslinks

    # the mass is in the segments / bonds
    volume = (
        (n_beads_per_chain_1) * n_chains_1
        + (n_beads_per_chain_2) * n_chains_2
        + (n_mono_beads_per_chain) * n_mono_chains
        + (n_crosslinks * n_beads_per_xlink)
        + (n_solvent_chains * n_beads_per_solvent_chain)
    ) / params.get_bead_density()
    box_l = np.cbrt(volume)

    # randomly sample chains
    universe_generator = MCUniverseGenerator(box_l, box_l, box_l)
    universe_generator.set_mean_squared_bead_distance(
        params.get("<b^2>").to(params.get("distance_units") ** 2).magnitude
    )

    if disable_primary_loops:
        universe_generator.config_primary_loop_probability(0.0)
    if disable_secondary_loops:
        universe_generator.config_secondary_loop_probability(0.0)

    if params.get_name() == "own-si":
        assert math.isclose(
            universe_generator.get_mean_bead_distance(),
            params.get("<b>").to(params.get("distance_units")).magnitude,
            rel_tol=0.1,
        )
    assert math.isclose(
        universe_generator.get_mean_squared_bead_distance(),
        params.get("<b^2>").to(params.get("distance_units") ** 2).magnitude,
        rel_tol=0.1,
    )
    universe_generator.config_nr_of_mc_steps(0)

    if n_beads_per_xlink <= 1:
        universe_generator.add_crosslinkers(
            nr_of_crosslinkers=n_crosslinks,
            crosslinker_type=crosslink_type,
            crosslinker_functionality=target_f,
        )
    else:
        if functionalize_discrete:
            functionality_per_bead = math.ceil(target_f / n_beads_per_xlink)
            functionalization_probability = (
                target_f / functionality_per_bead / (n_beads_per_xlink)
            )
            universe_generator.add_randomly_functionalized_strands(
                nr_of_strands=n_crosslinks,
                strand_length=[n_beads_per_xlink for _ in range(n_crosslinks)],
                functionalization_probability=functionalization_probability,
                crosslinker_functionality=functionality_per_bead,
            )
        else:
            universe_generator.add_randomly_functionalized_strands(
                nr_of_strands=n_crosslinks,
                strand_length=[n_beads_per_xlink for _ in range(n_crosslinks)],
                functionalization_probability=target_f / n_beads_per_xlink,
                crosslinker_functionality=1,
            )

    if n_mono_chains > 0:
        universe_generator.add_monofunctional_strands(
            n_mono_chains,
            np.repeat(n_mono_beads_per_chain, n_mono_chains).tolist(),
            monofunctional_chains_type,
        )
    universe_generator.add_strands(
        n_chains_1,
        np.repeat(n_beads_per_chain_1, n_chains_1).tolist(),
        strand_atom_type=normal_atom_type,
    )
    if n_chains_2 > 0 and n_beads_per_chain_2 > 0:
        universe_generator.add_strands(
            n_chains_2,
            np.repeat(n_beads_per_chain_2, n_chains_2).tolist(),
            strand_atom_type=normal_atom_type,
        )
    if n_solvent_chains > 0 and n_beads_per_solvent_chain > 0:
        universe_generator.add_solvent_chains(
            n_solvent_chains,
            n_beads_per_solvent_chain,
            solvent_atom_type=solvent_chain_type,
        )
    universe_generator.use_zscore_max_distance(
        z_score_std_mult,
        params.get("<b^2>").to(params.get("distance_units") ** 2).magnitude,
    )

    print_with_time(colored("Linking strands", "grey"))
    if target_p is not None:
        assert target_wsol is None
        universe_generator.link_strands_to_conversion(
            crosslinker_conversion=target_p,
        )
    else:
        assert target_wsol is not None
        universe_generator.link_strands_to_soluble_fraction(
            target_wsol,
        )

    if run_force_relaxation:
        universe_generator.relax_crosslinks()

    if remove_wsol:
        universe_generator.remove_soluble_fraction(True)

    return universe_generator


def is_active_network(
    params: Parameters,
    target_p: Union[float, None] = None,
    target_wsol: Union[float, None] = None,
    n_chains_crosslinks: int = 0,
    target_f: int = 4,
    n_solvent_chains: int = 0,
    n_beads_per_solvent_chain: int = 0,
    n_chains_1: int = 1000,
    n_beads_per_chain_1: int = 0,
    n_chains_2: int = 0,
    n_beads_per_chain_2: int = 0,
    n_mono_chains: int = 0,
    n_mono_beads_per_chain: int = 0,
    n_beads_per_xlink: int = 1,
    remove_wsol: bool = False,
    functionalize_discrete: bool = False,
    disable_primary_loops: bool = False,
    disable_secondary_loops: bool = False,
    z_score_std_mult: float = 3.0,
) -> bool:
    """
    Check whether the given network generation parameter choices may result in a network that
    has a modulus > 0.
    """
    try:
        universe_generator = prepare_structure_generation(
            params=params,
            target_p=target_p,
            target_wsol=target_wsol,
            n_chains_crosslinks=n_chains_crosslinks,
            target_f=target_f,
            n_solvent_chains=n_solvent_chains,
            n_beads_per_solvent_chain=n_beads_per_solvent_chain,
            n_chains_1=n_chains_1,
            n_beads_per_chain_1=n_beads_per_chain_1,
            n_chains_2=n_chains_2,
            n_beads_per_chain_2=n_beads_per_chain_2,
            n_mono_chains=n_mono_chains,
            n_mono_beads_per_chain=n_mono_beads_per_chain,
            n_beads_per_xlink=n_beads_per_xlink,
            remove_wsol=remove_wsol,
            functionalize_discrete=functionalize_discrete,
            disable_primary_loops=disable_primary_loops,
            disable_secondary_loops=disable_secondary_loops,
            z_score_std_mult=z_score_std_mult,
        )
    except Exception as e:
        print(f"Error during structure generation for network activity assessment: {e}")
        return False

    try:
        fb = universe_generator.get_force_balance2()
        assert isinstance(
            fb, MEHPForceBalance2
        ), "Expected ForceBalance2, got {}".format(type(fb))
        fb.run_force_relaxation()
    except Exception as e:
        print_with_time(
            colored(
                "Error during force relaxation of FB2 to investigate whether network is active: {}".format(
                    e
                ),
                "red",
            )
        )
        fb = universe_generator.get_force_balance()
        fb.run_force_relaxation()
    return fb.get_nr_of_active_springs() > 1


def generate_structure(
    params: Parameters,
    target_p: Union[float, None] = None,
    target_wsol: Union[float, None] = None,
    n_chains_crosslinks: int = 0,
    target_f: int = 4,
    n_solvent_chains: int = 0,
    n_beads_per_solvent_chain: int = 0,
    n_chains_1: int = 1000,
    n_beads_per_chain_1: int = 0,
    n_chains_2: int = 0,
    n_beads_per_chain_2: int = 0,
    n_mono_chains: int = 0,
    n_mono_beads_per_chain: int = 0,
    n_beads_per_xlink: int = 1,
    remove_wsol: bool = False,
    functionalize_discrete: bool = False,
    disable_primary_loops: bool = False,
    disable_secondary_loops: bool = False,
    z_score_std_mult: float = 3.0,
) -> Universe:
    """
    Generate a universe with the given network generation parameter choices.
    """
    universe_generator = prepare_structure_generation(
        params=params,
        target_p=target_p,
        target_wsol=target_wsol,
        n_chains_crosslinks=n_chains_crosslinks,
        target_f=target_f,
        n_solvent_chains=n_solvent_chains,
        n_beads_per_solvent_chain=n_beads_per_solvent_chain,
        n_chains_1=n_chains_1,
        n_beads_per_chain_1=n_beads_per_chain_1,
        n_chains_2=n_chains_2,
        n_beads_per_chain_2=n_beads_per_chain_2,
        n_mono_chains=n_mono_chains,
        n_mono_beads_per_chain=n_mono_beads_per_chain,
        n_beads_per_xlink=n_beads_per_xlink,
        remove_wsol=remove_wsol,
        functionalize_discrete=functionalize_discrete,
        disable_primary_loops=disable_primary_loops,
        disable_secondary_loops=disable_secondary_loops,
        z_score_std_mult=z_score_std_mult,
    )
    print_with_time(colored("Sampling beads", "grey"))
    universe = universe_generator.get_universe()

    # some final additional info
    universe.set_masses(
        {
            crosslink_type: params.get("Mw").magnitude,
            solvent_chain_type: params.get("Mw").magnitude,
            normal_atom_type: params.get("Mw").magnitude,
            monofunctional_chains_type: params.get("Mw").magnitude,
        }
    )

    # validation
    if target_p is not None and n_beads_per_xlink <= 1 and not remove_wsol:
        comp_p = compute_crosslinker_conversion(
            universe, crosslinker_type=crosslink_type, f=target_f
        )
        assert math.isclose(
            target_p,
            comp_p,
            rel_tol=0.05,
        ), "Wrong crosslinker conversion computed: expected {}, got {}".format(
            target_p, comp_p
        )
    # TODO: also check for w_sol instead

    bond_lengths = universe.compute_bond_lengths()
    print(
        "Bond lengths, squared mean: {}, mean: {} (median: {}, max: {}, min: {})".format(
            np.mean(np.square(bond_lengths)),
            np.mean(bond_lengths),
            np.median(bond_lengths),
            max(bond_lengths),
            min(bond_lengths),
        )
    )

    # check if bond vectors are normally distributed
    # bond_dir_distances = [a for bs in universe.compute_bond_vectors() for a in bs]
    # assert anderson(
    #     bond_dir_distances
    # ).fit_result.success, "Bond vectors not normally distributed according to Anderson"
    # b_is_normal = normaltest(bond_dir_distances)
    # assert (
    #     b_is_normal.pvalue > 0.05
    # ), "Bond vectors not normally distributed according to D'Agostino and Pearson"

    # assert math.isclose(params.get("<b>").magnitude, np.mean(bond_lengths), rel_tol=0.2)
    # assert math.isclose(
    #     params.get("<b^2>").magnitude, np.mean(np.square(bond_lengths)), rel_tol=0.2
    # )

    # end_to_end_distances = universe.compute_end_to_end_distances(
    #     crosslinker_type=crosslink_type, derive_image_flags=True
    # )
    # if n_mono_chains == 0:
    #     assert math.isclose(
    #         np.mean(np.square(end_to_end_distances)),
    #         n_beads_per_chain * params.get("<b^2>").magnitude,
    #         rel_tol=(
    #             0.05
    #             if (
    #                 (target_p is not None and target_p < 0.98)
    #                 or (target_wsol is not None and target_wsol < 0.05)
    #             )
    #             else 0.1
    #         ),
    #     )
    assert math.isclose(
        universe.get_nr_of_atoms() / universe.get_volume(),
        params.get_bead_density(),
    )

    return universe


def parse_args():
    import argparse

    parser = argparse.ArgumentParser(
        description="Generate a crosslinked polymer network."
    )
    parser.add_argument(
        "--polymer_name",
        type=str,
        default="PDMS",
        help="Name of the polymer to generate the network for",
        choices=get_supported_polymer_names(),
    )
    # Mutually exclusive group for target_p and target_wsol
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--target_p", type=float, help="Target crosslinker conversion (p)"
    )
    group.add_argument(
        "--target_wsol", type=float, help="Target soluble fraction (w_sol)"
    )
    parser.add_argument(
        "--n_chains_crosslinks", type=int, default=0, help="Number of crosslinks"
    )
    parser.add_argument(
        "--target_f", type=int, default=4, help="Crosslinker functionality"
    )
    parser.add_argument(
        "--n_solvent_chains", type=int, default=0, help="Number of solvent chains"
    )
    parser.add_argument(
        "--n_beads_per_solvent_chain",
        type=int,
        default=0,
        help="Beads per solvent chain",
    )
    parser.add_argument(
        "--n_chains_1", type=int, default=1000, help="Number of chains of type 1"
    )
    parser.add_argument(
        "--n_beads_per_chain_1",
        type=int,
        default=0,
        help="Beads per chain of type 1",
    )
    parser.add_argument(
        "--n_chains_2", type=int, default=0, help="Number of chains of type 2"
    )
    parser.add_argument(
        "--n_beads_per_chain_2",
        type=int,
        default=0,
        help="Beads per chain of type 2",
    )
    parser.add_argument(
        "--n_mono_chains", type=int, default=0, help="Number of monofunctional chains"
    )
    parser.add_argument(
        "--n_mono_beads_per_chain",
        type=int,
        default=0,
        help="Beads per monofunctional chain",
    )
    parser.add_argument(
        "--n_beads_per_xlink",
        type=int,
        default=1,
        help="Beads per crosslinker (xlink)",
    )
    parser.add_argument(
        "--remove_wsol",
        action="store_true",
        help="Remove soluble fraction after generation",
    )
    parser.add_argument(
        "--functionalize_discrete",
        action="store_false",
        help=(
            "Functionalize crosslinks discrete, i.e., each crosslinker bead is "
            "functionalized with a single functionality"
        ),
    )
    parser.add_argument(
        "--disable_primary_loops",
        action="store_true",
        help="Disable primary loops in the network",
    )
    parser.add_argument(
        "--disable_secondary_loops",
        action="store_true",
        help="Disable secondary loops in the network",
    )
    parser.add_argument(
        "--z_score_std_mult",
        type=float,
        default=3.0,
        help="Z-score standard deviation multiplier for max distance",
    )
    parser.add_argument(
        "--target_file",
        type=str,
        help="Output file for the generated network",
        required=True,
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    params = get_parameters_for_polymer(args.polymer_name)

    universe = generate_structure(
        params=params,
        target_p=args.target_p,
        target_wsol=args.target_wsol,
        n_chains_crosslinks=args.n_chains_crosslinks,
        target_f=args.target_f,
        n_solvent_chains=args.n_solvent_chains,
        n_beads_per_solvent_chain=args.n_beads_per_solvent_chain,
        n_chains_1=args.n_chains_1,
        n_beads_per_chain_1=args.n_beads_per_chain_1,
        n_chains_2=args.n_chains_2,
        n_beads_per_chain_2=args.n_beads_per_chain_2,
        n_mono_chains=args.n_mono_chains,
        n_mono_beads_per_chain=args.n_mono_beads_per_chain,
        n_beads_per_xlink=args.n_beads_per_xlink,
        remove_wsol=args.remove_wsol,
        functionalize_discrete=args.functionalize_discrete,
        disable_primary_loops=args.disable_primary_loops,
        disable_secondary_loops=args.disable_secondary_loops,
        z_score_std_mult=args.z_score_std_mult,
    )

    print_with_time(colored("Structure generation completed", "green"))

    writer = DataFileWriter(universe)
    writer.config_include_angles(True)
    writer.config_crosslinker_type(crosslink_type)
    writer.config_molecule_idx_for_swap(False)
    writer.config_move_into_box(True)
    writer.config_attempt_image_reset(True)
    writer.write_to_file(args.target_file)
    print_with_time(colored(f"Network written to {args.target_file}", "green"))
