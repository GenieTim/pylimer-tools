"""
Multi-seed simulation runner for Study 2: Finite Extensibility & Large-Strain Stiffening.

Compares linear Hookean (Gaussian) and non-Gaussian inverse-Langevin constitutive laws
across unimodal and bimodal PDMS-parameterized networks under uniaxial tension and
equibiaxial extension, exposes lock-up divergence and solver tolerance sensitivity,
implements a custom C1-regularized/breakable Python ``MEHPForceEvaluator``, and verifies
the fixed contour-partition constraint (``SlipMode.FIXED``) of
``NonlinearOptimizationSolver``.

This module is prefixed with an underscore so Sphinx-Gallery ignores it when building
documentation.
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
import copy
from enum import Enum
import math
from pathlib import Path
import sys
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

_EXAMPLES_DIR = str(Path(__file__).resolve().parents[1])
if _EXAMPLES_DIR not in sys.path:
    sys.path.insert(0, _EXAMPLES_DIR)

from _example_cache import (  # noqa: E402
    compute_ensemble_statistics,
    configure_process_threading,
    resolve_max_workers,
    worker_process_initializer,
)
from pylimer_tools.io.bead_spring_parameter_provider import (  # noqa: E402
    ParameterType,
    get_parameters_for_polymer,
)
import pylimer_tools_cpp as pt  # noqa: E402
from pylimer_tools_cpp import (  # noqa: E402
    Box,
    EntanglementMode,
    ExitReason,
    HookeanForceEvaluator,
    MCUniverseGenerator,
    MEHPForceEvaluator,
    MEHPLinearSolverConfig,
    MEHPNetwork,
    MEHPNonlinearSolverConfig,
    MEHPSimulator,
    NLoptAlgorithm,
    NonGaussianSpringForceEvaluator,
    SlipMode,
)


class NetworkTopology(Enum):
    """Precursor strand-length topology for end-linked networks."""

    UNIMODAL = "unimodal"
    BIMODAL = "bimodal"


class DeformationMode(Enum):
    """Volume-preserving macroscopic deformation protocol."""

    UNIAXIAL = "uniaxial"
    EQUIBIAXIAL = "equibiaxial"


class RegularizationMode(Enum):
    """Constitutive extension of the Padé inverse-Langevin law above x_c."""

    BOUNDED_SATURATION = "bounded_saturation"
    BOND_SCISSION = "bond_scission"


class SolverPreset(Enum):
    """Solver algorithm and tolerance configurations for lock-up benchmarking."""

    LBFGS_TIGHT = "lbfgs_tight"
    LBFGS_LOOSE = "lbfgs_loose"
    MMA_TIGHT = "mma_tight"
    MMA_LOOSE = "mma_loose"


SOLVER_PRESET_SPECS: Dict[SolverPreset, Tuple[NLoptAlgorithm, float, float, int]] = {
    SolverPreset.LBFGS_TIGHT: (NLoptAlgorithm.LD_LBFGS, 1e-12, 1e-9, 20000),
    SolverPreset.LBFGS_LOOSE: (NLoptAlgorithm.LD_LBFGS, 1e-5, 1e-4, 200),
    SolverPreset.MMA_TIGHT: (NLoptAlgorithm.LD_MMA, 1e-12, 1e-9, 20000),
    SolverPreset.MMA_LOOSE: (NLoptAlgorithm.LD_MMA, 1e-5, 1e-4, 200),
}


def pade_inverse_langevin(x: np.ndarray | float) -> np.ndarray:
    """Evaluate Cohen's Padé approximation to the inverse Langevin function."""
    x_arr = np.asarray(x, dtype=np.float64)
    x2 = x_arr * x_arr
    return (3.0 * x_arr) / ((1.0 - x2) * (1.0 + 0.5 * x2))


def pade_inverse_langevin_derivative(x: np.ndarray | float) -> np.ndarray:
    """Evaluate d/dx of Cohen's Padé inverse-Langevin approximation."""
    x_arr = np.asarray(x, dtype=np.float64)
    x2 = x_arr * x_arr
    denom = (1.0 - x2) * (1.0 + 0.5 * x2)
    return (3.0 * (1.0 + 0.5 * x2 + 1.5 * x2 * x2)) / (denom * denom)


def pade_antiderivative_energy(x: np.ndarray | float) -> np.ndarray:
    """Exact antiderivative int_0^x L_Pade^{-1}(t) dt = ln((1 + 0.5 x^2) / (1 - x^2))."""
    x_arr = np.asarray(x, dtype=np.float64)
    x2 = x_arr * x_arr
    return np.log((1.0 + 0.5 * x2) / (1.0 - x2))


class SoftenedLangevinForceEvaluator(MEHPForceEvaluator):
    """Custom C^1-force / C^2-energy regularized inverse-Langevin force evaluator.

    Below ``x_c`` (e.g. 0.88), the radial force matches
    :class:`~pylimer_tools_cpp.NonGaussianSpringForceEvaluator` exactly using Cohen's
    Padé approximation, paired with its exact analytical antiderivative so energy and
    gradient are variationally consistent to machine precision. Above ``x_c``, the
    force transitions smoothly to either:

    - :attr:`RegularizationMode.BOUNDED_SATURATION`: convex bounded plateau
      ``f(x) = f_c + df_max * tanh((f'_c / df_max) * (x - x_c))`` representing
      finite-extensibility / backbone-yield regularization.
    - :attr:`RegularizationMode.BOND_SCISSION`: smooth exponential bond-scission
      decay ``f(x) = f_c * (1 + alpha * (x - x_c)) * exp(-(x - x_c) / w)``
      representing sacrificial strand rupture and load transfer.
    """

    l: float
    one_over_l: float
    x_c: float
    regularization_mode: RegularizationMode
    saturation_factor: float
    scission_width: float
    e_c: float
    f_c: float
    fp_c: float
    alpha: float

    def __init__(
        self,
        l: float = 1.0,
        x_c: float = 0.88,
        regularization_mode: RegularizationMode = RegularizationMode.BOUNDED_SATURATION,
        saturation_factor: float = 0.5,
        scission_width: float = 0.06,
    ) -> None:
        super().__init__()
        if not math.isfinite(l) or l <= 0.0:
            raise ValueError("Kuhn segment length l must be finite and positive.")
        if not (0.1 < x_c < 0.99):
            raise ValueError(
                "Critical fractional extension x_c must be in (0.1, 0.99)."
            )
        if not math.isfinite(saturation_factor) or saturation_factor <= 0.0:
            raise ValueError("saturation_factor must be finite and positive.")
        if not math.isfinite(scission_width) or scission_width <= 0.0:
            raise ValueError("scission_width must be finite and positive.")

        self.l = float(l)
        self.one_over_l = 1.0 / self.l
        self.x_c = float(x_c)
        self.regularization_mode = regularization_mode
        self.saturation_factor = float(saturation_factor)
        self.scission_width = float(scission_width)

        self.e_c = float(pade_antiderivative_energy(self.x_c))
        self.f_c = float(pade_inverse_langevin(self.x_c))
        self.fp_c = float(pade_inverse_langevin_derivative(self.x_c))
        self.alpha = self.fp_c / self.f_c + 1.0 / self.scission_width

    def is_linear(self) -> bool:
        return False

    def get_kappa(self) -> float:
        return 3.0 / (self.l * self.l)

    def evaluate_dimensionless_curve(
        self, x: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Return dimensionless energy per segment e(x) and force f(x) = l * F_r(x)."""
        x_arr = np.asarray(x, dtype=np.float64)
        e = np.zeros_like(x_arr, dtype=np.float64)
        f = np.zeros_like(x_arr, dtype=np.float64)

        lo = (x_arr > 0.0) & (x_arr <= self.x_c)
        hi = x_arr > self.x_c

        if np.any(lo):
            xl = x_arr[lo]
            e[lo] = pade_antiderivative_energy(xl)
            f[lo] = pade_inverse_langevin(xl)

        if np.any(hi):
            u = x_arr[hi] - self.x_c
            if self.regularization_mode == RegularizationMode.BOUNDED_SATURATION:
                df_max = self.saturation_factor * self.f_c
                t = (self.fp_c / df_max) * u
                f[hi] = self.f_c + df_max * np.tanh(t)
                # log(cosh(t)) = t + log((1 + exp(-2*t)) / 2), stable for all t >= 0
                log_cosh = t + np.log(0.5 * (1.0 + np.exp(-2.0 * t)))
                e[hi] = (
                    self.e_c + self.f_c * u + ((df_max * df_max) / self.fp_c) * log_cosh
                )
            elif self.regularization_mode == RegularizationMode.BOND_SCISSION:
                exp_u = np.exp(-u / self.scission_width)
                f[hi] = self.f_c * (1.0 + self.alpha * u) * exp_u
                e[hi] = self.e_c + self.f_c * self.scission_width * (
                    (1.0 + self.alpha * self.scission_width)
                    - (1.0 + self.alpha * self.scission_width + self.alpha * u) * exp_u
                )
            else:
                raise ValueError(
                    f"Unsupported regularization mode: {self.regularization_mode}"
                )

        return e, f

    def evaluate_spring_energy(
        self, spring_vec: np.ndarray, contour_length: float
    ) -> float:
        vec = np.asarray(spring_vec, dtype=np.float64)
        r = float(np.linalg.norm(vec))
        if r <= 0.0:
            return 0.0
        x = (r * self.one_over_l) / contour_length
        e, _ = self.evaluate_dimensionless_curve(np.array([x], dtype=np.float64))
        return float(contour_length * e[0])

    def evaluate_spring_force(
        self, spring_vec: np.ndarray, contour_length: float
    ) -> np.ndarray:
        vec = np.asarray(spring_vec, dtype=np.float64)
        r = float(np.linalg.norm(vec))
        if r <= 0.0:
            return np.zeros(3, dtype=np.float64)
        x = (r * self.one_over_l) / contour_length
        _, f = self.evaluate_dimensionless_curve(np.array([x], dtype=np.float64))
        return ((self.one_over_l * f[0]) / r) * vec

    def evaluate_network_energy_and_gradient(
        self,
        net: MEHPNetwork,
        spring_vectors: np.ndarray,
        is_2d: bool,
        compute_gradient: bool,
    ) -> Tuple[float, Optional[np.ndarray]]:
        nr_springs = net.nr_of_springs
        vecs = (
            np.asarray(spring_vectors, dtype=np.float64).reshape(nr_springs, 3).copy()
        )
        if is_2d:
            vecs[:, 2] = 0.0
        r = np.linalg.norm(vecs, axis=1)
        ns = np.asarray(net.spring_contour_length, dtype=np.float64)
        x = (r * self.one_over_l) / ns

        e_dim, f_dim = self.evaluate_dimensionless_curve(x)
        energy = float(np.sum(ns * e_dim))
        if not compute_gradient:
            return (energy, None)

        f_over_r = np.zeros(nr_springs, dtype=np.float64)
        pos = r > 0.0
        f_over_r[pos] = (self.one_over_l * f_dim[pos]) / r[pos]
        forces = f_over_r[:, None] * vecs

        grad = np.zeros((net.nr_of_links, 3), dtype=np.float64)
        idx_a = np.asarray(net.spring_index_a, dtype=np.int64)
        idx_b = np.asarray(net.spring_index_b, dtype=np.int64)
        np.add.at(grad, idx_a, -forces)
        np.add.at(grad, idx_b, forces)
        if is_2d:
            grad[:, 2] = 0.0
        return (energy, grad.reshape(-1))


def _get_pdms_kuhn_constants() -> Tuple[float, float, float, float]:
    """Return (l_kuhn_nm, r02_nm2, bead_density_nm3, gamma_conv_mpa) for PDMS."""
    params = get_parameters_for_polymer("PDMS", parameter_type=ParameterType.KUHN)
    l_kuhn = float(params.get("<b>").to("nm").magnitude)
    r02 = float(params.get("<b^2>").to("nm^2").magnitude)
    rho_beads = float(params.get_bead_density())
    gamma_conv = float(params.get_gamma_conversion_factor().to("MPa").magnitude)
    return l_kuhn, r02, rho_beads, gamma_conv


def _strand_lengths_for_topology(
    topology: NetworkTopology, n_strands: int = 600
) -> List[int]:
    """Build strand bead-count list with matched mean length <N> ~ 25."""
    if topology == NetworkTopology.UNIMODAL:
        return [25] * n_strands
    if topology == NetworkTopology.BIMODAL:
        n_short = int(round(0.72 * n_strands))
        n_long = n_strands - n_short
        return [9] * n_short + [66] * n_long
    raise ValueError(f"Unsupported topology: {topology}")


def _build_pdms_network(
    topology: NetworkTopology,
    seed: int,
    n_strands: int = 600,
    conversion: float = 0.95,
) -> Tuple[MEHPNetwork, float, float]:
    """Generate a 3D f=4 phantom PDMS network via MCUniverseGenerator."""
    l_kuhn, r02, rho_beads, _ = _get_pdms_kuhn_constants()
    strand_lengths = _strand_lengths_for_topology(topology, n_strands=n_strands)
    n_xlinks = n_strands // 2
    total_beads = sum(strand_lengths) + n_xlinks
    vol = float(total_beads) / rho_beads
    box_l = vol ** (1.0 / 3.0)

    gen = MCUniverseGenerator(box_l, box_l, box_l)
    gen.set_seed(seed)
    gen.set_bead_distance(l_kuhn)
    gen.set_mean_squared_bead_distance(r02)
    gen.config_nr_of_mc_steps(0)
    gen.add_crosslinkers(n_xlinks, 4, 2, True)
    gen.add_strands(n_strands, strand_lengths, 1)
    gen.link_strands_to_conversion(conversion)
    net = gen.get_mehp_network(entanglement_mode=EntanglementMode.NONE)
    return net, box_l, vol


def _make_deformed_box(
    box_l: float, lam: float, mode: DeformationMode, axis: int
) -> Box:
    """Construct volume-preserving Box for uniaxial or equibiaxial deformation."""
    if mode == DeformationMode.UNIAXIAL:
        dims = [box_l / math.sqrt(lam)] * 3
        dims[axis] = box_l * lam
    elif mode == DeformationMode.EQUIBIAXIAL:
        dims = [box_l * lam] * 3
        dims[axis] = box_l / (lam * lam)
    else:
        raise ValueError(f"Unsupported deformation mode: {mode}")
    return Box(dims[0], dims[1], dims[2])


def _extract_stress_metrics(
    stress_tensor_mpa: np.ndarray,
    lam: float,
    mode: DeformationMode,
    axis: int,
) -> Tuple[float, float, float]:
    """Compute (sigma_true, sigma_eng, f_star) in MPa from virial stress tensor."""
    a1 = (axis + 1) % 3
    a2 = (axis + 2) % 3
    if mode == DeformationMode.UNIAXIAL:
        sigma_true = float(
            stress_tensor_mpa[axis, axis]
            - 0.5 * (stress_tensor_mpa[a1, a1] + stress_tensor_mpa[a2, a2])
        )
        sigma_eng = sigma_true / lam
        denom = lam - (1.0 / (lam * lam))
    elif mode == DeformationMode.EQUIBIAXIAL:
        sigma_true = float(
            0.5 * (stress_tensor_mpa[a1, a1] + stress_tensor_mpa[a2, a2])
            - stress_tensor_mpa[axis, axis]
        )
        sigma_eng = sigma_true / lam
        denom = lam - (lam ** (-5))
    else:
        raise ValueError(f"Unsupported deformation mode: {mode}")
    f_star = sigma_eng / denom
    return sigma_true, sigma_eng, f_star


def _make_nonlinear_config(preset: SolverPreset) -> MEHPNonlinearSolverConfig:
    """Create MEHPNonlinearSolverConfig from a SolverPreset enum."""
    alg, xtol, ftol, max_steps = SOLVER_PRESET_SPECS[preset]
    cfg = MEHPNonlinearSolverConfig()
    cfg.algorithm = alg
    cfg.xtol = xtol
    cfg.ftol = ftol
    cfg.max_nr_of_steps = max_steps
    cfg.slip_mode = SlipMode.FIXED
    return cfg


def _run_stepwise_curve(
    sim_eq: MEHPSimulator,
    box_l: float,
    lambdas: Sequence[float],
    mode: DeformationMode,
    l_kuhn: float,
    gamma_conv: float,
    solver_cfg: Any,
    average_axes: bool = True,
) -> Dict[str, np.ndarray]:
    """Deform ``sim_eq`` stepwise along ``lambdas`` and collect mechanical metrics."""
    n_pts = len(lambdas)
    axes = (0, 1, 2) if average_axes else (0,)
    n_axes = len(axes)

    sig_eng = np.zeros((n_pts, n_axes), dtype=np.float64)
    f_star = np.zeros((n_pts, n_axes), dtype=np.float64)
    max_x_all = np.zeros((n_pts, n_axes), dtype=np.float64)
    max_x_short = np.zeros((n_pts, n_axes), dtype=np.float64)
    max_x_long = np.zeros((n_pts, n_axes), dtype=np.float64)
    residuals = np.zeros((n_pts, n_axes), dtype=np.float64)
    iterations = np.zeros((n_pts, n_axes), dtype=np.float64)
    exit_codes = np.zeros((n_pts, n_axes), dtype=np.int64)

    contour_lengths = np.asarray(sim_eq.get_spring_contour_length(), dtype=np.float64)
    short_mask = contour_lengths < 20.0
    long_mask = contour_lengths > 40.0

    for ax_idx, axis in enumerate(axes):
        sim_cur = copy.copy(sim_eq)
        for li, lam in enumerate(lambdas):
            sim_cur.deform_to(_make_deformed_box(box_l, float(lam), mode, axis))
            res = sim_cur.run(solver_cfg)

            stress_mpa = (
                np.asarray(sim_cur.get_stress_tensor(), dtype=np.float64) * gamma_conv
            )
            _, se, fs = _extract_stress_metrics(stress_mpa, float(lam), mode, axis)
            sig_eng[li, ax_idx] = se
            f_star[li, ax_idx] = fs

            x_vals = np.asarray(sim_cur.get_spring_lengths(), dtype=np.float64) / (
                contour_lengths * l_kuhn
            )
            max_x_all[li, ax_idx] = float(np.max(x_vals))
            max_x_short[li, ax_idx] = (
                float(np.max(x_vals[short_mask]))
                if np.any(short_mask)
                else float(np.max(x_vals))
            )
            max_x_long[li, ax_idx] = (
                float(np.max(x_vals[long_mask]))
                if np.any(long_mask)
                else float(np.max(x_vals))
            )
            residuals[li, ax_idx] = float(res.residual)
            iterations[li, ax_idx] = float(res.nr_of_iterations)
            exit_codes[li, ax_idx] = int(res.exit_reason)

    return {
        "sigma_eng": np.mean(sig_eng, axis=1),
        "f_star": np.mean(f_star, axis=1),
        "max_x": np.max(max_x_all, axis=1),
        "max_x_short": np.max(max_x_short, axis=1),
        "max_x_long": np.max(max_x_long, axis=1),
        "residual": np.mean(residuals, axis=1),
        "iterations": np.mean(iterations, axis=1),
        "exit_code": exit_codes[:, 0],
    }


def _simulate_part1_seed(
    seed_idx: int,
    seed: int,
    uni_lams_uni: Sequence[float],
    uni_lams_bim: Sequence[float],
    biax_lams_uni: Sequence[float],
    biax_lams_bim: Sequence[float],
) -> Dict[str, Any]:
    """Run Part 1 (Hookean vs NonGaussian, Unimodal vs Bimodal) for one MC seed."""
    l_kuhn, r02, _, gamma_conv = _get_pdms_kuhn_constants()
    kappa_hook = 3.0 / (l_kuhn * l_kuhn)
    linear_cfg = MEHPLinearSolverConfig()
    nonlinear_cfg = _make_nonlinear_config(SolverPreset.LBFGS_TIGHT)

    out: Dict[str, Any] = {"seed_idx": seed_idx}
    specs = (
        ("uni", NetworkTopology.UNIMODAL, uni_lams_uni, biax_lams_uni),
        ("bim", NetworkTopology.BIMODAL, uni_lams_bim, biax_lams_bim),
    )
    for tag, topo, u_lams, b_lams in specs:
        net, box_l, vol = _build_pdms_network(topo, seed=seed)

        sim_hook = MEHPSimulator(net, HookeanForceEvaluator(kappa_hook), False)
        sim_hook.run(linear_cfg)
        g_eq_mpa = float(gamma_conv * np.sum(sim_hook.get_gamma_factors(r02)) / vol)
        w_active = float(sim_hook.get_active_weight_fraction())

        sim_ng = MEHPSimulator(net, NonGaussianSpringForceEvaluator(1.0, l_kuhn), False)
        sim_ng.run(nonlinear_cfg)

        out[f"{tag}_g_eq"] = g_eq_mpa
        out[f"{tag}_active"] = w_active

        for mode_tag, mode, lams in (
            ("uni", DeformationMode.UNIAXIAL, u_lams),
            ("biax", DeformationMode.EQUIBIAXIAL, b_lams),
        ):
            res_hook = _run_stepwise_curve(
                sim_hook,
                box_l,
                lams,
                mode,
                l_kuhn,
                gamma_conv,
                linear_cfg,
                average_axes=True,
            )
            res_ng = _run_stepwise_curve(
                sim_ng,
                box_l,
                lams,
                mode,
                l_kuhn,
                gamma_conv,
                nonlinear_cfg,
                average_axes=True,
            )
            out[f"{tag}_{mode_tag}_hook_sig"] = res_hook["sigma_eng"]
            out[f"{tag}_{mode_tag}_hook_fstar"] = res_hook["f_star"]
            out[f"{tag}_{mode_tag}_ng_sig"] = res_ng["sigma_eng"]
            out[f"{tag}_{mode_tag}_ng_fstar"] = res_ng["f_star"]
            out[f"{tag}_{mode_tag}_ng_max_x"] = res_ng["max_x"]

    return out


def _simulate_part2_reg_seed(
    seed_idx: int,
    seed: int,
    reg_lams: Sequence[float],
    x_c: float,
) -> Dict[str, Any]:
    """Run Part 2b (Hookean, Rigid Pade, Bounded Saturation, Bond Scission) for one seed."""
    l_kuhn, _, _, gamma_conv = _get_pdms_kuhn_constants()
    kappa_hook = 3.0 / (l_kuhn * l_kuhn)
    linear_cfg = MEHPLinearSolverConfig()
    nonlinear_cfg = _make_nonlinear_config(SolverPreset.LBFGS_TIGHT)

    net, box_l, _ = _build_pdms_network(NetworkTopology.BIMODAL, seed=seed)

    sim_hook = MEHPSimulator(net, HookeanForceEvaluator(kappa_hook), False)
    sim_hook.run(linear_cfg)
    res_hook = _run_stepwise_curve(
        sim_hook,
        box_l,
        reg_lams,
        DeformationMode.UNIAXIAL,
        l_kuhn,
        gamma_conv,
        linear_cfg,
        average_axes=True,
    )

    sim_pade = MEHPSimulator(net, NonGaussianSpringForceEvaluator(1.0, l_kuhn), False)
    sim_pade.run(nonlinear_cfg)
    res_pade = _run_stepwise_curve(
        sim_pade,
        box_l,
        reg_lams,
        DeformationMode.UNIAXIAL,
        l_kuhn,
        gamma_conv,
        nonlinear_cfg,
        average_axes=True,
    )

    ev_sat = SoftenedLangevinForceEvaluator(
        l=l_kuhn, x_c=x_c, regularization_mode=RegularizationMode.BOUNDED_SATURATION
    )
    sim_sat = MEHPSimulator(net, ev_sat, False)
    sim_sat.run(nonlinear_cfg)
    res_sat = _run_stepwise_curve(
        sim_sat,
        box_l,
        reg_lams,
        DeformationMode.UNIAXIAL,
        l_kuhn,
        gamma_conv,
        nonlinear_cfg,
        average_axes=True,
    )

    ev_sci = SoftenedLangevinForceEvaluator(
        l=l_kuhn, x_c=x_c, regularization_mode=RegularizationMode.BOND_SCISSION
    )
    sim_sci = MEHPSimulator(net, ev_sci, False)
    sim_sci.run(nonlinear_cfg)
    res_sci = _run_stepwise_curve(
        sim_sci,
        box_l,
        reg_lams,
        DeformationMode.UNIAXIAL,
        l_kuhn,
        gamma_conv,
        nonlinear_cfg,
        average_axes=True,
    )

    return {
        "seed_idx": seed_idx,
        "hook_sig": res_hook["sigma_eng"],
        "pade_sig": res_pade["sigma_eng"],
        "pade_max_x": res_pade["max_x"],
        "pade_exit": res_pade["exit_code"],
        "sat_sig": res_sat["sigma_eng"],
        "sat_max_x": res_sat["max_x"],
        "sat_exit": res_sat["exit_code"],
        "sci_sig": res_sci["sigma_eng"],
        "sci_max_x": res_sci["max_x"],
        "sci_exit": res_sci["exit_code"],
    }


def _run_solver_tolerance_benchmark(
    solver_lams: Sequence[float],
    seed: int = 101,
) -> Dict[str, Any]:
    """Compare LD_LBFGS vs LD_MMA and tight vs loose tolerances near bimodal lock-up."""
    l_kuhn, _, _, gamma_conv = _get_pdms_kuhn_constants()
    net, box_l, _ = _build_pdms_network(NetworkTopology.BIMODAL, seed=seed)

    results: Dict[str, Any] = {}
    for preset in SolverPreset:
        cfg = _make_nonlinear_config(preset)
        sim = MEHPSimulator(net, NonGaussianSpringForceEvaluator(1.0, l_kuhn), False)
        sim.run(cfg)
        curve = _run_stepwise_curve(
            sim,
            box_l,
            solver_lams,
            DeformationMode.UNIAXIAL,
            l_kuhn,
            gamma_conv,
            cfg,
            average_axes=False,
        )
        prefix = f"solver_{preset.value}"
        results[f"{prefix}_sig"] = curve["sigma_eng"]
        results[f"{prefix}_max_x"] = curve["max_x"]
        results[f"{prefix}_max_x_short"] = curve["max_x_short"]
        results[f"{prefix}_max_x_long"] = curve["max_x_long"]
        results[f"{prefix}_res"] = curve["residual"]
        results[f"{prefix}_iters"] = curve["iterations"]
        results[f"{prefix}_exit"] = curve["exit_code"]
    return results


def _verify_slip_mode_constraint() -> Dict[str, Any]:
    """Verify NonlinearOptimizationSolver supports SlipMode.FIXED and rejects SLIDING."""
    l_kuhn, r02, _, _ = _get_pdms_kuhn_constants()
    box_l = 18.0
    gen = MCUniverseGenerator(box_l, box_l, box_l)
    gen.set_seed(77)
    gen.set_bead_distance(l_kuhn)
    gen.set_mean_squared_bead_distance(r02)
    gen.config_nr_of_mc_steps(0)
    gen.add_crosslinkers(50, 4, 2, True)
    gen.add_strands(100, [25] * 100, 1)
    gen.link_strands_to_conversion(0.95)
    entanglements = gen.sample_entanglements(
        nr_of_entanglements_to_sample=60,
        upper_cutoff=2.6,
        seed="slip_check",
    )
    sim = gen.get_mehp_simulator(
        entanglements, entanglement_mode=EntanglementMode.MERGED_LINKS
    )
    sim.set_force_evaluator(NonGaussianSpringForceEvaluator(1.0, l_kuhn))

    cfg_default = MEHPNonlinearSolverConfig()
    default_is_fixed = cfg_default.slip_mode == SlipMode.FIXED
    res_fixed = sim.run(cfg_default)
    fixed_converged = res_fixed.exit_reason in (
        ExitReason.F_TOLERANCE,
        ExitReason.X_TOLERANCE,
    )

    cfg_sliding = MEHPNonlinearSolverConfig()
    cfg_sliding.slip_mode = SlipMode.SLIDING
    sliding_rejected = False
    sliding_error_msg = ""
    try:
        sim.run(cfg_sliding)
    except ValueError as exc:
        sliding_rejected = True
        sliding_error_msg = str(exc)

    return {
        "slip_default_is_fixed": np.array(default_is_fixed),
        "slip_fixed_converged": np.array(fixed_converged),
        "slip_sliding_rejected": np.array(sliding_rejected),
        "slip_sliding_error": np.array(sliding_error_msg),
    }


def _pack_ensemble_stats(
    target: Dict[str, Any], key_prefix: str, samples: np.ndarray
) -> None:
    """Store mean, std, and 95% Student's t CI into ``target``."""
    stats = compute_ensemble_statistics(samples, axis=-1)
    target[f"{key_prefix}_mean"] = stats["mean"]
    target[f"{key_prefix}_std"] = stats["std"]
    target[f"{key_prefix}_ci"] = stats["ci"]


def run_finite_extensibility_study(
    max_workers: Optional[int] = None,
) -> Dict[str, Any]:
    """Execute the full finite-extensibility & large-strain stiffening benchmark study."""
    effective_workers = resolve_max_workers(max_workers)
    configure_process_threading(workers=effective_workers, verbose=True)

    l_kuhn, r02, rho_beads, gamma_conv = _get_pdms_kuhn_constants()
    seeds = np.array([101, 102, 103, 104, 105, 106], dtype=np.int64)
    n_seeds = len(seeds)
    x_c = 0.88

    uni_lams_uni = np.array(
        [1.05, 1.2, 1.4, 1.6, 1.8, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0, 4.4, 4.8, 5.2],
        dtype=np.float64,
    )
    uni_lams_bim = np.array(
        [1.05, 1.2, 1.4, 1.6, 1.8, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0, 4.2],
        dtype=np.float64,
    )
    biax_lams_uni = np.array(
        [1.05, 1.2, 1.4, 1.6, 1.8, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0, 4.4, 4.8],
        dtype=np.float64,
    )
    biax_lams_bim = np.array(
        [1.05, 1.2, 1.4, 1.6, 1.8, 2.0, 2.4, 2.8, 3.2, 3.5, 3.8],
        dtype=np.float64,
    )
    solver_lams = np.array(
        [1.25, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.25, 4.5, 4.75, 5.0, 5.5, 6.0],
        dtype=np.float64,
    )
    reg_lams = np.array(
        [
            1.1,
            1.5,
            2.0,
            2.5,
            3.0,
            3.5,
            4.0,
            4.25,
            4.5,
            4.75,
            5.0,
            5.5,
            6.0,
            6.5,
            7.0,
            7.5,
        ],
        dtype=np.float64,
    )

    part1_results: List[Optional[Dict[str, Any]]] = [None] * n_seeds
    part2_results: List[Optional[Dict[str, Any]]] = [None] * n_seeds

    with ProcessPoolExecutor(
        max_workers=effective_workers,
        initializer=worker_process_initializer,
    ) as executor:
        p1_futures = {
            executor.submit(
                _simulate_part1_seed,
                idx,
                int(seed),
                uni_lams_uni,
                uni_lams_bim,
                biax_lams_uni,
                biax_lams_bim,
            ): idx
            for idx, seed in enumerate(seeds)
        }
        p2_futures = {
            executor.submit(
                _simulate_part2_reg_seed,
                idx,
                int(seed),
                reg_lams,
                x_c,
            ): idx
            for idx, seed in enumerate(seeds)
        }
        for fut in as_completed(p1_futures):
            res = fut.result()
            part1_results[res["seed_idx"]] = res
        for fut in as_completed(p2_futures):
            res = fut.result()
            part2_results[res["seed_idx"]] = res

    data: Dict[str, Any] = {
        "n_repeats": np.array(n_seeds),
        "seeds": seeds,
        "l_kuhn_nm": np.array(l_kuhn),
        "r02_nm2": np.array(r02),
        "rho_beads_nm3": np.array(rho_beads),
        "gamma_conv_mpa": np.array(gamma_conv),
        "x_c": np.array(x_c),
        "uni_lams_uni": uni_lams_uni,
        "uni_lams_bim": uni_lams_bim,
        "biax_lams_uni": biax_lams_uni,
        "biax_lams_bim": biax_lams_bim,
        "solver_lams": solver_lams,
        "reg_lams": reg_lams,
    }

    # Aggregate Part 1 ensembles
    for tag, u_lams, b_lams in (
        ("uni", uni_lams_uni, biax_lams_uni),
        ("bim", uni_lams_bim, biax_lams_bim),
    ):
        _pack_ensemble_stats(
            data,
            f"{tag}_g_eq",
            np.array([r[f"{tag}_g_eq"] for r in part1_results if r is not None]),
        )
        _pack_ensemble_stats(
            data,
            f"{tag}_active",
            np.array([r[f"{tag}_active"] for r in part1_results if r is not None]),
        )
        for mode_tag, lams in (("uni", u_lams), ("biax", b_lams)):
            for field in (
                "hook_sig",
                "hook_fstar",
                "ng_sig",
                "ng_fstar",
                "ng_max_x",
            ):
                key = f"{tag}_{mode_tag}_{field}"
                arr = np.column_stack([r[key] for r in part1_results if r is not None])
                assert arr.shape == (len(lams), n_seeds)
                _pack_ensemble_stats(data, key, arr)

    # Aggregate Part 2b ensembles
    for field in ("hook_sig", "sat_sig", "sat_max_x", "sci_sig", "sci_max_x"):
        arr = np.column_stack([r[field] for r in part2_results if r is not None])
        _pack_ensemble_stats(data, f"reg_{field}", arr)

    # Rigid Pade past lock-up: store raw 6-seed matrix + median/IQR + representative seed 101
    pade_sig_mat = np.column_stack(
        [r["pade_sig"] for r in part2_results if r is not None]
    )
    pade_max_x_mat = np.column_stack(
        [r["pade_max_x"] for r in part2_results if r is not None]
    )
    pade_exit_mat = np.column_stack(
        [r["pade_exit"] for r in part2_results if r is not None]
    )
    data["reg_pade_sig_all"] = pade_sig_mat
    data["reg_pade_sig_median"] = np.median(pade_sig_mat, axis=1)
    data["reg_pade_sig_q25"] = np.percentile(pade_sig_mat, 25.0, axis=1)
    data["reg_pade_sig_q75"] = np.percentile(pade_sig_mat, 75.0, axis=1)
    data["reg_pade_sig_rep"] = pade_sig_mat[:, 0]
    data["reg_pade_max_x_all"] = pade_max_x_mat
    data["reg_pade_max_x_rep"] = pade_max_x_mat[:, 0]
    data["reg_pade_max_x_max"] = np.max(pade_max_x_mat, axis=1)
    data["reg_pade_blowup_frac"] = np.mean(pade_max_x_mat >= 0.999, axis=1)
    data["reg_pade_exit_rep"] = pade_exit_mat[:, 0]
    data["reg_sat_exit_all"] = np.column_stack(
        [r["sat_exit"] for r in part2_results if r is not None]
    )
    data["reg_sci_exit_all"] = np.column_stack(
        [r["sci_exit"] for r in part2_results if r is not None]
    )

    # Constitutive single-spring force-extension curves for illustration
    x_curve = np.linspace(0.0, 1.25, 251)
    x_pade = x_curve[x_curve < 0.985]
    ev_sat_ref = SoftenedLangevinForceEvaluator(
        l=1.0, x_c=x_c, regularization_mode=RegularizationMode.BOUNDED_SATURATION
    )
    ev_sci_ref = SoftenedLangevinForceEvaluator(
        l=1.0, x_c=x_c, regularization_mode=RegularizationMode.BOND_SCISSION
    )
    _, f_sat_curve = ev_sat_ref.evaluate_dimensionless_curve(x_curve)
    _, f_sci_curve = ev_sci_ref.evaluate_dimensionless_curve(x_curve)
    data["constitutive_x"] = x_curve
    data["constitutive_x_pade"] = x_pade
    data["constitutive_f_hook"] = 3.0 * x_curve
    data["constitutive_f_pade"] = pade_inverse_langevin(x_pade)
    data["constitutive_f_sat"] = f_sat_curve
    data["constitutive_f_sci"] = f_sci_curve

    # Part 2a: Solver algorithm & tolerance benchmark near bimodal lock-up (seed 104)
    solver_seed = int(seeds[3])
    data["solver_seed"] = np.array(solver_seed)
    data.update(_run_solver_tolerance_benchmark(solver_lams, seed=solver_seed))

    # Part 2c: Fixed contour partition verification
    data.update(_verify_slip_mode_constraint())

    return data
