#!/usr/bin/env python3
"""Check that inputs which would run with a silent wrong answer are refused.

    python3 tools/input_guards.py --exe build/bin/yuzen --workdir /tmp/ig

Each case below is a shipped case with one thing changed. Without its
guard every one of them would run to completion -- a misspelled key keeping its
default, a one-sided Periodic or a CZ Pinning edge taken as Symmetry, a
Periodic edge on a rotating film (which is not periodic), a
pinned wall facing a film of full thickness failing only after many cut time
steps, with nothing in its messages about the wall. Each must stop at once
and say why; the
pinned_edgebump case is the same pinned wall set up correctly, which must
run. Then the convection schemes of [Numerics.Convection]: an unknown name
(for the concentration or for the height, the message listing the names
there are) and a Height scheme for CZ are refused, a known pair runs. Then
[Numerics.ManufacturedSolution]: a valid block runs, and a concentration
that leaves (0, 1) at either end or a film that touches zero is refused.
Then the ranges the manual promises and the defaults of what a file leaves
out: a step factor outside (0, 1) and counts below 1 are refused, the
disjoining pressure runs without its exponents, and a case with optional
blocks removed runs with those terms off. Last, what fuzzing found
(fuzz_found below): values that are not of their key's type, values with
which the equations mean nothing, a well-mixed run that breaks down, and a
fixed step whose solve does not converge; and a concentration solve that
does not converge, which must stop a run at a fixed step and cut the step
under step control (conc_unsolved). And a concentration that has left
[0, 1] on one side only, read from a field file, must stop the run as
diverged (conc_out_of_range); and the limit on the size of the grid holds
at either side of it (grid_size); a .inp string without its closing quote
is refused (inp_unterminated).

Only the Python standard library is used, so the test runs wherever ctest
does.

Copyright 2017-2026 Suguru Shiratori, Tokyo City University
SPDX-License-Identifier: Apache-2.0
"""
import argparse
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CASES = os.path.join(os.path.dirname(HERE), "cases")


def case(name):
    with open(os.path.join(CASES, name + ".toml"), encoding="utf-8") as f:
        return f.read()


def sub(s, pat, rep):
    s2, n = re.subn(pat, rep, s, count=1, flags=re.M)
    assert n == 1, pat
    return s2


def bc(s, edge, kind):
    return sub(s, r'(\[Boundary\.%sBC\]\n)Type = "[A-Za-z]+"' % edge, r'\g<1>Type = "%s"' % kind)


def listed(s, keys, keep_interval=False):
    """Output.FieldData with the given keys (and blocks after it) in place of
    its TimeInterval, unless keep_interval."""
    if not keep_interval:
        s = sub(s, r"^TimeInterval = 1e-1\n(NeutralFile)", r"\1")
    #  The keys right after the [Output.FieldData] header; a block given in
    #  them goes after the keys of [Output.FieldData]
    head, _, rest = s.partition("[Output.FieldData]\n")
    body, _, tail = rest.partition("\n[Output.GeneralHistory]")
    return head + "[Output.FieldData]\n" + body + "\n" + keys + "\n\n[Output.GeneralHistory]" + tail


def few_steps(s):
    return sub(s, r"^MaxSteps = .*$", "MaxSteps = 2")


EDGE_FILM = """[Initial.SurfaceShape]
TotalHeight = 10e-6
PrecursorRatio = 0.02

[Initial.SurfaceShape.EdgeBump]
Order = 1
CreateType = "add"
Base = 0.5
Height = 1.0
Location = "Xmax"
Alpha = 1.0
InnerLength = 0.05e-3
OuterLength = 0.15e-3

[Substrate]
TotalHeight = 10e-6

[Fluid]"""


def edge_film(s):
    return re.sub(r"^\[Initial\.SurfaceShape\]\n(?:.*\n)*?\[Fluid\]", EDGE_FILM, s,
                  count=1, flags=re.M)


#  (name, input text, whether it must run, text the output must contain)
def plan():
    lv = few_steps(case("leveling_1d_reg"))
    return [
        ("misspelled_key",
         sub(lv, r"^PrecursorRatio = 0\.5$", "PrecursorRatoi = 0.5"), False,
         "did you mean Initial.SurfaceShape.PrecursorRatio?"),
        ("key_in_wrong_block",
         sub(sub(lv, r"^PrecursorRatio = 0\.5\n", ""), r"^\[Initial\]$", "[Initial]\nPrecursorRatio = 0.5"),
         False, "this key is read in [Initial.SurfaceShape]"),
        ("misspelled_block",
         sub(lv, r"^\[Physics\.Evaporation\]$", "[Physics.Evaporaton]"), False,
         "did you mean [Physics.Evaporation]?"),
        ("future_format",
         sub(lv, r'^InputFormatVersion = .*$', 'InputFormatVersion = "2.0"'), False,
         "InputFormatVersion must be"),
        ("one_sided_periodic",
         bc(lv, "Xmin", "Periodic"), False,
         "Periodic must be set on both edges"),
        ("concz_pinning",
         bc(few_steps(case("concz_1d_reg")), "Xmax", "Pinning"), False,
         "the CZ model implements only Symmetry and Periodic"),
        ("pinned_full_film",
         bc(lv, "Xmax", "Pinning"), False,
         "Pinning holds the film thickness at zero"),
        ("pinned_edgebump",
         edge_film(bc(lv, "Xmax", "Pinning")), True,
         "Pinning Xmax: the initial film reaches"),
        #  The use EdgeBump is made for: a film flat at Base up to the
        #  meniscus, Height = Base, falling to zero at the pinned wall.
        #  MaxSteps is counted in 64 bits; 0 sets no limit, a negative count
        #  is refused.
        ("maxsteps_negative", sub(lv, r"^MaxSteps = .*$", "MaxSteps = -1"), False,
         "Numerics.MaxSteps must be 0 (no limit) or larger"),
        ("maxsteps_large", sub(lv, r"^MaxSteps = .*$", "MaxSteps = 50000000000"), True,
         "Program normally terminated"),
        ("pinned_edgebump_flat",
         sub(edge_film(bc(lv, "Xmax", "Pinning")), r"^Base = 0\.5$", "Base = 1.0"), True,
         "Pinning Xmax: the initial film reaches"),
        #  A rotating film is not periodic: the arm of the centrifugal force
        #  jumps by the length of the domain across the seam (VolErr drifts).
        #  Refused, but for a manufactured solution, whose source is formed
        #  for the arm of each cell (tools/mms.py runs every term on a
        #  periodic domain).
        ("periodic_centrifugal",
         spin(bc(bc(lv, "Xmin", "Periodic"), "Xmax", "Periodic")), False,
         "a rotating film is not periodic"),
        ("periodic_centrifugal_mms",
         spin(bc(bc(flat(lv), "Xmin", "Periodic"), "Xmax", "Periodic")) + mms(), True,
         "Program normally terminated"),
        #  The nozzle's liquid at its own concentration (ConcType =
        #  "supplied") takes Conc as a mass fraction: outside [0, 1] it is
        #  refused; with the default "local" Conc is not used, and runs.
        ("source_conc_above_one",
         sub(few_steps(case("nozzle_1d_fixed")), r"^Conc = 0\.2$", 'Conc = 1.5\nConcType = "supplied"'), False,
         "Physics.Source.Conc = 1.5: must be in [0, 1] with Physics.Source.ConcType = \"supplied\""),
        ("source_conc_negative",
         sub(few_steps(case("nozzle_1d_fixed")), r"^Conc = 0\.2$", 'Conc = -0.1\nConcType = "supplied"'), False,
         "Physics.Source.Conc = -0.1: must be in [0, 1]"),
        ("source_conc_local",
         sub(few_steps(case("nozzle_1d_fixed")), r"^Conc = 0\.2$", "Conc = 1.5"), True,
         "Program normally terminated"),
        ("source_conc_type",
         sub(few_steps(case("nozzle_1d_fixed")), r"^Conc = 0\.2$", 'Conc = 0.2\nConcType = "mixed"'), False,
         "Physics.Source.ConcType must be one of [`local',`supplied']"),
        #  A convection scheme Yuzen does not have is refused by name, and a
        #  scheme for the film height of CZ, which has no upwinded term, as
        #  a setting nothing reads; one it does have runs.
        ("unknown_scheme",
         lv + '\n[Numerics.Convection]\nConcentration = "WENO5"\n', False,
         "Numerics.Convection.Concentration = `WENO5' is not a scheme Yuzen has"),
        #  A value that differs from an accepted one only in case is refused
        #  too, with the spelling that would have been accepted.
        ("value_case",
         sub(lv, r'^(\[Numerics\.HeightEquation\.LinearSolver\]\n(?:[^\[\n].*\n)*?)PreconditionType = "ILU"',
             r'\1PreconditionType = "ilu"'), False,
         "Did you mean `ILU'? Values are case-sensitive."),
        ("scheme_case",
         lv + '\n[Numerics.Convection]\nConcentration = "koren"\n', False,
         "Did you mean `Koren'? Values are case-sensitive."),
        ("unknown_height_scheme",
         lv + '\n[Numerics.Convection]\nHeight = "minmod"\n', False,
         "Numerics.Convection.Height = `minmod' is not a scheme Yuzen has; it must be one of "
         "[`KK', `UTOPIA', `Central4', `Koren', `VanLeer', `QUICK', `Upwind1', `Upwind2', `Central2']"),
        ("concz_height_scheme",
         few_steps(case("concz_1d_reg")) + '\n[Numerics.Convection]\nHeight = "Koren"\n', False,
         "the CZ film-height equation has no upwinded convective term"),
        ("known_scheme",
         lv + '\n[Numerics.Convection]\nConcentration = "Koren"\nHeight = "VanLeer"\n', True,
         "Convection <conc>      = Koren"),
        #  The manufactured solution needs a film that never vanishes and a
        #  concentration strictly inside (0, 1), at both ends of its range:
        #  above 1 with the lower end fine, and touching 0 exactly
        #  (0.1 - 0.1 is 0 in floating point). The film without its relief
        #  (which the manufactured sources do not know), with a valid block,
        #  must run, so that each refusal is that check's alone.
        ("mms_valid",
         flat(lv) + mms(), True, "Program normally terminated"),
        ("mms_conc_above_one",
         flat(lv) + mms(ConcMean="0.95", ConcAmplitude="0.1"), False,
         "needs the concentration inside (0, 1)"),
        ("mms_conc_at_zero",
         flat(lv) + mms(ConcMean="0.1", ConcAmplitude="0.1"), False,
         "needs the concentration inside (0, 1)"),
        ("mms_film_touches_zero",
         flat(lv) + mms(Mean="2.0e-6", Amplitude="2.0e-6"), False,
         "needs Mean > Amplitude >= 0"),
        #  The ranges the manual promises (tools/doc_keys.py holds the two
        #  together): a domain of positive width, the step factor strictly
        #  inside (0, 1), and counts that must be at least 1.
        ("xmax_not_above_xmin",
         sub(lv, r"^Xmax = .*$", "Xmax = -0.5e-3"), False,
         "Grid.Xmax must exceed Grid.Xmin"),
        ("change_factor_one",
         sub(lv, r"^ChangeFactor = .*$", "ChangeFactor = 1.0"), False,
         "ChangeFactor must be between (0:1)"),
        ("change_factor_zero",
         sub(lv, r"^ChangeFactor = .*$", "ChangeFactor = 0.0"), False,
         "ChangeFactor must be between (0:1)"),
        ("keep_steps_zero",
         sub(lv, r"^KeepSteps = .*$", "KeepSteps = 0"), False,
         "KeepSteps must be 1 or larger"),
        ("lin_max_iterations_zero",
         sub(lv, r"^MaxIterations = 8000$", "MaxIterations = 0"), False,
         "Numerics.HeightEquation.LinearSolver.MaxIterations must be 1 or larger"),
        ("num_sub_region_zero",
         sub(lv, r"^SubRegions = .*$", "SubRegions = 0"), False,
         "SubRegions must be 1 or larger"),
        #  Enabling the disjoining pressure without its exponents takes the
        #  default pair (n, m) = (3, 2), as an absent block does (exponents
        #  of 0 would fail n > m > 1).
        ("disjoining_default_exponents",
         sub(sub(few_steps(case("stress_wetting")), r"^ExponentN = .*\n", ""), r"^ExponentM = .*\n", ""),
         True, "Program normally terminated"),
        #  An optional block that is left out reads as an empty one, and an
        #  absent Physics sub-block is a term that is off: no evaporation,
        #  no Marangoni term, no restart, no probes; none of its members is
        #  left indeterminate.
        ("absent_blocks",
         drop(drop(drop(drop(lv, "Physics.Evaporation"), "Physics.SolutoCapillarity"),
                   "Restart"), "Output.ProbeList"),
         True, " Ma = 0.000000e+00"),
        #  Field output at listed times (Output.FieldData.Times, LogSchedule)
        #  in place of TimeInterval: refused with it, out of order, at or
        #  below 0, all after EndTime, Times and a LogSchedule together, and
        #  ExactTimes with nothing listed
        ("times_with_interval", listed(lv, "Times = [0.05, 0.1]", keep_interval=True), False,
         "Output.FieldData.TimeInterval = 0.1: cannot be given with Output.FieldData.Times or LogSchedule"),
        ("times_unordered", listed(lv, "Times = [0.1, 0.05]"), False,
         "Output.FieldData.Times[1] = 0.05: must be after the time before it"),
        ("times_not_positive", listed(lv, "Times = [0.0, 0.1]"), False,
         "Output.FieldData.Times[0] = 0: must be positive"),
        ("times_after_end", listed(lv, "Times = [0.5, 1.0]"), False,
         "every listed time is after Numerics.EndTime"),
        ("times_and_log", listed(lv, "Times = [0.05]\n\n[Output.FieldData.LogSchedule]\nStart = 1e-3\nEnd = 0.1\nPerDecade = 5"),
         False, "give Times or a LogSchedule, not both"),
        ("exact_without_list", sub(lv, r"^\[Output\.FieldData\]$", "[Output.FieldData]\nExactTimes = true"), False,
         "Output.FieldData.ExactTimes = 1: applies to Output.FieldData.Times only"),
        ("log_start_zero", listed(lv, "\n[Output.FieldData.LogSchedule]\nStart = 0\nEnd = 0.1\nPerDecade = 5"), False,
         "Output.FieldData.LogSchedule.Start must be positive"),
        ("potential_concz", sub(few_steps(case("concz_1d_reg")), r"^\[Output\.FieldData\]$",
                                "[Output.FieldData]\nWritePotential = true"), False,
         "Output.FieldData.WritePotential = 1: is not implemented for Physics.ConcField.Model = \"CZ\""),
        ("times_empty", listed(lv, "Times = []"), False, "Output.FieldData.Times is empty"),
        ("energy_concz", few_steps(case("concz_1d_reg")) + "\n[Output.Energy]\nEnable = true\n", False,
         "Output.Energy.Enable = 1: is not implemented for Physics.ConcField.Model = \"CZ\""),
        ("slip_negative", lv + "\n[Physics.Slip]\nEnable = true\nLength = -1e-6\n", False,
         "Physics.Slip.Length = -1e-06: must not be negative"),
        ("slip_concz", few_steps(case("concz_1d_reg")) + "\n[Physics.Slip]\nEnable = true\nLength = 1e-6\n", False,
         "Physics.Slip is not implemented for Physics.ConcField.Model = \"CZ\""),
        ("slip_zero", lv + "\n[Physics.Slip]\nEnable = true\n", True,
         "Physics.Slip.Enable = true with Physics.Slip.Length = 0"),
        ("slip_valid", lv + "\n[Physics.Slip]\nEnable = true\nLength = 1e-6\n", True,
         "Lambda = 1.000000e-01"),
        ("times_valid", listed(lv, "Times = [1e-4, 0.05, 0.1, 0.5]\nExactTimes = true"), True,
         "FieldData = 3 listed times"),
        ("log_valid", listed(lv, "\n[Output.FieldData.LogSchedule]\nStart = 1e-4\nEnd = 0.1\nPerDecade = 5\nExactTimes = true"),
         True, "FieldData = 16 listed times"),
    ] + fuzz_found(lv)


def fuzz_found(lv):
    """What fuzzing the input found (yuzen/tests/fuzz, tools/proptest.py).

    Each of these must stop the run and say why: a value that is not of
    its key's type (taken for what strtod makes of it, "fast" would be 0
    and a NaN would pass every range check), a TOML date, a non-ASCII
    character in a key (undefined behavior in toml++), and values with
    which the equations mean nothing -- a viscosity, density, surface
    tension or height scale of zero or below, a time step or an output
    interval of zero, a film with no liquid in it, a shape that divides by
    a zero width -- which would end the run at once with exit code 0, run
    on with NaN, or never advance. So must a well-mixed run whose fixed
    time step breaks the film, which would otherwise go on to the end with
    NaN in its output and exit code 0.
    """
    two = few_steps(case("leveling_2d_reg"))
    return [
        ("negpow_exponent_one",
         sub(sub(lv, r'^SpatialFunc = "uniform"$', 'SpatialFunc = "negPow"'),
             r"^NegPowExponent = 0\.25$", "NegPowExponent = 1.0"), False,
         "Physics.Evaporation.NegPowExponent = 1: must be less than 1"),
        ("viscosity_zero", sub(lv, r"^BaseValue = 100\.0e-3$", "BaseValue = 0.0"), False,
         "Fluid.Viscosity.BaseValue = 0: must be positive"),
        ("surface_tension_negative", sub(lv, r"^BaseValue = 30\.0e-3$", "BaseValue = -30.0e-3"), False,
         "Fluid.SurfaceTension.BaseValue = -0.03: must be positive"),
        ("viscosity_overflow", sub(lv, r"^Exponent = 10\.0$", "Exponent = 2000.0"), False,
         "Fluid.Viscosity.Exponent = 2000: makes the viscosity overflow within c in [0, 1]"),
        ("density_zero", sub(lv, r"^Value = 1000\.0$", "Value = 0.0"), False,
         "Fluid.Density.Value = 0: must be positive"),
        ("height_scale_zero", sub(lv, r"^TotalHeight = 10e-6$", "TotalHeight = 0.0"), False,
         "Initial.SurfaceShape.TotalHeight = 0: must be positive"),
        ("precursor_negative", sub(lv, r"^PrecursorRatio = 0\.5$", "PrecursorRatio = -0.1"), False,
         "Initial.SurfaceShape.PrecursorRatio = -0.1: must be in [0, 1)"),
        ("initial_dt_zero", sub(lv, r"^InitialDt = .*$", "InitialDt = 0.0"), False,
         "Numerics.InitialDt = 0: must be positive"),
        ("end_time_negative", sub(lv, r"^EndTime = .*$", "EndTime = -1.0"), False,
         "Numerics.EndTime = -1: must not be negative"),
        ("output_interval_zero", sub(lv, r"^TimeInterval = 1e-1$", "TimeInterval = 0.0"), False,
         "Output.FieldData.TimeInterval = 0: must be positive"),
        #  a negative largest step: a run that never ended
        ("max_dt_factor_negative", sub(lv, r"^MaxDtFactor = .*$", "MaxDtFactor = -2147483648"), False,
         "Numerics.TimeStepControl.MaxDtFactor = -2.14748e+09: must be positive"),
        ("max_iter_nk_zero", sub(lv, r"^NewtonMaxIterations = .*$", "NewtonMaxIterations = 0"), False,
         "Numerics.HeightEquation.NewtonMaxIterations = 0: must be at least 1"),
        ("shape_width_zero", sub(lv, r"^Width = \[1\.0, 0\.5\]$", "Width = [0.0, 0.5]"), False,
         "Initial.SurfaceShape.EllipseBump[0].Width[0] = 0: must not be zero"),
        ("y_extent_2d", sub(two, r"^Ymax = .*$", "Ymax = -0.5e-3"), False,
         "Grid.Ymax = -0.0005: must exceed Grid.Ymin in a 2-D film"),
        ("no_liquid", sub(sub(lv, r"^PrecursorRatio = 0\.5$", "PrecursorRatio = 0.0"),
                          r"^\[Initial\.SurfaceShape\.EllipseBump\]\n(?:[^\[\n].*\n|\n)*", ""),
         False, "the initial film has no liquid in it"),
        ("not_a_number", sub(lv, r"^BaseValue = 100\.0e-3$", 'BaseValue = "fast"'), False,
         "Variable [Fluid.Viscosity.BaseValue] = `fast' is not a finite number"),
        ("nan_value", sub(lv, r"^Xmax = .*$", "Xmax = nan"), False,
         "Variable [Grid.Xmax] = `nan' is not a finite number"),
        ("fractional_count", sub(lv, r"^Nx = 300$", "Nx = 300.5"), False,
         "Variable [Grid.Nx] = `300.5' is not a whole number"),
        ("bool_word", sub(lv, r"^Overwrite = true$", 'Overwrite = "maybe"'), False,
         "Variable [Output.Overwrite] = `maybe' is not true or false"),
        ("empty_list", sub(lv, r"^Value = 1000\.0$", "Value = []"), False,
         "Variable [Fluid.Density.Value]: an empty value or list where one value is wanted"),
        ("toml_date", sub(lv, r"^EndTime = .*$", "EndTime = 1979-05-27"), False,
         "Variable [Numerics.EndTime]: a date or time is not a value Yuzen reads"),
        ("non_ascii_key", sub(lv, r"^\[Physics\.SolutoCapillarity\]$", "[Physics.Soluto\u011epillarity]"),
         False, "TOML parse error"),
        #  A fixed step far too large for the film stops the run, rather
        #  than running on to the end with exit code 0: its Newton
        #  correction stagnates far from convergence at the first step
        #  (found taken with a warning, the film then went negative and the
        #  divergence test stopped the run; that test is checked on a film
        #  read below zero, conc_out_of_range).
        ("wm_diverges", huge_fixed_step(lv), False, "stagnated far from convergence"),
        #  A fixed step whose linear solve cannot converge stops the run; it
        #  must not be taken as the failed iteration left it (the film it
        #  started from) with the run ending normally (proptest seed 3094).
        ("wm_fixed_unsolved", unsolvable_fixed_step(lv), False,
         "cannot be retried at a shorter one"),
    ] + conc_unsolved()


def conc_unsolved():
    """A concentration solve that does not converge (tools/proptest.py, the
    pre-release checks): the well-mixed and the ConcZ concentration
    equations took its answer as it stood, with the step control on or off,
    and the run ended normally. At a fixed step the run must stop and name
    the concentration solve; under step control the step must be cut and
    solved again: these inputs take no shorter step with the concentration
    solver of the case (allphys at 5e-6 s, concz at 1e-6 s), so the "dt-" of
    their output is the concentration's. (In 1-D the ILU preconditioner is
    exact, so the well-mixed cases are given none. One unpreconditioned
    iteration does not take the change of the well-mixed concentration of
    allphys over a step of 5e-6 s to a relative residual of 1e-8, and does
    over the steps cut from it, ten times; at 1e-12 no cut is short enough
    and the run stops on the shortest step.)"""
    wm = case("allphys_1d_reg")
    cz = case("concz_1d_reg")
    return [
        ("wm_conc_fixed_unsolved", fixed_step(conc_solver(few_steps(wm), 1, "1e-30", "None")), False,
         "Newton-Kantorovich: the step did not converge at the fixed step: the concentration solve"),
        ("cz_conc_fixed_unsolved", fixed_step(conc_solver(few_steps(cz), 1, "1e-30")), False,
         "ConcZ: the step did not converge at the fixed step: the concentration solve"),
        ("wm_conc_cut", steps_from(conc_solver(wm, 1, "1e-8", "None"), 20, "5e-6"), True,
         "dt- ("),
        ("cz_conc_cut", steps_from(conc_solver(cz, 2, "1e-10", "Jacobi"), 20, "1e-6"), True,
         "dt- ("),
    ]


def conc_solver(s, maxit, eps, precond=None):
    """The concentration's linear solver limited to maxit iterations at a
    tolerance eps, with another preconditioner if given."""
    s = sub(s, r"^(\[Numerics\.ConcEquation\.LinearSolver\]\n)MaxIterations = \d+\nEpsilon = .*$",
            r"\g<1>MaxIterations = %d\nEpsilon = %s" % (maxit, eps))
    if precond:
        s = sub(s, r'^(\[Numerics\.ConcEquation\.LinearSolver\]\n(?:[^\[\n].*\n)*?)PreconditionType = "ILU"',
                r'\g<1>PreconditionType = "%s"' % precond)
    return s


def fixed_step(s):
    return sub(s, r"^(\[Numerics\.TimeStepControl\]\n)Enable = true", r"\g<1>Enable = false")


def steps_from(s, n, dt):
    """n steps from a first step of dt."""
    return sub(sub(s, r"^MaxSteps = .*$", "MaxSteps = %d" % n), r"^InitialDt = .*$", "InitialDt = " + dt)


def unsolvable_fixed_step(s):
    """The step control off and a film-height solve that cannot meet its
    tolerance in the one iteration it is given, unpreconditioned: with ILU,
    exact in 1-D, one iteration meets the looser tolerance that the
    Eisenstat-Walker choice gives the early Newton passes (util/Opt.h)."""
    s = sub(s, r"^(\[Numerics\.TimeStepControl\]\n)Enable = true", r"\g<1>Enable = false")
    s = sub(s, r"^(\[Numerics\.HeightEquation\.LinearSolver\]\n(?:[^\[\n].*\n)*?)PreconditionType = .*$",
            r'\g<1>PreconditionType = "None"')
    return sub(s, r"^(\[Numerics\.HeightEquation\.LinearSolver\]\n)MaxIterations = \d+\nEpsilon = .*$",
               r"\g<1>MaxIterations = 1\nEpsilon = 1e-30")


def huge_fixed_step(s):
    """A fixed step of 10 s, a thousand times the capillary time of the
    film, with output intervals that do not cut it down."""
    s = re.sub(r"^TimeInterval = .*$", "TimeInterval = 1e3", s, flags=re.M)
    s = sub(sub(s, r"^InitialDt = .*$", "InitialDt = 10.0"), r"^EndTime = .*$", "EndTime = 1e3")
    return sub(s, r"^(\[Numerics\.TimeStepControl\]\n)Enable = true", r"\g<1>Enable = false")


def drop(s, block):
    """The case without the table [block] (its header and its keys)."""
    return sub(s, r"^\[%s\]\n(?:[^\[\n].*\n|\n)*" % re.escape(block), "")


def flat(s):
    """The case without its substrate relief."""
    return sub(s, r"^\[Substrate\.EllipseBump\]\n(?:[^\[\n].*\n|\n)*", "")


def spin(s):
    """The case with the centrifugal term on."""
    return sub(s, r"^(\[Physics\.Centrifugal\]\n)Enable = false", r"\g<1>Enable = true")


def mms(**kv):
    """A [Numerics.ManufacturedSolution] block, valid unless kv changes it."""
    p = {"Enable": "true", "Mean": "5.0e-6", "Amplitude": "1.5e-6", "Omega": "60.0",
         "Wavelength": "[0.00025, 0.00025]", "ConcMean": "0.5", "ConcAmplitude": "0.1"}
    p.update(kv)
    return "\n[Numerics.ManufacturedSolution]\n" + "".join("%s = %s\n" % i for i in p.items())


def conc_out_of_range(exe, workdir):
    """A concentration that has left [0, 1] by more than the slack (0.1)
    must stop the run as diverged, on either side alone: a run of two steps
    writes a field file, one column of it is set to 1.25 or to -0.25, and a
    run restarted from it (the field file alone, as a hand-made one is read)
    must stop at its first step. ConcZ both ways; the well-mixed model
    below 0 only, since it moves what lies above 1 to the cells around
    after every step. Both fill a cell below 0 from the cells around it,
    and the well-mixed step mixes the column back into range by itself:
    the test of divergence takes the lowest c of the start of the step and
    of the solve, before the fill. (Mutation testing: a test that stopped
    only when c left at both ends at once passed every other check.) And a
    well-mixed film read 1 um below zero in one column must stop the run as
    diverged, the film thickness having gone negative."""
    bad = 0
    for model, cs, val in (("CZ", "concz_1d_reg", 1.25), ("CZ", "concz_1d_reg", -0.25),
                           ("WM", "allphys_1d_reg", -0.25), ("WM", "allphys_1d_reg", "film")):
        if val == "film":
            name, col, val = "wm_film_below", "HH[um]", -1.0
        else:
            name, col = "%s_conc_%s" % (model.lower(), "above" if val > 1 else "below"), "CC[-]"
        d = os.path.join(workdir, name)
        os.makedirs(d)
        s = sub(case(cs), r"^JobName = .*$", 'JobName = "src"')
        s = sub(sub(s, r"^EndTime = .*$", "EndTime = 2e-3"), r"^InitialDt = .*$", "InitialDt = 1e-3")
        s = sub(s, r"^(\[Output\.FieldData\]\n)TimeInterval = .*$", r"\g<1>TimeInterval = 1e-3")
        outs = []
        for job, text in (("a", s), ("b", None)):
            if text is None:
                #  one column (the cell nearest x = 0, every layer) set to val
                lines = open(os.path.join(d, "src", "field00001.dat"), encoding="utf-8").read().split("\n")
                hdr = [h.strip() for h in lines[0].lstrip("#").split("\t")]
                ic = hdr.index(col)
                x0 = min({float(l.split("\t")[0]) for l in lines[1:] if l.strip()}, key=abs)
                for n, l in enumerate(lines[1:], 1):
                    t = l.split("\t")
                    if l.strip() and float(t[0]) == x0:
                        t[ic] = repr(val)
                        if col == "HH[um]":     # the free surface with it
                            t[hdr.index("HS[um]")] = repr(float(t[hdr.index("SS[um]")]) + val)
                        lines[n] = "\t".join(t)
                os.makedirs(os.path.join(d, "start"))
                with open(os.path.join(d, "start", "field00001.dat"), "w", encoding="utf-8") as f:
                    f.write("\n".join(lines))
                text = sub(sub(s, r"^JobName = .*$", 'JobName = "rst"'), r"^EndTime = .*$", "EndTime = 5e-3")
                text = sub(text, r"^\[Restart\]\n(?:[A-Za-z].*\n)*",
                           '[Restart]\nEnable = true\nASCII = true\nFileName = "start/field00001.dat"\n'
                           "RestartTime = 1e-3\nFileIndex = 2\n")
            with open(os.path.join(d, job + ".toml"), "w", encoding="utf-8") as f:
                f.write(text)
            try:
                r = subprocess.run([exe, job + ".toml"], cwd=d, capture_output=True, text=True, encoding="utf-8", errors="replace", timeout=300)
                outs.append((r.returncode, r.stdout + r.stderr))
            except subprocess.TimeoutExpired:
                outs.append((-1, "timed out"))
            if outs[-1][0] != 0:
                break
        expect = ("Diverged: the film thickness went negative" if col == "HH[um]"
                  else "Diverged: the concentration left [0,1]")
        if len(outs) < 2 or outs[0][0] != 0:
            print("FAIL %-20s the run that writes the field failed: %s" % (name, outs[0][1][-300:]))
            bad += 1
        elif outs[1][0] == 0 or expect not in outs[1][1]:
            print("FAIL %-20s %s = %g: %s, without saying: %s" % (
                name, "h" if col == "HH[um]" else "c", val, "ran" if outs[1][0] == 0 else "stopped (rc=%d)" % outs[1][0], expect))
            bad += 1
        else:
            print("ok   %-20s stopped" % name)
    return bad


def grid_size(exe, workdir):
    """The limit on the size of the grid, (Nx+4)(Ny+4)(Nz+4) <= 5e7, at
    either side of it (a 2-D well-mixed grid, Nz = 1): 3000 x 3330 must be
    refused, 3000 x 3320 accepted. Read with --check-input, which runs
    nothing: a grid let through would take gigabytes."""
    bad = 0
    two = case("leveling_2d_reg")
    for ny, ok in ((3330, False), (3320, True)):
        name = "grid_%dx%d" % (3000, ny)
        d = os.path.join(workdir, name)
        os.makedirs(d)
        s = sub(sub(two, r"^Nx = .*$", "Nx = 3000"), r"^Ny = .*$", "Ny = %d" % ny)
        with open(os.path.join(d, "a.toml"), "w", encoding="utf-8") as f:
            f.write(s)
        r = subprocess.run([exe, "--check-input", "a.toml"], cwd=d, capture_output=True,
                           text=True, encoding="utf-8", errors="replace", timeout=120)
        out = r.stdout + r.stderr
        expect = "(Grid.Nx+4)*(Grid.Ny+4)*(Grid.Nz+4) must not exceed 50000000"
        if (r.returncode == 0) != ok or (not ok and expect not in out):
            print("FAIL %-20s %s" % (name, "accepted" if r.returncode == 0 else "refused: " + out[-300:]))
            bad += 1
        else:
            print("ok   %-20s %s" % (name, "accepted" if ok else "refused"))
    return bad


def inp_unterminated(exe, workdir):
    """A .inp string without its closing quote must be refused, naming the
    key and the line it starts on. (It took the rest of the file as its
    value, and every key after it kept its default.)"""
    d = os.path.join(workdir, "inp_unterminated")
    os.makedirs(d)
    with open(os.path.join(d, "a.inp"), "w", encoding="utf-8") as f:
        f.write('SolverName = "Yuzen-1.0"\nJobName = "a\nblock Grid\n  Nx = 10\nend\n')
    r = subprocess.run([exe, "--check-input", "a.inp"], cwd=d, capture_output=True, text=True, encoding="utf-8", errors="replace", timeout=120)
    out = r.stdout + r.stderr
    expect = "a.inp:2: unterminated string for 'JobName'"
    if r.returncode == 0 or expect not in out:
        print("FAIL %-20s %s, without saying: %s" % ("inp_unterminated", "accepted" if r.returncode == 0
                                                     else "refused", expect))
        print("\n".join(out.strip().splitlines()[-4:]))
        return 1
    print("ok   %-20s refused" % "inp_unterminated")
    return 0


def thread_requests(exe, workdir):
    """More threads asked for (OMP_NUM_THREADS) than the run can use
    (util/ThreadCap.h: 300 cells, one thread per 1000 cells, at least 8) stop the run before its result
    directory is made, and say what to ask for; YUZEN_THREAD_CAP=0 lets it
    run on them; 8 runs."""
    bad = 0
    text = few_steps(case("leveling_1d_reg"))
    #  A build without OpenMP (AppleClang, or clang++ without libomp: CMake
    #  then builds serially) runs on one thread whatever is asked for, and
    #  has nothing to check here. It is recognized by the runtime it reports.
    d = os.path.join(workdir, "threads_probe")
    os.makedirs(d)
    with open(os.path.join(d, "probe.toml"), "w", encoding="utf-8") as f:
        f.write(re.sub(r"^JobName = .*$", 'JobName = "probe"', text, count=1, flags=re.M))
    r = subprocess.run([exe, "probe.toml"], cwd=d, capture_output=True, text=True, encoding="utf-8", errors="replace", timeout=300)
    if r.returncode == 0 and "OMP_MAX_THREADS" not in r.stdout + r.stderr:
        print("skip threads_*             built without OpenMP")
        return 0
    for name, env, must_run, expect in (
            #  YUZEN_THREAD_CAP named, so that the size term applies whatever
            #  the runtime (with libomp it is left out by default)
            ("threads_refused", {"OMP_NUM_THREADS": "16", "YUZEN_THREAD_CAP": "1000"}, False,
             "OMP_NUM_THREADS=16 is more threads than this run can use: 300 grid cells, so at most 8 threads"),
            ("threads_cap_off", {"OMP_NUM_THREADS": "16", "YUZEN_THREAD_CAP": "0"}, True, "OMP_MAX_THREADS : 16"),
            ("threads_at_limit", {"OMP_NUM_THREADS": "8"}, True, "OMP_MAX_THREADS : 8")):
        d = os.path.join(workdir, name)
        os.makedirs(d)
        t = re.sub(r"^JobName = .*$", 'JobName = "%s"' % name, text, count=1, flags=re.M)
        with open(os.path.join(d, name + ".toml"), "w", encoding="utf-8") as f:
            f.write(t)
        e = dict(os.environ, **env)
        if "YUZEN_THREAD_CAP" not in env:
            e.pop("YUZEN_THREAD_CAP", None)
        r = subprocess.run([exe, name + ".toml"], cwd=d, capture_output=True, text=True, encoding="utf-8", errors="replace", timeout=300, env=e)
        out = r.stdout + r.stderr
        ran = (r.returncode == 0)
        made = os.path.isdir(os.path.join(d, name))
        if ran != must_run or expect not in out or (not must_run and made):
            print("FAIL %-20s %s%s%s" % (name, "ran" if ran else "stopped (rc=%d)" % r.returncode,
                                         "" if expect in out else ", without saying: %s" % expect,
                                         ", the result directory made" if (not must_run and made) else ""))
            bad += 1
        else:
            print("ok   %-20s %s" % (name, "runs" if must_run else "refused"))
    return bad


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--exe", required=True)
    ap.add_argument("--workdir", required=True)
    args = ap.parse_args()
    exe = os.path.abspath(args.exe)
    shutil.rmtree(args.workdir, ignore_errors=True)

    bad = 0
    for name, text, must_run, expect in plan():
        d = os.path.join(args.workdir, name)
        os.makedirs(d)
        text = re.sub(r"^JobName = .*$", 'JobName = "%s"' % name, text, count=1, flags=re.M)
        with open(os.path.join(d, name + ".toml"), "w", encoding="utf-8") as f:
            f.write(text)
        try:
            r = subprocess.run([exe, name + ".toml"], cwd=d, capture_output=True,
                               text=True, encoding="utf-8", errors="replace", timeout=300)
            rc, out = r.returncode, r.stdout + r.stderr
        except subprocess.TimeoutExpired:
            rc, out = -1, "timed out"
        ran = (rc == 0)
        if ran != must_run or expect not in out:
            why = ("ran" if ran else "stopped (rc=%d)" % rc) + \
                  ("" if expect in out else ", without saying: %s" % expect)
            print("FAIL %-20s %s" % (name, why))
            print("\n".join(out.strip().splitlines()[-8:]))
            bad += 1
        else:
            print("ok   %-20s %s" % (name, "runs" if must_run else "refused"))
    bad += conc_out_of_range(exe, args.workdir)
    bad += grid_size(exe, args.workdir)
    bad += inp_unterminated(exe, args.workdir)
    bad += thread_requests(exe, args.workdir)
    print("%d failure(s)" % bad)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
