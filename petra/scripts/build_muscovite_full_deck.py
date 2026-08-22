#!/usr/bin/env python3
"""Render the tracked schema-v2 muscovite full-mechanism deck.

The generator keeps the isotope × direction × reservoir rule matrix explicit
without making the reviewed TOML hand-maintained. Its output is deterministic;
the focused test requires the tracked deck to be byte-identical.
"""

from __future__ import annotations

import argparse
from pathlib import Path


def render_deck() -> str:
    lines: list[str] = []

    def add(*values: str) -> None:
        lines.extend(values)

    add(
        "# E2 full muscovite mechanism: scheduled, species- and defect-resolved.",
        "# Published anchors and proxy barriers are documented in",
        "# docs/program/results/E2-muscovite-full-mechanism.md.",
        "[deck]",
        'name = "muscovite-full-mechanism"',
        "schema = 2",
        'units = "kcal/mol"',
        'comment = "mechanism comparison; proxy brackets are not computed kinetics"',
        "",
        "[structure]",
        'kind = "cell"',
    )
    for name, mass, charge in (
        ("K", None, 1.0),
        ("Ar40", 39.962383, None),
        ("Ar39", 38.964313, None),
        ("Ar36", 35.967545, None),
        ("OH", None, None),
        ("O", None, -2.0),
    ):
        add("[[structure.species]]", f'name = "{name}"')
        if mass is not None:
            add(f"mass = {mass}")
        if charge is not None:
            add(f"charge = {charge}")

    kinds = (
        (
            "K_gallery_A",
            "K",
            (
                ("K", "K"),
                ("vacant", "vacant"),
                ("Ar40", "Ar40"),
                ("Ar39", "Ar39"),
                ("Ar36", "Ar36"),
            ),
        ),
        (
            "K_gallery_B",
            "K",
            (
                ("K", "K"),
                ("vacant", "vacant"),
                ("Ar40", "Ar40"),
                ("Ar39", "Ar39"),
                ("Ar36", "Ar36"),
            ),
        ),
        ("OH_donor", "OH", (("OH", "OH"), ("water_removed", "vacant"))),
        ("OH_acceptor", "OH", (("OH", "OH"), ("O_residual", "O"))),
        ("Interface", "bonded", (("bonded", "vacant"), ("delaminated", "vacant"))),
        ("Defect_zone", "pristine", (("pristine", "vacant"), ("extended", "vacant"))),
        ("Octahedral_trap", "vacant", (("vacant", "vacant"), ("Ar39_trapped", "Ar39"))),
        ("Surface_gate", "closed", (("closed", "vacant"), ("open", "vacant"))),
    )
    for kind, initial, states in kinds:
        add("", "[[structure.kinds]]", f'name = "{kind}"', f'initial = "{initial}"')
        for state, occupant in states:
            add(
                "[[structure.kinds.states]]",
                f'name = "{state}"',
                f'occupant = "{occupant}"',
            )

    add(
        "",
        "[structure.cell]",
        "a = 5.19",
        "b = 9.02",
        "c = 20.05",
        "alpha = 90.0",
        "beta = 95.7",
        "gamma = 90.0",
    )
    sites = (
        ("K_gallery_A", (0.0, 0.0, 0.25)),
        ("K_gallery_B", (0.0, 0.0, 0.75)),
        ("OH_donor", (0.25, 0.25, 0.45)),
        ("OH_acceptor", (0.75, 0.75, 0.45)),
        ("Interface", (0.5, 0.5, 0.0)),
        ("Defect_zone", (0.1, 0.1, 0.5)),
        ("Octahedral_trap", (0.5, 0.0, 0.5)),
        ("Surface_gate", (0.9, 0.9, 0.5)),
    )
    for kind, frac in sites:
        add(
            "[[structure.cell.sites]]",
            f'kind = "{kind}"',
            f"frac = [{frac[0]}, {frac[1]}, {frac[2]}]",
        )

    def bond(i: int, j: int, dcell: tuple[int, int, int], label: str) -> None:
        add(
            "[[structure.cell.bonds]]",
            f"i = {i}",
            f"j = {j}",
            f"dcell = [{dcell[0]}, {dcell[1]}, {dcell[2]}]",
            f'label = "{label}"',
        )

    bond(0, 1, (0, 0, 0), "gallery_intra")
    bond(0, 1, (0, 0, -1), "gallery_inter")
    bond(2, 3, (0, 0, 0), "oh_pair")
    bond(2, 3, (0, 0, -1), "steric")
    bond(2, 3, (0, 0, 1), "steric")
    for gallery in (0, 1):
        bond(gallery, 3, (0, 0, 0), "local_structure")
        bond(gallery, 4, (0, 0, 0), "gallery_interface")
        bond(gallery, 5, (0, 0, 0), "local_defect")
        bond(gallery, 7, (0, 0, 0), "surface_gate")
    bond(4, 3, (0, 0, 0), "delam_driver")
    bond(4, 3, (0, 0, -1), "delam_driver")
    bond(6, 0, (0, 0, 0), "trap_gallery")

    add(
        "",
        "[structure.lattice]",
        "dims = [4, 4, 6]",
        'boundary = ["periodic", "periodic", "open"]',
    )

    def init(
        name: str,
        kind: str,
        state: str,
        probability: float,
        set_state: str,
        region: str | None = None,
    ) -> None:
        add(
            "[[structure.init]]",
            f'name = "{name}"',
            f'center = {{ kind = "{kind}", state = ["{state}"] }}',
        )
        if region:
            add(f"region = {region}")
        add(f"probability = {probability}", f'set = "{set_state}"')

    init("extended-defect-zones", "Defect_zone", "pristine", 0.2, "extended")
    for gallery in ("K_gallery_A", "K_gallery_B"):
        suffix = gallery[-1]
        init(f"vacancies-{suffix}", gallery, "K", 0.90, "vacant")
        init(f"ar40-{suffix}", gallery, "K", 0.55, "Ar40")
        init(f"ar39-{suffix}", gallery, "K", 0.55, "Ar39")
        init(f"ar36-{suffix}", gallery, "K", 0.55, "Ar36")
    init("recoil-octahedral-ar39", "Octahedral_trap", "vacant", 0.35, "Ar39_trapped")
    init(
        "lower-basal-surface",
        "Surface_gate",
        "closed",
        1.0,
        "open",
        "{ axis = 2, min = 0, max = 0 }",
    )
    init(
        "upper-basal-surface",
        "Surface_gate",
        "closed",
        1.0,
        "open",
        "{ axis = 2, min = 5, max = 5 }",
    )

    add(
        "",
        "[dynamics.thermo]",
        "temperature = 773.15",
        "",
        "# Published 57 kcal/mol bulk dehydroxylation anchor (Kodama & Brydon).",
        "[[dynamics.rules]]",
        'name = "dehydroxylate_pair"',
        'center = { kind = "OH_donor", state = ["OH"] }',
        'guards = [{ kind = "OH_acceptor", label = "oh_pair", state = ["OH"], min = 1 }]',
        "rate = { eyring = { dh = 57.0, ds = 0.0 } }",
        "[[dynamics.rules.modifiers]]",
        'select = { kind = "OH_acceptor", label = "steric", state = ["O_residual"] }',
        "by_count = { dea = [0.0, 3.0, 12.0] }",
        "[[dynamics.rules.effects]]",
        'target = "center"',
        'set = "water_removed"',
        "[[dynamics.rules.effects]]",
        'target = "neighbor"',
        'select = { kind = "OH_acceptor", label = "oh_pair", state = ["OH"] }',
        'set = "O_residual"',
        "",
        "# Explicit interface transition. 58 kcal/mol is a sensitivity proxy,",
        "# not a computed barrier; analysis brackets it by ±5 kcal/mol.",
        "[[dynamics.rules]]",
        'name = "delaminate_interface"',
        'center = { kind = "Interface", state = ["bonded"] }',
        'guards = [{ kind = "OH_acceptor", label = "delam_driver", state = ["O_residual"], min = 1 }]',
        "rate = { eyring = { dh = 58.0, ds = 0.0 } }",
        "[[dynamics.rules.effects]]",
        'target = "center"',
        'set = "delaminated"',
        "",
        "# Recoil 39Ar escape is likewise a phase-2 proxy pending phase-3 NEB.",
        "[[dynamics.rules]]",
        'name = "escape_octahedral_Ar39"',
        'center = { kind = "Octahedral_trap", state = ["Ar39_trapped"] }',
        'guards = [{ kind = "K_gallery_A", label = "trap_gallery", state = ["vacant"], min = 1 }]',
        "rate = { eyring = { dh = 64.0, ds = 0.0 } }",
        "[[dynamics.rules.effects]]",
        'target = "center"',
        'set = "vacant"',
        "[[dynamics.rules.effects]]",
        'target = "neighbor"',
        'select = { kind = "K_gallery_A", label = "trap_gallery", state = ["vacant"] }',
        'set = "Ar39"',
    )

    directions = (
        ("A_to_B_intra", "K_gallery_A", "K_gallery_B", "gallery_intra"),
        ("A_to_B_inter", "K_gallery_A", "K_gallery_B", "gallery_inter"),
        ("B_to_A_intra", "K_gallery_B", "K_gallery_A", "gallery_intra"),
        ("B_to_A_inter", "K_gallery_B", "K_gallery_A", "gallery_inter"),
    )
    for zone, barrier in (("pristine", 66.0), ("extended", 40.0)):
        add(
            "",
            f"# {zone} reservoir hops: {barrier:g} kcal/mol. Extended-zone value is",
            "# the Nteme-2023 six-order enhancement mapped to a proxy barrier.",
        )
        for isotope in ("Ar40", "Ar39", "Ar36"):
            for direction, center, neighbor, label in directions:
                add(
                    "[[dynamics.rules]]",
                    f'name = "hop_{zone}_{isotope}_{direction}"',
                    f'center = {{ kind = "{center}", state = ["{isotope}"] }}',
                    "guards = [",
                    f'  {{ kind = "{neighbor}", label = "{label}", state = ["vacant"], min = 1 }},',
                    f'  {{ kind = "{neighbor}", state = ["vacant"], min = 2 }},',
                    '  { kind = "OH_acceptor", label = "local_structure", state = ["O_residual"], min = 1 },',
                    f'  {{ kind = "Defect_zone", label = "local_defect", state = ["{zone}"], min = 1 }},',
                    "]",
                    f"rate = {{ eyring = {{ dh = {barrier}, ds = 0.0 }} }}",
                    "[[dynamics.rules.effects]]",
                    'target = "center"',
                    'set = "vacant"',
                    "[[dynamics.rules.effects]]",
                    'target = "neighbor"',
                    f'select = {{ kind = "{neighbor}", label = "{label}", state = ["vacant"] }}',
                    f'set = "{isotope}"',
                )

    for isotope in ("Ar40", "Ar39", "Ar36"):
        for gallery in ("K_gallery_A", "K_gallery_B"):
            suffix = gallery[-1]
            for mechanism, gate_kind, gate_label, gate_state in (
                ("surface", "Surface_gate", "surface_gate", "open"),
                ("delamination", "Interface", "gallery_interface", "delaminated"),
            ):
                add(
                    "[[dynamics.rules]]",
                    f'name = "release_{isotope}_{mechanism}_{suffix}"',
                    f'center = {{ kind = "{gallery}", state = ["{isotope}"] }}',
                    f'guards = [{{ kind = "{gate_kind}", label = "{gate_label}", state = ["{gate_state}"], min = 1 }}]',
                    "rate = { constant = 1.0 }",
                    "[[dynamics.rules.effects]]",
                    'target = "center"',
                    'set = "vacant"',
                )

    add(
        "",
        "[execution]",
        'strategy = "ctmc"',
    )
    for temperature in (773.15, 823.15, 873.15, 923.15, 973.15):
        add(
            "[[execution.schedule]]", f"temperature = {temperature}", "duration = 600.0"
        )
    add(
        "[execution.stop]",
        "steps = 250000",
        "[execution.ensemble]",
        "seed = 1998",
        "n_replicas = 1",
        'seed_policy = "increment"',
        "",
        "[observables]",
        "report_every = 1000",
        "[[observables.series]]",
        'kind = "state_counts"',
        "[[observables.series]]",
        'kind = "event_rates"',
    )
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    text = render_deck()
    if args.output:
        args.output.write_text(text, encoding="utf-8")
    else:
        print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
