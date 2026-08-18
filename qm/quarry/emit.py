"""Emit computed rates in the formats the KMC engines consume.

Two targets (HANDOFF.md §2):

- **Petra deck fragments** — TOML ``[[reactions]]`` entries with Eyring or
  Arrhenius parameters and ``by_count`` ΔEa modifier tables, each number
  annotated with provenance (method, geometry hash, date). The primary
  target; emitted in kJ/mol and saying so is the caller's job via the
  deck's ``[deck] units`` (see ``splice_into_deck``).
- **Legacy ``data.rxn``** — the dissertation model's input, regenerated
  verbatim-format for the Phase 3 comparison against the original
  parameterization.

The reaction *mechanism* (center/guards/effects) is deck knowledge that
quarry passes through untouched; quarry owns the numbers.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from quarry.store import Provenance

# --- minimal TOML writer for the deck-fragment shapes we emit -------------
# stdlib tomllib is read-only; a general writer dependency is not worth it
# for the small closed set of shapes below (petra's deck style is inline
# tables for selectors and dotted arrays-of-tables for reactions).


def _toml_value(v: object) -> str:
    match v:
        case bool():
            return "true" if v else "false"
        case int():
            return str(v)
        case float():
            return repr(v)
        case str():
            return '"' + v.replace("\\", "\\\\").replace('"', '\\"') + '"'
        case list():
            return "[" + ", ".join(_toml_value(x) for x in v) + "]"
        case dict():
            inner = ", ".join(f"{k} = {_toml_value(x)}" for k, x in v.items())
            return "{ " + inner + " }"
        case _:
            raise TypeError(f"cannot render {type(v).__name__} as TOML")


@dataclass(frozen=True)
class EyringRate:
    """ΔH‡ and ΔS‡ in deck units (kJ/mol and kJ mol^-1 K^-1 for us)."""

    dh: float
    ds: float

    def spec(self) -> dict:
        return {"eyring": {"dh": self.dh, "ds": self.ds}}


@dataclass(frozen=True)
class ArrheniusRate:
    prefactor: float
    ea: float

    def spec(self) -> dict:
        return {"arrhenius": {"prefactor": self.prefactor, "ea": self.ea}}


@dataclass(frozen=True)
class ByCountModifier:
    """A ``[[reactions.modifiers]]`` entry: ΔEa table by neighbor count."""

    select: dict
    dea: list[float]
    provenance: Provenance | None = None


@dataclass
class ReactionEmit:
    """One ``[[reactions]]`` entry: mechanism passed through, rate computed."""

    name: str
    center: dict
    rate: EyringRate | ArrheniusRate
    guards: list[dict] = field(default_factory=list)
    produces: list[str] = field(default_factory=list)
    modifiers: list[ByCountModifier] = field(default_factory=list)
    effects: list[dict] = field(default_factory=list)
    provenance: Provenance | None = None


def _provenance_comment(p: Provenance) -> str:
    return (
        f"# provenance: method={p.method} engine={p.engine}"
        f" geometry=sha256:{p.geometry_hash[:12]} date={p.date}"
    )


def petra_fragment(reactions: list[ReactionEmit]) -> str:
    """Render ``[[reactions]]`` entries as a TOML fragment.

    Energies are emitted as given — the surrounding deck's
    ``[deck] units`` declares them; quarry's convention is kJ/mol
    (``splice_into_deck`` refuses templates that say otherwise).
    """
    out: list[str] = []
    for r in reactions:
        if r.provenance:
            out.append(_provenance_comment(r.provenance))
        out.append("[[reactions]]")
        out.append(f"name = {_toml_value(r.name)}")
        out.append(f"center = {_toml_value(r.center)}")
        if r.guards:
            out.append(f"guards = {_toml_value(r.guards)}")
        if r.produces:
            out.append(f"produces = {_toml_value(r.produces)}")
        out.append(f"rate = {_toml_value(r.rate.spec())}")
        for m in r.modifiers:
            out.append("")
            if m.provenance:
                out.append(_provenance_comment(m.provenance))
            out.append("[[reactions.modifiers]]")
            out.append(f"select = {_toml_value(m.select)}")
            out.append(f"by_count = {_toml_value({'dea': m.dea})}")
        for e in r.effects:
            out.append("")
            out.append("[[reactions.effects]]")
            for k, v in e.items():
                out.append(f"{k} = {_toml_value(v)}")
        out.append("")
    return "\n".join(out)


REACTIONS_MARKER = "# {QUARRY_REACTIONS}"


def splice_into_deck(template: str, reactions: list[ReactionEmit]) -> str:
    """Replace the marker line in a template deck with emitted reactions.

    The template must declare ``units = "kJ/mol"`` — quarry emits kJ/mol
    numbers and refuses a deck that would silently reinterpret them.
    """
    if REACTIONS_MARKER not in template:
        raise ValueError(f"template deck has no '{REACTIONS_MARKER}' marker line")
    if 'units = "kJ/mol"' not in template:
        raise ValueError('template deck must declare units = "kJ/mol"')
    return template.replace(REACTIONS_MARKER, petra_fragment(reactions))


# --- legacy data.rxn ------------------------------------------------------


@dataclass(frozen=True)
class RxnBlock:
    """One legacy reaction block: (k, dE-kcal/mol) pairs by environment.

    ``to_state = None`` renders a desorption block (single state column).
    """

    from_state: int
    comment: str
    pairs: list[tuple[float, float]]
    to_state: int | None = None


@dataclass(frozen=True)
class AdsorptionBlock:
    """Legacy adsorption block: prefactor + chemical-potential term."""

    from_state: int
    to_state: int
    comment: str
    prefactor: float
    mu: float


def _pair_lines(pairs: list[tuple[float, float]], per_line: int) -> list[str]:
    lines = []
    for i in range(0, len(pairs), per_line):
        chunk = pairs[i : i + per_line]
        lines.append("\t".join(f"{k:g} {de:g}" for k, de in chunk))
    return lines


def data_rxn(
    temperature: float,
    dmu_si: float,
    dmu_al: float,
    hydrolysis: list[RxnBlock],
    adsorption: list[AdsorptionBlock],
    desorption: list[RxnBlock],
    *,
    pairs_per_line: int = 10,
    provenance: Provenance | None = None,
) -> str:
    """Render a complete legacy ``data.rxn`` (format of legacy/cpp-model)."""
    out: list[str] = []
    out.append(f"{temperature:g}   # temperature in kelvin")
    out.append(f"{dmu_si:g}    # delta chemical potential for si")
    out.append(f"{dmu_al:g}    # delta chemical potential for al")
    if provenance:
        out.append(_provenance_comment(provenance))
    out.append("# hydrolysis reactions " + "*" * 41)
    out.append("# hydrolysis rxn data in (k+ dE) pairs")
    out.append("# k in relative units, dE in kcal / mole")
    for b in hydrolysis:
        if b.to_state is None:
            raise ValueError(f"hydrolysis block '{b.comment}' needs a to_state")
        out.append(f"{b.from_state}\t{b.to_state}\t# {b.comment}")
        out.append(str(len(b.pairs)))
        out.extend(_pair_lines(b.pairs, pairs_per_line))
        out.append("")
    out.append("# adsorption reactions " + "*" * 41)
    for a in adsorption:
        out.append(f"{a.from_state}\t{a.to_state}\t# {a.comment}")
        out.append(f"{a.prefactor:g}   {a.mu:g}")
        out.append("")
    out.append("# desorption reactions " + "*" * 41)
    for b in desorption:
        if b.to_state is not None:
            raise ValueError(f"desorption block '{b.comment}' takes no to_state")
        out.append(f"{b.from_state}\t\t# {b.comment}")
        out.append(str(len(b.pairs)))
        out.extend(_pair_lines(b.pairs, pairs_per_line))
        out.append("")
    return "\n".join(out)
