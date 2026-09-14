#!/usr/bin/env python3
"""Independent raw-to-derived verifier for the A9b sensitivity campaign.

This module intentionally does not import the campaign analyzer. It separately
regenerates decks, replays likelihood/event/lineage/state evidence, and
recomputes every derived byte.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import importlib.util
import json
import math
import random
import re
import sys
import tempfile
from collections import defaultdict
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import tomllib

try:
    from a9b_artifacts import (
        canonical_json_bytes,
        read_jsonl,
        remove_file_durable,
        sha256_file,
        write_json_atomic,
    )
except ModuleNotFoundError:
    _artifact_spec = importlib.util.spec_from_file_location(
        "a9b_artifacts", Path(__file__).with_name("a9b_artifacts.py")
    )
    assert _artifact_spec and _artifact_spec.loader
    _artifact_module = importlib.util.module_from_spec(_artifact_spec)
    _artifact_spec.loader.exec_module(_artifact_module)
    canonical_json_bytes = _artifact_module.canonical_json_bytes
    read_jsonl = _artifact_module.read_jsonl
    remove_file_durable = _artifact_module.remove_file_durable
    sha256_file = _artifact_module.sha256_file
    write_json_atomic = _artifact_module.write_json_atomic

SEEDS = (90401, 90403, 90407, 90409, 90413, 90419, 90421, 90427)
FAMILIES = (
    "siloxane-neutral",
    "sioal-si-neutral",
    "sioal-al-neutral",
    "connectivity-ladder",
    "al-o-al-analogue",
    "adsorption",
    "cation-desorption",
)
PERTURBATIONS = ("ea-minus-3", "ea-plus-3", "prefactor-x0.1", "prefactor-x10")
TARGETS = {
    "siloxane-neutral": ("R0-sio-si-hydrolysis", "R1-sio-si-condensation"),
    "sioal-si-neutral": (
        "R2-sioal2-si-hydrolysis",
        "R3-sioal2-si-condensation",
        "R10-sial-hydrolysis",
        "R11-sial-condensation",
    ),
    "sioal-al-neutral": (
        "R4-sioal2-al-hydrolysis",
        "R5a-sioal2-al-condensation",
        "R5b-sioal2-al-condensation",
        "R12-albr-hydrolysis",
        "R13-albr-condensation",
    ),
    "connectivity-ladder": (
        "R6-second-stage-hydrolysis",
        "R7-second-stage-condensation",
        "R8a-sioha-hydrolysis",
        "R8b-sioha-hydrolysis",
        "R9-sioha-condensation",
    ),
    "al-o-al-analogue": ("R14-alohal-hydrolysis", "R15-alohal-condensation"),
    "adsorption": ("adsorb-al", "adsorb-si"),
    "cation-desorption": ("desorb-al", "desorb-si"),
}
REACTION_BLOCK = re.compile(
    r"(?ms)^\[\[reactions\]\]\n(?P<body>.*?)(?=^\[\[reactions\]\]\n|\Z)"
)
RATE = re.compile(
    r"rate = \{ arrhenius = \{ prefactor = (?P<prefactor>[^,]+), ea = (?P<ea>[^}]+) \} \}"
)
AVOGADRO = 6.02214076e23
A2_TO_M2 = 1e-20
WEIGHT_GATE = 0.1
CONTRIBUTION_GATE = 2.0
CONTRIBUTION_GATE_ROLE = "analysis_precision_ranking_gate_not_sampling_contract"
STATIONARITY_LIMIT = 0.05
BOOTSTRAPS = 2000
EQUIVALENCE = (-0.1, 0.1)
TOP_K = 3
PACKAGE_FILES = {
    "metadata.json",
    "initial.pgif.json",
    "events.jsonl",
    "likelihood.jsonl",
    "checkpoints.jsonl",
    "final.pgif.json",
    "receipt.json",
}
PETRA = Path(__file__).resolve().parents[1]
REPO = PETRA.parent
ORIGINAL_BASE = PETRA / "examples" / "kaolinite-approx.toml"
ORIGINAL_CONTRACT = PETRA / "examples" / "kaolinite-reservoirs.toml"
ORIGINAL_PATH_REPORT = (
    REPO / "docs/program/results/a9b-mechanism-reachability/path-report.json"
)
ORIGINAL_SNAPSHOT = ORIGINAL_PATH_REPORT.with_name("initial-snapshot.pgif.json")


def _scenarios() -> list[tuple[str, str | None, str | None]]:
    rows = [("nominal", None, None)]
    rows.extend(
        (f"{family}__{perturbation}", family, perturbation)
        for family in FAMILIES
        for perturbation in PERTURBATIONS
    )
    if len(rows) != 29 or len({row[0] for row in rows}) != 29:
        raise AssertionError("independent scenario table drifted")
    return rows


def _load_reservoir_module() -> Any:
    path = Path(__file__).with_name("reservoir_origin_contract.py")
    spec = importlib.util.spec_from_file_location("a9b_verifier_reservoir", path)
    if not spec or not spec.loader:
        raise RuntimeError("cannot load reservoir contract implementation")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _independent_perturb(nominal: str, family: str, perturbation: str) -> str:
    wanted = set(TARGETS[family])
    seen: set[str] = set()

    def transform(match: re.Match[str]) -> str:
        block = match.group(0)
        name_match = re.search(r'(?m)^name = "([^"]+)"$', block)
        if name_match is None or name_match.group(1) not in wanted:
            return block
        name = name_match.group(1)
        rate_match = RATE.search(block)
        if rate_match is None:
            raise ValueError(f"{name}: missing canonical Arrhenius rate")
        prefactor = float(rate_match.group("prefactor"))
        barrier = float(rate_match.group("ea"))
        if perturbation == "ea-minus-3":
            barrier -= 3.0
        elif perturbation == "ea-plus-3":
            barrier += 3.0
        elif perturbation == "prefactor-x0.1":
            prefactor *= 0.1
        elif perturbation == "prefactor-x10":
            prefactor *= 10.0
        else:
            raise ValueError("unknown perturbation")
        seen.add(name)
        replacement = (
            "rate = { arrhenius = { "
            f"prefactor = {prefactor:.12g}, ea = {barrier:.12g}"
            " } }"
        )
        return block[: rate_match.start()] + replacement + block[rate_match.end() :]

    generated = REACTION_BLOCK.sub(transform, nominal)
    if seen != wanted:
        raise ValueError(f"independent perturbation missed reactions for {family}")
    before = tomllib.loads(nominal)
    after = tomllib.loads(generated)
    for left, right in zip(before["reactions"], after["reactions"], strict=True):
        name = left["name"]
        left_copy = dict(left)
        right_copy = dict(right)
        left_rate = left_copy.pop("rate")
        right_rate = right_copy.pop("rate")
        if left_copy != right_copy:
            raise ValueError(f"scenario changed non-rate semantics for {name}")
        if name not in wanted and left_rate != right_rate:
            raise ValueError(f"scenario changed unrelated rate for {name}")
    return generated


def _verify_and_regenerate_decks(root: Path, manifest: Mapping[str, Any]) -> None:
    inputs = manifest["inputs"]
    source_bindings = (
        (ORIGINAL_BASE, root / inputs["base_deck"], inputs["base_deck_sha256"]),
        (
            ORIGINAL_CONTRACT,
            root / inputs["reservoir_contract"],
            inputs["reservoir_contract_sha256"],
        ),
        (
            ORIGINAL_PATH_REPORT,
            root / inputs["path_report"],
            inputs["path_report_sha256"],
        ),
        (
            ORIGINAL_SNAPSHOT,
            root / inputs["lineage_snapshot"],
            inputs["lineage_snapshot_sha256"],
        ),
    )
    for original, bundled, digest in source_bindings:
        if sha256_file(original) != digest or sha256_file(bundled) != digest:
            raise ValueError(f"source/bundle hash mismatch for {original.name}")
    ph4_deck = root / inputs["ph4_deck"]
    ph4_evidence = root / inputs["ph4_evidence"]
    if (
        sha256_file(ph4_deck) != inputs["ph4_deck_sha256"]
        or sha256_file(ph4_evidence) != inputs["ph4_evidence_sha256"]
    ):
        raise ValueError("pH4 deck/evidence manifest hash mismatch")
    report = json.loads((root / inputs["path_report"]).read_text(encoding="utf-8"))
    targets = report.get("targets", [])
    target_sites = [row.get("site") for row in targets]
    snapshot = _pgif(root / inputs["lineage_snapshot"])
    snapshot_names = snapshot["states"]
    snapshot_states = snapshot["state_data"]
    snapshot_frozen = snapshot["frozen"]
    if (
        report.get("schema") != "a9b-mechanism-reachability-v2"
        or report.get("temperature_kelvin") != 298.0
        or report.get("target_count") != 320
        or report.get("reachable_count") != 260
        or report.get("unreachable_count") != 60
        or report.get("snapshot_sha256") != inputs["lineage_snapshot_sha256"]
        or snapshot["step"] != 0
        or snapshot["time"] != 0.0
        or len(targets) != 320
        or len(set(target_sites)) != 320
        or sum(row.get("reachable") is True for row in targets) != 260
        or sum(row.get("reachable") is False for row in targets) != 60
        or any(row.get("kind") not in {"Si", "Al"} for row in targets)
        or any(
            type(row.get("site")) is not int
            or not 0 <= row["site"] < len(snapshot_states)
            or snapshot_frozen[row["site"]]
            or row.get("origin") != "original-lattice"
            or row.get("initial_state") != snapshot_names[snapshot_states[row["site"]]]
            or not row["initial_state"].startswith(f"{row['kind']}.")
            or (
                row.get("reachable") is True
                and (
                    not isinstance(row.get("legal_path"), dict)
                    or type(row.get("release_event_index")) is not int
                    or row.get("missing_predicates") != []
                )
            )
            or (
                row.get("reachable") is False
                and (
                    row.get("legal_path") is not None
                    or row.get("release_event_index") is not None
                    or not isinstance(row.get("structural_upper_bound"), dict)
                    or not row.get("missing_predicates")
                )
            )
            for row in targets
        )
    ):
        raise ValueError("independent reachability census/path binding failed")
    with tempfile.TemporaryDirectory() as directory:
        temporary = Path(directory)
        regenerated_deck = temporary / "kaolinite-ph4.toml"
        regenerated_evidence = temporary / "kaolinite-ph4-evidence.json"
        reservoir = _load_reservoir_module()
        reservoir.materialize(
            root / inputs["base_deck"],
            root / inputs["reservoir_contract"],
            4,
            regenerated_deck,
            regenerated_evidence,
        )
        if regenerated_deck.read_bytes() != (root / inputs["ph4_deck"]).read_bytes():
            raise ValueError("independently regenerated pH4 deck differs")
        if (
            regenerated_evidence.read_bytes()
            != (root / inputs["ph4_evidence"]).read_bytes()
        ):
            raise ValueError("independently regenerated pH4 evidence differs")
    nominal = (root / inputs["ph4_deck"]).read_text(encoding="utf-8")
    records = manifest["scenarios"]
    if [record["name"] for record in records] != [row[0] for row in _scenarios()]:
        raise ValueError("manifest scenario order/identity mismatch")
    for record, (name, family, perturbation) in zip(records, _scenarios(), strict=True):
        expected = (
            nominal
            if family is None
            else _independent_perturb(nominal, family, perturbation)
        ).encode()
        path = root / record["deck"]
        if path.read_bytes() != expected or sha256_file(path) != record["deck_sha256"]:
            raise ValueError(f"independently regenerated deck mismatch for {name}")


def _number(value: object, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError(f"{label} must be numeric")
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"{label} must be finite")
    return number


def _same(left: float, right: float) -> bool:
    return math.isclose(left, right, rel_tol=2e-12, abs_tol=1e-18)


def _pgif(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if set(payload) != {"pgif", "meta", "nodes", "edges"} or payload.get("pgif") != 1:
        raise ValueError("not canonical PGIF v1")
    meta = payload["meta"]
    if (
        set(meta) != {"directed", "kind", "petra", "producer"}
        or meta["directed"] is not False
        or meta["kind"] != "kmc-lattice"
        or meta["producer"] != "petra"
        or set(meta["petra"])
        != {"deck", "state_types", "states", "step", "temperature", "time"}
    ):
        raise ValueError("PGIF graph metadata is not canonical Petra metadata")
    petra = meta["petra"]
    states = petra["states"]
    state_types = petra["state_types"]
    step = petra["step"]
    time = _number(petra["time"], "PGIF time")
    _number(petra["temperature"], "PGIF temperature")
    if (
        not isinstance(petra["deck"], str)
        or not petra["deck"]
        or not isinstance(states, list)
        or not states
        or any(not isinstance(value, str) for value in states)
        or len(states) != len(set(states))
        or not isinstance(state_types, list)
        or len(state_types) != len(states)
        or any(not isinstance(value, str) for value in state_types)
        or type(step) is not int
        or step < 0
        or time < 0
    ):
        raise ValueError("PGIF Petra metadata values are malformed")

    nodes = payload["nodes"]
    if set(nodes) != {"count", "columns"} or type(nodes["count"]) is not int:
        raise ValueError("PGIF node table is malformed")
    count = nodes["count"]
    columns = nodes["columns"]
    expected_node_types = {
        "frozen": "bool",
        "kind": "categorical",
        "state": "categorical",
        "type": "categorical",
        "x": "f32",
        "y": "f32",
        "z": "f32",
    }
    if count < 0 or set(columns) != set(expected_node_types):
        raise ValueError("PGIF node columns are not canonical")
    for name, expected_type in expected_node_types.items():
        column = columns[name]
        expected_keys = (
            {"type", "data", "dict"}
            if expected_type == "categorical"
            else {"type", "data"}
        )
        if (
            not isinstance(column, dict)
            or set(column) != expected_keys
            or column.get("type") != expected_type
            or not isinstance(column.get("data"), list)
            or len(column["data"]) != count
        ):
            raise ValueError(f"PGIF {name} column is malformed")
        if expected_type == "categorical":
            dictionary = column["dict"]
            if (
                not isinstance(dictionary, list)
                or not dictionary
                or any(not isinstance(value, str) for value in dictionary)
                or len(dictionary) != len(set(dictionary))
                or any(
                    type(value) is not int or not 0 <= value < len(dictionary)
                    for value in column["data"]
                )
            ):
                raise ValueError(f"PGIF {name} categorical values are malformed")
    if (
        columns["state"]["dict"] != states
        or any(type(value) is not bool for value in columns["frozen"]["data"])
        or any(
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(float(value))
            for name in ("x", "y", "z")
            for value in columns[name]["data"]
        )
    ):
        raise ValueError("PGIF state/frozen/coordinate values are malformed")
    state_data = columns["state"]["data"]
    type_dict = columns["type"]["dict"]
    kind_dict = columns["kind"]["dict"]
    for site, state_id in enumerate(state_data):
        state_name = states[state_id]
        if (
            state_types[state_id] not in type_dict
            or columns["type"]["data"][site] != type_dict.index(state_types[state_id])
            or "." not in state_name
            or state_name.partition(".")[0] not in kind_dict
            or columns["kind"]["data"][site]
            != kind_dict.index(state_name.partition(".")[0])
        ):
            raise ValueError(f"PGIF type/kind disagrees with state at node {site}")

    edges = payload["edges"]
    if (
        set(edges) != {"count", "src", "dst", "columns"}
        or type(edges["count"]) is not int
    ):
        raise ValueError("PGIF edge table is malformed")
    edge_count = edges["count"]
    seam = edges["columns"].get("seam") if isinstance(edges["columns"], dict) else None
    if (
        edge_count < 0
        or set(edges["columns"]) != {"seam"}
        or not isinstance(seam, dict)
        or set(seam) != {"type", "data"}
        or seam.get("type") != "bool"
        or any(
            not isinstance(values, list) or len(values) != edge_count
            for values in (edges["src"], edges["dst"], seam["data"])
        )
        or any(type(value) is not bool for value in seam["data"])
        or any(
            type(value) is not int or not 0 <= value < count
            for values in (edges["src"], edges["dst"])
            for value in values
        )
        or any(
            source == destination
            for source, destination in zip(edges["src"], edges["dst"], strict=True)
        )
    ):
        raise ValueError("PGIF edges are malformed")
    return {
        "payload": payload,
        "states": states,
        "state_data": state_data,
        "frozen": columns["frozen"]["data"],
        "step": step,
        "time": time,
    }


def _pgif_static(document: Mapping[str, Any]) -> dict[str, Any]:
    payload = document["payload"]
    petra = payload["meta"]["petra"]
    columns = payload["nodes"]["columns"]
    return {
        "meta": {
            **payload["meta"],
            "petra": {
                key: value
                for key, value in petra.items()
                if key not in {"step", "time"}
            },
        },
        "node_count": payload["nodes"]["count"],
        "static_node_columns": {
            name: column
            for name, column in columns.items()
            if name not in {"state", "type"}
        },
        "state_dict_and_schema": {
            key: value for key, value in columns["state"].items() if key != "data"
        },
        "type_dict_and_schema": {
            key: value for key, value in columns["type"].items() if key != "data"
        },
        "edges": payload["edges"],
    }


def _check_package_hashes(
    root: Path,
    scenario: str,
    seed: int,
    manifest: Mapping[str, Any],
) -> tuple[Path, dict[str, Any]]:
    external_path = root / "receipts" / scenario / f"seed-{seed}.json"
    external = json.loads(external_path.read_text(encoding="utf-8"))
    scenario_row = next(row for row in manifest["scenarios"] if row["name"] == scenario)
    expected_package = Path("runs") / scenario / f"seed-{seed}"
    if (
        external.get("schema") != "a9b-job-receipt-v2"
        or external.get("scenario") != scenario
        or external.get("seed") != seed
        or Path(external.get("package", "")) != expected_package
        or external.get("deck_sha256") != scenario_row["deck_sha256"]
        or external.get("runner_sha256") != manifest.get("runner_sha256")
        or sha256_file(root / external["log"]) != external.get("log_sha256")
    ):
        raise ValueError("external run receipt identity/input mismatch")
    package = root / external["package"]
    if {path.name for path in package.iterdir()} != PACKAGE_FILES:
        raise ValueError("replica package file set mismatch")
    if set(external["package_files_sha256"]) != PACKAGE_FILES:
        raise ValueError("external receipt hash inventory mismatch")
    for name, digest in external["package_files_sha256"].items():
        if sha256_file(package / name) != digest:
            raise ValueError(f"external receipt hash mismatch: {name}")
    internal = json.loads((package / "receipt.json").read_text(encoding="utf-8"))
    if (
        internal.get("schema") != "a9b-replica-receipt-v2"
        or internal.get("seed") != seed
        or set(internal.get("files_sha256", {})) != PACKAGE_FILES - {"receipt.json"}
    ):
        raise ValueError("internal receipt identity/inventory mismatch")
    for name, digest in internal["files_sha256"].items():
        if sha256_file(package / name) != digest:
            raise ValueError(f"internal receipt hash mismatch: {name}")
    return package, external


def _audit_replica(root: Path, scenario: str, seed: int) -> dict[str, Any]:
    manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
    package, external = _check_package_hashes(root, scenario, seed, manifest)
    meta = json.loads((package / "metadata.json").read_text(encoding="utf-8"))
    if meta.get("schema") != "a9b-replica-package-v2" or meta.get("seed") != seed:
        raise ValueError("replica metadata mismatch")
    if meta["census"] != {
        "whole_finite_deck_original_centers": 320,
        "reachable_original_centers": 260,
        "topology_no_go_original_centers": 60,
        "whole_deck_estimand_includes_topology_no_go": True,
    }:
        raise ValueError("replica census contract mismatch")
    scenario_row = next(row for row in manifest["scenarios"] if row["name"] == scenario)
    inputs = manifest["inputs"]
    contract = manifest["run_contract"]
    bound = meta.get("inputs", {})
    command = external.get("command")
    expected_argv = [
        str((root / scenario_row["deck"]).resolve()),
        str((root / inputs["path_report"]).resolve()),
        str((root / inputs["lineage_snapshot"]).resolve()),
        str(package.resolve()),
        "--seed",
        str(seed),
        "--horizon",
        str(contract["horizon_s"]),
        "--checkpoints",
        str(contract["checkpoints"]),
        "--bias-factor",
        str(contract["bias_factor"]),
        "--max-events",
        str(contract["max_events"]),
    ]
    if (
        meta.get("survey_tier_platform_test") is not True
        or meta.get("not_production") is not True
        or meta.get("temperature_k") != 298.0
        or meta.get("horizon_s") != contract["horizon_s"]
        or meta.get("fixed_checkpoint_count") != contract["checkpoints"]
        or meta.get("bias_factor") != contract["bias_factor"]
        or meta.get("max_events") != contract["max_events"]
        or bound.get("deck_sha256") != scenario_row["deck_sha256"]
        or bound.get("path_report_sha256") != inputs["path_report_sha256"]
        or bound.get("lineage_snapshot_sha256") != inputs["lineage_snapshot_sha256"]
        or bound.get("runner_executable_sha256") != external.get("runner_sha256")
        or not isinstance(command, list)
        or len(command) != 18
        or command[:3] != ["nice", "-n", "10"]
        or command[3] != bound.get("runner_executable")
        or meta.get("argv") != expected_argv
        or command[4:] != expected_argv
    ):
        raise ValueError("replica metadata/command/input binding mismatch")
    reactions = [row["name"] for row in meta["reaction_mapping"]]
    bound_reactions = [
        row["name"]
        for row in tomllib.loads(
            (root / scenario_row["deck"]).read_text(encoding="utf-8")
        )["reactions"]
    ]
    if reactions != bound_reactions:
        raise ValueError("reaction mapping differs from bound scenario deck")
    bias = []
    for expected_id, row in enumerate(meta["reaction_mapping"]):
        factor = _number(row["bias_factor"], "reaction bias factor")
        expected_factor = (
            contract["bias_factor"]
            if "hydrolysis" in bound_reactions[expected_id]
            or bound_reactions[expected_id].startswith("desorb-")
            else 1.0
        )
        if (
            row.get("id") != expected_id
            or factor <= 0
            or not _same(factor, expected_factor)
        ):
            raise ValueError("reaction mapping/bias must match the bound deck")
        bias.append(factor)

    likelihood = read_jsonl(package / "likelihood.jsonl")
    if (
        not likelihood
        or likelihood[0].get("schema") != "a9b-likelihood-stream-v2"
        or likelihood[0].get("seed") != seed
    ):
        raise ValueError("likelihood header mismatch")
    segments = likelihood[1:]
    if len(segments) != meta.get("likelihood_segment_count"):
        raise ValueError("likelihood segment count mismatch")
    total_log = 0.0
    last_time = 0.0
    fired_segments: dict[int, int] = {}
    for expected_index, segment in enumerate(segments):
        if segment.get("segment") != expected_index:
            raise ValueError("likelihood segment identity mismatch")
        start = _number(segment["start_time_s"], "segment start")
        end = _number(segment["end_time_s"], "segment end")
        physical = _number(segment["physical_total_rate_s-1"], "physical rate")
        tilted = _number(segment["biased_total_rate_s-1"], "biased rate")
        if not _same(start, last_time) or end < start or physical < 0 or tilted < 0:
            raise ValueError("likelihood time/rate sequence invalid")
        fired = segment["fired_reaction_id"]
        event_penalty = 0.0
        if fired is not None:
            if physical <= 0 or tilted <= 0:
                raise ValueError(
                    "fired likelihood segments require positive physical and biased rates"
                )
            if type(fired) is not int or not 0 <= fired < len(reactions):
                raise ValueError("likelihood fired reaction invalid")
            factor = _number(segment["fired_bias_factor"], "fired factor")
            if factor <= 0 or not _same(factor, bias[fired]):
                raise ValueError("likelihood factor does not match reaction mapping")
            event_penalty = math.log(factor)
            fired_segments[expected_index] = fired
        elif segment["fired_bias_factor"] is not None:
            raise ValueError("deadline segment contains event factor")
        independently_computed = (tilted - physical) * (end - start) - event_penalty
        if not _same(independently_computed, segment["log_likelihood_increment"]):
            raise ValueError("likelihood increment tampering detected")
        total_log += independently_computed
        if not _same(total_log, segment["cumulative_log_likelihood"]):
            raise ValueError("likelihood cumulative tampering detected")
        last_time = end
    if not _same(total_log, meta["final_log_likelihood"]):
        raise ValueError("final likelihood summary mismatch")

    initial_pgif = _pgif(package / "initial.pgif.json")
    final_pgif = _pgif(package / "final.pgif.json")
    state_names = initial_pgif["states"]
    initial = initial_pgif["state_data"]
    final = final_pgif["state_data"]
    final_step = final_pgif["step"]
    final_time = final_pgif["time"]
    if (
        initial_pgif["step"] != 0
        or initial_pgif["time"] != 0.0
        or _pgif_static(initial_pgif) != _pgif_static(final_pgif)
    ):
        raise ValueError("PGIF initial/final canonical identity mismatch")
    source_pgif = _pgif(root / inputs["lineage_snapshot"])
    if initial_pgif["payload"] != source_pgif["payload"]:
        raise ValueError("initial PGIF differs from full bound lineage source")
    events = read_jsonl(package / "events.jsonl")
    if (
        not events
        or events[0].get("schema") != "a9b-event-stream-v2"
        or events[0].get("seed") != seed
        or events[0]["states"] != state_names
        or events[0]["reactions"] != reactions
        or len(events) - 1 != meta.get("event_count")
    ):
        raise ValueError("event stream header/count mismatch")
    state_occupants = events[0]["state_occupants"]
    if len(state_occupants) != len(state_names):
        raise ValueError("event state occupant table malformed")
    path_report = json.loads((root / inputs["path_report"]).read_text())
    ledger = {row["site"]: row["kind"] for row in path_report["targets"]}
    state = list(initial)
    releases = {"Al": 0, "Si": 0}
    used_segments: set[int] = set()
    prior_event_time = 0.0
    checkpoints = read_jsonl(package / "checkpoints.jsonl")
    fixed_count = meta["fixed_checkpoint_count"]
    horizon = _number(meta["horizon_s"], "horizon")
    if (
        not checkpoints
        or len(checkpoints) != meta.get("recorded_checkpoint_count")
        or (meta["complete"] and len(checkpoints) != fixed_count)
    ):
        raise ValueError("checkpoint stream count contract invalid")
    event_position = 0
    prior_checkpoint_time = -1.0
    kinds = {name.partition(".")[0] for name in state_names}
    immutable_kind_totals = {
        kind: sum(
            count
            for name, count in checkpoints[0]["state_counts"].items()
            if name.startswith(f"{kind}.")
        )
        for kind in kinds
    }
    for checkpoint_number, checkpoint in enumerate(checkpoints):
        checkpoint_time = _number(checkpoint["actual_time_s"], "checkpoint time")
        target_time = _number(checkpoint["target_time_s"], "checkpoint target")
        role = checkpoint.get("checkpoint_role")
        if role == "fixed_physical_time":
            expected = horizon * checkpoint_number / (fixed_count - 1)
            if checkpoint_number >= fixed_count or not _same(target_time, expected):
                raise ValueError("fixed checkpoint cadence mismatch")
        elif role == "terminal_censor_snapshot":
            if meta["complete"] or checkpoint_number != len(checkpoints) - 1:
                raise ValueError("terminal censor checkpoint invalid")
        else:
            raise ValueError("checkpoint role invalid")
        if (
            checkpoint.get("checkpoint") != checkpoint_number
            or not _same(checkpoint_time, target_time)
            or checkpoint_time <= prior_checkpoint_time
            or checkpoint.get("kind_totals") != immutable_kind_totals
        ):
            raise ValueError("checkpoint identity/time/kind census invalid")
        expected_log = 0.0
        has_boundary = checkpoint_time == 0.0
        for segment in segments:
            end = _number(segment["end_time_s"], "segment end")
            if end <= checkpoint_time or _same(end, checkpoint_time):
                expected_log = _number(
                    segment["cumulative_log_likelihood"], "segment cumulative"
                )
            has_boundary |= _same(end, checkpoint_time)
        if not has_boundary or not _same(
            expected_log, checkpoint["cumulative_log_likelihood"]
        ):
            raise ValueError("checkpoint likelihood boundary mismatch")
        while event_position + 1 < len(events):
            event = events[event_position + 1]
            event_time = _number(event["time_s"], "event time")
            if event_time > checkpoint_time:
                break
            if (
                event.get("event") != event_position
                or event.get("step") != event_position + 1
                or event_time < prior_event_time
            ):
                raise ValueError("event identity/timestamp sequence invalid")
            reaction_id = event["reaction_id"]
            segment_id = event["likelihood_segment"]
            if (
                reactions[reaction_id] != event["reaction"]
                or fired_segments.get(segment_id) != reaction_id
                or segment_id in used_segments
                or not _same(segments[segment_id]["end_time_s"], event_time)
            ):
                raise ValueError("event/likelihood/reaction mapping mismatch")
            used_segments.add(segment_id)
            center = event["center_site"]
            if type(center) is not int or not 0 <= center < len(state):
                raise ValueError("event center site invalid")
            center_before = state[center]
            for change in event["changes"]:
                site = change["site"]
                old = change["old_state_id"]
                new = change["new_state_id"]
                if (
                    state[site] != old
                    or state_names[old] != change["old_state"]
                    or state_names[new] != change["new_state"]
                ):
                    raise ValueError("event/PGIF state tampering detected")
                state[site] = new
            lineage_release = None
            if event["reaction"] in ("desorb-al", "desorb-si"):
                species = event["reaction"].split("-", 1)[1].capitalize()
                if state_occupants[center_before] != species or center not in {
                    change["site"] for change in event["changes"]
                }:
                    raise ValueError("desorption center/occupant semantics invalid")
                if ledger.get(center) == species:
                    lineage_release = species
                    del ledger[center]
                    releases[species] += 1
            if event["lattice_release_kind"] != lineage_release:
                raise ValueError("lineage annotation tampering detected")
            prior_event_time = event_time
            event_position += 1
        counts = {name: 0 for name in state_names}
        for state_id in state:
            counts[state_names[state_id]] += 1
        if (
            counts != checkpoint["state_counts"]
            or checkpoint.get("step") != event_position
        ):
            raise ValueError("checkpoint state-count/step tampering detected")
        if releases != checkpoint["original_lattice_releases"]:
            raise ValueError("checkpoint lineage-count tampering detected")
        prior_checkpoint_time = checkpoint_time
    if event_position != len(events) - 1 or state != final:
        raise ValueError("final PGIF does not match event replay")
    if set(fired_segments) != used_segments:
        raise ValueError("fired likelihood segments are missing events")
    if final_step != len(events) - 1 or not _same(final_time, meta["final_time_s"]):
        raise ValueError("final PGIF timestamp/step mismatch")
    if (
        not _same(last_time, final_time)
        or not _same(prior_checkpoint_time, final_time)
        or meta.get("remaining_original_ion_ledger") != len(ledger)
    ):
        raise ValueError("terminal likelihood/checkpoint/ledger mismatch")
    if meta["complete"]:
        if meta["stop_reason"] != "physical_horizon" or not _same(final_time, horizon):
            raise ValueError("complete run did not reach physical horizon")
    elif meta["stop_reason"] == "physical_horizon" or final_time > horizon:
        raise ValueError("censored run stop contract invalid")
    if releases != meta["original_lattice_releases"]:
        raise ValueError("lineage summary mismatch")

    times = [_number(row["actual_time_s"], "checkpoint time") for row in checkpoints]
    areas = [_number(row["geometric_area_a2"], "area") for row in checkpoints]
    trapezoids = []
    for left, right, area_left, area_right in zip(times, times[1:], areas, areas[1:]):
        if right <= left or area_left <= 0 or area_right <= 0:
            raise ValueError("area-time samples invalid")
        trapezoids.append((right - left) * (area_left + area_right) / 2.0)
    exposure = math.fsum(trapezoids)
    denominator = AVOGADRO * exposure * A2_TO_M2
    rates = {species: releases[species] / denominator for species in ("Si", "Al")}
    reachable_denominator = denominator * 260 / 320
    reachable = {
        species: releases[species] / reachable_denominator for species in ("Si", "Al")
    }
    fixed_rows = [
        row for row in checkpoints if row["checkpoint_role"] == "fixed_physical_time"
    ]
    checkpoint_states = [row["state_counts"] for row in fixed_rows]
    checkpoint_log_weights = [
        _number(row["cumulative_log_likelihood"], "fixed checkpoint log weight")
        for row in fixed_rows
    ]
    checkpoint_times = [
        _number(row["actual_time_s"], "fixed checkpoint time") for row in fixed_rows
    ]
    return {
        "scenario": scenario,
        "seed": seed,
        "complete": meta["complete"],
        "stop_reason": meta["stop_reason"],
        "log_weight": total_log,
        "release_counts": releases,
        "whole_deck_rates": rates,
        "reachable_census_rates": reachable,
        "area_time_a2_s": exposure,
        "checkpoint_times_s": checkpoint_times,
        "checkpoint_state_counts": checkpoint_states,
        "checkpoint_log_weights": checkpoint_log_weights,
        "stationarity": _stationarity(checkpoint_states),
        "event_count": len(events) - 1,
        "segment_count": len(segments),
    }


def _log_total(values: Sequence[float]) -> float:
    maximum = max(values)
    return maximum + math.log(sum(math.exp(value - maximum) for value in values))


def _ess(values: Sequence[float]) -> float:
    if all(value == values[0] for value in values):
        return float(len(values))
    reference = max(values)
    normalized = [math.exp(value - reference) for value in values]
    normalized_total = sum(normalized)
    normalized_square_total = sum(weight * weight for weight in normalized)
    answer = normalized_total * normalized_total / normalized_square_total
    return min(float(len(values)), answer)


def _log_mean_unscaled(
    log_weights: Sequence[float], samples: Sequence[float]
) -> float | None:
    terms = [
        log_weight + math.log(sample)
        for log_weight, sample in zip(log_weights, samples, strict=True)
        if sample > 0
    ]
    return None if not terms else _log_total(terms) - math.log(len(samples))


def _linear_positive(log_value: float) -> tuple[float | None, str]:
    try:
        value = math.exp(log_value)
    except OverflowError:
        return None, "censored_numerical_overflow"
    if not math.isfinite(value):
        return None, "censored_numerical_overflow"
    if value == 0.0:
        return None, "censored_numerical_underflow"
    return value, "represented"


def _state_population_point(
    log_weights: Sequence[float], samples: Sequence[float]
) -> dict[str, Any]:
    log_mean = _log_mean_unscaled(log_weights, samples)
    if log_mean is None:
        return {"value": 0.0, "log_mean": None, "status": "represented"}
    value, status = _linear_positive(log_mean)
    return {"value": value, "log_mean": log_mean, "status": status}


def _percentile_log_band(
    log_weights: Sequence[float], samples: Sequence[float], seed: int
) -> list[float | None]:
    generator = random.Random(seed)
    size = len(samples)
    boot = []
    for _ in range(BOOTSTRAPS):
        indices = [generator.randrange(size) for _ in range(size)]
        boot.append(
            _log_mean_unscaled(
                [log_weights[index] for index in indices],
                [samples[index] for index in indices],
            )
        )
    boot.sort(key=lambda value: -math.inf if value is None else value)
    return [
        boot[round((BOOTSTRAPS - 1) * 0.025)],
        boot[round((BOOTSTRAPS - 1) * 0.975)],
    ]


def _estimate(
    log_weights: Sequence[float], samples: Sequence[float], seed: int
) -> dict[str, Any]:
    if any(not math.isfinite(value) for value in (*log_weights, *samples)):
        raise ValueError("estimator inputs must be finite")
    if any(sample < 0 for sample in samples):
        raise ValueError("estimator samples must be nonnegative")
    weight_ess = _ess(log_weights)
    weight_fraction = weight_ess / len(log_weights)
    result = {
        "estimator": "unnormalized_importance_sampling_log_domain",
        "replicas": len(samples),
        "weight_ess": weight_ess,
        "weight_ess_fraction": weight_fraction,
        "contribution_ess": None,
        "contribution_ess_gate_role": CONTRIBUTION_GATE_ROLE,
        "log_mean": None,
        "log_ci95": None,
        "mean": None,
        "ci95": None,
        "upper_95": None,
        "poisson_bound": "unavailable",
        "rankable": False,
    }
    if weight_fraction < WEIGHT_GATE:
        return {**result, "status": "censored_low_weight_ess"}
    if not any(sample > 0 for sample in samples):
        return {**result, "status": "censored_zero_biased_observations"}
    contribution_logs = [
        weight + math.log(sample)
        for weight, sample in zip(log_weights, samples, strict=True)
        if sample > 0
    ]
    contribution_ess = _ess(contribution_logs)
    result["contribution_ess"] = contribution_ess
    if contribution_ess < CONTRIBUTION_GATE:
        return {**result, "status": "censored_low_contribution_ess"}
    log_mean = _log_mean_unscaled(log_weights, samples)
    assert log_mean is not None
    log_interval = _percentile_log_band(log_weights, samples, seed)
    result["log_mean"] = log_mean
    result["log_ci95"] = log_interval
    mean, mean_status = _linear_positive(log_mean)
    interval: list[float] = []
    interval_statuses = []
    for endpoint in log_interval:
        if endpoint is None:
            interval.append(0.0)
            continue
        value, status = _linear_positive(endpoint)
        interval_statuses.append(status)
        if value is not None:
            interval.append(value)
    statuses = [mean_status, *interval_statuses]
    if "censored_numerical_overflow" in statuses:
        return {**result, "status": "censored_numerical_overflow"}
    if "censored_numerical_underflow" in statuses:
        return {**result, "status": "censored_numerical_underflow"}
    assert mean is not None and len(interval) == 2
    return {
        **result,
        "status": "estimated",
        "mean": mean,
        "ci95": interval,
        "rankable": True,
    }


def _enforce_scientific_gates(
    estimate: dict[str, Any], *, complete: bool, stationary: bool
) -> dict[str, Any]:
    if not complete:
        return {
            **estimate,
            "status": "censored_incomplete_replicas",
            "mean": None,
            "ci95": None,
            "upper_95": None,
            "rankable": False,
        }
    if estimate["rankable"] and not stationary:
        return {
            **estimate,
            "status": "censored_nonstationary_state_distribution",
            "mean": None,
            "ci95": None,
            "upper_95": None,
            "rankable": False,
        }
    return estimate


def _stationarity(
    evolution: Sequence[Mapping[str, float]],
    fixed_kind_totals: Mapping[str, float] | None = None,
) -> dict[str, Any]:
    if len(evolution) < 5:
        return {"passed": False, "status": "insufficient_checkpoints", "kinds": {}}
    by_kind: dict[str, list[str]] = defaultdict(list)
    for state in sorted(evolution[0]):
        by_kind[state.split(".", 1)[0]].append(state)
    output = {}
    all_pass = True
    for kind, states in by_kind.items():
        if fixed_kind_totals is None:
            totals = [sum(row[state] for state in states) for row in evolution]
        else:
            totals = [fixed_kind_totals[kind]] * len(evolution)
        if any(total <= 0 for total in totals):
            output[kind] = {"status": "invalid_zero_kind_total", "states": {}}
            all_pass = False
            continue
        state_rows = {}
        kind_pass = True
        for state in states:
            fractions = [
                row[state] / total for row, total in zip(evolution, totals, strict=True)
            ]
            tail = fractions[-5:]
            span = max(tail) - min(tail)
            endpoint = tail[-1] - tail[0]
            passed = span <= STATIONARITY_LIMIT and abs(endpoint) <= STATIONARITY_LIMIT
            kind_pass &= passed
            state_rows[state] = {
                "tail_fraction_range": span,
                "tail_endpoint_change": endpoint,
                "status": "stationary" if passed else "nonstationary",
            }
        output[kind] = {
            "status": "stationary" if kind_pass else "nonstationary",
            "states": state_rows,
        }
        all_pass &= kind_pass
    return {
        "passed": all_pass,
        "status": "stationary" if all_pass else "nonstationary",
        "tolerance_fraction": STATIONARITY_LIMIT,
        "kinds": output,
    }


def _seed(*parts: str) -> int:
    return int.from_bytes(hashlib.sha256("\0".join(parts).encode()).digest()[:8], "big")


def _paired(
    nominal: Mapping[int, tuple[float, float]],
    changed: Mapping[int, tuple[float, float]],
    nominal_rankable: bool,
    changed_rankable: bool,
    seed: int,
) -> dict[str, Any]:
    if set(nominal) != set(changed):
        return {
            "status": "incomplete_pairs",
            "same_seed_linear_deltas": [],
            "mean_linear_delta": None,
            "mean_delta_sign": None,
            "log_abs_mean_delta_mol_m2_s": None,
            "delta_log10": None,
            "delta_log10_ci95": None,
        }
    seeds = sorted(nominal)

    def log_contribution(pair: tuple[float, float]) -> float | None:
        log_weight = _number(pair[0], "paired log weight")
        sample = _number(pair[1], "paired sample")
        if sample < 0:
            raise ValueError("paired samples must be nonnegative")
        return None if sample == 0 else log_weight + math.log(sample)

    def signed_log_delta(
        left: float | None, right: float | None
    ) -> tuple[int, float | None]:
        if left is None and right is None:
            return 0, None
        if left is None:
            return 1, right
        if right is None:
            return -1, left
        if left == right:
            return 0, None
        high, low = (right, left) if right > left else (left, right)
        return (1 if right > left else -1), high + math.log(-math.expm1(low - high))

    def linear_if_representable(sign: int, log_abs: float | None) -> float | None:
        if sign == 0:
            return 0.0
        assert log_abs is not None
        try:
            magnitude = math.exp(log_abs)
        except OverflowError:
            return None
        return sign * magnitude if math.isfinite(magnitude) else None

    left_logs = [log_contribution(nominal[value]) for value in seeds]
    right_logs = [log_contribution(changed[value]) for value in seeds]
    signed_deltas = [
        signed_log_delta(left, right)
        for left, right in zip(left_logs, right_logs, strict=True)
    ]
    delta_logs = [value for sign, value in signed_deltas if sign and value is not None]
    if not delta_logs:
        mean_delta = 0.0
        mean_sign = 0
        mean_log_abs = None
    else:
        scale = max(delta_logs)
        scaled_sum = sum(
            sign * math.exp(log_abs - scale)
            for sign, log_abs in signed_deltas
            if sign and log_abs is not None
        )
        if scaled_sum == 0.0:
            mean_delta = 0.0
            mean_sign = 0
            mean_log_abs = None
        else:
            mean_sign = 1 if scaled_sum > 0 else -1
            mean_log_abs = scale + math.log(abs(scaled_sum)) - math.log(len(seeds))
            mean_delta = linear_if_representable(mean_sign, mean_log_abs)
    result = {
        "status": "unranked_censored",
        "same_seed_linear_deltas": [
            {
                "seed": value,
                "delta_mol_m2_s": linear_if_representable(sign, log_abs),
                "delta_sign": sign,
                "log_abs_delta_mol_m2_s": log_abs,
            }
            for value, (sign, log_abs) in zip(seeds, signed_deltas, strict=True)
        ],
        "mean_linear_delta": mean_delta,
        "mean_delta_sign": mean_sign,
        "log_abs_mean_delta_mol_m2_s": mean_log_abs,
        "delta_log10": None,
        "delta_log10_ci95": None,
    }

    def log_mean(
        values: Mapping[int, tuple[float, float]], chosen: Sequence[int]
    ) -> float | None:
        terms = [log_contribution(values[value]) for value in chosen]
        positive = [term for term in terms if term is not None]
        return None if not positive else _log_total(positive) - math.log(len(chosen))

    nominal_log_mean = log_mean(nominal, seeds)
    changed_log_mean = log_mean(changed, seeds)
    if (
        not nominal_rankable
        or not changed_rankable
        or nominal_log_mean is None
        or changed_log_mean is None
    ):
        return result
    result["status"] = "rankable"
    result["delta_log10"] = (changed_log_mean - nominal_log_mean) / math.log(10.0)
    generator = random.Random(seed)
    boot = []
    for _ in range(BOOTSTRAPS):
        chosen = [seeds[generator.randrange(len(seeds))] for _ in seeds]
        left = log_mean(nominal, chosen)
        right = log_mean(changed, chosen)
        if left is not None and right is not None:
            boot.append((right - left) / math.log(10.0))
    if len(boot) == BOOTSTRAPS:
        boot.sort()
        result["delta_log10_ci95"] = [
            boot[round((BOOTSTRAPS - 1) * 0.025)],
            boot[round((BOOTSTRAPS - 1) * 0.975)],
        ]
    return result


def _irrelevant(interval: Sequence[float] | None) -> bool:
    return bool(
        interval is not None
        and len(interval) == 2
        and EQUIVALENCE[0] <= interval[0] <= interval[1] <= EQUIVALENCE[1]
    )


def _lab(species_rates: Mapping[str, float | None]) -> dict[str, Any]:
    rates = {
        str(ph): 10**-11.31 * (10 ** (-ph)) ** 0.777 + 10**-13.18 for ph in (3, 4, 5)
    }
    gaps = {
        "ph3_bound_not_simulated": None,
        "ph4_simulated": None,
        "ph5_bound_not_simulated": None,
    }
    si_rate = species_rates.get("Si")
    al_rate = species_rates.get("Al")
    formula_rate = None
    if (
        isinstance(si_rate, (int, float))
        and not isinstance(si_rate, bool)
        and isinstance(al_rate, (int, float))
        and not isinstance(al_rate, bool)
        and si_rate > 0
        and al_rate > 0
        and math.isfinite(si_rate)
        and math.isfinite(al_rate)
    ):
        formula_rate = min(si_rate / 2.0, al_rate / 2.0)
        gaps = {
            "ph3_bound_not_simulated": math.log10(formula_rate)
            - math.log10(rates["3"]),
            "ph4_simulated": math.log10(formula_rate) - math.log10(rates["4"]),
            "ph5_bound_not_simulated": math.log10(formula_rate)
            - math.log10(rates["5"]),
        }
    return {
        "simulated_ph": 4,
        "comparison_bound_ph": [3, 5],
        "ph3_ph5_role": "laboratory_comparison_bounds_only_not_simulated",
        "laboratory_reference": "Palandri-Kharaka acid_plus_neutral_298K_approximation",
        "laboratory_formula_unit_rates_mol_m2_s": rates,
        "laboratory_rate_basis": "kaolinite_formula_units_mol_m2_s",
        "kaolinite_stoichiometry": "Al2Si2O5(OH)4",
        "formula_unit_conversion": "min(Si_mol_flux/2, Al_mol_flux/2)",
        "simulated_species_rates_mol_m2_s": {
            "Si": si_rate,
            "Al": al_rate,
        },
        "simulated_kaolinite_formula_unit_rate_mol_m2_s": formula_rate,
        "formula_unit_rate_status": (
            "represented"
            if formula_rate is not None
            else "censored_species_rate_unavailable"
        ),
        "log10_gaps_vs_laboratory": gaps,
    }


def independently_derive(root: Path, manifest: Mapping[str, Any]) -> dict[str, Any]:
    smoke = bool(manifest["smoke"])
    expected_seeds = SEEDS[:1] if smoke else SEEDS
    identities = [
        (scenario, seed) for scenario, _, _ in _scenarios() for seed in expected_seeds
    ]
    declared = [(row["scenario"], row["seed"]) for row in manifest["required_runs"]]
    if len(declared) != len(set(declared)):
        raise ValueError("duplicate run identities")
    if declared != identities:
        raise ValueError("missing or reordered run identities")
    raw = [_audit_replica(root, scenario, seed) for scenario, seed in identities]
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in raw:
        grouped[row["scenario"]].append(row)
    scenarios: dict[str, dict[str, Any]] = {}
    replica_output = []
    paired_contributions: dict[str, dict[int, tuple[float, float]]] = {}
    for scenario, family, perturbation in _scenarios():
        rows = sorted(grouped[scenario], key=lambda row: row["seed"])
        logs = [row["log_weight"] for row in rows]
        all_replicas_complete = all(row["complete"] for row in rows)
        cadence = min(len(row["checkpoint_state_counts"]) for row in rows)
        cadence_times = rows[0]["checkpoint_times_s"][:cadence]
        if (
            cadence == 0
            or len(cadence_times) != cadence
            or (
                all_replicas_complete
                and any(len(row["checkpoint_state_counts"]) != cadence for row in rows)
            )
            or any(
                len(row["checkpoint_times_s"]) < cadence
                or len(row["checkpoint_log_weights"]) < cadence
                or any(
                    not _same(left, right)
                    for left, right in zip(
                        row["checkpoint_times_s"][:cadence],
                        cadence_times,
                        strict=True,
                    )
                )
                for row in rows
            )
        ):
            raise ValueError("common fixed checkpoint cadence mismatch")
        evolution_points = []
        for checkpoint in range(cadence):
            names = rows[0]["checkpoint_state_counts"][checkpoint]
            checkpoint_logs = [
                row["checkpoint_log_weights"][checkpoint] for row in rows
            ]
            if checkpoint == 0:
                initial_counts = dict(names)
                if any(not _same(weight, 0.0) for weight in checkpoint_logs) or any(
                    row["checkpoint_state_counts"][0] != initial_counts for row in rows
                ):
                    raise ValueError(
                        "time-zero checkpoint weights/states do not reproduce initial census"
                    )
                evolution_points.append(
                    {
                        name: {
                            "value": float(count),
                            "log_mean": math.log(count) if count > 0 else None,
                            "status": "represented",
                        }
                        for name, count in initial_counts.items()
                    }
                )
                continue
            evolution_points.append(
                {
                    name: _state_population_point(
                        checkpoint_logs,
                        [
                            row["checkpoint_state_counts"][checkpoint][name]
                            for row in rows
                        ],
                    )
                    for name in names
                }
            )
        evolution_counts = [
            {name: point["value"] for name, point in points.items()}
            for points in evolution_points
        ]
        evolution_statuses = {
            point["status"]
            for points in evolution_points
            for point in points.values()
            if point["status"] != "represented"
        }
        evolution = [
            {
                "time_s": instant,
                "state_counts": counts,
                "log_state_counts": {
                    name: point["log_mean"] for name, point in points.items()
                },
                "numerical_status": (
                    "represented"
                    if not evolution_statuses
                    else next(iter(evolution_statuses))
                    if len(evolution_statuses) == 1
                    else "censored_numerical_range"
                ),
            }
            for instant, counts, points in zip(
                cadence_times, evolution_counts, evolution_points, strict=True
            )
        ]
        fixed_kind_totals: dict[str, float] = defaultdict(float)
        for state, count in rows[0]["checkpoint_state_counts"][0].items():
            fixed_kind_totals[state.partition(".")[0]] += count
        if evolution_statuses:
            stationarity = {
                "passed": False,
                "status": "censored_numerical_state_evolution",
                "numerical_statuses": sorted(evolution_statuses),
                "tolerance_fraction": STATIONARITY_LIMIT,
                "kinds": {},
            }
        else:
            stationarity = _stationarity(evolution_counts, fixed_kind_totals)
        all_replica_stationarity_passed = all(
            row["stationarity"]["passed"] for row in rows
        )
        stationarity_passed = bool(
            stationarity["passed"] and all_replica_stationarity_passed
        )
        species_rows = {}
        for species in ("Si", "Al"):
            estimate = _estimate(
                logs,
                [row["whole_deck_rates"][species] for row in rows],
                _seed(scenario, species),
            )
            species_rows[species] = _enforce_scientific_gates(
                estimate,
                complete=all_replicas_complete,
                stationary=stationarity_passed,
            )
        combined_samples = [
            row["whole_deck_rates"]["Si"] + row["whole_deck_rates"]["Al"]
            for row in rows
        ]
        combined = _estimate(logs, combined_samples, _seed(scenario, "combined"))
        combined = _enforce_scientific_gates(
            combined,
            complete=all_replicas_complete,
            stationary=stationarity_passed,
        )
        paired_contributions[scenario] = {
            row["seed"]: (
                row["log_weight"],
                row["whole_deck_rates"]["Si"] + row["whole_deck_rates"]["Al"],
            )
            for row in rows
        }
        si_mean = species_rows["Si"]["mean"]
        al_mean = species_rows["Al"]["mean"]
        ratio = (
            si_mean / al_mean
            if isinstance(si_mean, float) and isinstance(al_mean, float) and al_mean > 0
            else None
        )
        scenarios[scenario] = {
            "family": family,
            "perturbation": perturbation,
            "species": species_rows,
            "combined_si_al": combined,
            "stoichiometry_si_to_al_from_ensemble_means": ratio,
            "state_population_evolution_unnormalized_is": evolution,
            "stationarity": stationarity,
            "all_replica_stationarity_passed": all_replica_stationarity_passed,
            "all_replicas_complete": all_replicas_complete,
            "stop_reasons": sorted({row["stop_reason"] for row in rows}),
        }
        for row in rows:
            replica_output.append(
                {
                    "scenario": row["scenario"],
                    "seed": row["seed"],
                    "log_weight": row["log_weight"],
                    "original_lattice_release_counts": row["release_counts"],
                    "whole_finite_deck_rates_mol_m2_s": row["whole_deck_rates"],
                    "reachable_census_rates_mol_m2_s": row["reachable_census_rates"],
                    "area_time_a2_s": row["area_time_a2_s"],
                    "event_count": row["event_count"],
                    "likelihood_segment_count": row["segment_count"],
                    "stationarity": row["stationarity"],
                    "stop_reason": row["stop_reason"],
                }
            )
    families = []
    for family in FAMILIES:
        responses = {}
        for perturbation in PERTURBATIONS:
            scenario = f"{family}__{perturbation}"
            responses[perturbation] = _paired(
                paired_contributions["nominal"],
                paired_contributions[scenario],
                bool(scenarios["nominal"]["combined_si_al"]["rankable"]),
                bool(scenarios[scenario]["combined_si_al"]["rankable"]),
                _seed(family, perturbation, "response"),
            )
        rankable = all(
            response["status"] == "rankable"
            and response["delta_log10_ci95"] is not None
            for response in responses.values()
        )
        score = (
            max(abs(response["delta_log10"]) for response in responses.values())
            if rankable
            else None
        )
        irrelevant = rankable and all(
            _irrelevant(response["delta_log10_ci95"]) for response in responses.values()
        )
        families.append(
            {
                "family": family,
                "rankable": rankable,
                "score_max_abs_delta_log10": score,
                "demonstrated_irrelevant": irrelevant,
                "responses": responses,
            }
        )
    ordered = sorted(
        (row for row in families if row["rankable"]),
        key=lambda row: (-row["score_max_abs_delta_log10"], row["family"]),
    )
    for rank, row in enumerate(ordered, start=1):
        row["ordinal_rank_among_rankable"] = rank
    irrelevant_names = [
        row["family"] for row in families if row["demonstrated_irrelevant"]
    ]
    all_complete = all(row["complete"] for row in raw)
    all_stationary = all(
        row["stationarity"]["passed"] and row["all_replica_stationarity_passed"]
        for row in scenarios.values()
    )
    if smoke:
        acceptance = "not_evaluated_smoke"
        ranking_status = "not_evaluated_smoke"
    elif all_complete and all_stationary:
        if not ordered:
            acceptance = "accepted_typed_all_censored_no_go"
            ranking_status = "no_rankable_families_no_go"
        elif len(ordered) == len(FAMILIES):
            acceptance = "accepted_typed_outcomes"
            ranking_status = "complete"
        else:
            acceptance = "accepted_typed_partial_ranking"
            ranking_status = "partial_rankable_families"
    else:
        acceptance = "rejected_incomplete_or_nonstationary"
        ranking_status = "not_scientifically_accepted"
    nominal_species = scenarios["nominal"]["species"]
    lab_species_rates = {
        species: estimate["mean"] if estimate["rankable"] else None
        for species, estimate in nominal_species.items()
    }
    return {
        "schema": "a9b-sensitivity-analysis-v2",
        "survey_tier_platform_test": True,
        "not_production": True,
        "artifact_verified": True,
        "scientific_acceptance": acceptance,
        "ordinal_ranking_complete": len(ordered) == len(FAMILIES),
        "ordinal_ranking_status": ranking_status,
        "simulated_scope": {
            "temperature_k": 298.0,
            "ph": 4,
            "ph3_ph5": "laboratory_comparison_bounds_only_not_simulated",
        },
        "estimator_contract": {
            "name": "unnormalized_importance_sampling",
            "aggregation": "log_domain",
            "self_normalized": False,
            "weight_ess_fraction_min": WEIGHT_GATE,
            "positive_species_contribution_ess_min": CONTRIBUTION_GATE,
            "positive_species_contribution_ess_role": CONTRIBUTION_GATE_ROLE,
            "zero_biased_observation_poisson_bound": "unavailable",
        },
        "census": manifest["estimand"],
        "run_count": len(raw),
        "expected_run_count": 29 if smoke else 232,
        "per_replica": sorted(
            replica_output, key=lambda row: (row["scenario"], row["seed"])
        ),
        "scenarios": scenarios,
        "family_ranking": families,
        "top_k_preregistered": TOP_K,
        "top3_among_rankable": [row["family"] for row in ordered[:TOP_K]],
        "families_demonstrated_irrelevant": irrelevant_names,
        "irrelevance_statement": (
            "none demonstrated irrelevant"
            if not irrelevant_names
            else "only families whose complete response intervals lie inside the documented equivalence band"
        ),
        "equivalence_band_log10": list(EQUIVALENCE),
        "laboratory_comparison": _lab(lab_species_rates),
    }


def _pending_analysis(finalized: dict[str, Any]) -> dict[str, Any]:
    """Independently construct the only report the analyzer may publish."""
    pending = copy.deepcopy(finalized)
    pending["preliminary_conclusions"] = {
        "status": "unverified_preliminary",
        "accepted_conclusions_withheld_until": "independent_verification_receipt",
    }
    pending.update(
        {
            "artifact_verified": False,
            "verification_status": "pending_independent_verification",
            "conclusions_status": "unverified_preliminary",
            "scientific_acceptance": "unverified_pending_independent_verification",
            "ordinal_ranking_complete": False,
            "ordinal_ranking_status": "unverified_pending_independent_verification",
            "top3_among_rankable": [],
            "families_demonstrated_irrelevant": [],
            "irrelevance_statement": "unverified pending independent verification",
        }
    )
    for scenario in pending["scenarios"].values():
        scenario["verification_status"] = "unverified_preliminary"
        scenario["combined_si_al"]["verification_status"] = "unverified_preliminary"
        for estimate in scenario["species"].values():
            estimate["verification_status"] = "unverified_preliminary"
    for family in pending["family_ranking"]:
        family["verification_status"] = "unverified_preliminary"
        family["demonstrated_irrelevant"] = False
        if "ordinal_rank_among_rankable" in family:
            family.pop("ordinal_rank_among_rankable")
        for response in family["responses"].values():
            response["verification_status"] = "unverified_preliminary"
    pending["laboratory_comparison"]["verification_status"] = "unverified_preliminary"
    return pending


def _validate_manifest_contract(manifest: Mapping[str, Any]) -> None:
    smoke = manifest.get("smoke")
    if type(smoke) is not bool:
        raise ValueError("campaign smoke flag must be boolean")
    expected_seeds = list(SEEDS[:1] if smoke else SEEDS)
    expected_run_count = 29 if smoke else 232
    expected_checkpoints = 3 if smoke else 21
    exact_values = {
        "schema": "a9b-sensitivity-campaign-v2",
        "smoke": smoke,
        "survey_tier_platform_test": True,
        "not_production": True,
        "temperature_k": 298.0,
        "simulated_ph": 4,
        "ph3_ph5_scope": "laboratory_comparison_bounds_only_not_simulated",
        "seeds": expected_seeds,
        "production_fixed_seeds": list(SEEDS),
        "scenario_count": 29,
        "replicas_per_scenario": len(expected_seeds),
        "required_run_count": expected_run_count,
        "top_k_preregistered": TOP_K,
        "equivalence_band_log10": list(EQUIVALENCE),
        "weight_ess_fraction_min": WEIGHT_GATE,
        "species_contribution_ess_min": CONTRIBUTION_GATE,
        "species_contribution_ess_role": CONTRIBUTION_GATE_ROLE,
        "state_population_stationarity": {
            "required": "all_states_within_every_kind_for_ensemble_and_each_replica",
            "tail_checkpoint_count": 5,
            "max_fractional_span_and_endpoint_change": STATIONARITY_LIMIT,
        },
        "estimand": {
            "whole_finite_deck_centers": 320,
            "reachable_centers": 260,
            "topology_no_go_centers": 60,
            "primary": "whole_finite_deck_including_topology_no_go_in_denominator_and_area",
            "secondary": "reachable_center_census",
        },
        "run_contract": {
            "horizon_s": 1e-12 if smoke else 1e-6,
            "checkpoints": expected_checkpoints,
            "bias_factor": 1e6,
            "max_events": 10_000 if smoke else 1_000_000,
        },
    }
    for key, expected in exact_values.items():
        if manifest.get(key) != expected:
            raise ValueError(f"campaign manifest contract mismatch in {key}")
    required_keys = set(exact_values) | {
        "status",
        "inputs",
        "scenarios",
        "required_runs",
    }
    optional_keys = {"workers", "runner_sha256", "completed_run_count"}
    if set(manifest) != required_keys | (set(manifest) & optional_keys):
        raise ValueError("campaign manifest key set is not canonical")
    inputs = manifest.get("inputs")
    expected_paths = {
        "base_deck": "inputs/kaolinite-approx.toml",
        "reservoir_contract": "inputs/kaolinite-reservoirs.toml",
        "ph4_deck": "inputs/kaolinite-ph4.toml",
        "ph4_evidence": "inputs/kaolinite-ph4-evidence.json",
        "path_report": "inputs/path-report.json",
        "lineage_snapshot": "inputs/initial-snapshot.pgif.json",
    }
    expected_input_keys = {
        key for name in expected_paths for key in (name, f"{name}_sha256")
    }
    if not isinstance(inputs, dict) or set(inputs) != expected_input_keys:
        raise ValueError("campaign input manifest keys mismatch")
    for key, expected_path in expected_paths.items():
        if inputs.get(key) != expected_path:
            raise ValueError(f"campaign input path mismatch in {key}")


def verify_campaign(root: Path, *, write_verification: bool) -> dict[str, Any]:
    root = root.resolve()
    verification_path = root / "verification.json"
    if write_verification:
        remove_file_durable(verification_path)
    manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
    if (
        manifest.get("schema") != "a9b-sensitivity-campaign-v2"
        or manifest.get("status") != "raw_complete"
    ):
        raise ValueError("campaign manifest is not raw-complete A9b v2")
    _validate_manifest_contract(manifest)
    _verify_and_regenerate_decks(root, manifest)
    independently_derived = independently_derive(root, manifest)
    expected_bytes = canonical_json_bytes(_pending_analysis(independently_derived))
    actual_path = root / "analysis.json"
    actual_bytes = actual_path.read_bytes()
    if actual_bytes != expected_bytes:
        raise ValueError("independently regenerated analysis bytes differ")
    finalized_analysis = copy.deepcopy(independently_derived)
    finalized_analysis["verification_status"] = "independently_verified"
    finalized_analysis["conclusions_status"] = "independently_verified"
    finalized_bytes = canonical_json_bytes(finalized_analysis)
    verification = {
        "schema": "a9b-sensitivity-independent-verification-v3",
        "status": "verified",
        "independent_from_analyzer": True,
        "scenario_decks_regenerated": 29,
        "raw_replicas_replayed": independently_derived["run_count"],
        "likelihood_increments_recomputed": True,
        "lineage_recomputed": True,
        "state_stationarity_recomputed": True,
        "estimator_and_ranking_recomputed": True,
        "derived_bytes_equal": True,
        "analysis_sha256": hashlib.sha256(actual_bytes).hexdigest(),
        "finalized_analysis_sha256": hashlib.sha256(finalized_bytes).hexdigest(),
        "artifact_verified": finalized_analysis["artifact_verified"],
        "scientific_acceptance": finalized_analysis["scientific_acceptance"],
        "ordinal_ranking_complete": finalized_analysis["ordinal_ranking_complete"],
        "ordinal_ranking_status": finalized_analysis["ordinal_ranking_status"],
        "finalized_analysis": finalized_analysis,
    }
    if write_verification:
        write_json_atomic(verification_path, verification)
    return verification


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("campaign_root", type=Path)
    parser.add_argument("--write-verification", action="store_true")
    args = parser.parse_args(argv)
    result = verify_campaign(
        args.campaign_root, write_verification=args.write_verification
    )
    print(
        json.dumps(
            {
                key: value
                for key, value in result.items()
                if key != "finalized_analysis"
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
