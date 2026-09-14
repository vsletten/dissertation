#!/usr/bin/env python3
"""Prove that the caller is in the bounded A9b user unit."""

from __future__ import annotations

import json
import os
import subprocess
import sys
from collections.abc import Mapping
from pathlib import Path
from typing import Any

BUSCTL = Path("/usr/bin/busctl")
SYSTEMD_DESTINATION = "org.freedesktop.systemd1"
SYSTEMD_MANAGER_PATH = "/org/freedesktop/systemd1"
SYSTEMD_MANAGER_INTERFACE = "org.freedesktop.systemd1.Manager"
SYSTEMD_UNIT_INTERFACE = "org.freedesktop.systemd1.Unit"
SYSTEMD_SERVICE_INTERFACE = "org.freedesktop.systemd1.Service"


def validate_live_proof(
    proof: Mapping[str, object], expected_unit: str, maximum_seconds: int
) -> int:
    """Validate exact unit membership and return its bounded runtime in usec."""
    if proof.get("unit_id") != expected_unit:
        raise ValueError("live unit id does not match the intended unit")
    if proof.get("active_state") != "active":
        raise ValueError("intended unit is not active")
    control_group = proof.get("control_group")
    process_cgroup = proof.get("process_cgroup")
    if (
        not isinstance(control_group, str)
        or not control_group.startswith("/")
        or control_group == "/"
        or not isinstance(process_cgroup, str)
        or not (
            process_cgroup == control_group
            or process_cgroup.startswith(f"{control_group}/")
        )
    ):
        raise ValueError("current process is not in the intended unit control group")
    runtime_max_usec = proof.get("runtime_max_usec")
    maximum_usec = maximum_seconds * 1_000_000
    if (
        isinstance(runtime_max_usec, bool)
        or not isinstance(runtime_max_usec, int)
        or runtime_max_usec <= 0
        or runtime_max_usec > maximum_usec
    ):
        raise ValueError("live RuntimeMaxUSec is absent or exceeds the ceiling")
    return runtime_max_usec


def _busctl_json(*arguments: str) -> dict[str, Any]:
    runtime_directory = f"/run/user/{os.getuid()}"
    environment = {
        "DBUS_SESSION_BUS_ADDRESS": f"unix:path={runtime_directory}/bus",
        "XDG_RUNTIME_DIR": runtime_directory,
    }
    completed = subprocess.run(
        [str(BUSCTL), "--user", "--json=short", "--no-pager", *arguments],
        check=False,
        capture_output=True,
        text=True,
        env=environment,
        timeout=10,
    )
    if completed.returncode != 0:
        raise RuntimeError("cannot query the user systemd manager")
    payload = json.loads(completed.stdout)
    if not isinstance(payload, dict) or "data" not in payload:
        raise ValueError("malformed systemd D-Bus response")
    return payload


def _property(object_path: str, interface: str, name: str) -> object:
    payload = _busctl_json(
        "get-property",
        SYSTEMD_DESTINATION,
        object_path,
        interface,
        name,
    )
    return payload["data"]


def _process_cgroup(path: Path = Path("/proc/self/cgroup")) -> str:
    unified = []
    for line in path.read_text(encoding="utf-8").splitlines():
        hierarchy, controllers, cgroup = line.split(":", 2)
        if hierarchy == "0" and controllers == "":
            unified.append(cgroup)
    if len(unified) != 1 or not unified[0].startswith("/"):
        raise ValueError("cannot identify the current unified control group")
    return unified[0]


def process_niceness(path: Path = Path("/proc/self/stat")) -> int:
    """Read Linux proc stat field 19 without trusting process environment."""
    stat = path.read_text(encoding="utf-8").strip()
    command_end = stat.rfind(")")
    if command_end < 0:
        raise ValueError("cannot parse process stat command field")
    fields = stat[command_end + 1 :].split()
    if len(fields) <= 16:
        raise ValueError("process stat is missing niceness field 19")
    niceness = int(fields[16])
    if not -20 <= niceness <= 19:
        raise ValueError("process niceness is outside the Linux range")
    return niceness


def require_minimum_niceness(observed: object, minimum: int) -> int:
    if isinstance(observed, bool) or not isinstance(observed, int):
        raise TypeError("observed niceness must be an integer")
    if isinstance(minimum, bool) or not isinstance(minimum, int):
        raise TypeError("minimum niceness must be an integer")
    if observed < minimum:
        raise ValueError(
            f"observed process niceness {observed} is below required {minimum}"
        )
    return observed


def collect_live_proof(expected_unit: str) -> dict[str, object]:
    unit_lookup = _busctl_json(
        "call",
        SYSTEMD_DESTINATION,
        SYSTEMD_MANAGER_PATH,
        SYSTEMD_MANAGER_INTERFACE,
        "GetUnit",
        "s",
        expected_unit,
    )
    object_paths = unit_lookup["data"]
    if (
        not isinstance(object_paths, list)
        or len(object_paths) != 1
        or not isinstance(object_paths[0], str)
    ):
        raise ValueError("systemd returned an invalid unit object path")
    object_path = object_paths[0]
    invocation_bytes = _property(object_path, SYSTEMD_UNIT_INTERFACE, "InvocationID")
    if (
        not isinstance(invocation_bytes, list)
        or len(invocation_bytes) != 16
        or any(
            isinstance(value, bool)
            or not isinstance(value, int)
            or not 0 <= value <= 255
            for value in invocation_bytes
        )
    ):
        raise ValueError("systemd returned an invalid invocation id")
    return {
        "unit_id": _property(object_path, SYSTEMD_UNIT_INTERFACE, "Id"),
        "active_state": _property(object_path, SYSTEMD_UNIT_INTERFACE, "ActiveState"),
        "control_group": _property(
            object_path, SYSTEMD_SERVICE_INTERFACE, "ControlGroup"
        ),
        "runtime_max_usec": _property(
            object_path, SYSTEMD_SERVICE_INTERFACE, "RuntimeMaxUSec"
        ),
        "invocation_id": "".join(f"{value:02x}" for value in invocation_bytes),
        "process_cgroup": _process_cgroup(),
        "process_niceness": process_niceness(),
    }


def main(argv: list[str]) -> int:
    if len(argv) not in {3, 4}:
        print(
            f"usage: {argv[0]} UNIT MAXIMUM_SECONDS [MINIMUM_NICENESS]",
            file=sys.stderr,
        )
        return 64
    expected_unit = argv[1]
    try:
        maximum_seconds = int(argv[2])
        if maximum_seconds <= 0:
            raise ValueError("maximum must be positive")
        proof = collect_live_proof(expected_unit)
        runtime_max_usec = validate_live_proof(proof, expected_unit, maximum_seconds)
        invocation_id = proof["invocation_id"]
        niceness = require_minimum_niceness(proof["process_niceness"], -20)
        if argv[3:]:
            niceness = require_minimum_niceness(niceness, int(argv[3]))
    except (
        OSError,
        RuntimeError,
        TypeError,
        ValueError,
        json.JSONDecodeError,
    ) as error:
        print(f"bounded systemd proof failed: {error}", file=sys.stderr)
        return 77
    print(runtime_max_usec, invocation_id, niceness)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
