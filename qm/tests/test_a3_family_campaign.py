"""Unit gates for the bounded A3 family campaign supervisor."""

from __future__ import annotations

import json

import pytest

from scripts import a3_family_campaign as campaign


def test_wait_for_gpu_retries_until_available():
    responses = iter([(False, "busy pid=7"), (False, "still busy"), (True, "free")])
    clock = iter([0.0, 0.0, 1.0, 2.0])
    sleeps: list[float] = []
    heartbeats: list[str] = []

    message = campaign.wait_for_gpu(
        timeout_seconds=10,
        poll_seconds=1,
        heartbeat=heartbeats.append,
        lane_probe=lambda: next(responses),
        monotonic=lambda: next(clock),
        sleep=sleeps.append,
    )

    assert message == "free"
    assert heartbeats == ["busy pid=7", "still busy", "free"]
    assert sleeps == [1, 1]


def test_wait_for_gpu_times_out_with_last_owner():
    clock = iter([0.0, 0.0, 2.0])
    with pytest.raises(TimeoutError, match="live owner"):
        campaign.wait_for_gpu(
            timeout_seconds=1,
            poll_seconds=1,
            heartbeat=lambda _message: None,
            lane_probe=lambda: (False, "live owner"),
            monotonic=lambda: next(clock),
            sleep=lambda _seconds: None,
        )


def test_validate_result_binds_identity_finite_values_and_hashes(tmp_path):
    result = tmp_path / "results.json"
    store = tmp_path / "store.sqlite"
    result.write_text(
        json.dumps(
            {
                "family": "oss",
                "state": "neutral",
                "n_intact": 2,
                "dG_kj": 123.4,
                "dH_kj": 111.0,
                "ts_imaginary_cm": 98.0,
                "route": "proton-neb",
            }
        )
    )
    store.write_bytes(b"sqlite-evidence")

    record = campaign.validate_result(result, family="oss", state="neutral", n_intact=2)

    assert record["dG_kj"] == 123.4
    assert record["route"] == "proton-neb"
    assert record["result_sha256"] == campaign.sha256_path(result)
    assert record["store_sha256"] == campaign.sha256_path(store)


@pytest.mark.parametrize(
    ("change", "message"),
    [
        ({"n_intact": 3}, "identity mismatch"),
        ({"dG_kj": float("nan")}, "invalid dG_kj"),
        ({"route": ""}, "no route provenance"),
    ],
)
def test_validate_result_rejects_false_green_payloads(tmp_path, change, message):
    payload = {
        "family": "oss",
        "state": "neutral",
        "n_intact": 2,
        "dG_kj": 123.4,
        "dH_kj": 111.0,
        "ts_imaginary_cm": 98.0,
        "route": "direct",
    }
    payload.update(change)
    result = tmp_path / "results.json"
    result.write_text(json.dumps(payload))
    (tmp_path / "store.sqlite").write_bytes(b"sqlite-evidence")

    with pytest.raises(RuntimeError, match=message):
        campaign.validate_result(result, family="oss", state="neutral", n_intact=2)


def test_atomic_json_never_emits_nan(tmp_path):
    path = tmp_path / "receipt.json"
    with pytest.raises(ValueError):
        campaign.atomic_json(path, {"bad": float("nan")})
    assert not path.exists()
