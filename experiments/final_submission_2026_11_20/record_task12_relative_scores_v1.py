"""Bind the displayed Task12 table to its exact accepted submission and email."""

from __future__ import annotations

import hashlib
import json
import re
from datetime import datetime, timezone
from decimal import Decimal
from pathlib import Path


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def main() -> None:
    folder = Path("/private/tmp/ghist-final-submissions-20261120/task12")
    repo = Path(__file__).resolve().parent / "task12"
    ax_path = folder / "leaderboard_current.ax.txt"
    ax_raw = ax_path.read_bytes()
    ax = ax_raw.decode()
    url = "https://www.synapse.org/Synapse:syn74803798/wiki/640983"
    assert url in ax.splitlines()[0]
    heading = "Relatedness in a multifamily cohort with close inbreeding"
    assert f"heading {heading}, Value: 2" in ax
    assert ax.count(" table\n") == 1
    table = ax.split(" table\n")[1].split("container pagination navigation")[0]
    blocks = re.split(r"(?m)^\t+\d+ row\n", table)[1:]
    assert re.findall(r"\d+ text (.+)", blocks[0]) == [
        "createdBy", "id", "name", "createdOn", "meioses_count", "relation"]
    entries = []
    for block in blocks[1:]:
        users = re.findall(r"link Description: @([^,]+), Value: synapse.org/Profile:(\d+)", block)
        assert len(users) == 1
        texts = re.findall(r"\d+ text (.+)", block)
        assert len(texts) == 5
        sid, name, created, meioses, relation = texts
        entries.append({"submitter": users[0][0], "submitter_profile_id": users[0][1],
                        "submission_id": sid, "submission_name": name, "created_on_display": created,
                        "meioses_count": float(meioses), "relation": float(relation),
                        "meioses_count_display_text": meioses, "relation_display_text": relation})
    assert len(entries) == 2 and re.findall(r"\d+ text (\d+) total rows /", ax) == ["2"]
    assert "button (disabled) Go to next page" in ax
    assert {entry["submission_id"] for entry in entries} == {"9785373", "9783965"}
    own = next(entry for entry in entries if entry["submission_id"] == "9785373")
    assert own["submitter"] == "ancestrum"
    email_path = folder / "score_email_evidence_9785373.json"
    email_raw = email_path.read_bytes()
    email = json.loads(email_raw)
    assert email["submission_id"] == "9785373" and email["evaluation_queue_id"] == "9619678"
    assert email["task"] == 12 and email["status"] == "SCORED"
    message = next(item for item in email["records"] if item["gmail_message_id"] == email["status_email_id"])
    assert "(id: 9785373) is scored" in message["body_text"]
    email_text = dict(re.findall(r"(?m)^(meioses_count|relation) : ([0-9.eE+-]+)$", message["body_text"]))
    assert set(email_text) == {"meioses_count", "relation"}
    receipt_path = folder / "submission_receipt.json"
    receipt_raw = receipt_path.read_bytes()
    receipt = json.loads(receipt_raw)
    assert receipt["submission_id"] == "9785373" and receipt["evaluation_queue_id"] == "9619678"
    assert receipt["score_evidence"]["sha256"] == sha(email_raw) and receipt["scores"] == email["scores"]
    comparisons = {}
    for metric in ("meioses_count", "relation"):
        ours = Decimal(own[f"{metric}_display_text"])
        assert Decimal(email_text[metric]) == Decimal(str(email["scores"][metric]))
        assert Decimal(email_text[metric]).quantize(Decimal(".001")) == ours
        others = [entry for entry in entries if entry["submission_id"] != "9785373"]
        best = max(Decimal(entry[f"{metric}_display_text"]) for entry in others)
        comparisons[metric] = {
            "own_displayed_score": float(ours), "own_email_score": email["scores"][metric],
            "own_email_score_text": email_text[metric],
            "rank": 1 + sum(Decimal(entry[f"{metric}_display_text"]) > ours for entry in others),
            "rank_denominator": len(entries),
            "tie_count_including_own": sum(Decimal(entry[f"{metric}_display_text"]) == ours for entry in entries),
            "best_other_displayed_score": float(best),
            "best_other_submitters": [entry["submitter"] for entry in others if Decimal(entry[f"{metric}_display_text"]) == best],
            "displayed_gap_to_best_other": float(best - ours),
        }
    record = {
        "schema_version": 2, "task": 12, "own_submission_id": "9785373", "evaluation_queue_id": "9619678",
        "comparison_status": "VERIFIED_BY_OWN_SUBMISSION_ID", "mapping_verified": True, "own_listed": True,
        "uploaded_entity_id": receipt["uploaded_entity_id"], "uploaded_entity_version": receipt["uploaded_entity_version"],
        "candidate_sha256": receipt["candidate_sha256"], "source_url": url,
        "observed_at_utc": datetime.fromtimestamp(ax_path.stat().st_mtime, timezone.utc).isoformat(),
        "observation_time_basis": "Filesystem save time of root-provided fresh authenticated full AX",
        "recorded_at_utc": datetime.now(timezone.utc).isoformat(), "entries": entries, "comparisons": comparisons,
        "visible_total_rows": 2, "distinct_displayed_submitters": 2, "other_displayed_competitors": 1,
        "own_scores_from_email": email["scores"], "own_scores_from_email_text": email_text,
        "email_score_time_utc": email["status_at_utc"], "email_message_id": email["status_email_id"],
        "metric": "Challenge truth-class-size-weighted Jaccard, separately for relation and meioses_count; higher is better.",
        "precision_note": "Ranks and gaps use displayed rounded scores only; exact own email scores are separately retained. No hidden competitor precision is inferred.",
        "page_heading_as_displayed": heading,
        "page_title_mismatch": "Official wiki640983 describes the multifamily close-inbreeding Task11, but its table contains the independently scored Task12 submission9785373 and bgrel submission9783965. Wiki640984 instead contains Task11 submission9785361. The official dashboard titles and table associations appear swapped.",
        "mapping_verification_scope": "Our displayed submission ID exactly matches the Task12 Relatedness4 queue9619678 scoring email and uploaded candidate receipt; both email metrics round to displayed scores. Competitor is in this same published table. Underlying table scope metadata was not independently retrieved.",
        "scope": "Current accepted task-specific displayed table only; no overall challenge rank or final closed ranking asserted.",
        "table_evidence": {"path": str(ax_path), "sha256": sha(ax_raw)},
        "email_evidence": {"path": str(email_path), "sha256": sha(email_raw)},
        "submission_receipt_evidence": {"path": str(receipt_path), "sha256": sha(receipt_raw)},
        "model_or_prediction_changes": False, "truth_access": False,
        "source_sha256": sha(Path(__file__).read_bytes()),
    }
    payload = (json.dumps(record, indent=2, sort_keys=True) + "\n").encode()
    for destination in (folder, repo):
        target = destination / "relative_score_record_current.json"
        assert not target.exists(), f"Preserve prior record: {target}"
        target.write_bytes(payload)
        assert target.read_bytes() == payload
    if (repo / ax_path.name).exists():
        assert (repo / ax_path.name).read_bytes() == ax_raw
    else:
        (repo / ax_path.name).write_bytes(ax_raw)
    assert ax_path.read_bytes() == ax_raw and email_path.read_bytes() == email_raw
    assert receipt_path.read_bytes() == receipt_raw
    print(json.dumps({"comparisons": comparisons, "record_sha256": sha(payload)}))


if __name__ == "__main__":
    main()
