"""Read saved authenticated scoreboard AX and email evidence, without model access."""

from __future__ import annotations

import hashlib
import json
import re
from datetime import datetime, timezone
from decimal import Decimal
from pathlib import Path


ROOT = Path("/private/tmp/ghist-final-submissions-20261120")
REPO = Path(__file__).resolve().parent
TASKS = {
    5: ("9785369", "9619656", "640748"),
    6: ("9785370", "9619680", "640744"),
    7: ("9785371", "9619658", "640749"),
    8: ("9785372", "9619662", "640745"),
}


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def main() -> None:
    for task, (own_id, queue_id, wiki_id) in TASKS.items():
        folder = ROOT / f"task{task}"
        ax_path = folder / "leaderboard_current.ax.txt"
        ax_raw = ax_path.read_bytes()
        ax = ax_raw.decode()
        source_url = f"https://www.synapse.org/Synapse:syn74803798/wiki/{wiki_id}"
        assert source_url in ax.splitlines()[0]
        assert ax.count(" table\n") == 1
        table = ax.split(" table\n", 1)[1].split("container pagination navigation", 1)[0]
        blocks = re.split(r"(?m)^\t+\d+ row\n", table)[1:]
        assert re.findall(r"\d+ text (.+)", blocks[0]) == ["createdBy", "id", "name", "createdOn", "f1"]
        rows = []
        for block in blocks[1:]:
            submitters = re.findall(r"link Description: @([^,]+), Value: synapse.org/Profile:(\d+)", block)
            assert len(submitters) == 1
            texts = re.findall(r"\d+ text (.+)", block)
            assert len(texts) == 4
            submission_id, name, created_on, displayed = texts
            assert submission_id.isdecimal() and Decimal(displayed).is_finite()
            rows.append({"submitter": submitters[0][0], "submitter_profile_id": submitters[0][1],
                         "submission_id": submission_id, "submission_name": name,
                         "created_on_display": created_on, "f1": float(displayed),
                         "f1_display_text": displayed})
        total_matches = re.findall(r"\d+ text (\d+) total rows /", ax)
        assert len(total_matches) == 1 and len(rows) == int(total_matches[0]) == 6
        assert "button (disabled) Go to next page" in ax
        assert len({row["submission_id"] for row in rows}) == len(rows)
        assert len({row["submitter"] for row in rows}) == len(rows)
        scores = [Decimal(row["f1_display_text"]) for row in rows]
        assert scores == sorted(scores, reverse=True)
        own_matches = [row for row in rows if row["submission_id"] == own_id]
        assert len(own_matches) == 1 and own_matches[0]["submitter"] == "ancestrum"
        own = own_matches[0]
        own_display = Decimal(own["f1_display_text"])
        email_path = folder / f"score_email_evidence_{own_id}.json"
        email_raw = email_path.read_bytes()
        email = json.loads(email_raw)
        assert email["submission_id"] == own_id and email["evaluation_queue_id"] == queue_id
        assert email["task"] == task and email["status"] == "SCORED"
        messages = [r for r in email["records"] if r["gmail_message_id"] == email["status_email_id"]]
        assert len(messages) == 1 and f"(id: {own_id}) is scored" in messages[0]["body_text"]
        raw_scores = dict(re.findall(r"(?m)^(recall|size|f1) : ([0-9.eE+-]+)$", messages[0]["body_text"]))
        assert set(raw_scores) == {"f1", "recall", "size"}
        assert all(Decimal(raw_scores[key]) == Decimal(str(value)) for key, value in email["scores"].items())
        assert Decimal(raw_scores["f1"]).quantize(Decimal(".001")) == own_display
        receipt_path = folder / "submission_receipt.json"
        receipt_raw = receipt_path.read_bytes()
        receipt = json.loads(receipt_raw)
        assert receipt["submission_id"] == own_id and receipt["evaluation_queue_id"] == queue_id
        assert receipt["score_evidence"]["sha256"] == sha(email_raw)
        assert receipt["scores"] == email["scores"]
        better = [row for row in rows if Decimal(row["f1_display_text"]) > own_display]
        equal = [row for row in rows if Decimal(row["f1_display_text"]) == own_display]
        worse = [row for row in rows if Decimal(row["f1_display_text"]) < own_display]
        rank = 1 + len(better)
        best = max(scores)
        nearest_better = min((Decimal(row["f1_display_text"]) for row in better), default=None)
        nearest_worse = max((Decimal(row["f1_display_text"]) for row in worse), default=None)
        comparison = {
            "own_displayed_score": float(own_display), "own_email_score": email["scores"]["f1"],
            "own_email_score_text": raw_scores["f1"],
            "rank": rank, "rank_denominator": len(rows),
            "rank_span_including_displayed_ties": [rank, rank + len(equal) - 1],
            "tie_count_including_own": len(equal),
            "tied_other_submitters": [row["submitter"] for row in equal if row["submission_id"] != own_id],
            "better_other_count": len(better), "equal_other_count": len(equal) - 1,
            "worse_other_count": len(worse), "tied_for_last": len(equal) > 1 and not worse,
            "best_other_displayed_score": float(best),
            "best_other_submitters": [row["submitter"] for row in rows if Decimal(row["f1_display_text"]) == best],
            "displayed_gap_to_best_other": float(best - own_display),
            "nearest_better_displayed_score": float(nearest_better) if nearest_better is not None else None,
            "displayed_gap_to_nearest_better": float(nearest_better - own_display) if nearest_better is not None else None,
            "nearest_worse_displayed_score": float(nearest_worse) if nearest_worse is not None else None,
            "displayed_lead_over_nearest_worse": float(own_display - nearest_worse) if nearest_worse is not None else None,
        }
        record = {
            "schema_version": 2, "task": task, "own_submission_id": own_id,
            "evaluation_queue_id": queue_id, "uploaded_entity_id": receipt["uploaded_entity_id"],
            "uploaded_entity_version": receipt["uploaded_entity_version"],
            "candidate_sha256": receipt["candidate_sha256"], "own_listed": True,
            "comparison_status": "VERIFIED_BY_OWN_SUBMISSION_ID", "mapping_verified": True,
            "source_url": source_url,
            "observed_at_utc": datetime.fromtimestamp(ax_path.stat().st_mtime, timezone.utc).isoformat(),
            "observation_time_basis": "Filesystem save time of root-provided authenticated full AX capture, approximately 21:52 UTC",
            "recorded_at_utc": datetime.now(timezone.utc).isoformat(),
            "visible_total_rows": len(rows), "distinct_displayed_submitters": len(rows),
            "other_displayed_competitors": len(rows) - 1, "entries": rows,
            "metric": "F1; higher is better. Other exact own email metrics are retained without cross-task combination.",
            "comparisons": {"f1": comparison}, "own_scores_from_email": email["scores"],
            "own_scores_from_email_text": raw_scores,
            "email_score_time_utc": email["status_at_utc"], "email_message_id": email["status_email_id"],
            "table_evidence": {"path": str(ax_path), "sha256": sha(ax_raw)},
            "email_evidence": {"path": str(email_path), "sha256": sha(email_raw)},
            "submission_receipt_evidence": {"path": str(receipt_path), "sha256": sha(receipt_raw)},
            "rank_basis": "Competition rank = one plus number of displayed entries strictly higher; all six rows present, each a distinct submitter",
            "precision_note": "Ranks, apparent ties, and gaps use rounded displayed F1 only; exact own email score is separately bound. Hidden competitor precision is unavailable.",
            "scope": "Current displayed task-specific final dashboard only; no overall challenge rank or final closed ranking asserted.",
            "model_or_prediction_changes": False, "truth_access": False,
            "source_sha256": sha(Path(__file__).read_bytes()),
        }
        payload = (json.dumps(record, indent=2, sort_keys=True) + "\n").encode()
        for destination in (folder, REPO / f"task{task}"):
            destination.mkdir(parents=True, exist_ok=True)
            target = destination / "relative_score_record_current.json"
            assert not target.exists(), f"Preserve existing record: {target}"
            target.write_bytes(payload)
            assert target.read_bytes() == payload
        assert ax_path.read_bytes() == ax_raw and email_path.read_bytes() == email_raw
        assert receipt_path.read_bytes() == receipt_raw
        print(json.dumps({"task": task, "comparison": comparison, "record_sha256": sha(payload)}))


if __name__ == "__main__":
    main()
