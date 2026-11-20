"""Record Tasks1--4 posted RRMSE ranks from saved authenticated evidence."""

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
    1: ("9785374", "9619646", "640737", 6),
    2: ("9785375", "9619648", "640747", 4),
    3: ("9785376", "9619650", "640742", 3),
    4: ("9785377", "9619652", "640750", 3),
}


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def main() -> None:
    for task, (sid, queue, wiki, total) in TASKS.items():
        folder = ROOT / f"task{task}"
        ax_path = folder / "leaderboard_current.ax.txt"
        ax_raw = ax_path.read_bytes()
        ax = ax_raw.decode()
        url = f"https://www.synapse.org/Synapse:syn74803798/wiki/{wiki}"
        assert url in ax.splitlines()[0] and ax.count(" table\n") == 1
        table = ax.split(" table\n")[1].split("container pagination navigation")[0]
        blocks = re.split(r"(?m)^\t+\d+ row\n", table)[1:]
        assert re.findall(r"\d+ text (.+)", blocks[0]) == ["createdBy", "id", "name", "createdOn", "RRMSE"]
        entries = []
        for block in blocks[1:]:
            user = re.findall(r"link Description: @([^,]+), Value: synapse.org/Profile:(\d+)", block)
            texts = re.findall(r"\d+ text (.+)", block)
            assert len(user) == 1 and len(texts) == 4
            row_sid, name, created, score = texts
            assert row_sid.isdecimal() and Decimal(score).is_finite() and Decimal(score) >= 0
            entries.append({"submitter": user[0][0], "submitter_profile_id": user[0][1],
                            "submission_id": row_sid, "submission_name": name,
                            "created_on_display": created, "RRMSE": float(score), "RRMSE_display_text": score})
        assert len(entries) == total and re.findall(r"\d+ text (\d+) total rows /", ax) == [str(total)]
        assert "button (disabled) Go to next page" in ax
        assert len({row["submission_id"] for row in entries}) == len({row["submitter"] for row in entries}) == total
        own_rows = [row for row in entries if row["submission_id"] == sid]
        assert len(own_rows) == 1 and own_rows[0]["submitter"] == "ancestrum"
        own = own_rows[0]
        ours = Decimal(own["RRMSE_display_text"])
        email_path = folder / f"score_email_evidence_{sid}.json"
        email_raw = email_path.read_bytes()
        email = json.loads(email_raw)
        assert email["submission_id"] == sid and email["evaluation_queue_id"] == queue
        assert email["task"] == task and email["status"] == "SCORED"
        messages = [record for record in email["records"] if record["gmail_message_id"] == email["status_email_id"]]
        assert len(messages) == 1 and f"(id: {sid}) is scored" in messages[0]["body_text"]
        exact_matches = re.findall(r"(?m)^RRMSE : ([0-9.eE+-]+)$", messages[0]["body_text"])
        assert len(exact_matches) == 1 and set(email["scores"]) == {"RRMSE"}
        exact_text = exact_matches[0]
        assert Decimal(exact_text) == Decimal(str(email["scores"]["RRMSE"]))
        assert Decimal(exact_text).quantize(Decimal(".001")) == ours
        receipt_path = folder / "submission_receipt.json"
        receipt_raw = receipt_path.read_bytes()
        receipt = json.loads(receipt_raw)
        assert receipt["submission_id"] == sid and receipt["evaluation_queue_id"] == queue
        assert receipt["status"] == "SCORED" and receipt["scores"] == email["scores"]
        assert receipt["score_evidence"]["sha256"] == sha(email_raw)
        candidate = Path(receipt["candidate_temporary_path"])
        assert sha(candidate.read_bytes()) == receipt["candidate_sha256"]
        others = [row for row in entries if row["submission_id"] != sid]
        better = [row for row in others if Decimal(row["RRMSE_display_text"]) < ours]
        tied = [row for row in others if Decimal(row["RRMSE_display_text"]) == ours]
        worse = [row for row in others if Decimal(row["RRMSE_display_text"]) > ours]
        best = min(Decimal(row["RRMSE_display_text"]) for row in others)
        next_better = max((Decimal(row["RRMSE_display_text"]) for row in better), default=None)
        next_worse = min((Decimal(row["RRMSE_display_text"]) for row in worse), default=None)
        rank = 1 + len(better)
        comparison = {
            "own_displayed_score": float(ours), "own_email_score": email["scores"]["RRMSE"],
            "own_email_score_text": exact_text, "rank": rank, "rank_denominator": total,
            "rank_span_including_displayed_ties": [rank, rank + len(tied)],
            "tie_count_including_own": 1 + len(tied), "tied_other_submitters": [r["submitter"] for r in tied],
            "better_other_count": len(better), "equal_other_count": len(tied), "worse_other_count": len(worse),
            "best_other_displayed_score": float(best),
            "best_other_submitters": [r["submitter"] for r in others if Decimal(r["RRMSE_display_text"]) == best],
            "displayed_gap_to_best_other": float(ours - best),
            "nearest_better_displayed_score": float(next_better) if next_better is not None else None,
            "displayed_gap_to_nearest_better": float(ours - next_better) if next_better is not None else None,
            "nearest_worse_displayed_score": float(next_worse) if next_worse is not None else None,
            "displayed_advantage_over_nearest_worse": float(next_worse - ours) if next_worse is not None else None,
        }
        record = {
            "schema_version": 2, "task": task, "own_submission_id": sid, "evaluation_queue_id": queue,
            "comparison_status": "VERIFIED_BY_OWN_SUBMISSION_ID", "mapping_verified": True, "own_listed": True,
            "uploaded_entity_id": receipt["uploaded_entity_id"], "uploaded_entity_version": receipt["uploaded_entity_version"],
            "candidate_sha256": receipt["candidate_sha256"], "source_url": url,
            "observed_at_utc": datetime.fromtimestamp(ax_path.stat().st_mtime, timezone.utc).isoformat(),
            "observation_time_basis": "Filesystem save time of root-provided fresh authenticated full AX capture; all corrected own IDs now visible",
            "recorded_at_utc": datetime.now(timezone.utc).isoformat(), "entries": entries,
            "visible_total_rows": total, "distinct_displayed_submitters": total, "other_displayed_competitors": total - 1,
            "metric": "RRMSE; lower is better.", "comparisons": {"RRMSE": comparison},
            "own_scores_from_email": email["scores"], "own_scores_from_email_text": {"RRMSE": exact_text},
            "email_score_time_utc": email["status_at_utc"], "email_message_id": email["status_email_id"],
            "rank_basis": "Competition rank = one plus count of displayed RRMSE values strictly lower; all displayed rows parsed, each a distinct submitter",
            "precision_note": "Ranks, apparent ties, and gaps use rounded displayed scores; exact own email RRMSE is separately retained. Hidden competitor precision is unavailable.",
            "scope": "Current task-specific final dashboard only; no overall challenge rank or final closed ranking asserted.",
            "table_evidence": {"path": str(ax_path), "sha256": sha(ax_raw)},
            "email_evidence": {"path": str(email_path), "sha256": sha(email_raw)},
            "submission_receipt_evidence": {"path": str(receipt_path), "sha256": sha(receipt_raw)},
            "model_or_prediction_changes": False, "truth_access": False,
            "source_sha256": sha(Path(__file__).read_bytes()),
        }
        payload = (json.dumps(record, indent=2, sort_keys=True) + "\n").encode()
        for destination in (folder, REPO / f"task{task}"):
            target = destination / "relative_score_record_current.json"
            assert not target.exists(), f"Preserve prior record: {target}"
            target.write_bytes(payload)
            assert target.read_bytes() == payload
        repo_ax = REPO / f"task{task}" / ax_path.name
        if repo_ax.exists() and repo_ax.read_bytes() != ax_raw:
            previous = repo_ax.read_bytes()
            archive = repo_ax.with_name(f"leaderboard_prior_{sha(previous)[:16]}.ax.txt")
            assert not archive.exists() or archive.read_bytes() == previous
            archive.write_bytes(previous)
        repo_ax.write_bytes(ax_raw)
        assert ax_path.read_bytes() == ax_raw and email_path.read_bytes() == email_raw
        assert receipt_path.read_bytes() == receipt_raw
        print(json.dumps({"task": task, "observed_at_utc": record["observed_at_utc"],
                          "comparison": comparison, "record_sha256": sha(payload)}))


if __name__ == "__main__":
    main()
