"""Read-only bounded monitor of this isolated experiment's retained outputs."""
import datetime
import json
import sys
import time
from collections import Counter
from pathlib import Path

root = Path(sys.argv[1]).resolve()
last_count, reported_fit, reported_audit = -16, False, False
while datetime.datetime.now(datetime.UTC) < datetime.datetime(2026, 11, 20, 22, 20, tzinfo=datetime.UTC):
    paths = sorted((root / "results_crossed_v1").glob("*/result.json"))
    records = [json.loads(path.read_text()) for path in paths]
    failed = [r["cohort"] for r in records if not r["completed"]]
    if len(records) >= last_count + 16 or failed:
        print(json.dumps({"utc": datetime.datetime.now(datetime.UTC).isoformat(), "completed": len(records),
            "splits": dict(Counter(r["cohort"]["split"] for r in records)), "failed": failed}), flush=True)
        last_count = len(records)
    if failed:
        raise SystemExit(1)
    fit = root / "crossed_context_fit_v1/report.json"
    if fit.exists() and not reported_fit:
        report = json.loads(fit.read_text())
        print(json.dumps({"fit_complete": True, "reports": report["reports"]}), flush=True)
        reported_fit = True
    audit = root / "crossed_development_audit_v1/audit.json"
    if audit.exists() and not reported_audit:
        data = json.loads(audit.read_text())
        print(json.dumps({"audit_complete": True, "fresh_stress_launch_allowed": data["fresh_stress_launch_allowed"],
                          "summaries": data["summaries"]}), flush=True)
        reported_audit = True
    archive = root / "crossed_development_compact_v1.tar.gz"
    manifest = root / "crossed_all_raw_metadata_v1/manifest.json"
    if archive.exists() and manifest.exists():
        print(json.dumps({"export_started": True, "archive": str(archive)}), flush=True)
        break
    time.sleep(45)
else:
    print(json.dumps({"hard_monitor_end_reached": True}), flush=True)
