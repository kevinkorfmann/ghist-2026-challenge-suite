"""Prepare new batch files only; preserve older runs and source snapshots."""
from pathlib import Path

for name in ["run_development", "fit", "confirm", "summarize"]:
    text = Path(name + "_v3.sbatch").read_text().replace("_v3", "_v4")
    if name == "run_development":
        text = text.replace('--baseline "$TASK_ROOT/results_v2/', '--baseline "$TASK_ROOT/results_v3/')
    path = Path(name + "_v4.sbatch")
    with path.open("x") as handle:
        handle.write(text)
text = Path("fit_v4.sbatch").read_text()
text = text[:text.rindex('"$PROJECT_ROOT/.venv/bin/python"')]
text += '"$PROJECT_ROOT/.venv/bin/python" "$TASK_ROOT/phase_development_ablation.py" "$TASK_ROOT"\n'
with Path("ablation_v4.sbatch").open("x") as handle:
    handle.write(text)
