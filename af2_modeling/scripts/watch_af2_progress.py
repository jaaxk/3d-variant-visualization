#!/usr/bin/env python3
"""Run the cluster's official AlphaFold2 (run-alphafold.py) as a subprocess
and show live progress -- an elapsed-time heartbeat plus timestamped
stage-change lines -- instead of a silent multi-hour/multi-day black box.

AlphaFold2's own CLI has no built-in progress bar, and at the scale we're
attempting here (7,207 residues -- well beyond any published AF2/AF-Multimer
benchmark), whether it finishes in a reasonable time at all is genuinely
unknown. This wrapper doesn't modify the read-only cluster install; it just
launches it and tails its log.

Stage markers below are exact log strings confirmed by reading
run_alphafold-msa-only.py (a modified copy of the same base predict_structure()
used inside the container, with the model-inference block commented out --
the comments preserve the exact log format the real script emits, e.g.
'Running model %s on %s' and 'Total JAX model %s on %s predict time'). Marker
matching is deliberately a case-insensitive substring search rather than an
exact-format parse, since we haven't observed this cluster's real AF2 log
output live yet -- calibrate/extend this list after the smoke test.

Usage:
    python watch_af2_progress.py --log-file run.log -- \\
        /share/apps/af2/pyenv/bin/python3 /share/apps/af2/2.3.2-20231225/run-alphafold.py \\
        --fasta_paths=... --output_dir=... --max_template_date=... ...
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

from tqdm import tqdm

STAGE_MARKERS = [
    "predicting ",
    "running model ",
    "total jax model",
    "launching subprocess",
    "jackhmmer",
    "hhblits",
    "hhsearch",
    "hmmsearch",
    "relax",
    "started ",
    "finished ",
    "error",
    "traceback",
    "oom",
    "out of memory",
]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--log-file", required=True)
    ap.add_argument("command", nargs=argparse.REMAINDER)
    args = ap.parse_args()

    command = args.command
    if command and command[0] == "--":
        command = command[1:]
    if not command:
        raise SystemExit("no command given after '--'")

    log_path = Path(args.log_file)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"[watch_af2_progress] command: {' '.join(command)}")
    print(f"[watch_af2_progress] logging to: {log_path}")
    sys.stdout.flush()

    start = time.monotonic()
    with log_path.open("w") as log_f:
        proc = subprocess.Popen(command, stdout=log_f, stderr=subprocess.STDOUT)

        last_seen_line = ""
        with log_path.open("r") as tail_f, tqdm(
            desc="elapsed", unit="s", bar_format="{desc}: {n:.0f}s [{elapsed}]"
        ) as bar:
            last_update = start
            while proc.poll() is None:
                line = tail_f.readline()
                if line:
                    lowered = line.lower()
                    if line.strip() != last_seen_line and any(m in lowered for m in STAGE_MARKERS):
                        elapsed = time.monotonic() - start
                        tqdm.write(f"[+{elapsed / 60:6.1f} min] {line.strip()}")
                        last_seen_line = line.strip()
                else:
                    time.sleep(2)
                now = time.monotonic()
                if now - last_update >= 5:
                    bar.n = now - start
                    bar.refresh()
                    last_update = now
            # Drain any remaining buffered lines after the process exits.
            for line in tail_f:
                lowered = line.lower()
                if line.strip() != last_seen_line and any(m in lowered for m in STAGE_MARKERS):
                    elapsed = time.monotonic() - start
                    tqdm.write(f"[+{elapsed / 60:6.1f} min] {line.strip()}")
            bar.n = time.monotonic() - start
            bar.refresh()

    returncode = proc.returncode
    total_elapsed = time.monotonic() - start
    print(f"[watch_af2_progress] process exited with code {returncode} "
          f"after {total_elapsed / 60:.1f} minutes ({total_elapsed / 3600:.2f} hours)")
    sys.exit(returncode)


if __name__ == "__main__":
    main()
