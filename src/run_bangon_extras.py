"""After TH173 100×40 finishes, run the one extra check that actually tightens the headline.

Full SBDB-covariance clones of TH173, 20 clones × 10 Myr. That tests whether
ignoring element correlations can restore Q. Everything else stays a limitation.
"""
from __future__ import annotations

import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TARGET = ROOT / "results" / "fullscale" / "th173_100x40" / "2005_TH173.json"
LOG = ROOT / "results" / "limitation_suite" / "bangon_extras.log"


def main() -> None:
    LOG.parent.mkdir(parents=True, exist_ok=True)
    print(f"waiting for {TARGET}", flush=True)
    while not TARGET.exists():
        time.sleep(30)
    print("TH173 100×40 complete; starting covariance Q-test (20 clones, 10 Myr, 1 process)", flush=True)
    cmd = [
        sys.executable,
        str(ROOT / "src" / "run_covariance_campaign.py"),
        "--campaign",
        "cov_th173_q10",
        "--kernel-dir",
        r"C:\Users\Hp\Downloads\Zenith",
    ]
    with LOG.open("a", encoding="utf-8") as log:
        p = subprocess.Popen(cmd, cwd=str(ROOT), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        assert p.stdout is not None
        for line in p.stdout:
            print(line, end="", flush=True)
            log.write(line)
        rc = p.wait()
    print(f"extras exit {rc}", flush=True)


if __name__ == "__main__":
    main()
