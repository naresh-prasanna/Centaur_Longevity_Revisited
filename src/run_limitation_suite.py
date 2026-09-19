"""Master orchestrator for limitation-suite computations."""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
LOG = ROOT / "results" / "limitation_suite"
KERNEL = r"C:\Users\Hp\Downloads\Zenith"


def run_step(name: str, cmd: list[str], log_file: Path) -> int:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    header = f"\n{'='*60}\n{time.strftime('%Y-%m-%d %H:%M:%S')} START {name}\n{'='*60}\n"
    print(header, flush=True)
    with log_file.open("a", encoding="utf-8") as log:
        log.write(header)
        log.flush()
        p = subprocess.Popen(
            cmd,
            cwd=str(ROOT),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        assert p.stdout is not None
        for line in p.stdout:
            print(line, end="", flush=True)
            log.write(line)
        rc = p.wait()
        footer = f"\nEXIT {name} rc={rc}\n"
        log.write(footer)
        print(footer, flush=True)
    return rc


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--kernel-dir", default=KERNEL)
    ap.add_argument("--skip-fast", action="store_true")
    ap.add_argument("--skip-heavy", action="store_true")
    args = ap.parse_args()
    py = sys.executable

    status_path = LOG / "status.json"
    status = {"started": time.time(), "steps": {}}

    def record(name: str, rc: int) -> None:
        status["steps"][name] = {"rc": rc, "ts": time.time()}
        status_path.write_text(json.dumps(status, indent=2), encoding="utf-8")

    fast = [
        ("reclass_v2_json", [py, str(SRC / "reclass_v2.py")]),
        ("reclass_fullscale_v1", [py, str(SRC / "reclass_fullscale_v1.py")]),
        ("derive_22myr_cut", [py, str(SRC / "derive_22myr_cut.py")]),
    ]

    heavy = [
        (
            "th173_100x40",
            [py, str(SRC / "run_fullscale.py"), "--phase", "th173_100x40", "--workers", str(args.workers), "--kernel-dir", args.kernel_dir],
        ),
        (
            "d20_census",
            [py, str(SRC / "run_fullscale.py"), "--phase", "d20", "--workers", str(args.workers), "--kernel-dir", args.kernel_dir],
        ),
        (
            "mercurius_100",
            [py, str(SRC / "run_th173_integrator_check.py"), "--n-clones", "100", "--t-max-myr", "10", "--kernel-dir", args.kernel_dir],
        ),
        (
            "exp_c_clones",
            [py, str(SRC / "run_exp_c_clones.py"), "--workers", str(args.workers), "--kernel-dir", args.kernel_dir],
        ),
        (
            "covariance_all",
            [py, str(SRC / "run_covariance_campaign.py"), "--campaign", "all", "--kernel-dir", args.kernel_dir],
        ),
        (
            "arc_tiers_10",
            [
                py,
                str(SRC / "run_pilot.py"),
                "--production",
                "--n-tiers",
                "10",
                "--n-clones",
                "20",
                "--t-max-myr",
                "25",
                "--kernel-dir",
                args.kernel_dir,
            ],
        ),
    ]

    steps = []
    if not args.skip_fast:
        steps.extend(fast)
    if not args.skip_heavy:
        steps.extend(heavy)

    for name, cmd in steps:
        rc = run_step(name, cmd, LOG / f"{name}.log")
        record(name, rc)

    # Post-pass reclass including new campaigns
    if not args.skip_fast:
        rc = run_step("reclass_fullscale_post", [py, str(SRC / "reclass_fullscale_v1.py")], LOG / "reclass_fullscale_post.log")
        record("reclass_fullscale_post", rc)

    status["finished"] = time.time()
    status_path.write_text(json.dumps(status, indent=2), encoding="utf-8")
    print(f"\nSuite complete. Status: {status_path}")


if __name__ == "__main__":
    main()
