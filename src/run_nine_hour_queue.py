"""Sequential 9-hour compute queue. One Python job at a time, workers=1."""
from __future__ import annotations

import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
LOG = ROOT / "results" / "limitation_suite"
KERNEL = r"C:\Users\Hp\Downloads\Zenith"
DEADLINE_S = 9 * 3600


def remaining(deadline: float) -> float:
    return deadline - time.time()


def run_step(name: str, cmd: list[str], deadline: float) -> int:
    left = remaining(deadline)
    if left < 120:
        print(f"SKIP {name}: {left:.0f}s left in 9h window", flush=True)
        return -2
    LOG.mkdir(parents=True, exist_ok=True)
    log_file = LOG / f"{name}.log"
    header = f"\n{'='*60}\n{time.strftime('%Y-%m-%d %H:%M:%S')} START {name}  left={left/3600:.2f}h\n{'='*60}\n"
    print(header, flush=True)
    with log_file.open("a", encoding="utf-8") as log:
        log.write(header)
        log.flush()
        p = subprocess.Popen(cmd, cwd=str(ROOT), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        assert p.stdout is not None
        for line in p.stdout:
            print(line, end="", flush=True)
            log.write(line)
            if remaining(deadline) <= 0:
                print(f"\nDEADLINE: stopping {name} pid={p.pid}", flush=True)
                p.terminate()
                try:
                    p.wait(timeout=30)
                except subprocess.TimeoutExpired:
                    p.kill()
                log.write("\nSTOPPED by 9h deadline\n")
                return -1
        rc = p.wait()
        log.write(f"\nEXIT {name} rc={rc}\n")
        print(f"\nEXIT {name} rc={rc}\n", flush=True)
        return rc


def main() -> None:
    py = sys.executable
    deadline = time.time() + DEADLINE_S
    status_path = LOG / "nine_hour_status.json"
    LOG.mkdir(parents=True, exist_ok=True)

    # Preserve published 20-clone Mercurius result before n=100 overwrite
    src20 = ROOT / "results" / "th173_integrator_check.json"
    bak20 = ROOT / "results" / "th173_integrator_check_n20.json"
    if src20.exists() and not bak20.exists():
        shutil.copy2(src20, bak20)

    steps = [
        (
            "th173_100x40",
            [py, str(SRC / "run_fullscale.py"), "--phase", "th173_100x40", "--workers", "1", "--kernel-dir", KERNEL],
        ),
        (
            "mercurius_100",
            [py, str(SRC / "run_th173_integrator_check.py"), "--n-clones", "100", "--t-max-myr", "10", "--kernel-dir", KERNEL],
        ),
        (
            "cov_th173",
            [py, str(SRC / "run_covariance_campaign.py"), "--campaign", "cov_th173", "--kernel-dir", KERNEL],
        ),
        (
            "cov_qp112",
            [py, str(SRC / "run_covariance_campaign.py"), "--campaign", "cov_qp112", "--kernel-dir", KERNEL],
        ),
        (
            "exp_c_clones",
            [py, str(SRC / "run_exp_c_clones.py"), "--workers", "1", "--kernel-dir", KERNEL],
        ),
        (
            "d20_census",
            [py, str(SRC / "run_fullscale.py"), "--phase", "d20", "--workers", "1", "--kernel-dir", KERNEL],
        ),
        (
            "reclass_fullscale_post",
            [py, str(SRC / "reclass_fullscale_v1.py")],
        ),
    ]

    status = {"started": time.time(), "deadline": deadline, "steps": {}}
    status_path.write_text(json.dumps(status, indent=2), encoding="utf-8")
    for name, cmd in steps:
        rc = run_step(name, cmd, deadline)
        status["steps"][name] = {"rc": rc, "ts": time.time(), "hours_left": remaining(deadline) / 3600}
        status_path.write_text(json.dumps(status, indent=2), encoding="utf-8")
        if remaining(deadline) <= 0:
            print("9h window closed.", flush=True)
            break
    status["finished"] = time.time()
    status_path.write_text(json.dumps(status, indent=2), encoding="utf-8")
    print(f"Queue done. Status: {status_path}", flush=True)


if __name__ == "__main__":
    main()
