#!/usr/bin/env python3
"""perf_metrics.py — turn run.sh's raw `perf stat` dumps into a report table.

Parses results/stage{0..3}.perf (the `perf stat -r 3` output) and prints, per
stage, the derived metrics that actually tell the cache-optimization story:

    IPC            instructions / cycles            (higher = fewer stalls)
    cache-miss %   cache-misses / cache-references  (lower  = better locality)
    L1-miss %      L1-dcache-load-misses / loads
    LLC-miss %     LLC-load-misses / LLC-loads
    MPKI (LLC)     LLC-load-misses / instr * 1000   (misses per 1k insns)
    dTLB-miss %    dTLB-load-misses / dTLB-loads

Then a comparison block vs stage 0 (speedup + % change per metric). If a
matching stage{N}.topdown file exists, its memory-bound % is shown too.

Usage: perf_metrics.py [results_dir]   (default: ./results)
Output goes to stdout and is also written to <results_dir>/metrics.md.
"""
import os
import re
import sys

# Raw counters we pull straight out of the perf stat dump.
COUNTERS = [
    "cycles", "instructions",
    "cache-references", "cache-misses",
    "L1-dcache-loads", "L1-dcache-load-misses",
    "LLC-loads", "LLC-load-misses",
    "dTLB-loads", "dTLB-load-misses",
]

STAGE_LABELS = {
    0: "naive AoS",
    1: "SoA",
    2: "SoA + skewing",
    3: "SoA + skewing + SIMD",
}


def norm_event(tok):
    """Strip a hybrid-PMU prefix (cpu_core/cycles/ -> cycles) and any :u suffix."""
    tok = tok.strip().rstrip(":u")
    if "/" in tok:
        parts = [p for p in tok.split("/") if p]
        tok = parts[-1] if parts else tok
    return tok


def parse_perf(path):
    """Return (counters dict, elapsed_seconds or None) from one perf stat file."""
    counters = {}
    elapsed = None
    counter_re = re.compile(r"^([\d.,]+)\s+(\S+)")
    # With `-r N` the line is "<mean> +- <stddev> seconds time elapsed"; grab the
    # leading mean, not the stddev that sits just before the words.
    elapsed_re = re.compile(r"^([\d.,]+)\s+(?:\+-\s+[\d.,]+\s+)?seconds time elapsed")
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            m = elapsed_re.search(s)
            if m:
                elapsed = float(m.group(1).replace(",", ""))
                continue
            m = counter_re.match(s)
            if not m:
                continue
            name = norm_event(m.group(2))
            if name not in COUNTERS:
                continue
            try:
                val = float(m.group(1).replace(",", ""))
            except ValueError:
                continue
            counters[name] = counters.get(name, 0.0) + val  # sum split PMUs
    return counters, elapsed


def memory_bound_pct(results_dir, stage):
    """Best-effort: pull a memory-bound % from stage{N}.topdown if present."""
    path = os.path.join(results_dir, f"stage{stage}.topdown")
    if not os.path.exists(path):
        return None
    try:
        with open(path) as fh:
            for line in fh:
                low = line.lower()
                if "memory bound" in low or "tma_memory_bound" in low:
                    m = re.search(r"([\d.]+)\s*%", line)
                    if m:
                        return float(m.group(1))
    except OSError:
        return None
    return None


def derive(c, elapsed, mem_bound):
    """Compute the derived metrics from raw counters. Missing -> None."""
    def ratio(a, b, scale=1.0):
        if c.get(a) and c.get(b):
            return scale * c[a] / c[b]
        return None
    return {
        "time_s":   elapsed,
        "IPC":      ratio("instructions", "cycles"),
        "cache%":   ratio("cache-misses", "cache-references", 100.0),
        "L1%":      ratio("L1-dcache-load-misses", "L1-dcache-loads", 100.0),
        "LLC%":     ratio("LLC-load-misses", "LLC-loads", 100.0),
        "MPKI":     ratio("LLC-load-misses", "instructions", 1000.0),
        "dTLB%":    ratio("dTLB-load-misses", "dTLB-loads", 100.0),
        "mem-bnd%": mem_bound,
    }


COLS = ["time_s", "IPC", "cache%", "L1%", "LLC%", "MPKI", "dTLB%", "mem-bnd%"]


def fmt(v):
    if v is None:
        return "n/a"
    return f"{v:.3f}" if v < 100 else f"{v:,.0f}"


def main():
    results_dir = sys.argv[1] if len(sys.argv) > 1 else "results"
    stages = {}
    for stage in range(4):
        path = os.path.join(results_dir, f"stage{stage}.perf")
        if not os.path.exists(path):
            continue
        counters, elapsed = parse_perf(path)
        if not counters and elapsed is None:
            continue
        stages[stage] = derive(counters, elapsed,
                               memory_bound_pct(results_dir, stage))

    if not stages:
        print("  No stage*.perf files found in", results_dir)
        return 0

    out = []
    out.append("## Derived perf metrics")
    out.append("")
    header = "| stage | layout | " + " | ".join(COLS) + " |"
    sep = "|" + "---|" * (len(COLS) + 2)
    out.append(header)
    out.append(sep)
    for stage in sorted(stages):
        row = stages[stage]
        cells = " | ".join(fmt(row[c]) for c in COLS)
        out.append(f"| {stage} | {STAGE_LABELS.get(stage, '')} | {cells} |")

    # Comparison vs stage 0.
    if 0 in stages:
        base = stages[0]
        for stage in sorted(stages):
            if stage == 0:
                continue
            cur = stages[stage]
            out.append("")
            out.append(f"### Stage {stage} ({STAGE_LABELS.get(stage,'')}) vs stage 0")
            if base["time_s"] and cur["time_s"]:
                spd = base["time_s"] / cur["time_s"]
                out.append(f"- speedup: **{spd:.2f}x**  "
                           f"({base['time_s']:.3f}s -> {cur['time_s']:.3f}s)")
            for c in ["IPC", "cache%", "L1%", "LLC%", "MPKI", "dTLB%", "mem-bnd%"]:
                b, v = base[c], cur[c]
                if b is None or v is None or b == 0:
                    continue
                delta = 100.0 * (v - b) / b
                arrow = "↑" if delta > 0 else "↓"
                out.append(f"- {c}: {fmt(b)} -> {fmt(v)}  ({arrow} {abs(delta):.0f}%)")

    text = "\n".join(out)
    print(text)
    try:
        with open(os.path.join(results_dir, "metrics.md"), "w") as fh:
            fh.write(text + "\n")
        print(f"\n  (also written to {os.path.join(results_dir, 'metrics.md')})")
    except OSError:
        pass
    return 0


if __name__ == "__main__":
    sys.exit(main())
