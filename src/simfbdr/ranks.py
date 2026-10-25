#!/usr/bin/env python3
"""Rank SBC for one analysis: the fraction of posterior draws below the truth, per rate.

Ranks, not point estimates: (est-true)/true divides by a drawn rate that can be near
zero and reads as a large bias even when the posterior is centred.

Usage:  simfbdr ranks <batch>/runs/<name>
"""
import os
import sys
import numpy as np

from . import manifest

ESS_MIN = 100.0
BURN = 0.25


def ess(x):
    """Single-chain effective sample size, Geyer initial positive sequence."""
    n = len(x)
    x = x - x.mean()
    if np.var(x) == 0:
        return 0.0
    f = np.fft.rfft(x, 2 * n)
    ac = np.fft.irfft(f * np.conjugate(f))[:n].real
    ac /= ac[0]
    t, s = 1, 0.0
    while t + 1 < n:
        p = ac[t] + ac[t + 1]
        if p <= 0:
            break
        s += p
        t += 2
    return n / (1.0 + 2.0 * s)


def collect(run):
    """Per-replicate ranks for one run: (rate names, replicate numbers, rank matrix, counts).

    A replicate is kept only where every rate's chain reached ESS_MIN, so two runs keep
    different replicates and have to be aligned on the ones both kept before differencing.
    """
    batch = os.path.dirname(os.path.dirname(run))
    tv = np.genfromtxt(os.path.join(batch, "true_vals.tsv"), delimiter="\t", names=True)
    ni = sum(1 for c in tv.dtype.names if c.startswith("lambda"))

    # batches with no true_vals_order key predate the fix and are oldest-first
    youngest = manifest.flat(batch).get("true_vals_order") == "youngest_first"
    idx = list(range(1, ni + 1)) if youngest else list(range(ni, 0, -1))
    rates = [f"{r}[{i}]" for r in ("lambda", "mu", "psi") for i in range(1, ni + 1)]
    tcols = [f"{r}{i}" for r in ("lambda", "mu", "psi") for i in idx]

    # a chain still being written is not a posterior, so require the full sample count
    conf = manifest.flat(run)
    if "GENS" not in conf:
        sys.exit(f"{run} has no manifest: run simfbdr analyze on it first")
    want = int(conf["GENS"]) // int(conf["PRINTGEN"]) + 1

    out = os.path.join(run, "output")
    ranks, reps, seen, partial = [], [], 0, 0
    for fn in sorted(os.listdir(out)):
        if not fn.endswith(".log") or "_" not in fn:
            continue
        seen += 1
        path = os.path.join(out, fn)
        with open(path) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
        try:
            cols = [hdr.index(c) for c in rates]
        except ValueError:
            continue
        try:
            a = np.loadtxt(path, delimiter="\t", skiprows=1, usecols=cols, ndmin=2)
        except Exception:
            continue
        if a.shape[0] != want:
            partial += 1
            continue
        a = a[int(a.shape[0] * BURN):]
        if min(ess(a[:, j]) for j in range(a.shape[1])) < ESS_MIN:
            continue
        reps.append(int(fn.rsplit("_", 1)[1][:-4]))
        ranks.append((a < np.array([tv[c][reps[-1] - 1] for c in tcols])).mean(axis=0))

    return rates, reps, np.array(ranks), (seen, partial, youngest)


def main(argv=None):
    run = (argv or sys.argv[1:])[0].rstrip("/")
    rates, reps, R, (seen, partial, youngest) = collect(run)
    kept = len(reps)

    print(f"\n{run}  truth {'youngest' if youngest else 'oldest'}-first"
          f"  kept {kept} of {seen}  ({partial} still running, {seen - partial - kept} low ESS)")
    if kept < 2:
        return
    se = 1.0 / np.sqrt(12 * kept)
    crit = 1.36 / np.sqrt(kept)
    # the mean alone cannot see a width mismatch, which is what a wrong prior scale gives
    print(f"mean rank 0.5 (SE {se:.3f}), sd 0.289, tails 0.20, KS below {crit:.3f}\n")
    print(f"{'':10} {'mean':>7} {'z':>7} {'sd':>7} {'tails':>7} {'KS':>7}")
    for j, name in enumerate(rates):
        c = np.sort(R[:, j])
        ks = max(np.abs(np.arange(1, kept + 1) / kept - c).max(),
                 np.abs(c - np.arange(kept) / kept).max())
        tail = ((R[:, j] < 0.1) | (R[:, j] > 0.9)).mean()
        m = R[:, j].mean()
        print(f"{name:10} {m:7.3f} {(m - 0.5) / se:+7.2f} {R[:, j].std(ddof=1):7.3f}"
              f" {tail:7.3f} {ks:7.3f}{'  <--' if ks > crit else ''}")
    print("\nrank below 0.5 means the posterior sits above the truth")


if __name__ == "__main__":
    main()
