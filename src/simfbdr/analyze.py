"""Analyze a run: every replicate of its batch not yet analyzed.

A run chooses only what to read and how to analyze it. It inherits the batch's hyperpriors and
has no way to set them, so the generating priors and the analysis's cannot drift. Its manifest
records the batch config it ran against, and a config edited after its output exists is
refused rather than mixed.

A replicate whose log is complete is skipped, so a rerun fills gaps, and a run under a batch
that has been extended analyzes only the new replicates.
"""
import fcntl
import os
import random
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor

from . import config, manifest, tools
from .batch import die, replicates, warn_if_dirty

TEMPLATES = {"fbdrp": "infer-fbdrp.Rev",      # dnFBDRP: accounts for unsampled lineages
             "bds": "infer-bds.Rev",         # dnBDS: complete lineage sampling, PyRate's assumption
             "fbdsp": "infer-fbdsp.Rev"}     # dnFBDSP: the range process on a sampled budding tree

# what a template is handed for a setting the batch left out
UNSET = {"PRESENT": "0", "GAMMA_L_RATE": "-1", "GAMMA_M_RATE": "-1", "GAMMA_M_SHAPE": "-1",
         "TURN_A": "0", "TURN_B": "0", "TURN_LMEAN": "0", "TURN_LSD": "0",
         "GAMMA_RATE_MAX": "20"}
PRIORS = ["GAMMA_L_SHAPE", "GAMMA_L_RATE", "GAMMA_M_SHAPE", "GAMMA_M_RATE", "GAMMA_P_SHAPE",
          "TURN_A", "TURN_B", "TURN_LMEAN", "TURN_LSD",
          "GAMMA_P_RATE", "GAMMA_RATE_MAX", "AGE_MIN", "AGE_MAX", "PRESENT"]


def lines_in(path):
    with open(path, "rb") as fh:
        return sum(1 for _ in fh)


def max_first_appearance(taxa):
    """The oldest bin floor across taxa: PyRate's max FA, and the tightest lower bound the
    data puts on the origin. Only used when the origin prior is exponential."""
    oldest = {}
    with open(taxa) as fh:
        next(fh, None)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            try:
                oldest[f[0]] = max(oldest.get(f[0], 0.0), float(f[1]))
            except (IndexError, ValueError):
                continue
    x = max(oldest.values(), default=0.0)
    return f"{x if x > 0 else 1:.6f}"


def midpoints_on_a_shift(batch, occurrences, n):
    """Replicates in which some bin's midpoint is a rate-shift time. Such an age is in neither
    interval, so rb finds no starting state. It is every replicate when a grid as wide as the
    rate intervals is offset by half a bin."""
    hit = []
    for rep in range(1, n + 1):
        with open(os.path.join(batch, "times", f"times_{rep}.tsv")) as fh:
            shifts = [float(t) for t in fh.read().split()]
        with open(os.path.join(batch, occurrences, f"taxa_{rep}.tsv")) as fh:
            next(fh, None)
            ages = {(float(f[1]) + float(f[2])) / 2 for f in (l.split("\t") for l in fh)
                    if f[1] != f[2]}
        if any(abs(a - t) < 1e-9 for a in ages for t in shifts):
            hit.append(rep)
    return hit


def treat_ages(taxa, out, how, seed):
    """Write an occurrences table with every occurrence at an exact age: its bin's midpoint, or
    a uniform draw within the bin. This is what an analysis without age ranges is given."""
    rng = random.Random(seed)
    with open(taxa) as src, open(out, "w") as dst:
        dst.write(next(src))
        for line in src:
            f = line.rstrip("\n").split("\t")
            lo, hi = float(f[1]), float(f[2])
            if lo != hi:
                age = (lo + hi) / 2 if how == "midpoint" else rng.uniform(lo, hi)
                f[1] = f[2] = f"{age:.6f}"
            dst.write("\t".join(f) + "\n")


BURN = 0.25      # what summarize.R and ranks.py drop from the head of a log


def pool_draws(logs, out):
    """One log standing for several draws' chains, the length of one of them, so whatever
    reads a replicate's log reads the pooled posterior. Its head is the first draw's burn-in,
    which every reader drops. The rest takes each sample after burn-in from the next draw in
    turn, so each draw counts the same and no stretch of the log is one draw's alone: a
    single-chain effective sample size would read a block per draw as autocorrelation."""
    chains = []
    for path in logs:
        with open(path) as fh:
            header, *rows = fh.readlines()
        chains.append(rows)
    burn = int(len(chains[0]) * BURN)
    kept = [chains[i % len(chains)][i] for i in range(burn, len(chains[0]))]
    with open(out, "w") as fh:
        fh.writelines([header] + chains[0][:burn] + kept)


def analyze(run, force=False, jobs=1, slots=None, slot_dir=None):
    run = os.path.normpath(run)
    batch = os.path.dirname(os.path.dirname(run))
    if os.path.basename(os.path.dirname(run)) != "runs":
        die("analyze takes a run, <batch>/runs/<name>")
    bcfg, rcfg = os.path.join(batch, "config.yml"), os.path.join(run, "config.yml")
    for f in (bcfg, rcfg):
        if not os.path.exists(f):
            die(f"no {f}")
    b = config.load(bcfg, "batch")
    r = config.load(rcfg, "run")
    b.pop("_uses")
    stated = r.pop("_uses")
    # which rb: the one the config states, else the one that analyzed this run before
    before = manifest.find(run)[1] or {}
    tools.choose_rb(stated=manifest.used(stated, "rb"),
                    recorded=manifest.used(before.get("uses"), "rb"),
                    config=rcfg, run=run, force=force)
    tools.require(stated, rcfg)
    bhash, rhash = manifest.digest(bcfg), manifest.digest(rcfg)
    _, bdoc = manifest.find(batch)
    if bdoc is None or not replicates(batch):
        die(f"{batch} has no replicates yet: simfbdr extend {batch} <how many>")
    simulated_from = manifest.flat(batch).get("config_hash")
    if simulated_from != bhash:
        die(f"{bcfg} changed since its replicates were simulated ({simulated_from} -> {bhash}).\n"
            f"       Analyzing them under it would describe data that was not drawn that way. Make a new batch.")
    n = replicates(batch)

    # rho is the generating parameter and the one the analysis assumes. A batch that sets it wins, and a run
    # that disagrees is refused rather than silently assuming a rho the data was not drawn under.
    if "RHO" in b and float(b["RHO"]) != float(r["RHO"]):
        die(f"{rcfg} assumes an extant sampling probability of {r['RHO']}, and {batch} was "
            f"generated under {b['RHO']}.")
    if r["MODEL"] != "bds" and r["DROP_PRESENT_ONLY"] == "true":
        die(f"{rcfg} sets inference.drop_present_only_taxa, which only the bds model reads.")
    if r["MODEL"] == "bds" and r["RANGE_START"] == "divergence":
        die(f"{rcfg} asks for divergence times, which dnBDS does not offer: under complete "
            f"lineage sampling a birth is the species' own origination.")
    mix_k = int(b.get("MIX_K") or 1)
    if r["PSI_MIX"] == "true" and r["MODEL"] not in ("fbdrp", "fbdsp"):
        die(f"{rcfg} asks for the psi mixture, which only dnFBDRP and dnFBDSP have.")
    if r["PSI_MIX"] == "true" and mix_k < 2:
        die(f"{rcfg} asks for the psi mixture, and {batch} was simulated without one: the "
            f"analysis takes its components and their prior from the batch.")
    if r["CHARS"] == "true" and r["MODEL"] != "fbdsp":
        die(f"{rcfg} asks for the characters, which need the tree process: model fbdsp.")
    if r["CHARS"] == "true" and not b.get("N_CHARS"):
        die(f"{rcfg} asks for the characters, and {batch} simulated none.")
    if r["ANA"] == "true" and r["MODEL"] != "fbdsp":
        die(f"{rcfg} asks for anagenesis, which only the tree process has: model fbdsp.")
    if r["ANA"] == "true" and b.get("ANAGENESIS") != "true":
        die(f"{rcfg} asks for anagenesis, and {batch} was simulated without it: the analysis "
            f"takes the anagenetic rate's prior from the batch.")
    if r["MODEL"] == "fbdsp" and r["ORIGIN_PRIOR"] == "none":
        die(f"{rcfg} supplies no origin, which dnFBDSP needs.")
    if r["MODEL"] == "fbdsp" and r["RANGE_START"] == "divergence":
        die(f"{rcfg} asks for divergence times, which dnFBDSP does not report: on a tree a "
            f"range's start is its own origination.")
    template = tools.script(TEMPLATES[r["MODEL"]])
    occurrences = f"occurrences_{r['REPORTING']}"
    if not os.path.isdir(os.path.join(batch, occurrences)):
        die(f"the batch has no {occurrences}")
    sky = "true" if int(b["NINTERVALS"]) > 1 else "false"
    if r["AGES"] == "midpoint":
        hit = midpoints_on_a_shift(batch, occurrences, n)
        if hit:
            die(f"in {len(hit)} of {n} replicates of {occurrences} a bin's midpoint is a rate-shift "
                f"time, which is in neither interval, so those chains cannot start. Read a grid "
                f"whose midpoints miss the shifts.")

    path, _ = manifest.find(run)
    if path is not None and not force:
        was = manifest.flat(run)
        if was.get("batch_config_hash") != bhash:
            die(f"{run} ran against batch config {was.get('batch_config_hash')}, now {bhash}. "
                f"Its output would mix. Delete output/ or pass --force.")
        if was.get("config_hash") != rhash:
            die(f"{rcfg} changed since its output was written ({was.get('config_hash')} -> {rhash}). "
                f"Delete output/ or pass --force.")

    out, aux = os.path.join(run, "output"), os.path.join(run, "aux")
    for d in (out, aux, os.path.join(run, "results")):
        os.makedirs(d, exist_ok=True)
    print(f"{os.path.basename(run)}: batch={os.path.basename(batch)} reads={r['REPORTING']} "
          f"ages {r['AGES']} reporting {r['INFER']} condition {r['COND']} rho={r['RHO']} reps={n}", flush=True)
    warn_if_dirty()
    tools.announce()

    # a finished replicate's log has this many lines; a rerun fills gaps rather than redo work
    full = int(r["GENS"]) // int(r["PRINTGEN"]) + 2
    failures = os.path.join(run, "failures.log")
    open(failures, "w").close()
    absolute = os.path.abspath
    rb, env = tools.rb(), tools.environment()
    slots = int(slots or os.environ.get("SLOTS") or 0)
    slot_dir = slot_dir or os.environ.get("SLOT_DIR") or os.path.join(os.path.dirname(batch), ".slots")
    if slots:
        os.makedirs(slot_dir, exist_ok=True)

    draws = int(r["AGE_DRAWS"]) if r["AGES"] == "drawn" else 1

    def chain(rep, draw):
        """One chain: the replicate as reported, or one draw of its ages. True if it finished."""
        # a draw beyond the first keeps its files apart, and its log is pooled, not read
        work = os.path.join(aux, f"draw_{draw}") if draws > 1 else aux
        to = work if draws > 1 else out
        os.makedirs(work, exist_ok=True)
        log = os.path.join(to, f"rep_{rep}.log")
        if os.path.exists(log) and lines_in(log) == full:
            return True
        # one self-contained file per replicate: the settings as Rev variables, then the template
        v = {**UNSET, **b, **r}
        taxa = os.path.join(batch, occurrences, f"taxa_{rep}.tsv")
        treated = []
        if r["AGES"] != "reported":
            src, taxa = taxa, absolute(os.path.join(work, f"taxa_{rep}.tsv"))
            treat_ages(src, taxa, r["AGES"], f"{r['SEED']}:{rep}:{draw}")
            treated = [f'TAXA <- "{taxa}"']
        # seeded, so the same replicate analyzed twice gives the same chain. Each replicate has
        # its own seed, or their Monte Carlo errors would move together, and so has each draw.
        head = [f"seed({int(r['SEED']) + rep + 1000003 * (draw - 1)})",
                f'rep <- "{rep}"', f'ORIGIN_PRIOR <- "{r["ORIGIN_PRIOR"]}"',
                f"MAXFA <- {max_first_appearance(taxa)}",
                f"SKYLINE <- {sky}", f'REPORTING <- "{r["INFER"]}"', f'COND <- "{r["COND"]}"',
                f"RHO <- {r['RHO']}", f"GENS <- {r['GENS']}", f"PRINTGEN <- {r['PRINTGEN']}",
                f'BATCHDIR <- "{absolute(batch)}"', f'OCCURRENCES <- "{occurrences}"',
                f'OUTDIR <- "{absolute(to)}"']
        head += [f"{k} <- {v[k]}" for k in PRIORS]
        # which extinction prior the batch drew from: its own gammas, or a turnover shared by all
        turnover = "beta" if b.get("TURN_A") else "lognormal" if b.get("TURN_LMEAN") else "none"
        head += [f'TURNOVER <- "{turnover}"']
        # the coexistence factor belongs to the published model, so a config that leaves it
        # out gets the process's own default: dnFBDRP weights each birth, dnBDS does not
        coex = r.get("COEXISTENCE") or ("false" if r["MODEL"] == "bds" else "true")
        head += [f"DROP_PRESENT_ONLY <- {r['DROP_PRESENT_ONLY']}",
                 f"LATENTS <- {r['LATENTS']}", f"COEXISTENCE <- {coex}",
                 f'RANGE_START <- "{r["RANGE_START"]}"']
        # only when the config asks for it, so a template running against an rb that has no
        # such argument never sees it
        if r.get("REQUIRE_PARENT"):
            head += [f"REQUIRE_PARENT <- {r['REQUIRE_PARENT']}"]
        # the psi mixture and the characters take their hyperpriors from the batch, which drew them
        head += [f"PSI_MIX <- {r['PSI_MIX']}", f"MIX_K <- {mix_k}",
                 f"MIX_ALPHA_RATE <- {b.get('MIX_ALPHA_RATE', 1)}",
                 f'MIX_MODEL <- "{b.get("MIX_MODEL") or "gamma"}"', f"MIX_W_CONC <- {b.get('MIX_W_CONC') or 1}",
                 f"CHARS <- {r['CHARS']}",
                 f"ANAGENESIS <- {r['ANA']}"]
        if r["CHARS"] == "true":
            n_chars = int(b["N_CHARS"])
            states = str(b.get("CHAR_STATES", "2")).split(",")
            states = states * n_chars if len(states) == 1 else states
            charfile = absolute(os.path.join(batch, "characters", f"characters_{rep}.nex"))
            head += [f'CHARFILE <- "{charfile}"', f"NSTATES <- v({', '.join(states)})",
                     f"EPS_A <- {b.get('EPS_A', 1)}", f"EPS_B <- {b.get('EPS_B', 9)}"]
        head += [f'BINDIR <- "{tools.SCRIPTS}"'] + treated
        head += [f'source("{template}")']
        rev = os.path.join(work, f"run_{rep}.Rev")
        with open(rev, "w") as fh:
            fh.write("\n".join(head) + "\n")
        # slots cap chains across every run sharing the directory, so concurrent runs fill the
        # cores between them: each chain holds one slot lock for as long as rb runs
        held = None
        while slots and held is None:
            for s in range(1, slots + 1):
                fh = open(os.path.join(slot_dir, f"slot_{s}"), "w")
                try:
                    fcntl.flock(fh, fcntl.LOCK_EX | fcntl.LOCK_NB)
                    held = fh
                    break
                except OSError:
                    fh.close()
            if held is None:
                time.sleep(10)
        try:
            with open(os.path.join(work, f"rb_{rep}.out"), "w") as so:
                subprocess.run([rb, absolute(rev)], stdin=subprocess.DEVNULL, stdout=so,
                               stderr=subprocess.STDOUT, env=env)
        finally:
            if held is not None:
                held.close()
        if os.path.exists(log) and lines_in(log) == full:
            return True
        why = "no log written"
        for line in open(os.path.join(work, f"rb_{rep}.out"), errors="replace"):
            if "error" in line.lower():
                why = line.strip()
                break
        with open(failures, "a") as fh:
            fh.write(f"{rep} :: {why}" + (f" (draw {draw})" if draws > 1 else "") + "\n")
        return False

    def one(rep):
        log = os.path.join(out, f"rep_{rep}.log")
        if os.path.exists(log) and lines_in(log) == full:
            return
        done = [chain(rep, k) for k in range(1, draws + 1)]
        if draws > 1 and all(done):
            pool_draws([os.path.join(aux, f"draw_{k}", f"rep_{rep}.log")
                        for k in range(1, draws + 1)], log)

    with ThreadPoolExecutor(max_workers=max(1, int(jobs))) as pool:
        list(pool.map(one, range(1, n + 1)))

    _, rb_md5, rb_describe = tools.rb_identity()
    manifest.write_run(run, {"batch_dir": batch, "replicates": n, "occurrences": occurrences,
                             "inputs": {"config.yml": rhash,
                                        "batch": {"config.yml": bhash,
                                                  occurrences: manifest.digest(os.path.join(batch, occurrences))},
                                        "simfbdr": {template.name: manifest.digest(template)}},
                             "simfbdr_commit": tools.version(),
                             "rb_md5": rb_md5, "rb_commit": rb_describe})
    done = sum(1 for f in os.listdir(out) if f.startswith("rep_") and f.endswith(".log"))
    print(f"logs: {done}/{n}, {lines_in(failures)} failed")
    print(f"manifest: {os.path.join(run, 'manifest.yml')}")
    print(f"now run: simfbdr summarize {run}")
    return lines_in(failures) == 0
