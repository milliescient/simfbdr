"""The simfbdr command.

    simfbdr extend <batch> <n>          add n replicates to a batch, or give a new one its first
    simfbdr analyze <batch>/runs/<name> analyze every replicate of the batch not yet analyzed
    simfbdr status <batch or run> ...   does it still match its manifest
    simfbdr summarize <run> ...         SBC ranks, coverage and diagnostics for a run
    simfbdr ranks <run>                 mean rank per parameter
    simfbdr config default [batch|run]  every setting a config can hold, with its default
    simfbdr script <name> [args]        run one of the other scripts simfbdr carries
    simfbdr path <name>                 where one of those scripts is, for something that sources it
    simfbdr doctor                      the tools this machine would run, and whether they are there
    simfbdr uses [name]                 what it would run, as name@version, or one version
    simfbdr manifest flat|get <dir>     a batch's or a run's settings and record, for a script
    simfbdr --version                   the commit of this code, -dirty if it has uncommitted changes

A batch is a directory holding a config.yml, anywhere on disk, and a run is a directory
under its runs/ holding one too.
"""
import argparse
import os
import subprocess
import sys

from . import analyze, batch, config, manifest, tools


def cmd_extend(a):
    d = os.path.abspath(a.batch)
    if "/runs/" in d + "/":
        raise SystemExit(f"error: {a.batch} is a run; extend its batch")
    batch.extend(d if os.path.isabs(a.batch) else a.batch, a.n, a.allow_sim_change)


def cmd_analyze(a):
    run = a.run.rstrip("/")
    if not analyze.analyze(run, force=a.force, jobs=a.jobs, slots=a.slots, slot_dir=a.slot_dir):
        raise SystemExit(1)


def cmd_status(a):
    raise SystemExit(1 if sum(manifest.status(p.rstrip("/")) for p in a.paths) else 0)


def run_script(name, args):
    path = tools.script(name)
    how = {".R": [tools.rscript()], ".py": [sys.executable], ".sh": ["bash"]}.get(path.suffix)
    if how is None:
        raise SystemExit(f"simfbdr: {name} is not something to run; simfbdr path {name} says where it is")
    raise SystemExit(subprocess.run([*how, str(path), *args], env=tools.environment()).returncode)


def cmd_config(a):
    for kind in ([a.arg] if a.arg in ("batch", "run") else ["batch", "run"]):
        sys.stdout.write(open(config.defaults_file(kind)).read() + "\n")


def cmd_doctor(a):
    import shutil
    print(f"simfbdr     {tools.version()}   {tools.home()}")
    ok = True
    for label, at, rule in tools.found():
        ok = ok and not at.startswith("NOT FOUND")
        print(f"{label:<10} {at}  ({rule})")
    _, _, describe = tools.rb_identity()
    print(f"           rb was built from {describe}")
    builds = tools.rb_builds()
    if len(builds) > 1:
        print("           a run can be analyzed with any of these, by stating it under uses in its config:")
        for path, commit in builds:
            print(f"             rb@{commit}  {path}")
    r = subprocess.run([tools.rscript(), "-e",
                        'for (p in c("paleobuddy","posterior","SBC")) cat(sprintf("%-10s %s\\n", p, '
                        'if (requireNamespace(p, quietly=TRUE)) as.character(packageVersion(p)) else "NOT INSTALLED"))'],
                       capture_output=True, text=True)
    print(r.stdout.rstrip() or r.stderr.strip()[-200:])
    print(f"tools file {tools.TOOLS_FILE}" + ("" if tools.TOOLS_FILE.exists() else "  (none; $RBIN, $RSCRIPT or PATH are used)"))
    raise SystemExit(0 if ok and "NOT INSTALLED" not in r.stdout else 1)


def before_the_rename(argv):
    """A notebook runs a step against a checkout of the commit it states, named by $SIMFBDR.
    A commit from before the rename has the package as simfbd, and that is what runs."""
    at = os.environ.get("SIMFBDR") or ""
    src = os.path.join(at, "src")
    if not os.path.isfile(os.path.join(src, "simfbd", "cli.py")) or os.path.isdir(os.path.join(src, "simfbdr")):
        return
    env = dict(os.environ, PYTHONPATH=os.pathsep.join(filter(None, [src, os.environ.get("PYTHONPATH")])))
    args = sys.argv[1:] if argv is None else list(argv)
    raise SystemExit(subprocess.run([sys.executable, "-m", "simfbd", *args], env=env).returncode)


def main(argv=None):
    before_the_rename(argv)
    ap = argparse.ArgumentParser(prog="simfbdr", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--version", action="store_true")
    ap.add_argument("--home", action="store_true", help="where this code is")
    sub = ap.add_subparsers(dest="cmd")

    p = sub.add_parser("extend", help="add replicates to a batch")
    p.add_argument("batch")
    p.add_argument("n", type=int)
    p.add_argument("--allow-sim-change", action="store_true",
                   help="add to a batch whose earlier replicates a different simulator made")
    p = sub.add_parser("analyze", help="analyze a run")
    p.add_argument("run")
    p.add_argument("--force", action="store_true", help="go on although a config has changed since the output was written")
    p.add_argument("--jobs", type=int, default=1, help="replicates analyzed in parallel")
    p.add_argument("--slots", type=int, help="cap chains across every run sharing --slot-dir")
    p.add_argument("--slot-dir")
    p = sub.add_parser("status", help="does a batch or a run still match its manifest")
    p.add_argument("paths", nargs="+")
    p = sub.add_parser("summarize", help="SBC ranks, coverage and diagnostics")
    p.add_argument("runs", nargs="+")
    p.add_argument("--jobs", type=int, help="logs read in parallel (default: every core)")
    p = sub.add_parser("ranks", help="mean rank per parameter")
    p.add_argument("run")
    p = sub.add_parser("config", help="the settings a config can hold")
    p.add_argument("what", choices=["default"])
    p.add_argument("arg", nargs="?")
    p = sub.add_parser("script", help="run one of the other scripts simfbdr carries")
    p.add_argument("name")
    p.add_argument("args", nargs=argparse.REMAINDER)
    p = sub.add_parser("path", help="where one of those scripts is")
    p.add_argument("name")
    sub.add_parser("doctor", help="the tools this machine would run")
    p = sub.add_parser("uses", help="what simfbdr would run, as name@version")
    p.add_argument("name", nargs="?")
    p = sub.add_parser("manifest", help="a batch's or a run's settings and record, for a script")
    p.add_argument("rest", nargs=argparse.REMAINDER)

    a = ap.parse_args(argv)
    if a.version:
        print(tools.version())
    elif a.home:
        print(tools.home())
    elif a.cmd == "extend":
        cmd_extend(a)
    elif a.cmd == "analyze":
        cmd_analyze(a)
    elif a.cmd == "status":
        cmd_status(a)
    elif a.cmd == "summarize":
        if a.jobs:
            os.environ["SUMMARIZE_JOBS"] = str(a.jobs)
        run_script("summarize.R", a.runs)
    elif a.cmd == "ranks":
        from . import ranks
        ranks.main([a.run])
    elif a.cmd == "config":
        cmd_config(a)
    elif a.cmd == "script":
        run_script(a.name, a.args)
    elif a.cmd == "path":
        print(tools.script(a.name))
    elif a.cmd == "doctor":
        cmd_doctor(a)
    elif a.cmd == "uses":
        print(tools.version_of(a.name) if a.name else "\n".join(tools.uses(tools.KNOWN)))
    elif a.cmd == "manifest":
        if len(a.rest) < 2 or a.rest[0] not in ("flat", "get"):
            raise SystemExit("usage: simfbdr manifest flat <dir> | simfbdr manifest get <dir> <name>")
        manifest.main(a.rest)
    else:
        ap.print_help()
