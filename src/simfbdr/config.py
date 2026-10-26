"""A batch's or a run's config.yml: what was asked for, in the names settings.py gives.

    timeline:
      rate_intervals: 8
      interval_width: 1.5
    origin_age:
      min: 15
      max: 21

A name that is not a setting is refused, with the nearest one that is, so a typo cannot
quietly leave a default in place. A config may also say which versions it is to be run with,

    uses: [rb@28e5b2, paleobuddy@1.1.0]

and then it is refused where they are not those. The defaults are written down once, in defaults/batch.yml
and defaults/run.yml beside this file, which list every setting there is. A batch does not say how many replicates it has: `simfbdr
extend` adds them and the manifest counts them.
"""
import difflib
import os

import yaml

from . import settings as S

HERE = os.path.dirname(os.path.abspath(__file__))
KINDS = {"batch": (S.BATCH_SETTINGS, S.BATCH_REQUIRED), "run": (S.RUN_SETTINGS, S.RUN_REQUIRED)}


def defaults_file(kind):
    return os.path.join(HERE, "defaults", f"{kind}.yml")


def defaults(kind):
    """{internal name: value} for the settings that have one when a config leaves them out."""
    by_name = {name: key for key, name, _ in KINDS[kind][0]}
    doc = yaml.safe_load(open(defaults_file(kind))) or {}
    return {by_name[name]: _word(v) for name, v in _flatten(doc) if v is not None}


def _flatten(doc, prefix=""):
    for k, v in (doc or {}).items():
        name = f"{prefix}{k}"
        if isinstance(v, dict):
            yield from _flatten(v, name + ".")
        else:
            yield name, v


def _word(v):
    if isinstance(v, list):
        return ",".join(_word(x) for x in v)
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, float) and v == int(v) and abs(v) < 1e15:
        return repr(v)                  # 1.0 stays 1.0, as it was written
    return str(v)


def load(path, kind, with_defaults=True):
    """{internal name: value as the scripts take it}, validated. Without the defaults it is
    only what the config itself says."""
    table, required = KINDS[kind]
    try:
        doc = yaml.safe_load(open(path)) or {}
    except yaml.YAMLError as e:
        raise SystemExit(f"simfbdr: {path} is not valid YAML\n  {e}")
    if not isinstance(doc, dict):
        raise SystemExit(f"simfbdr: {path} must be settings under their section names")
    by_name = {name: key for key, name, _ in table}
    out = defaults(kind) if with_defaults else {}
    stated = doc.pop("uses", None) or []
    out["_uses"] = [stated] if isinstance(stated, str) else [str(u) for u in stated]
    for name, v in _flatten(doc):
        if name not in by_name:
            other = {n: k for k, n, _ in KINDS["run" if kind == "batch" else "batch"][0]}
            if name in other:
                raise SystemExit(f"simfbdr: {path} sets {name}, which belongs to the "
                                 f"{'run' if kind == 'batch' else 'batch'}, not the {kind}")
            near = difflib.get_close_matches(name, by_name, n=1)
            raise SystemExit(f"simfbdr: {path} sets {name}, which is not a setting"
                             + (f". Did you mean {near[0]}?" if near else "")
                             + f"\n  simfbdr config default {kind} lists them")
        if v is None:
            continue
        key = by_name[name]
        out[key] = _word(v)
        if key in S.CHOICES and out[key] not in S.CHOICES[key]:
            raise SystemExit(f"simfbdr: {path}: {name} is {out[key]}; it is one of "
                             f"{', '.join(S.CHOICES[key])}")
    # the invisible mixture has exactly three components, so it needs no count
    if out.get("MIX_MODEL") == "invisible":
        if out.setdefault("MIX_K", "3") != "3":
            raise SystemExit(f"simfbdr: {path}: psi_mixture.model invisible has 3 components")
    name_of = {key: name for key, name, _ in table}
    if with_defaults is False:
        out.pop("_uses")
    missing = [name_of[k] for k in required if k not in out]
    if missing:
        raise SystemExit(f"simfbdr: {path} does not set {', '.join(missing)}")
    return out

