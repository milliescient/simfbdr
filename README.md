# simfbdr

Simulates batches of fossil records and analyzes them in RevBayes, for simulation-based
calibration of the fossilized birth-death range models.

## Install

    pip install -e /path/to/simfbdr        # or: uv tool install /path/to/simfbdr

That puts `simfbdr` on the path. It runs R and RevBayes, which are not Python's to install.
Each is found from `$RBIN` or `$RSCRIPT`, else `~/.config/simfbdr/tools.yml`, else `PATH`:

```yaml
rb:                               # one build, or a list of them
  - /opt/revbayes/build/rb
  - /opt/revbayes/build-older/rb
rscript: /opt/R/bin/Rscript
```

`simfbdr doctor` says which this machine would use, which rule supplied each, and whether the R
packages are there. Every `extend` and `analyze` prints the same before it starts. `simfbdr uses`
prints what would run as `name@version`, and `simfbdr uses rb` one version.

A config can say which versions it is to be run with, and is then refused where they are not
those, however the tools were found:

```yaml
uses: [rb@28e5b2]                 # in a run's config: any of simfbdr, rb
uses: [paleobuddy@1.1.0, R@4.5.3] # in a batch's: any of simfbdr, R, paleobuddy
```

### Which rb analyzes a run

A machine can hold several builds of rb, and each is known by the commit its binary carries,
which `simfbdr doctor` lists. `simfbdr analyze` picks one for the run, by the first of these that
applies:

1. the version the run's config states under `uses`. If no build here is that one, the run is
   refused and the builds that are here are listed.
2. the build that analyzed the run's earlier replicates, read from its manifest, so a run that
   is topped up after its batch is extended keeps one rb throughout. If that build is gone the
   run is refused. `--force` goes on with the default and the manifest then names the new one.
3. the first build the tools file lists, else `rb` on `PATH`.

`$RBIN` overrides all three for one command, and is still held to what the config states.

The run's manifest records which build it was, as `uses: ["simfbdr@...", "rb@..."]`.

## Layout

A **batch** owns the generating process and every hyperprior. One simulation pass emits
all reporting models from one record, so the runs under it share a timeline and one
`true_vals.tsv`, and comparisons between them are paired on identical truth.

A **run** under a batch's `runs/` chooses only what to read and how to analyze it. It inherits
the batch's hyperpriors and has no way to set them.

```
a-batch/
  config.yml         the generating process and every hyperprior, written by you
  manifest.yml       the record of what was simulated from it, written by simfbdr
  occurrences_complete/  occurrences_firstlast/  occurrences_presence/
  times/  true_ranges/  true_vals.tsv  sim_list.RData  seeds.RData
  true_components/  characters/   with a psi mixture, with characters
  runs/fbdrp_firstlast/
    config.yml       which occurrences to read, and how to analyze them
    manifest.yml     + what it read and wrote, and which code and which rb ran
    output/  aux/  results/
```

A batch is a directory anywhere on disk. A manifest does not record where it is, and a run
finds its batch as the directory two above it.

```
src/simfbdr/          the package
  cli.py             the simfbdr command
  settings.py        every setting: the name it is written under and what it means
  defaults/          batch.yml and run.yml: every setting, with the value it takes when left out
  config.py          reads and checks a config.yml
  batch.py  analyze.py  simulate replicates, analyze a run
  manifest.py        writes, reads and checks a manifest.yml
  ranks.py           mean rank per parameter
  scripts/           what simfbdr runs or carries: the simulator sim.R, the Rev templates,
                     and summarize.R
```

## Running

```sh
simfbdr extend sims/batch1 200                # simulate the batch's first 200 replicates
simfbdr analyze sims/batch1/runs/fbdrp_firstlast --jobs 6  # analyze every replicate not yet analyzed
simfbdr summarize sims/batch1/runs/fbdrp_firstlast
simfbdr extend sims/batch1 800                # later: 800 more, and the batch has 1000
simfbdr status sims/batch1 sims/batch1/runs/fbdrp_firstlast
```

## Config and manifest

You write the config and simfbdr writes the manifest. A setting is written once, in the
config. The manifest holds the record: the hash of the config it was made from, how many
replicates there are, which code ran, the hashes of what it read and wrote, and what the
filter refused.

```yaml
# a batch's config.yml
timeline:
  rate_intervals: 3
  interval_width: 3
origin_age:
  min: 9
  max: 15
rate_prior:
  lambda: {shape: 1.1}
  mu: {shape: 1.1}
  psi: {shape: 1.5, rate: 1.1}
record:
  occurrence_bins: {width: 3, oldest_edge: 30}
simulation_limits:
  max_expected_lineages: 250
```

```yaml
# a run's config.yml
inference:
  model: fbdrp
  occurrences: complete
  reporting_model: complete
  condition: time
  rho: 1.0
mcmc:
  generations: 200000
  sample_every: 1000
```

Every setting, with what it means and its default, is in `defaults/batch.yml` and
`defaults/run.yml` in the package. `simfbdr config default batch` prints one, which can be
copied to start a config. A name that is not a setting is refused, with the nearest one that
is, and so is a run setting in a batch's config or the reverse. Times are in whatever unit
the rates are per.

A batch's config does not say how many replicates it has. `simfbdr extend` adds them and the
manifest counts them. `extend` refuses if the config or the simulator has changed since the
earlier replicates, which are never redrawn or rewritten. An extension simulates only what is
missing, appends rows to `true_vals.tsv` and `nums.tsv`, and adds its per-replicate files
beside the others, so the hashes of those files change.

### Checking a batch or a run

A manifest records what a notebook's lock records, in the same words. `uses` names the code
that ran as `name@version`, and `inputs` and `outputs` give each path with its hash:

```yaml
uses: ["simfbdr@f0c6fe5", "rb@28e5b2"]    # a batch's names R and paleobuddy too
inputs:
  config.yml: d26610a2e386               # this run's config
  batch:
    config.yml: 60319e066c3d             # the batch's
    occurrences_complete: 5d415912e147
  simfbdr:
    infer-fbdrp.Rev: 12eb464f40e7        # the template simfbdr carried
outputs:
  output: 0b1b35e37934                   # written by this run
```

Everything a batch or a run read is under `inputs`, and everything is hashed one way. The hash
is `digest()` in `manifest.py`: SHA-1 over one line per file, its relative name and the SHA-1
of its bytes, files in sorted order, first 12 hex digits. It is the hash sloop computes, so a
page that pins one of these paths, a config included, records the same value.

`simfbdr status` recomputes the hashes and reports each path as `current`, `changed` or
`missing`. For a run it also checks the edge: what the run read must be what the batch says
it wrote, and a run analyzed before its batch was extended is reported `stale`. It exits
non-zero if anything is changed, missing or stale. A script of simfbdr's own that has changed
since is reported as `differs` and not counted.

A batch's manifest keeps the output hashes of every extension, and only the last can still
be checked, since an extension changes the files.

`simfbdr --version` prints the commit of the code, ending in `-dirty` when it has uncommitted
changes.

## The rate prior

The rates are drawn from PyRate's prior, an independent gamma for each rate in each interval,
and analyzed under the same one.

```yaml
rate_prior:
  lambda: {shape: 1.1}               # speciation, PyRate's -pL
  mu: {shape: 1.1}                   # extinction, PyRate's -pM
  psi: {shape: 1.5, rate: 1.1}       # fossil sampling, PyRate's -pP
```

A `lambda` or `mu` `rate` left out is PyRate's hyperprior, which it uses whenever the
shift times are fixed: one rate parameter shared by the intervals, Uniform(0, 20), whose upper
bound is `rate_parameter_max`. The simulator draws it per replicate and the analysis estimates
it, as `gamma_rate_lambda` and `gamma_rate_mu`, which are in `true_vals.tsv` and ranked with
the rates. A stated `rate` is fixed.
The `psi` `rate` is always stated. PyRate's hyperprior on it, which `-pP 1.5 0`
asks for, is not built.

A rate parameter near zero is a speciation rate with no ceiling, so a batch that leaves one
out wants `simulation_limits.max_expected_lineages`, and its `prior_filter` says how much of the
uniform that refused.

Lognormal rates, a random walk or a stationary process through time, the dnFBDSP model with
anagenetic speciation, and batches made before configs were YAML are on the `all-priors`
branch.

## The two ceilings on a simulation

`simulation_limits.max_expected_lineages` is a ceiling on the lineages a draw is expected to
bear. A deep origin needs it, because a draw in the upper tail of the speciation rate
otherwise never finishes simulating.

`simulation_limits.max_waiting_time_draws` is the second ceiling, on the work one simulation
may do. It counts the waiting times the birth-death simulator draws, across every tree it
tries, and refuses the draw when the count is passed. The count follows from the random
stream, so the same seeds refuse the same draws on any machine and a batch reproduces from
`seeds.RData`. A tree of n lineages takes about 3n draws, so keep it several times
`max_expected_lineages`: 100000 cuts a runaway off after roughly 20 seconds on the machine it was
measured on, 200000 after roughly 60. The simulator prints how many draws were refused and
the most any finished simulation used.

The first ceiling refuses a draw before it is simulated, so it truncates the prior the
simulator draws from while the analysis keeps the untruncated one, and the SBC ranks are
non-uniform by that much. Every batch records `prior_filter`, the fraction of draws the
ceiling refused, in its manifest. A truncation of fraction e shifts a mean rank by roughly
e/2, against a standard error of 0.289/sqrt(n). At n = 10000 that is 0.0029, so e has to stay
near or below 0.5% to be invisible.

A floor and ceiling on the number of reported taxa, `accepted_replicates`, is a function of
the record and does not bias the ranks.

## Occurrences and reporting model

- `inference.occurrences` = `complete` | `firstlast` | `presence`: which occurrences the
  analysis reads. All come from the same simulated history. `presence` reports each bin a taxon
  was sampled in once, so a bin left out truly held nothing.
- `inference.reporting_model` = `complete` | `firstlast` | `presence`: how the model treats them.
  Under `presence` a reported count is a lower bound, and `mvSingletonRJSwitch` samples whether
  a taxon reported once has one appearance age or two.

A run that matches them validates the model. A run that does not is a misspecification
analysis.

## Sampling heterogeneity

A batch with `psi_mixture.components: K` gives each species a relative sampling rate drawn at its
birth, independently of its parent: one of the K rates of a discretized gamma of shape alpha,
with equal weights, so psi stays the mean rate. alpha ~ Exponential(`alpha_prior_rate`) is drawn
per replicate and is a truth like the rates; each species' component is in
`true_components/`. A run with `inference.psi_mixture: true` assumes the mixture and estimates
alpha (dnFBDRP only); `false` analyzes the same data under a single psi.

## Morphospeciation

A batch with `characters.count` simulates characters that change only at speciation, down the
complete tree: the first species' states are uniform, and at each budding every character of
the new species changes with probability epsilon to a uniformly drawn other state. epsilon ~
Beta(`epsilon_prior`) is drawn per replicate. The characters of every reported taxon are in
`characters/characters_<rep>.nex`. A run with `inference.model: fbdsp` and
`inference.characters: true` adds them under `dnPhyloMorphospeciation` on the dnFBDSP tree and
estimates epsilon.

A batch with neither draws no extra random number, so it is the same batch it was before
either existed.

## Dating resolution

A batch can report its record on other occurrence grids beside its own, so that analyses differing
only in how finely the fossils are dated are paired on the same fossils:

```yaml
record:
  occurrence_bins:
    width: 1.5
    oldest_edge: 36
    also_widths: [0, 0.75, 3, 6]     # 0 is the exact ages
    also_offsets: [0, 0.25]          # fractions of a bin
```

Each width at each offset is written as `occurrences_<model>_w<width>`, with `_o<offset>` when
the offset is not 0, and a run reads one with `occurrences: complete_w3_o0.25`. These grids
stop at the present and at `origin_age.max`, so an offset grid has no bin past the origin
prior. The batch's own grid is written as it always was.

How an analysis treats the ages it reads is the run's choice, `inference.occurrence_ages`:
`reported` hands rb the bins, `midpoint` puts each occurrence at its bin's middle, and `drawn`
at a uniform draw within it. `inference.age_draws: 10` makes ten such draws of each replicate
and runs a chain on each, which is PyRate's practice, and pools them. `output/rep_<n>.log` is
then the length of one chain, and after burn-in takes each sample from the next draw in turn,
so each draw counts the same. The chains themselves are under `aux/draw_<k>/`. A midpoint
run is refused where a midpoint is a rate-shift time, since rb cannot start there.

## What is enforced

The batch's `config.yml` is the only definition of the hyperpriors. simfbdr hands it to the
simulator and writes the same values into the per-replicate `aux/run_<rep>.Rev` the Rev
template reads.

- A missing setting is fatal. The simulator has no defaults for what a batch must say.
- A run cannot set what the batch owns. A batch setting in a run's config is refused.
- The simulator reports what it drew with. `sim_params.tsv` is compared against the config.
- An edited config is refused. Each manifest records the hash of the config that produced its
  data, and a run records the batch hash it ran against. `simfbdr analyze --force` overrides
  for a run.

`settings.py` has one line a setting: the name it is written under, what it means, and the
name the simulator and the Rev templates know it by.

## What is not reproducible

The same seeds give the same batch, byte for byte. An analysis reproduces on one machine with
one binary, since each replicate's chain is seeded with `mcmc.seed` plus the replicate's
number. Across machines or builds of rb a chain may differ, and its summaries agree to within
Monte Carlo error.

## Scope

`inference.model` is `fbdrp`, `bds` or `fbdsp`, the range process on a sampled budding tree.
PyRate itself is not run. `bds` under
`condition: sampling` is its likelihood, and the script that ran the program on a batch is on
the `all-priors` branch.
