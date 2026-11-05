"""Every setting a batch or a run has: the name it is written under, and what it means.

A config is what you ask for, in these names, and a manifest is the record of what was done. The first column is the name the simulator and the Rev templates know a setting by,
which is how it reaches them and nothing a person writes or reads.
"""

# internal name, the name it is written under, what it means
BATCH_SETTINGS = [
    ("NINTERVALS", "timeline.rate_intervals",
     "piecewise-constant rate intervals"),
    ("INTERVAL_WIDTH", "timeline.interval_width",
     "width of each interval; the oldest runs on to the origin"),
    ("PRESENT", "timeline.present",
     "the process and its record end here, not at zero; younger occurrences are discarded"),
    ("AGE_MIN", "origin_age.min", "origin ~ Uniform(min, max)"),
    ("AGE_MAX", "origin_age.max", ""),
    ("GAMMA_L_SHAPE", "rate_prior.lambda.shape",
     "the speciation rate. PyRate's prior: an independent gamma draw per interval"),
    ("GAMMA_L_RATE", "rate_prior.lambda.rate",
     "left out, it is drawn per replicate from Uniform(0, rate_parameter_max) and estimated"),
    ("GAMMA_M_SHAPE", "rate_prior.mu.shape", "the extinction rate, unless turnover is set"),
    ("GAMMA_M_RATE", "rate_prior.mu.rate", ""),
    ("TURN_A", "rate_prior.turnover.beta.alpha",
     "set instead of mu: one turnover for every interval, mu = turnover * lambda, "
     "turnover ~ Beta(alpha, beta), always below 1"),
    ("TURN_B", "rate_prior.turnover.beta.beta", ""),
    ("TURN_LMEAN", "rate_prior.turnover.lognormal.log_mean",
     "or turnover ~ LogNormal(log_mean, log_sd), which can exceed 1, so a clade can decline"),
    ("TURN_LSD", "rate_prior.turnover.lognormal.log_sd", ""),
    ("GAMMA_P_SHAPE", "rate_prior.psi.shape", "the fossil sampling rate"),
    ("GAMMA_P_RATE", "rate_prior.psi.rate", "fixed; PyRate's is 1.1"),
    ("GAMMA_RATE_MAX", "rate_prior.rate_parameter_max",
     "upper bound of the uniform a lambda or mu rate parameter left out is drawn from"),

    ("RHO", "record.rho",
     "the extant sampling probability: the chance a lineage alive at the present is seen there"),
    ("BIN_WIDTH", "record.occurrence_bins.width",
     "occurrence ages are reported to bins of this width, edges at every multiple of it without end"),
    ("BIN_WIDTHS", "record.occurrence_bins.also_widths",
     "the same record on other grids too, as occurrences_<model>_w<width>; 0 is the exact ages"),
    ("JITTER_MIN", "record.occurrence_intervals.min_width",
     "set with max_width, the complete record is also written as occurrences_complete_jitter, each "
     "fossil reported over its own interval of width Uniform(min_width, max_width) holding its age "
     "at a uniform position, so intervals overlap"),
    ("JITTER_MAX", "record.occurrence_intervals.max_width", ""),
    ("BIN_OFFSETS", "record.occurrence_bins.also_offsets",
     "each of those shifted by this fraction of a bin as well, as ..._o<offset>"),

    ("MIN_TAXA", "accepted_replicates.min_taxa",
     "a replicate outside these bounds is redrawn"),
    ("MAX_TAXA", "accepted_replicates.max_taxa", ""),
    ("TAXA_COUNT", "accepted_replicates.taxa_counted",
     "fossil: fossil-sampled taxa only; record: survivors seen at the present too"),
    ("MIN_FOSSIL_TAXA", "accepted_replicates.min_fossil_taxa", ""),
    ("ANAGENESIS", "tree.anagenesis",
     "true: species also arise by anagenesis, at a rate lambda_a drawn per replicate from the "
     "speciation prior and constant in time"),
    ("SYM_A", "tree.symmetric_prior.alpha",
     "beta ~ Beta(alpha, beta), drawn per replicate: each speciation's chance to be symmetric, the "
     "parent ending and both daughters new species. Unset, every speciation is budding"),
    ("SYM_B", "tree.symmetric_prior.beta", ""),
    ("N_EXTANT", "tree.extant_species",
     "set, trees are drawn backward from exactly this many extant species and the origin is an "
     "output rather than a draw"),

    ("MIX_MODEL", "psi_mixture.model",
     "gamma: K rates of a discretized gamma of shape alpha, with equal weights; invisible: three "
     "components, one that never samples and a low and a high rate, with free weights"),
    ("MIX_K", "psi_mixture.components",
     "each species draws one of this many relative sampling rates at its birth; unset or 1, "
     "every species samples at psi. invisible has 3"),
    ("MIX_ALPHA_RATE", "psi_mixture.alpha_prior_rate",
     "gamma: alpha ~ Exponential(this rate), drawn per replicate"),
    ("MIX_W_CONC", "psi_mixture.weights_prior",
     "invisible: the weights ~ Dirichlet(this, this, this), and the ratio of the low rate to "
     "the high ~ Uniform(0, 1), drawn per replicate"),
    ("N_CHARS", "characters.count",
     "morphological characters that change only at speciation; unset, none are simulated"),
    ("CHAR_STATES", "characters.states",
     "states of each character, one number for all or a list with one per character"),
    ("EPS_A", "characters.epsilon_prior.alpha",
     "epsilon ~ Beta(alpha, beta), drawn per replicate: each character's chance to change at a speciation"),
    ("EPS_B", "characters.epsilon_prior.beta", ""),
    ("EPSA_A", "characters.epsilon_a_prior.alpha",
     "epsilon_a ~ Beta(alpha, beta), drawn per replicate: each character's chance to change at an "
     "anagenetic speciation, epsilon then being budding's alone. Unset, anagenesis changes with epsilon"),
    ("EPSA_B", "characters.epsilon_a_prior.beta", ""),
    ("CHAR_VARIABLE", "characters.variable_only",
     "keep only characters that vary among the reported taxa, as a published matrix does: each "
     "invariant one is drawn again along the same history until it varies, so the count is kept"),
    ("RET_A", "characters.retention_prior.alpha",
     "keep only characters that are parsimony informative among the reported taxa, and each varying "
     "only by autapomorphies with probability retention ~ Beta(alpha, beta), drawn per replicate. A "
     "character dropped is drawn again along the same history, as for variable_only"),
    ("RET_B", "characters.retention_prior.beta", ""),
    ("CHAR_MODEL", "characters.model",
     "morphospeciation: each character changes only at speciation, with probability epsilon; mk: instead "
     "along each lineage in time, under Mk at a rate drawn per replicate"),
    ("MKR_SHAPE", "characters.rate_prior.shape",
     "mk: the rate ~ Gamma(shape, rate), in expected changes per character per unit time"),
    ("MKR_RATE", "characters.rate_prior.rate", ""),

    ("MAX_LINEAGES", "simulation_limits.max_expected_lineages",
     "a draw expected to bear more lineages than this is refused unsimulated, which truncates the prior"),
    ("SIM_MAX_DRAWS", "simulation_limits.max_waiting_time_draws",
     "a draw whose simulation needs more waiting times than this is redrawn"),
]

RUN_SETTINGS = [
    ("MODEL", "inference.model", "fbdrp: dnFBDRP; bds: dnBDS; fbdsp: dnFBDSP, the tree process"),
    ("REPORTING", "inference.occurrences", "the batch's occurrences_<this> directory"),
    ("INFER", "inference.reporting_model",
     "complete: every occurrence is reported; firstlast: each taxon's oldest and youngest; "
     "presence: each reported interval holds at least its count"),
    ("PSI_MIX", "inference.psi_mixture",
     "assume the batch's psi mixture, estimating its hyperparameters and each species' "
     "component (fbdrp or fbdsp); false analyzes the same data under a single psi"),
    ("ANA", "inference.anagenesis",
     "estimate the anagenetic rate lambda_a; false analyzes the data as pure budding. fbdsp only"),
    ("CHARS", "inference.characters",
     "analyze the batch's characters on the tree, estimating epsilon; fbdsp only"),
    ("EPSILON_A", "inference.epsilon_a",
     "estimate a separate chance to change at anagenesis, with the batch's epsilon_a prior or else "
     "epsilon's; needs characters and anagenesis"),
    ("CODING", "inference.character_coding",
     "all: the characters are every one simulated; variable: each is conditioned on varying, "
     "which a batch with characters.variable_only needs; informative: each is conditioned on "
     "being parsimony informative"),
    ("RETENTION", "inference.retention",
     "refused: it was dnPhyloMorphospeciation's own correction for keeping a character varying only "
     "by autapomorphies, and neither character model here has it"),
    ("CHAR_INFER", "inference.character_model",
     "morphospeciation: each character changes at a speciation with probability epsilon, the "
     "speciations on a branch being the latents fnTreeAssembly carries as its length; mk: the same "
     "characters along each lineage in time at a strict clock rate, whichever model the batch simulated"),
    ("MKI_SHAPE", "inference.character_rate_prior.shape",
     "mk: the rate's prior, Gamma(shape, rate). Left out, the batch's characters.rate_prior"),
    ("MKI_RATE", "inference.character_rate_prior.rate", ""),
    ("AGES", "inference.occurrence_ages",
     "reported: the bins as they are; midpoint: each at its bin's middle; drawn: each uniform in its bin"),
    ("AGE_DRAWS", "inference.age_draws",
     "when drawn: this many draws, a chain each, pooled into one posterior per replicate"),
    ("COND", "inference.condition", "what the process is conditioned on, as RevBayes names it"),
    ("RHO", "inference.rho", "the extant sampling probability the analysis assumes"),
    ("ORIGIN_PRIOR", "inference.origin_prior", ""),
    ("LATENTS", "inference.monitor_latents",
     "false drops the per-taxon b, d, tau_1 and tau_K from the log, which a rates-only "
     "comparison never reads and which costs 63 MB a replicate on a thousand ranges"),
    ("RANGE_START", "inference.range_start",
     "which start of a range the log carries: origination, the species' own, which both "
     "processes report and which a cross-model comparison needs; or divergence, where it "
     "attaches to the nearest sampled range, which only dnFBDRP has"),
    ("REQUIRE_PARENT", "inference.require_parent",
     "reject a birth with no lineage alive to have produced it. Left out it follows "
     "coexistence; turn it off where ranges were left out of the analysis"),
    ("COEXISTENCE", "inference.coexistence",
     "weight each birth by the lineages alive to produce it. Left out, each process keeps "
     "its own published default: dnFBDRP weights them, dnBDS does not"),
    ("DROP_PRESENT_ONLY", "inference.drop_present_only_taxa", "bds only"),
    ("GENS", "mcmc.generations", ""),
    ("PRINTGEN", "mcmc.sample_every", ""),
    ("SEED", "mcmc.seed", "each replicate's chain is seeded with this plus its number"),
]

# what a manifest records beyond settings: the name a script asks for, and where it sits
BATCH_RECORD = {
    "kind": "kind", "batch": "batch", "config_hash": "config_hash",
    "true_vals_order": "true_vals_order", "NREPS": "replicates",
    "filtered": "prior_filter.refused", "filtered_of": "prior_filter.draws",
    "filtered_frac": "prior_filter.fraction",
}
RUN_RECORD = {
    "kind": "kind", "run": "run", "batch": "batch",
    "batch_config_hash": "batch_config_hash", "config_hash": "config_hash",
    "true_vals_order": "true_vals_order", "NREPS": "replicates", "occurrences": "occurrences",
    "rb_md5": "rb_md5",
}


# what a batch cannot be simulated without, and a run cannot be analyzed without
BATCH_REQUIRED = ["NINTERVALS", "INTERVAL_WIDTH", "GAMMA_L_SHAPE", "GAMMA_P_SHAPE",
                  "GAMMA_P_RATE", "AGE_MIN", "AGE_MAX", "BIN_WIDTH"]
RUN_REQUIRED = ["MODEL", "REPORTING", "INFER", "COND", "RHO", "GENS", "PRINTGEN"]
# The value a setting takes when a config leaves it out is in defaults/batch.yml and
# defaults/run.yml, not here.

# a setting that takes one of a few words
CHOICES = {"RANGE_START": ("origination", "divergence"), "MODEL": ("fbdrp", "bds", "fbdsp"),
           "INFER": ("complete", "firstlast", "presence"),
           "PSI_MIX": ("true", "false"), "CHARS": ("true", "false"),
           "CHAR_VARIABLE": ("true", "false"), "CODING": ("all", "variable", "informative"),
           "CHAR_MODEL": ("morphospeciation", "mk"), "CHAR_INFER": ("morphospeciation", "mk"), "EPSILON_A": ("true", "false"),
           "RETENTION": ("true", "false"),
           "MIX_MODEL": ("gamma", "invisible"),
           "ANAGENESIS": ("true", "false"), "ANA": ("true", "false"),
           "COND": ("time", "sampling", "survival", "nTaxa"),
           "TAXA_COUNT": ("fossil", "record"),
           "ORIGIN_PRIOR": ("uniform", "exponential", "none"),
           "AGES": ("reported", "midpoint", "drawn")}
