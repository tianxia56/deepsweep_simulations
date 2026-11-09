# DeepSweep simulations: CMS component statistics for neutral and sweep replicates

This repository records exactly how the simulated training data were made,
and lets anyone remake them on a Linux machine. One command builds the tools
from pinned sources, simulates neutral and selected-sweep replicates with
cosi2, computes nine CMS component statistics per site, normalizes them
against neutral reference bins, and writes zipped per-replicate tables.

```bash
git clone <this repo> && cd deepsweep_simulations
./run.sh config/test_100.json        # 100 sims of each kind: 19 min on 4 cores, measured
./run.sh                             # config/config.json: production scale
```

`run.sh` builds the environment the first time (24 min measured), verifies it,
then runs every stage to the end. Running the same command again resumes an
interrupted run. Requirements: Linux x86-64, conda, git, internet access for
the first build.

---

## 1. What is computed

For every polymorphic site of every replicate, in this column order:

| column | definition |
|---|---|
| `sim_batch_no`, `sim_id`, `pos` | batch ID, replicate number within the batch, position (bp) in the simulated region |
| `mean_fst` | Weir–Hill FST from haplotype allele counts (as in `broadinstitute/cms` `freqs_stats.c`), averaged over target-vs-outgroup pairs; a pair fixed for the same allele is excluded |
| `deldaf` | derived allele frequency (DAF) in the target minus the mean DAF of the outgroups; signed |
| `daf`, `maf` | DAF and minor allele frequency in the target population |
| `norm_nsl`, `norm_ihh12` | selscan v3.0.0 `--nsl` / `--ihh12 --maf 0.05`, normalized |
| `norm_ihs_binom`, `norm_delihh_binom` | hapbin `ihsbin --binom --minmaf 0.05`; iHS and ΔiHH = iHH₁ − iHH₀ (signed), normalized |
| `max_xpehh` | hapbin `xpehhbin --minmaf 0`, one run per outgroup; row-wise maximum of the per-pair normalized scores (joined on variant ID) |
| `deri_gen`, `sel_gen`, `s` | sweep replicates only: derived-allele age and selection onset (generations), selection coefficient, as sampled by cosi2 |
| `selpos` | sweep replicates only: 1 at the causal site, 0 elsewhere |

Normalization: z = (x − mean_k) / sd_k, with one bin per exact derived-allele
count k in the target (k = 0 … n_hap). The bins come from neutral replicates
of the same model, region length, target population and tool flags. A bin
with fewer than `min_sites_per_bin` neutral sites gives NaN. Derived alleles
are coded 1 throughout: cosi2 writes ancestral = 1, so every TPED is flipped
before use.

Outgroups are the model's other populations. The number of XP-EHH pairs
therefore varies (1 for cousins, 2 for gravel, 3 for the others), and so does
the neutral baseline of `max_xpehh`.

## 2. Pipeline

| stage | script | input → output |
|---|---|---|
| 1 plan | `run.py` | config → `plan/*.tsv`: every reference set, neutral batch and sweep batch with its seeds; seed ranges are checked to be disjoint and below 2³² |
| 2 reference bins | `pipeline/01_make_neutral_bins.py` | `n_sims` neutral cosi2 replicates per (model, region, target) → per-replicate arrays (`reference_pools/`) → bin tables + `MANIFEST.txt` (`reference_bins/<combo>/`) |
| 3 neutral corpus | `pipeline/run_neutral_batch.py` | neutral cosi2 replicates, scored with stage 4's own code against stage 2's bins → `neutral/<model>/*_cms_stats_all_*.zip` and seed tables |
| 4 sweep simulation | `pipeline/simulate_sweep_batch.py` | model `.par` + sweep line → TPED zip, `par_inputs` (sampled sweep parameters, causal DAF), drop list |
| 5 sweep scoring | `pipeline/02_recompute_batch.py`, `check_output.py`, `make_metafile.py` | TPED zip → `sweeps/stats/<model>/*_cms_stats_recomputed.zip`; independent re-check; one causal row per replicate (`sweeps/metafiles/`) |

Shared formulas, tool flags, binning and input guards are in
`pipeline/cms_stats.py`. Each batch stops with an error, and writes nothing,
if any of these fail:
- tool and TPED allele frequencies disagree;
- polarization looks inverted;
- a variant ID repeats;
- a statistic is never produced;
- the declared target population has no causal site at the midpoint, or its
  causal DAF is outside [0.05, 0.95];
- a replicate is missing a population;
- any replicate fails.

cosi2 command, identical for neutral and sweep replicates (a sweep replicate
is the model `.par` with one sweep line appended):

```
COSI_NEWSIM=1 COSI_MAXATTEMPTS=5000 coalescent -p <par> -v -g --genmapRandomRegions \
    --nTriesPerSim 5000 --drop-singletons .25 --tped <prefix> -n 1 -M -r <seed>
```

Seeds:

| stage | seed for replicate i | on retry |
|---|---|---|
| reference bins | `seed0 + i` | `+ k·10⁸` |
| neutral corpus | `seed0 + i·16 + attempt`, attempt 0–15 | next attempt |
| sweeps | `seed0 + i·16 + attempt`, attempt 1–15 | next attempt |

An attempt is retried after 120 s, and a replicate is given up after 600 s.
The seed actually used is recorded for every replicate.

## 3. Configuration

One JSON file sets everything. Lower the `n_sims` values to test; the rest
can stay as it is.

| key | meaning |
|---|---|
| `output_dir` | where outputs go (relative paths are inside the repo) |
| `workers` | parallel processes |
| `models[]` | `name` (a file in `models/`), `pop_ids` (populations scored), `target_pops`, `region_mb` (3 and/or 5) |
| `reference_bins.n_sims` | neutral replicates per (model, region, target) for the bins; production 7,000 |
| `reference_bins.seed_base`, `seed_step` | reference set c uses `seed_base + c·seed_step` |
| `neutral.n_sims_per_model_per_mb` | neutral corpus size per model and region length, split evenly over target pops; production 25,000 |
| `neutral.max_batch`, `seed_base`, `seed_block` | batch size cap, and the seed block given to each batch |
| `sweeps.designs[]` | `models`, `region_mb`, `target_pops` (`"all"` or a list), `age` (onset prior U(0, age)), `daf_window` (cosi2's final-frequency condition), `n_sims` per (model, region, target) |
| `sweeps.template` | the cosi2 sweep line; `{pop}`, `{age}`, `{daf_window}` are filled from each design |
| `sweeps.batch_table`, `replay_sims` | instead of designs: an explicit batch list, and exact (sim_id, seed) pairs to regenerate |
| `min_sites_per_bin` | 5,000 in production; lower values are accepted only with `"test": true` |
| `keep_sweep_tpeds` | keep the sweep TPED zips after scoring |

The configs in `config/`:
- `config.json`: the default. Production-scale counts, all six models, and
  one sweep design (U(0, 5000) onset, final frequency 0.05–0.95).
- `test_100.json`: one model (gravel), one target population, 3 Mb, 100
  replicates per stage, `min_sites_per_bin` 50. Measured on 4 cores: reference
  bins 315 s, neutral 322 s, sweeps 503 s; 148 MB of output.
- `production_2026oct.json`: the settings that produced the published
  corpus; see section 4.

## 4. Provenance of the 2026oct corpus

`config/production_2026oct.json` regenerates the published reference bins,
neutral corpus and seeded sweep batches. Its plan is identical to the
production plan for all 42 reference sets (seeds, populations, `.par`) and all
330 neutral batches (IDs, sizes, seeds, file names).

| component | count | regenerable here |
|---|---|---|
| reference bins | 42 sets × 7,000 replicates (vitti_gamma 6,996–6,999 pooled) | yes; 1,212 replicates were reseeded after a timeout, and their seeds are in `provenance/reference_reseeded_2026oct.tsv` |
| neutral corpus | 300,000 replicates, 330 batches | yes; batches are identical to `provenance/batches.tsv` |
| sweeps, seeded | 43,032 replicates kept from 162 batches | yes, exactly: `replay_sims` runs each kept replicate once at its recorded seed (`provenance/sweep_sims_kept_2026oct.tsv`) |
| sweeps, archived | 966,968 replicates from 1,018 batches | **no**: these were simulated before seeds were recorded; their TPEDs were rescored with stage 5 (`provenance/archived_batches_2026oct.tsv`) |

Total sweep replicates: 43,032 + 966,968 = 1,010,000.

Files in `provenance/`:

| file | rows | content |
|---|---|---|
| `combos.tsv` | 42 | production reference sets: model, `.par`, region (Mb), target pop, pops scored, seed0 |
| `batches.tsv` | 330 | production neutral-corpus plan |
| `reference_reseeded_2026oct.tsv` | 1,212 | reference replicates whose first seed timed out: seed formula vs seed used |
| `reference_pooled_2026oct.tsv` | 42 | replicates pooled per reference set, and the IDs not pooled |
| `sweep_batches_2026oct.tsv` | 162 | seeded sweep batches: model, batch, `.par`, target, seed0, sweep line, sims planned and kept |
| `sweep_sims_kept_2026oct.tsv` | 43,032 | kept seeded sweep replicates: seed, attempt, causal DAF, deri_gen, sel_gen, s |
| `archived_batches_2026oct.tsv` | 1,018 | unseeded archived batches and how many replicates each kept |
| `runplan_v2.csv` | 89 | the replacement-batch plan (DAF windows, ages, expected yields) |

Seeded batches were oversampled and then pruned to fill per-DAF-bin quotas,
so a batch kept fewer replicates than it simulated. The kept list is what
defines the corpus.

The production XP-EHH bin tables were not written by
`01_make_neutral_bins.py` itself. They were rebuilt on 2026-10-16, when
target-monomorphic sites were kept (k = 0 and k = n_hap), from pools
regenerated at the same seeds. The rebuilt pools match what `01` writes now:
`tests/reproduce_production.sh` finds them identical for the replicates it
checks.

One caveat for the reference and neutral stages: which attempt succeeds
depends on the 120 s timeout. On a slower machine a replicate that finished
in production could time out and be reseeded. The seed used is always
recorded, so any such difference shows up in the outputs.

## 5. Verification

- `env/check_env.sh` runs before every `run.sh` run. It checks that:
  - every tool resolves inside the env, links, and exits 0;
  - selscan matches v3.0.0 by md5;
  - the built cosi2 reproduces two production replicates byte for byte
    (`env/reference_md5.txt`).
- `tests/reproduce_production.sh`, about 2 min, regenerates production data
  with this repo and compares it:
  - a production sweep replicate: its TPEDs, sampled parameters and every
    statistic at every site;
  - three reference-pool replicates.

  The reference data is in `tests/reference/`.
- The empirical pipeline was checked against the 2026oct production files
  for YRI chromosome 22 (with its outgroups CEU, CHB and BEB). This needs the
  production intermediates (about 1.4 TB), so it is not shipped as a test.
  Identical to production:
  - the downloaded v5a VCF and panel (byte for byte), and the Zenodo pyrho
    hg19 maps (every file);
  - the converted maps (4 populations × 22 autosomes), the gap mask and the
    chr22 ancestral alleles;
  - the annotated TPED, `.hap` and `.map` of all four populations;
  - the raw statistics, the per-pair XP-EHH tables and the exact-k arrays;
  - the YRI genome-wide bin tables, rebuilt from the production raw arrays;
  - the normalized chr22 table, the genome-wide YRI table, and its
    `maf > 0.01` version (34,303,372 rows → 13,988,717).

  The YRI chr22 run took 14.5 min on 4 cores.

## 6. Tools

| tool | version / source | built by `env/setup_env.sh` |
|---|---|---|
| cosi2 (`coalescent`) | `broadinstitute/cosi2` commit `5da717c` (branch `is-251008-drop-imposs-traj`; reports "2.3.0rc4") | from source |
| hapbin (`ihsbin`, `xpehhbin`) | `evotools/hapbin` commit `3babf68` | from source |
| selscan | `szpiech/selscan` tag v3.0.0, `bin/linux/selscan` | prebuilt binary, md5-checked |
| Python | 3.9.23, numpy 2.0.2, pandas 2.3.1 | conda (`env/environment.yml`) |
| vcftools | 0.1.16, empirical pipeline only | conda (bioconda) |

Changing any tool flag invalidates every bin table built with the old flags.
Each `MANIFEST.txt` records the flags and versions used.

## 7. Empirical data: 1000 Genomes, GRCh37

The model is applied to 1000 Genomes data. The empirical pipeline computes
the same nine statistics with the same tools, flags and formulas
(`pipeline/cms_stats.py`), with two deliberate differences: the
normalization reference is genome-wide empirical, built per population at
its own n_hap, and the site list is the target plus its three outgroups.

```bash
bash empirical/download_inputs.sh inputs 22            # public inputs, chr22 only
./run.sh --empirical config/empirical_test_chr22.json  # YRI, chr22
bash empirical/download_inputs.sh                      # all autosomes
./run.sh --empirical config/empirical.json             # 22 populations x 22 autosomes
```

Inputs (every one is publicly available):

| input | source | build |
|---|---|---|
| genotypes | 1000 Genomes phase 3 release 20130502, `ALL.chr<N>.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz`, autosomes 1–22, and `integrated_call_samples_v3.20130502.ALL.panel` [1]. The v5a files are on the NCBI mirror (https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/); the EBI directory now lists v5b, and the differences were not assessed | GRCh37 |
| genetic maps | pyrho per-population maps for the 26 phase 3 populations (https://github.com/popgenmethods/pyrho#human-recombination-maps; Zenodo 10.5281/zenodo.11437540, `hg19_maps.tar.gz`) [2, 3]; each site's cM is linearly interpolated from its own population's map | GRCh37 (hg19) |
| ancestral alleles | `site.ancestral_state` of the unified tree sequences of 1000 Genomes, HGDP, SGDP and ancient genomes (Zenodo 10.5281/zenodo.5512994) [4], processed by E. Gilbertson's code (https://github.com/erin-n-gilbertson/3DGenome-diversity/blob/main/activeNotebooks/AncestralSequence_trees.ipynb), which also converted the coordinates between GRCh38 and GRCh37. Expected as `chr<N>_ErinG_AA_hg19.tsv.gz` with columns `CHR POS ErinG_ARG_AA hg38_position hg38_REF`. At a VCF site with no ARG ancestral allele, the global 1000 Genomes major allele (from the VCF `AF` field) is used instead (`make_unified_aa_map.py`) | trees GRCh38; tables used carry GRCh37 positions (with the GRCh38 position alongside) |
| excluded rsIDs | 11,963,907 historical rsIDs that have been merged into other rsIDs and are no longer used, from the NCBI dbSNP `RsMergeArch` table [5]. Records carrying them are removed during VCF conversion. The dbSNP build of the table was not recorded. | |

The simulations use a different map: cosi2 draws each replicate's region from
`models/test_recom.recom`. The pyrho maps are used only for the empirical
data.

Pipeline (`run_empirical.py`; each stage resumes where it stopped):

| stage | script | output |
|---|---|---|
| maps | `empirical/convert_pyrho_maps.py` | pyrho HapMap-format maps rewritten as `pos chr cM` (values copied as text) |
| gaps | `empirical/build_gap_mask.py` | intervals > 1 Mb with no map data in any of YRI, CEU, CHB and BEB (19 intervals, 97.4 Mb; `provenance/gaps_hg19.tsv`), dropped identically for every population |
| ancestral | `empirical/make_unified_aa_map.py` | one ancestral allele per VCF position |
| tped | `empirical/vcf_to_tped.py` | `vcftools --keep <pop> --min-alleles 2 --max-alleles 2`, no frequency filter; drops excluded rsIDs, duplicate IDs, duplicate positions (the first record is kept) and records with any missing or non-0/1 genotype; renames non-`rs` IDs to `rspos<pos>` |
| annotate | `empirical/annotate_tped.py` | derived = 1 against the ancestral allele (sites without one are dropped); cM linearly interpolated from the population's own map; sites before or after the map or inside a gap dropped; `.` IDs renamed `rspos<pos>` |
| stats | `empirical/emp_stats.py` | raw statistics per (target, chromosome); one-population statistics on the target's full site list, XP-EHH on each pair's shared sites (`harmonize.py`), site list applied last |
| bins | `empirical/build_bins.py` | exact-k tables per target from all 22 autosomes; bins with < 5,000 sites are NaN |
| normalize | `empirical/normalize.py` | z-scores; `max_xpehh` = maximum of the three per-pair z-scores, joined on variant ID |
| tables | `empirical/stage_tables.py`, `maf_filter.py` | `<POP>.empirical.2026oct.tsv.gz`, and a version keeping `maf > 0.01` (strict) |

Output columns: `chr pos ID mean_fst deldaf daf maf norm_nsl norm_ihh12
norm_ihs_binom norm_delihh_binom max_xpehh`, defined as in section 1. The
genome-wide reference is not neutral: it contains whatever selection the
genome carries.

Memory: production used 64 GB per stats task and 25 GB per annotate task
(chr1–2 are the largest), so on a PC set `workers` to what memory allows.

**Populations.** 22 of the 26 phase 3 populations are scored:

| superpopulation | populations |
|---|---|
| AFR | ACB, ASW, ESN, GWD, LWK, MSL, YRI |
| EUR | CEU, FIN, GBR, IBS, TSI |
| EAS | CDX, CHB, CHS, JPT, KHV |
| SAS | BEB, GIH, ITU, PJL, STU |

The four admixed American (AMR) populations, CLM, MXL, PEL and PUR, are
excluded. Their recent three-way admixture includes Indigenous American
ancestry that none of the simulated demographic models represents, so the
model was not trained to identify sweeps in such populations.

Each target is compared with three outgroups: the core populations (YRI,
CEU, CHB, BEB) outside its own superpopulation. A site is kept when the
target or any of its outgroups is polymorphic there.

### References

1. The 1000 Genomes Project Consortium. A global reference for human genetic variation. *Nature* 526, 68–74 (2015). doi:10.1038/nature15393
2. Spence, J.P. & Song, Y.S. Inference and analysis of population-specific fine-scale recombination maps across 26 diverse human populations. *Science Advances* 5, eaaw9206 (2019). doi:10.1126/sciadv.aaw9206
3. Kamm, J.A., Spence, J.P., Chan, J. & Song, Y.S. Two-locus likelihoods under variable population size and fine-scale recombination rate estimation. *Genetics* 203, 1381–1399 (2016). doi:10.1534/genetics.115.184820
4. Wohns, A.W. et al. A unified genealogy of modern and ancient genomes. *Science* 375, eabi8264 (2022). doi:10.1126/science.abi8264. Data: Zenodo, doi:10.5281/zenodo.5512994
5. NCBI dbSNP, `RsMergeArch` table (merged rsIDs), https://www.ncbi.nlm.nih.gov/snp/

## 8. Layout

```
run.sh, run.py        one-command driver (run_empirical.py for the empirical pipeline)
config/               JSON configs
models/               the 12 cosi2 demographic models (6 models x 3/5 Mb) + test_recom.recom (recombination map)
pipeline/             statistic and simulation code
empirical/            1000 Genomes pipeline (run by run_empirical.py) and input download
env/                  environment build, activation, check
provenance/           production plans, seeds and kept-replicate lists
tests/                reproduction test and its reference data
runs/<run_name>/      outputs (not tracked): plan/, reference_bins/, reference_pools/,
                      neutral/, sweeps/{tped,stats,metafiles}/, logs/, OUTPUTS.txt
```
