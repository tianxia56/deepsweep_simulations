#!/usr/bin/env python3
"""
Shared helpers for the CMS component-statistic pipelines.

Conventions (see CLAUDE.md for the reasoning behind each)
---------------------------------------------------------
TPED coding : derived allele is '1'. cosi2 emits ancestral=1 / derived=0, so
              every TPED is flipped before anything else touches it.
FST         : allele-count variance-component estimator from
              broadinstitute/cms  cms/cms/model/freqs_stats.c ("Weir-Hill").
              No within-individual term. A pair fixed for the same allele in
              both pops is NaN and left out of mean_fst (CMS convention);
              mean_fst is NaN where no pair is informative.
deldaf      : D_target - mean(D_outgroups)                            [CMS]
delihh      : iHH_1 - iHH_0, SIGNED (deliberate deviation from the CMS
              |iHH_1 - iHH_0|).
iHS/delihh  : production is ihsbin --binom only. The 'default' flavour exists
              for the finished archive comparison and is never written in
              production (ihs_flavours(include_default=False)).
XP-EHH      : binned on TARGET DAF; target-monomorphic sites are KEPT
              (XPEHH_KEEP_TARGET_MONOMORPHIC = True, bins k = 0..n_hap).
norm        : z = (x - mean_k) / std_k, one bin per EXACT derived allele
              count k (count_from_daf); bins specific to n_hap. Bins with
              fewer than MIN_SITES_PER_BIN neutral sites are NaN.
"""

import os
import subprocess

import numpy as np
import pandas as pd

# Normalization bins are keyed on the EXACT derived allele count k, not on
# equal-width frequency intervals. Equal-width bins left large residual
# frequency dependence INSIDE each bin - see CLAUDE.md, Normalization.
MIN_SITES_PER_BIN = 5000     # thin bins -> NaN rather than an unreliable sigma
# Smoke tests only (workflow/plan.py refuses it unless the config has
# "test": true): a few dozen neutral sims cannot fill 5,000 sites per bin.
# Every bin MANIFEST records the value used ("min_per_bin").
if os.environ.get("DEEPSWEEP_MIN_SITES_PER_BIN"):
    MIN_SITES_PER_BIN = int(os.environ["DEEPSWEEP_MIN_SITES_PER_BIN"])

# Frequency filters, always passed explicitly - never rely on a tool default.
# Each value reproduces the tool's effective default as of 2026-10 (verified
# byte-for-byte against the gravel_3mb_pop1 pool):
#   selscan  --maf 0.05     MAF filter on core SNPs; v3.0.0 keeps MAF == 0.05
#   ihsbin   --minmaf 0.05  gates on maxEHH (= allele freq at the focus);
#                           main-ihs.cpp:45
#   xpehhbin --minmaf 0     main-xpehh.cpp:47 sets 0.00 although --help says
#                           0.05. A nonzero value requires BOTH pops within
#                           [m, 1-m] (ehhfinder-impl.hpp:250-254), which would
#                           drop sweeps near fixation in the target.
SELSCAN_FLAGS = ["--maf", "0.05"]
IHSBIN_FLAGS = ["--minmaf", "0.05"]
XPEHHBIN_FLAGS = ["--minmaf", "0"]

# XP-EHH keeps sites monomorphic in the TARGET (polymorphic only in the
# outgroup), binned at k = 0 / k = n_hap on target DAF. Decision 2026-10-16:
# matches the empirical site-list rule, which keeps every site polymorphic in at
# least one core population, target-monomorphic included. Bin tables must
# have been built the same way: their MANIFEST carries this name.
XPEHH_KEEP_TARGET_MONOMORPHIC = True

# ==============================================================================
# INPUT GUARDS
# ==============================================================================
# Five assumptions the pipeline relied on without checking (code review
# 2026-10-28). Each is now a hard error, and each has a deliberate failure
# case in diag/review_guard_tests.py showing it fires.
class InputError(ValueError):
    """A hard input/consistency failure. The caller must abort and write
    nothing: every one of these means a statistic would be silently wrong
    rather than absent."""


# Share of the symmetric DAF window [0.05, 0.95) lying below 0.5. cosi2's
# spectrum has a large rare-derived excess, so this is well above 0.5 when
# derived = 1 and maps to 1 - f_low when the TPED is inverted.
# Calibration (2026-10-28). Pooled per batch, all 1,180 sweep + 330 neutral
# batches: minimum 0.6427. Per sim, 5,574 sims across six models: minimum
# 0.6314. Per sim inside the two WORST-polarized batches in the corpus
# (vitti_gamma 43015040246 and 330180030242), which the six-model sample did
# not cover: minimum 0.5864 -- that is the real floor. An inverted TPED maps
# f_low -> 1 - f_low, so it reaches at most 0.4136. Nothing real lies in
# [0.414, 0.586] and the threshold sits at its midpoint, 0.086 from either
# side. Re-measure this if the corpus ever gains a demography with a flatter
# spectrum than vitti_gamma's.
POLARIZATION_MIN = 0.50


def polarization_f_low(daf):
    d = np.asarray(daf, dtype=float)
    d = d[(d >= 0.05) & (d < 0.95)]
    return float((d < 0.5).mean()) if len(d) else float("nan")


def check_polarization(daf, label, min_f_low=POLARIZATION_MIN):
    """Is this TPED polarized derived = 1?

    02 never calls flip_tped; it assumes the archive is already flipped. The
    only other guard is check_target_label comparing the causal DAF against
    the metafile, which an inverted TPED passes whenever the causal DAF is
    exactly 0.5 -- 6,795 of the corpus's 1,010,000 sims (0.67%).
    """
    f = polarization_f_low(daf)
    if not np.isfinite(f):
        return f
    if f < min_f_low:
        raise InputError(
            f"{label}: site-frequency spectrum looks INVERTED "
            f"(f_low={f:.4f} < {min_f_low}). Share of DAF in [0.05, 0.95) "
            f"below 0.5; every real sim measured is >= 0.586. The TPED is "
            f"probably ancestral = 1, i.e. flip_tped was not applied.")
    return f


def check_unique_ids(ids, label):
    """Every statistic is joined on variant ID, so a duplicate multiplies
    rows silently. This is what inflated the empirical chr12 by 12,152 rows
    (241 sites sharing the rsID '.')."""
    a = np.asarray(ids, dtype=str)
    u, c = np.unique(a, return_counts=True)
    dups = u[c > 1]
    if len(dups):
        raise InputError(f"{label}: {len(dups)} duplicate variant ID(s) "
                         f"(e.g. {list(dups[:5])}) among {len(a)} rows")
    return len(a)


def safe_merge(left, right, on, how, label):
    """A left join that refuses to change the row count."""
    out = left.merge(right, on=on, how=how)
    if len(out) != len(left):
        raise InputError(f"{label}: join on '{on}' changed the row count "
                         f"{len(left)} -> {len(out)}; the right side has "
                         f"duplicate keys")
    return out


def check_tool_freq(tool_ids, tool_daf, tped_ids, tped_daf, n_hap, label):
    """Does the tool's own frequency column agree with the TPED's DAF?

    Closes two gaps at once.

    (1) The bin tables are keyed on k = rint(tool_freq * n_hap), because
        01_make_neutral_bins stores the tool's printed frequency in the pool,
        while 02_recompute_batch normalizes on k = rint(tped_daf * n_hap).
        Two routes to the same bin; nothing checked they agree.

    (2) hapbin reports the frequency of whatever .hap it was handed. If a
        population's .hap were paired with the target's .map, this is the
        only place the mix-up would show -- sfs_table's position assertion
        validates `parsed`, not the hap files.
    """
    k_tool = count_from_daf(np.asarray(tool_daf, dtype=float), n_hap)
    lut = dict(zip(np.asarray(tped_ids, dtype=str),
                   count_from_daf(np.asarray(tped_daf, dtype=float), n_hap)))
    tids = np.asarray(tool_ids, dtype=str)
    k_tped = np.array([lut.get(i, -1) for i in tids])
    unknown = int((k_tped < 0).sum())
    if unknown:
        raise InputError(f"{label}: {unknown} of {len(tids)} output rows carry "
                         f"a variant ID absent from the target TPED")
    bad = k_tool != k_tped
    if bad.any():
        i = int(np.flatnonzero(bad)[0])
        raise InputError(
            f"{label}: derived count from the tool's frequency disagrees with "
            f"the TPED at {int(bad.sum())} of {len(tids)} sites "
            f"(first: ID {tids[i]}, tool k={k_tool[i]}, tped k={k_tped[i]}, "
            f"max |diff| {int(np.abs(k_tool - k_tped).max())}). The bin tables "
            f"are keyed on the tool's frequency and the batch on the TPED's, "
            f"so these sites would be scored against the wrong reference; a "
            f"wholesale disagreement means the .hap and .map come from "
            f"different populations.")
    return len(tids)


def tool_versions():
    """Version strings for the MANIFEST, each with where it came from."""
    def first_line(cmd, pat):
        try:
            out = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
            for line in (out.stdout + out.stderr).splitlines():
                if pat in line:
                    return line.strip()
        except OSError as e:
            return f"unavailable ({e})"
        return "not reported"
    import shutil
    hb_src = os.environ.get("HAPBIN_SRC", "/home/tx56/hapbin")
    try:
        commit = subprocess.run(["git", "-C", hb_src, "rev-parse", "--short", "HEAD"],
                                capture_output=True, text=True, timeout=30).stdout.strip()
    except OSError:
        commit = ""
    return {
        "selscan": f"{first_line(['selscan', '--help'], 'selscan v')} ({shutil.which('selscan')})",
        "hapbin": f"ihsbin --version: {first_line(['ihsbin', '--version'], 'Version')}; "
                  f"source tree {hb_src} at commit {commit or 'unknown'} "
                  f"(binaries in {os.path.dirname(shutil.which('ihsbin') or '?')}; "
                  f"build provenance not verified)",
    }


# iHS / delihh flavours: label -> extra ihsbin flags. Production is binom
# only; "default" existed for the archive comparison and stays available via
# ihs_flavours(include_default=True) (--ihs-default on 01/02).
ALL_IHS_FLAVOURS = {
    "default": [],
    "binom":   ["--binom"],
}


def ihs_flavours(include_default=False):
    return {k: v for k, v in ALL_IHS_FLAVOURS.items()
            if include_default or k != "default"}


# ==============================================================================
# TPED
# ==============================================================================
def parse_tped(path):
    """-> (ids, positions, daf, n_hap). ids are TPED col 2 (hapbin's 'ID')."""
    ids, positions, daf = [], [], []
    n_hap = 0
    with open(path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 5:
                continue
            alleles = parts[4:]
            if n_hap == 0:
                n_hap = len(alleles)
            ids.append(parts[1])
            positions.append(int(float(parts[3])))
            daf.append(alleles.count("1") / len(alleles))
    return (np.array(ids, dtype=str),
            np.array(positions, dtype=np.int64),
            np.array(daf, dtype=np.float64),
            n_hap)


def write_hap_map(tped, hap_out, map_out=None):
    """Split a TPED into hapbin's .hap (cols 5+) and .map (cols 1-4)."""
    fmap = open(map_out, "w") if map_out else None
    try:
        with open(tped) as fin, open(hap_out, "w") as fhap:
            for line in fin:
                parts = line.split()
                if len(parts) < 5:
                    continue
                if fmap:
                    fmap.write(" ".join(parts[:4]) + "\n")
                fhap.write(" ".join(parts[4:]) + "\n")
    finally:
        if fmap:
            fmap.close()


def parse_and_write_hap_map(tped, hap_out, map_out=None):
    """One pass over a TPED: return (ids, positions, daf, n_hap) AND write
    hapbin's .hap (cols 5+) and optional .map (cols 1-4).

    parse_tped() and write_hap_map() each read the whole file and split every
    line; doing both separately read each TPED twice. Per simulation that was
    8 full reads at 3 populations and 10 at 4 - about 50 GB per gravel batch
    and 173 GB per medina-munoz batch, against 18.8/69 GB of actual data.
    Fusing them removes one read per population per simulation.

    The two remaining reads of the target TPED are selscan's own (--nsl and
    --ihh12) and cannot be removed without changing its input format.
    """
    ids, positions, daf = [], [], []
    n_hap = 0
    fmap = open(map_out, "w") if map_out else None
    try:
        with open(tped) as fin, open(hap_out, "w") as fhap:
            for line in fin:
                parts = line.split()
                if len(parts) < 5:
                    continue
                alleles = parts[4:]
                if n_hap == 0:
                    n_hap = len(alleles)
                ids.append(parts[1])
                positions.append(int(float(parts[3])))
                daf.append(alleles.count("1") / len(alleles))
                if fmap:
                    fmap.write(" ".join(parts[:4]) + "\n")
                fhap.write(" ".join(alleles) + "\n")
    finally:
        if fmap:
            fmap.close()
    return (np.array(ids, dtype=str),
            np.array(positions, dtype=np.int64),
            np.array(daf, dtype=np.float64),
            n_hap)


def flip_tped(src, dst):
    """cosi2 writes ancestral=1 / derived=0; swap so derived is 1."""
    swap = {"0": "1", "1": "0"}
    with open(src) as fin, open(dst, "w") as fout:
        for line in fin:
            parts = line.split()
            if len(parts) < 5:
                continue
            fout.write(" ".join(parts[:4]) + " " +
                       " ".join(swap.get(a, a) for a in parts[4:]) + "\n")


# ==============================================================================
# FST
# ==============================================================================
def weir_hill_fst(p1, n1, p2, n2):
    """n1, n2 are ALLELE counts (haplotypes), not individuals."""
    nt = n1 + n2
    pmean = (n1 * p1 + n2 * p2) / nt
    nc = nt - (n1 * n1 + n2 * n2) / nt
    msp = n1 * (p1 - pmean) ** 2 + n2 * (p2 - pmean) ** 2
    msg = (n1 * p1 * (1 - p1) + n2 * p2 * (1 - p2)) / (n1 - 1 + n2 - 1)

    # Both pops fixed for the same allele: msp = msg = 0, so den = 0 -> NaN,
    # and the pair is left out of mean_fst. This is the CMS convention
    # (cms/model/freqs_stats.c skips it: "if the same allele is fixed in both
    # populations, we skip it"). Until 2026-10-17 these were forced to 0.0,
    # which averaged a spurious 0 into mean_fst at ~15% of sites.
    num, den = msp - msg, msp + (nc - 1) * msg
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(den != 0, num / den, np.nan)


def sfs_table(parsed, target_pop, alt_pops):
    """parsed: {pop_id: (ids, pos, daf, n_hap)} -> DataFrame of SFS stats."""
    ids_t, pos_t, daf_t, nhap_t = parsed[target_pop]

    for p in alt_pops:
        if not np.array_equal(parsed[p][1], pos_t):
            raise ValueError(f"position mismatch: pop {target_pop} vs pop {p}")

    alt_daf = np.column_stack([parsed[p][2] for p in alt_pops])
    fst = np.column_stack([weir_hill_fst(daf_t, nhap_t, parsed[p][2], parsed[p][3])
                           for p in alt_pops])
    # mean over informative pairs only; NaN where no pair is informative
    n_inf = np.isfinite(fst).sum(axis=1)
    mean_fst = np.full(len(fst), np.nan)
    ok = n_inf > 0
    mean_fst[ok] = np.nansum(fst[ok], axis=1) / n_inf[ok]

    return pd.DataFrame({
        "ID": ids_t,
        "pos": pos_t,
        "daf": daf_t,
        "maf": np.minimum(daf_t, 1.0 - daf_t),
        "deldaf": daf_t - alt_daf.mean(axis=1),
        "mean_fst": mean_fst,
    })


# ==============================================================================
# TOOL OUTPUT READERS
# ==============================================================================
def _read_selscan_v3(path, value_col):
    """selscan v3 output: header 'chr id pos freq ... <value_col>'.

    Parse by header name and refuse anything else. v2.0.0 output is
    headerless with no chr column; reading it positionally would silently
    shift every column, and v2.0.0's iHH12 is one-sided (see CLAUDE.md).
    """
    df = pd.read_csv(path, sep="\t", header=0, dtype={"id": str})
    need = {"id", "pos", "freq", value_col}
    if list(df.columns[:4]) != ["chr", "id", "pos", "freq"] or not need <= set(df.columns):
        raise ValueError(f"{path}: not selscan v3 output (header {list(df.columns)})")
    df = df.rename(columns={"id": "ID", "freq": "daf"})
    return df.dropna(subset=["daf", value_col])[["ID", "pos", "daf", value_col]]


def read_nsl(path):
    """selscan v3 --nsl: chr id pos freq ihh1 ihh0 nsl (with header)."""
    return _read_selscan_v3(path, "nsl")


def read_ihh12(path):
    """selscan v3 --ihh12: chr id pos freq ihh12 (with header)."""
    return _read_selscan_v3(path, "ihh12")


def read_ihs_hapbin(path):
    """
    hapbin ihsbin: Index ID Freq iHH_0 iHH_1 iHS 'Std iHS'
    delihh is derived here, SIGNED. hapbin's own 'Std iHS' is ignored: we
    normalise ourselves against our own neutral bins.
    """
    df = pd.read_csv(path, sep="\t", header=0, dtype={"ID": str})
    df = df.rename(columns={"Freq": "daf", "iHS": "ihs"})
    df["delihh"] = df["iHH_1"] - df["iHH_0"]
    return df.dropna(subset=["daf", "ihs", "delihh"])[["ID", "daf", "ihs", "delihh"]]


def read_xpehh(path):
    """hapbin xpehhbin: ID ... Freq ... XPEHH."""
    df = pd.read_csv(path, sep="\t", header=0, dtype={"ID": str})
    df = df.rename(columns={"Freq": "daf", "XPEHH": "xpehh"})
    missing = {"ID", "daf", "xpehh"} - set(df.columns)
    if missing:
        raise ValueError(f"{path}: missing columns {sorted(missing)}")
    df["daf"] = pd.to_numeric(df["daf"], errors="coerce")
    df["xpehh"] = pd.to_numeric(df["xpehh"], errors="coerce")
    return df.dropna(subset=["daf", "xpehh"])[["ID", "daf", "xpehh"]]


# ==============================================================================
# TOOL RUNNERS
# ==============================================================================
def xpehh_on_target_daf(xp, target_ids, target_daf,
                        keep_monomorphic=None):
    """Replace xpehhbin's POOLED Freq with the TARGET pop's DAF, joined on ID.

    Target-monomorphic sites (DAF 0 or 1) are KEPT by default
    (XPEHH_KEEP_TARGET_MONOMORPHIC, decision 2026-10-16) and land in bins
    k = 0 and k = n_hap. With keep_monomorphic=False they are dropped, which
    is how the pre-2026-10-16 bins were built.
    Returns (df[ID, daf, xpehh], n_dropped_monomorphic).
    """
    if keep_monomorphic is None:
        keep_monomorphic = XPEHH_KEEP_TARGET_MONOMORPHIC
    t = pd.DataFrame({"ID": pd.Series(target_ids, dtype=str),
                      "daf_t": np.asarray(target_daf, dtype=float)})
    # A duplicate ID on either side turns this into a Cartesian expansion,
    # which no later step would notice.
    check_unique_ids(t["ID"], "xpehh_on_target_daf: target TPED IDs")
    check_unique_ids(xp["ID"], "xpehh_on_target_daf: xpehhbin output IDs")
    d = safe_merge(xp.drop(columns=["daf"]), t, "ID", "left",
                   "xpehh_on_target_daf")
    if d["daf_t"].isna().any():
        raise ValueError(f"{d['daf_t'].isna().sum()} XP-EHH sites have no "
                         f"target DAF (ID join failed)")
    mono = (d["daf_t"] <= 0) | (d["daf_t"] >= 1)
    if keep_monomorphic:
        mono = pd.Series(False, index=d.index)
    d = d[~mono].rename(columns={"daf_t": "daf"})
    return d[["ID", "daf", "xpehh"]], int(mono.sum())


def check_tools(tools):
    """Verify each tool is not just present but actually RUNS.

    `command -v` plus an executable-bit test does not catch a dynamic-linker
    failure, which is how a whole batch once came out missing four statistics
    while every Slurm job reported success.
    """
    import shutil
    for t in tools:
        if shutil.which(t) is None:
            raise RuntimeError(f"{t} not on PATH")
        try:
            subprocess.run([t, "--help"], capture_output=True, timeout=60)
        except OSError as e:
            raise RuntimeError(
                f"{t} found but not runnable ({e}). "
                f"Check LD_LIBRARY_PATH for hapbin's shared library.")


def run_selscan(tped, out_prefix, mode, threads=1):
    subprocess.run(["selscan", f"--{mode}", "--tped", tped,
                    "--out", out_prefix, "--threads", str(threads),
                    *SELSCAN_FLAGS],
                   check=True, capture_output=True)


def run_ihsbin(hap, mapf, out, extra_flags=()):
    subprocess.run(["ihsbin", "--hap", hap, "--map", mapf, "--out", out,
                    *IHSBIN_FLAGS, *extra_flags],
                   check=True, capture_output=True)


def run_xpehhbin(hapA, hapB, mapf, out):
    subprocess.run(["xpehhbin", "--hapA", hapA, "--hapB", hapB,
                    "--map", mapf, "--out", out, *XPEHHBIN_FLAGS],
                   check=True, capture_output=True)


# ==============================================================================
# BINNING AND NORMALIZATION
# ==============================================================================
def count_from_daf(daf, n_hap):
    """DAF -> exact derived allele count k. Every DAF is k/n_hap by
    construction, so this is a lossless relabelling, not a binning decision."""
    k = np.rint(np.asarray(daf, dtype=float) * n_hap)
    return np.clip(k, 0, n_hap).astype(int)


def build_bins(daf, values, n_hap, out_csv=None):
    """Per-count reference table: one row per derived count k.

    Keyed on k rather than a frequency interval, so a bin holds exactly one
    frequency and there is no within-bin frequency dependence left to correct.
    """
    df = pd.DataFrame({"daf": daf, "v": values}).dropna()
    df = df[(df["daf"] >= 0) & (df["daf"] <= 1)]
    if df.empty:
        return None
    df["k"] = count_from_daf(df["daf"].to_numpy(), n_hap)

    stats = (df.groupby("k")["v"].agg(["mean", "std", "count"]).reset_index())
    stats.loc[stats["count"] < MIN_SITES_PER_BIN, ["mean", "std"]] = np.nan
    stats.insert(1, "n_hap", n_hap)

    if out_csv:
        stats.to_csv(out_csv, index=False)
    return stats


def _bins_to_lut(bins_df, n_hap):
    """Reference table -> (mean, sd) arrays indexed by k, length n_hap + 1."""
    if "k" not in bins_df.columns:
        raise ValueError("bin table is not keyed on k; it predates exact-count "
                         "binning and must be rebuilt")
    stored = bins_df["n_hap"].unique()
    if len(stored) != 1 or int(stored[0]) != int(n_hap):
        raise ValueError(f"bin table built for n_hap={list(stored)} but the data "
                         f"has n_hap={n_hap}; counts are not comparable")
    m = np.full(n_hap + 1, np.nan)
    sd = np.full(n_hap + 1, np.nan)
    k = bins_df["k"].to_numpy(dtype=int)
    ok = (k >= 0) & (k <= n_hap)
    m[k[ok]] = bins_df["mean"].to_numpy()[ok]
    sd[k[ok]] = bins_df["std"].to_numpy()[ok]
    return m, sd


def normalize(values, dafs, bins_df, n_hap):
    """z-score each value against the reference for its own derived count."""
    lut_m, lut_s = _bins_to_lut(bins_df, n_hap)
    idx = count_from_daf(dafs, n_hap)
    m, s = lut_m[idx], lut_s[idx]
    with np.errstate(divide="ignore", invalid="ignore"):
        z = (np.asarray(values, dtype=float) - m) / s
    return np.where(np.isfinite(z) & (s > 0), z, np.nan)


# ==============================================================================
# BATCH ID
# ==============================================================================
def decode_batch_id(batch_id):
    """batch_id[0] = target pop; batch_id[1] = region length in Mb."""
    batch_id = str(batch_id)
    target_pop = batch_id[0]
    region_mb = int(batch_id[1])
    return target_pop, region_mb, region_mb * 1_000_000 // 2


def bins_dir_for(root, model, region_mb, target_pop):
    return os.path.join(root, f"{model}_{region_mb}mb_pop{target_pop}")
