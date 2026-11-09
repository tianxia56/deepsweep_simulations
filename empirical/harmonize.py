#!/usr/bin/env python3
"""Stream-harmonise a (target, comparison) population pair for xpehhbin.

One pass over each of the four input files, no accumulation: both maps are
sorted ascending on physical position with no duplicates (verified), so the
intersection is a merge-join. The old pairwise_harmonize.py held every kept
hap line of a chromosome in a Python dict, which is ~10 GB of str objects for
chr2; this holds one line per population.

Target DAF is computed in the same pass rather than in a second sweep over the
tped. hap and tped genotype fields were verified identical, and the pipeline
that produced them codes 1 = derived (unified_annotation_pipeline.py:95), the
same convention as the simulations.

Writes:
  out_map     joint map, target's lines, intersection only
  out_hap_a   target haplotypes, intersection only
  out_hap_b   comparison haplotypes, intersection only
  out_sites   npz: pos (int64), k (int32, target derived count), n_hap (int)
              row i corresponds to xpehhbin's Index i
"""
import sys
import numpy as np


def fields(map_line):
    p = map_line.split()
    return int(p[3])


def harmonise(map_a, hap_a, map_b, hap_b, out_map, out_hap_a, out_hap_b,
              out_sites):
    pos_out = []
    k_out = []
    n_hap_a = None

    with open(map_a) as fma, open(hap_a) as fha, \
         open(map_b) as fmb, open(hap_b) as fhb, \
         open(out_map, "w") as oma, open(out_hap_a, "w") as oha, \
         open(out_hap_b, "w") as ohb:

        la, ha = fma.readline(), fha.readline()
        lb, hb = fmb.readline(), fhb.readline()
        prev = -1

        while la and lb:
            pa, pb = fields(la), fields(lb)
            if pa < pb:
                la, ha = fma.readline(), fha.readline()
            elif pb < pa:
                lb, hb = fmb.readline(), fhb.readline()
            else:
                if pa <= prev:
                    raise ValueError(f"{map_a}: positions not strictly "
                                     f"increasing at {pa}")
                prev = pa
                if n_hap_a is None:
                    n_hap_a = len(ha.split())
                oma.write(la)
                oha.write(ha)
                ohb.write(hb)
                pos_out.append(pa)
                k_out.append(ha.count("1"))
                la, ha = fma.readline(), fha.readline()
                lb, hb = fmb.readline(), fhb.readline()

    if not pos_out:
        sys.exit(2)

    # an open handle, not a path: numpy appends .npz to a path that lacks it
    with open(out_sites, "wb") as fh:
        np.savez_compressed(fh,
                            pos=np.asarray(pos_out, dtype=np.int64),
                            k=np.asarray(k_out, dtype=np.int32),
                            n_hap=np.int32(n_hap_a))


if __name__ == "__main__":
    harmonise(*sys.argv[1:9])
