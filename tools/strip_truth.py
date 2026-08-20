#!/usr/bin/env python3
"""Drop the MCParticle collection from an edm4hep simulation file.

ddsim has no option to skip writing truth: DDSim's outputConfig exposes only
forceDD4HEP / forceEDM4HEP / forceLCIO / useRNTuple / userOutputPlugin, and
SIM.part.minimalKineticEnergy only *shrinks* the collection. So this rewrites the
file with the MCParticle collection removed.

Why bother: the MCParticle collection, not the detector hits, is what makes these
files large. Measured with TBranch::GetZipBytes("*"), which accounts for
99.7-100% of the file on disk:

    pT>16 MeV incoherent-pair sample, 26.5 MB/event
        MC truth (particles + relations)   99.37%
        all detector hits                   0.63%

    unfiltered incoherent-pair sample, 223.3 MB/event
        MC truth                           61.51%
        detector hits (+contributions)     38.49%

Caveat: the SimTrackerHit / CaloHitContribution -> MCParticle references are kept
as stored indices but no longer resolve, so anything that walks from a hit back to
a particle stops working. Use this for occupancy / overlay-input files where only
the hits matter, not for truth-matching studies.

Run inside an environment with podio + edm4hep available:

  python3 strip_truth.py in.edm4hep.root out.edm4hep.root
"""
import gc
import os
import sys

from podio import root_io


def _rewrite(inpath, outpath, drop):
    """Do the copy in its own scope so the Writer is destroyed -- and therefore
    flushed and closed -- before the caller stats the output. podio's Writer has
    no finish() in all builds, so scope exit is the portable way to close it."""
    reader = root_io.Reader(inpath)
    writer = root_io.Writer(outpath)
    categories = list(reader.categories)
    n_written = 0
    dropped_names = set()
    for cat in categories:
        for frame in reader.get(cat):
            present = [c for c in drop if c in frame.getAvailableCollections()]
            dropped_names.update(present)
            # keep everything except the truth collections
            keep = [c for c in frame.getAvailableCollections() if c not in present]
            writer.write_frame(frame, cat, keep)
            n_written += 1
    return categories, n_written, dropped_names


def main(inpath, outpath, drop=("MCParticles", "MCParticle")):
    categories, n_written, dropped_names = _rewrite(inpath, outpath, drop)
    gc.collect()
    print(f"{inpath} -> {outpath}")
    print(f"  frames written : {n_written} across {len(categories)} category(ies)")
    print(f"  dropped        : {sorted(dropped_names) or 'nothing (no MCParticle collection found)'}")
    print(f"  input size     : {os.path.getsize(inpath)/1e6:,.2f} MB")
    # The output is only flushed when podio's Writer is finalised at interpreter
    # exit, so its size cannot be read reliably from inside this process --
    # stat it after the command returns.


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
