#!/usr/bin/env python3
"""Convert a GuineaPig bunch-crossing ("event") output into a HepMC3 ASCII
(Asciiv3) event file.

Inputs (all paths are optional except the main GuineaPig output file):
  --out       Main GuineaPig output file (e.g. muon_pairs_10tev_event1.out).
              Used to read the two beam energies (energy.1 / energy.2).
  --photons   photon.dat produced when store_photons/write_photons is set.
              Each line is "energy vx vy" with vx, vy the photon direction
              slopes (px/pz, py/pz) and a positive/negative energy sign used
              to indicate which beam the photon came from.
  --pairs     pairs.dat or pairs0.dat produced when track_pairs>0 or
              store_pairs>1. Each line is "energy vx vy vz [x y z]" with
              vx, vy, vz the normalised momentum direction (px/E, py/E, pz/E)
              of a pair-production lepton.
  --event-number   GuineaPig bunch-crossing number to record (default: 1).
  --sub-events     Split the outgoing particles into this many random HepMC
                   sub-events (default: 1).
  --output    Path to write the resulting HepMC3 ASCII file to.

The incoming muon beams are recorded as the two incoming particles of a
single vertex at the origin; every stored photon and pair-production lepton
is attached to that vertex as an outgoing (final-state) particle.  This is
not a full simulation of the underlying physics process, but it gives a
standard HepMC record of the particles GuineaPig produced for one bunch
crossing so it can be read by downstream HEP tools.
"""

import argparse
import math
import os
import random
import re
import sys

MUON_MASS = 0.1056583715  # GeV
ELECTRON_MASS = 0.00051099895  # GeV
MUON_PDGID = 13
ELECTRON_PDGID = 11
PHOTON_PDGID = 22

# Which sign of the energy field in pairs0.dat/pairs.dat corresponds to the
# electron.  GuineaPig guarantees the two members of a pair have opposite
# signs but the absolute species assignment is a convention; flip this single
# flag if the downstream convention differs.
POSITIVE_ENERGY_IS_ELECTRON = True

# Default transverse-momentum cut applied to the charged pair-production
# leptons, in GeV.  GuineaPig itself applies no pt cut: store_pair() cuts on
# ENERGY (pair_ecut, default 5 MeV) and the pt>20 MeV / theta>0.15 test in
# background.c only feeds the printed anzahl_* counters, never the stored
# output.  So without this the HepMC keeps pairs down to ~1 keV of pt.
#
# For reference, the MAIA inside-beam-pipe value is 17 MeV: below that a
# charged particle's helix in the solenoid field keeps it inside the beam
# pipe, so it never reaches the detector.  The default here is 15 MeV,
# slightly looser so that the beam-pipe threshold can still be applied
# downstream without having to regenerate the sample.
#
# How much this cut costs depends entirely on WHICH pair file it is applied to:
#   pairs0.dat (production time) : pT>15 MeV keeps ~1% of pairs
#   pairs.dat  (tracked)         : pT>15 MeV keeps ~41% of pairs
# On the tracked sample this is therefore a consequential physics choice, not a
# cheap cleanup -- it removes most of the pairs that actually reach the detector.
# It is also NOT a useful way to shrink the downstream ddsim output: that file is
# dominated by the MCParticle collection (>99% of bytes for a cut sample), not by
# hits, so SIM.part.minimalKineticEnergy is the lever for file size.  Pass 0 to
# disable this cut.
DEFAULT_PT_MIN = 0.015


def parse_beam_energies(out_path):
    """Read energy.1 / energy.2 (GeV) from a GuineaPig output file."""
    energy1 = energy2 = None
    if out_path and os.path.isfile(out_path):
        with open(out_path) as f:
            text = f.read()
        m1 = re.search(r"energy\.1\s*=\s*([-0-9.eE+]+)", text)
        m2 = re.search(r"energy\.2\s*=\s*([-0-9.eE+]+)", text)
        if m1:
            energy1 = float(m1.group(1))
        if m2:
            energy2 = float(m2.group(1))
    return energy1, energy2


def read_photons(path):
    """Yield (energy, vx, vy) tuples from a photon.dat file."""
    if not path or not os.path.isfile(path):
        return
    with open(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 3:
                continue
            yield float(parts[0]), float(parts[1]), float(parts[2])


def read_pairs(path):
    """Yield (energy, vx, vy, vz) tuples from a pairs.dat/pairs0.dat file."""
    if not path or not os.path.isfile(path):
        return
    with open(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 4:
                continue
            yield (float(parts[0]), float(parts[1]),
                   float(parts[2]), float(parts[3]))


def photon_four_vector(energy, vx, vy):
    e = abs(energy)
    pz = e / math.sqrt(1.0 + vx * vx + vy * vy)
    px = vx * pz
    py = vy * pz
    return px, py, pz, e


def pair_four_vector(energy, vx, vy, vz):
    e = abs(energy)
    px = vx * e
    py = vy * e
    pz = vz * e
    return px, py, pz, e


def build_particles(photon_path, pairs_path, pt_min=DEFAULT_PT_MIN):
    """Return (particles, n_pairs_cut).

    particles is a list of (pdgid, px, py, pz, e, mass) outgoing particles.
    The pt cut is applied only to the charged pair leptons: the cut models
    containment in the beam pipe by the solenoid field, which does not bend
    photons, so beamstrahlung photons are always kept.
    """
    particles = []
    n_cut = 0
    for energy, vx, vy in read_photons(photon_path):
        px, py, pz, e = photon_four_vector(energy, vx, vy)
        particles.append((PHOTON_PDGID, px, py, pz, e, 0.0))

    for energy, vx, vy, vz in read_pairs(pairs_path):
        px, py, pz, e = pair_four_vector(energy, vx, vy, vz)
        if math.hypot(px, py) < pt_min:
            n_cut += 1
            continue
        # GuineaPig records the charge in the SIGN of the energy: store_full_pair()
        # in background.c negates the second member of every pair ("e2 = -e2"), so
        # the two members always carry opposite signs.  Index parity is NOT usable:
        # pair_ecut and pair_ratio drop members individually inside store_pair(),
        # which desynchronises the alternation (measured agreement with the true
        # sign on a real pairs0.dat: 49.5%, i.e. a coin flip).
        pdgid = ELECTRON_PDGID if (energy > 0.0) == POSITIVE_ENERGY_IS_ELECTRON \
            else -ELECTRON_PDGID
        particles.append((pdgid, px, py, pz, e, ELECTRON_MASS))

    return particles, n_cut


def split_particles(particles, n_sub_events, rng):
    """Shuffle particles and split them into near-equal random chunks."""
    shuffled = list(particles)
    rng.shuffle(shuffled)
    n_particles = len(shuffled)
    base_size, remainder = divmod(n_particles, n_sub_events)
    chunks = []
    start = 0
    for i in range(n_sub_events):
        size = base_size + (1 if i < remainder else 0)
        chunks.append(shuffled[start:start + size])
        start += size
    return chunks


def write_hepmc(output_path, first_event_number, energy1, energy2, particle_groups):
    energy1 = energy1 if energy1 is not None else 0.0
    energy2 = energy2 if energy2 is not None else 0.0

    beam1 = (MUON_PDGID, 0.0, 0.0, energy1, energy1, MUON_MASS)
    beam2 = (-MUON_PDGID, 0.0, 0.0, -energy2, energy2, MUON_MASS)

    with open(output_path, "w") as f:
        f.write("HepMC::Version 3.02.05\n")
        f.write("HepMC::Asciiv3-START_EVENT_LISTING\n")
        for event_offset, particles in enumerate(particle_groups):
            event_number = first_event_number + event_offset
            n_out = len(particles)
            n_particles = 2 + n_out
            n_vertices = 1
            vertex_id = -1

            f.write("E %d %d %d\n" % (event_number, n_vertices, n_particles))
            f.write("U GEV MM\n")

            barcode = 1
            beam_ids = []
            for pdgid, px, py, pz, e, m in (beam1, beam2):
                f.write(
                    "P %d 0 %d %.9g %.9g %.9g %.9g %.9g 4\n"
                    % (barcode, pdgid, px, py, pz, e, m)
                )
                beam_ids.append(barcode)
                barcode += 1

            f.write(
                "V %d 0 [%s]\n" % (vertex_id, ",".join(str(i) for i in beam_ids))
            )

            for pdgid, px, py, pz, e, m in particles:
                f.write(
                    "P %d %d %d %.9g %.9g %.9g %.9g %.9g 1\n"
                    % (barcode, vertex_id, pdgid, px, py, pz, e, m)
                )
                barcode += 1

        f.write("HepMC::Asciiv3-END_EVENT_LISTING\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", help="GuineaPig main output file")
    parser.add_argument("--photons", help="photon.dat file")
    parser.add_argument("--pairs", help="pairs.dat or pairs0.dat file")
    parser.add_argument("--event-number", type=int, default=1)
    parser.add_argument(
        "--sub-events", type=int, default=1, metavar="N",
        help="split the outgoing particles into N random HepMC sub-events "
             "(default: %(default)s)")
    parser.add_argument(
        "--pt-min", type=float, default=DEFAULT_PT_MIN, metavar="GEV",
        help="transverse-momentum cut in GeV applied to the charged pair "
             "leptons (default: %(default)s; MAIA inside-beam-pipe value is "
             "0.017). Pass 0 to disable.")
    parser.add_argument("--output", required=True, help="HepMC file to write")
    args = parser.parse_args()
    if args.sub_events < 1:
        parser.error("--sub-events must be a positive integer")

    energy1, energy2 = parse_beam_energies(args.out)
    particles, n_cut = build_particles(args.photons, args.pairs, args.pt_min)
    particle_groups = split_particles(
        particles, args.sub_events, random.Random(args.event_number)
    )
    first_event_number = (args.event_number - 1) * args.sub_events + 1
    write_hepmc(
        args.output, first_event_number, energy1, energy2, particle_groups
    )

    print(
        "Wrote %d outgoing particle(s) across %d HepMC event(s) to %s "
        "(pt cut %g GeV removed %d pair lepton(s))"
        % (
            len(particles), args.sub_events, args.output,
            args.pt_min, n_cut
        )
    )


if __name__ == "__main__":
    sys.exit(main())
