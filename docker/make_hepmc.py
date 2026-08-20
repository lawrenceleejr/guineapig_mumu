#!/usr/bin/env python3
"""Convert a GuineaPig bunch-crossing ("event") output into a HepMC2 ASCII
event file.

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
  --event-number   HepMC event number to record (default: 1).
  --output    Path to write the resulting HepMC2 ASCII file to.

The incoming muon beams are recorded as the two incoming particles of a
single vertex at the origin; every stored photon and pair-production lepton
is attached to that vertex as an outgoing (final-state) particle.  This is
not a full simulation of the underlying physics process, but it gives a
standard HepMC record of the particles GuineaPig produced for one bunch
crossing so it can be read by downstream HEP tools.
"""

import argparse
import math
import random
import os
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

# Default minimum total ENERGY for an outgoing particle, in GeV.  Applied to the
# pair leptons and to the beamstrahlung photons alike.
#
# Be aware of what this does and does not buy you.  GuineaPig's own pair_ecut
# (default 5 MeV) already discards low-energy pairs inside store_pair(), before
# tracking, so on a tracked pairs.dat almost nothing survives below a few MeV.
# Measured on one bunch crossing (560,566 tracked pair leptons, min E 531 keV,
# median E 68 MeV):
#     E < 0.5 MeV :       0  (0.000%)
#     E < 1   MeV :       4  (0.001%)
#     E < 2   MeV :     425  (0.076%)
#     E < 5   MeV :   9,207  (1.642%)
#     E < 10  MeV :  82,949  (14.797%)
# So 2 MeV removes <0.1% of the pairs and will not reduce simulation cost
# measurably.  It does remove essentially all of the beamstrahlung photons from
# photon.dat, which are ~0.1 keV and irrelevant to detector hits anyway.
# For a real reduction see --pt-min and the note there.
DEFAULT_E_MIN = 0.002


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


def build_particles(photon_path, pairs_path, pt_min=DEFAULT_PT_MIN,
                    e_min=DEFAULT_E_MIN):
    """Return (particles, n_pt_cut, n_e_cut).

    particles is a list of (pdgid, px, py, pz, e, mass) outgoing particles.
    The pt cut is applied only to the charged pair leptons: it models containment
    in the beam pipe by the solenoid field, which does not bend photons.  The
    energy cut is applied to every outgoing particle.
    """
    particles = []
    n_cut = 0
    n_ecut = 0
    for energy, vx, vy in read_photons(photon_path):
        px, py, pz, e = photon_four_vector(energy, vx, vy)
        if e < e_min:
            n_ecut += 1
            continue
        particles.append((PHOTON_PDGID, px, py, pz, e, 0.0))

    for energy, vx, vy, vz in read_pairs(pairs_path):
        px, py, pz, e = pair_four_vector(energy, vx, vy, vz)
        if e < e_min:
            n_ecut += 1
            continue
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

    return particles, n_cut, n_ecut


def write_hepmc(output_path, event_number, energy1, energy2, subevents):
    """Write one HepMC2 ASCII file containing one event per entry in subevents.

    Splitting a bunch crossing into N sub-events is what makes the downstream
    Geant4 stage tractable: ddsim processes the events one at a time, so peak
    memory scales with the size of a sub-event rather than the whole crossing.
    A full 10 TeV crossing of tracked pairs (~562k primaries) needs >8 GB and
    over an hour as a single event; in N sub-events it is N cheap events instead.
    The sub-events are a random partition, so each is an unbiased 1/N sample and
    the union is exactly the original crossing -- pair leptons do not interact
    with one another, so summing the hits of the N events reconstructs the
    crossing exactly. The user is expected to merge them downstream.

    The incoming beam muons are repeated in every sub-event because each event
    needs its own vertex; they carry status 4 and are not simulated.
    """
    energy1 = energy1 if energy1 is not None else 0.0
    energy2 = energy2 if energy2 is not None else 0.0

    beam1 = (MUON_PDGID, 0.0, 0.0, energy1, energy1, MUON_MASS)
    beam2 = (-MUON_PDGID, 0.0, 0.0, -energy2, energy2, MUON_MASS)

    with open(output_path, "w") as f:
        f.write("HepMC::Version 2.06.09\n")
        f.write("HepMC::IO_GenEvent-START_EVENT_LISTING\n")

        for j, particles in enumerate(subevents):
            n_particles = 2 + len(particles)
            f.write(
                "E %d -1 -1.0 -1.0 -1.0 0 0 %d 0 0 0 0\n"
                % (event_number + j, 1)
            )
            f.write("U GEV MM\n")
            f.write("V -1 0 0.0 0.0 0.0 0.0 0 %d 0\n" % n_particles)

            barcode = 1
            for pdgid, px, py, pz, e, m in (beam1, beam2):
                f.write(
                    "P %d %d %.9g %.9g %.9g %.9g %.9g 4 0.0 0.0 -1 0\n"
                    % (barcode, pdgid, px, py, pz, e, m)
                )
                barcode += 1

            for pdgid, px, py, pz, e, m in particles:
                f.write(
                    "P %d %d %.9g %.9g %.9g %.9g %.9g 1 0.0 0.0 0 0\n"
                    % (barcode, pdgid, px, py, pz, e, m)
                )
                barcode += 1

        f.write("HepMC::IO_GenEvent-END_EVENT_LISTING\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", help="GuineaPig main output file")
    parser.add_argument("--photons", help="photon.dat file")
    parser.add_argument("--pairs", help="pairs.dat or pairs0.dat file")
    parser.add_argument("--event-number", type=int, default=1)
    parser.add_argument(
        "--pt-min", type=float, default=DEFAULT_PT_MIN, metavar="GEV",
        help="transverse-momentum cut in GeV applied to the charged pair "
             "leptons (default: %(default)s; MAIA inside-beam-pipe value is "
             "0.017). Pass 0 to disable.")
    parser.add_argument(
        "--e-min", type=float, default=DEFAULT_E_MIN, metavar="GEV",
        help="minimum total energy in GeV for any outgoing particle "
             "(default: %(default)s). Pass 0 to disable. NOTE: on a tracked "
             "pairs.dat this removes <0.1%% of the pairs -- see the comment in "
             "this file. It is not a compute-reduction lever.")
    parser.add_argument(
        "--n-subevents", type=int, default=1, metavar="N",
        help="randomly split the bunch crossing into N sub-events written as N "
             "events in the one output file, so the simulation stage treats "
             "them as separate events and peak memory scales with 1/N. The "
             "union of the N events is exactly the original crossing; merge "
             "the hits downstream. (default: %(default)s)")
    parser.add_argument(
        "--seed", type=int, default=0,
        help="RNG seed for the sub-event partition (default: %(default)s)")
    parser.add_argument("--output", required=True, help="HepMC file to write")
    args = parser.parse_args()

    if args.n_subevents < 1:
        parser.error("--n-subevents must be >= 1")
    if args.n_subevents == 1:
        print("WARNING: with --n-subevents 1 this writes a single-event HepMC2 "
              "file. DD4hep's hepmc4 input action mis-handles that case: it "
              "reports EOF, writes an event with no MCParticles and still exits "
              "0. Use --n-subevents >= 2, or convert to HepMC3 ASCII before "
              "ddsim.", file=sys.stderr)

    energy1, energy2 = parse_beam_energies(args.out)
    particles, n_cut, n_ecut = build_particles(
        args.photons, args.pairs, args.pt_min, args.e_min)

    # Random partition, not a contiguous slice: pairs.dat is written by walking a
    # linked list, so its order can correlate with when a particle was stored
    # during tracking and a contiguous slice would not be an unbiased sample.
    random.seed(args.seed)
    random.shuffle(particles)
    n = args.n_subevents
    subevents = [particles[i::n] for i in range(n)]

    write_hepmc(args.output, args.event_number, energy1, energy2, subevents)

    sizes = [len(se) for se in subevents]
    print(
        "Wrote %d outgoing particle(s) to %s in %d sub-event(s) of %d-%d "
        "(pt cut %g GeV removed %d, energy cut %g GeV removed %d)"
        % (len(particles), args.output, n, min(sizes), max(sizes),
           args.pt_min, n_cut, args.e_min, n_ecut)
    )


if __name__ == "__main__":
    sys.exit(main())
